"""
Enhanced doublet detection for single-cell RNA-seq data.

This module provides comprehensive functions for identifying potential doublet cells
using an algorithmic method (Scrublet) and a flexible heuristic approach based on
mutually exclusive lineage marker co-expression.
"""

import gc
import logging
import random
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
import seaborn as sns
from anndata import AnnData
from upsetplot import plot as upset_plot

from ..utils.marker_manager import get_marker_manager
from .config import DoubletConfig, MarkerConfig

log = logging.getLogger(__name__)

# --- Use constants for column names for easier maintenance ---
HEURISTIC_PRED_COL = "heuristic_predicted"
FINAL_PRED_COL = "predicted_doublet"
LINEAGE_PRED_KEY = "lineage_predictions"

__all__ = [
    "generate_doublet_rates",
    "create_custom_marker_dict",
    "run_heuristic_analysis",
    "analyze_lineage_coexpression",
    "predict_doublets",
]


# --- Helper Functions ---
def _create_doublet_marker_config_from_manager(
    adata: AnnData, cfg: DoubletConfig
) -> Dict[str, MarkerConfig]:
    """
    Create MarkerConfig objects using the marker manager and the main DoubletConfig.

    Args:
        adata: The AnnData object to intersect markers with.
        cfg: The main DoubletConfig object, providing context like species, tissue,
             and default evaluation parameters.

    Returns:
        A dictionary mapping lineage names to MarkerConfig objects.
    """
    try:
        manager = get_marker_manager(
            species=cfg.marker_species,
            tissue=cfg.marker_tissue,
            states=None,  # States can be added in future versions if needed
            case_sensitive=False,
        )
        manager.intersect_with(adata)
        markers_dict = manager.get_doublet_lineage_markers()

        if not markers_dict:
            log.warning(
                "No tissue-specific lineage markers found, falling back to base markers."
            )
            base_manager = get_marker_manager(
                species=cfg.marker_species, case_sensitive=False
            )
            base_manager.intersect_with(adata)
            markers_dict = base_manager.get_doublet_lineage_markers()

    except Exception as e:
        log.error(
            f"Failed to load markers via marker_manager: {e}. Cannot create heuristic configs."
        )
        return {}

    marker_configs = {}
    for lineage, genes in markers_dict.items():
        if genes:
            marker_configs[lineage] = MarkerConfig(
                genes=genes,
                expression_threshold=cfg.default_expression_threshold,
                min_genes_required=cfg.default_min_genes_required,
                use_raw=cfg.default_use_raw,
            )
    log.info(
        f"Auto-generated {len(marker_configs)} marker configurations for doublet detection."
    )
    return marker_configs


def _evaluate_lineage_expression(
    adata: AnnData, lineage_name: str, marker_config: MarkerConfig
) -> pd.Series:
    """
    Evaluate expression of markers for a specific lineage based on a MarkerConfig.

    Args:
        adata: AnnData object containing expression data.
        lineage_name: Name of the lineage being evaluated.
        marker_config: Configuration object defining marker evaluation parameters.

    Returns:
        A boolean pandas Series indicating cells positive for this lineage.
    """
    source_adata = (
        adata.raw.to_adata() if marker_config.use_raw and adata.raw else adata
    )
    var_names_upper = {name.upper(): name for name in source_adata.var_names}

    if marker_config.is_regex:
        pattern = marker_config.genes
        matching_genes = source_adata.var_names.str.contains(
            pattern, regex=True, na=False, case=False
        )
        valid_genes = source_adata.var_names[matching_genes].tolist()
    else:
        marker_genes_upper = [g.upper() for g in marker_config.genes]
        valid_genes_upper = [g for g in marker_genes_upper if g in var_names_upper]
        valid_genes = [var_names_upper[g] for g in valid_genes_upper]

    if not valid_genes:
        log.warning(f"No valid genes found for lineage '{lineage_name}'")
        return pd.Series(False, index=adata.obs_names)

    expr_data = source_adata[:, valid_genes].X
    if hasattr(expr_data, "toarray"):
        expr_data = expr_data.toarray()

    expr_binary = expr_data > marker_config.expression_threshold
    genes_expressed = expr_binary.sum(axis=1)
    lineage_positive = genes_expressed >= marker_config.min_genes_required
    lineage_positive_series = pd.Series(lineage_positive, index=source_adata.obs_names)

    positive_count = lineage_positive_series.sum()
    log.info(
        f"Lineage '{lineage_name}': {positive_count} cells positive "
        f"({positive_count / len(lineage_positive_series) * 100:.2f}%)"
    )
    return lineage_positive_series


def _merge_doublet_predictions(
    adata: AnnData,
    algorithm_col: str,
    heuristic_col: str,
    strategy: str = "weighted",
    algorithm_weight: float = 0.7,
    random_state: int = 61,
) -> pd.Series:
    """
    Merge algorithmic and heuristic doublet predictions using different strategies.

    Args:
        adata: AnnData object containing prediction results
        algorithm_col: Column name for algorithmic predictions
        heuristic_col: Column name for heuristic predictions
        strategy: Merge strategy ("union", "weighted, "intersection", "algorithm_priority", "heuristic_priority")
        algorithm_weight: Weight for algorithmic predictions (0-1)

    Returns:
        Boolean series with merged predictions
    """
    algo_pred = adata.obs[algorithm_col].fillna(False).astype(bool)
    heur_pred = adata.obs[heuristic_col].fillna(False).astype(bool)

    if strategy == "union":
        merged = algo_pred | heur_pred
    elif strategy == "intersection":
        merged = algo_pred & heur_pred
    elif strategy == "algorithm_priority":
        merged = algo_pred
    elif strategy == "heuristic_priority":
        merged = heur_pred
    elif strategy == "weighted":
        # Set random seed for reproducible results
        random.seed(random_state)
        agree = algo_pred == heur_pred
        disagree = ~agree

        # Start with points of agreement
        merged = pd.Series(False, index=adata.obs_names)
        merged[agree] = algo_pred[agree]  # Use agreement value

        # For disagreements, use algorithm with probability algorithm_weight
        if np.any(disagree):
            disagree_cells_indices = np.where(disagree)[0]
            disagree_cells_obs_names = adata.obs_names[disagree_cells_indices]

            score_col = algorithm_col.replace("_predicted", "_score")
            if score_col in adata.obs:
                disagree_scores = adata.obs.loc[disagree_cells_obs_names, score_col]

                n_disagree = len(disagree_cells_obs_names)
                n_algo = int(n_disagree * algorithm_weight)

                top_cells = disagree_scores.nlargest(n_algo).index
                other_cells = disagree_cells_obs_names.difference(top_cells)

                merged[top_cells] = algo_pred[top_cells]
                merged[other_cells] = heur_pred[other_cells]
            else:
                log.warning(
                    f"Doublet score column '{score_col}' not found. Using random assignment for 'weighted' strategy."
                )
                n_disagree = disagree.sum()
                n_algo = int(n_disagree * algorithm_weight)

                disagree_obs_names_list = list(adata.obs_names[disagree])
                algo_cells = random.sample(disagree_obs_names_list, n_algo)
                other_cells = set(disagree_obs_names_list) - set(algo_cells)

                merged.loc[algo_cells] = algo_pred.loc[algo_cells]
                merged.loc[list(other_cells)] = heur_pred.loc[list(other_cells)]
    else:
        raise ValueError(f"Unknown merge strategy: {strategy}")

    log.info(
        f"Merged predictions using '{strategy}' strategy. Final count: {merged.sum()}"
    )
    return merged


def _run_scrublet(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run Scrublet algorithm for doublet detection on a single AnnData view.
    Returns (scores, predicted) arrays.
    """
    rate = config.expected_doublet_rate
    current_rate = rate.get(sample_name, 0.1) if isinstance(rate, dict) else rate
    if current_rate is None:  # Handle case where rate is not provided for a sample
        log.warning(
            f"No doublet rate provided for sample '{sample_name}', using default of 0.1."
        )
        current_rate = 0.1

    actual_n_pcs = min(config.n_pcs, adata_view.n_obs - 1, adata_view.n_vars - 1)

    try:
        scrub = scr.Scrublet(adata_view.X, expected_doublet_rate=current_rate)
        scores, _ = scrub.scrub_doublets(n_prin_comps=actual_n_pcs, verbose=False)
        predicted = scrub.call_doublets(verbose=False)

        doublet_count = sum(predicted)
        doublet_rate = doublet_count / len(predicted)
        log.info(
            f"  Found {doublet_count} potential doublets via Scrublet ({doublet_rate:.2%})"
        )

        if config.plot_umap:
            try:
                scrub.set_embedding(
                    "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
                )
                fig, ax = scrub.plot_embedding("UMAP", order_points=True)
                if config.save_dir:
                    save_path = (
                        Path(config.save_dir) / f"{sample_name}_doublets_umap.png"
                    )
                    save_path.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(save_path, dpi=300, bbox_inches="tight")
                if config.show_plots:
                    plt.show()
                else:
                    plt.close(fig)
            except Exception as e:
                log.warning(f"Could not generate UMAP for sample {sample_name}: {e}")

        return scores, predicted

    except Exception as e:
        log.error(f"Scrublet failed for sample {sample_name}: {e}")
        return None, None
    finally:
        gc.collect()


def _run_heuristic(
    adata: AnnData, cfg: DoubletConfig
) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Runs the full heuristic doublet detection workflow.

    This helper function loads markers, identifies co-expressing cells,
    and applies ignore rules. It returns both the final boolean prediction
    and the detailed lineage-by-cell matrix.

    Args:
        adata: AnnData object.
        cfg: The DoubletConfig object.

    Returns:
        A tuple containing:
        - A boolean pd.Series of heuristic doublet predictions.
        - A pd.DataFrame (`lineage_df`) with detailed predictions for each lineage.
    """
    marker_configs = cfg.marker_configs
    if marker_configs is None:
        marker_configs = _create_doublet_marker_config_from_manager(adata, cfg)

    if not marker_configs:
        log.warning(
            "Heuristics enabled, but no marker configurations were found. Skipping."
        )
        empty_df = pd.DataFrame(index=adata.obs_names)
        return pd.Series(False, index=adata.obs_names), empty_df

    lineage_results = {}
    for name, marker_config in marker_configs.items():
        # If a marker_config from a manual dict is missing a parameter, it will use the default from DoubletConfig
        final_mc = MarkerConfig(
            genes=marker_config.genes,
            expression_threshold=getattr(
                marker_config, "expression_threshold", cfg.default_expression_threshold
            ),
            min_genes_required=getattr(
                marker_config, "min_genes_required", cfg.default_min_genes_required
            ),
            use_raw=getattr(marker_config, "use_raw", cfg.default_use_raw),
        )
        lineage_results[name] = _evaluate_lineage_expression(adata, name, final_mc)

    lineage_df = pd.DataFrame(lineage_results, index=adata.obs_names)

    # Add prevalence filtering for lineages
    # Only consider lineages that are present in a minimum percentage of cells
    min_prevalence = cfg.min_lineage_prevalence
    lineage_prevalence = lineage_df.mean()
    valid_lineages = lineage_prevalence[
        lineage_prevalence >= min_prevalence
    ].index.tolist()

    if len(valid_lineages) < 2:
        log.warning(
            "Insufficient lineages with minimum prevalence. Heuristic analysis may be unreliable."
        )
        # Fall back to using all lineages if needed
        if len(lineage_df.columns) >= 2:
            valid_lineages = lineage_df.columns.tolist()

    # Focus on valid lineages only
    lineage_df = lineage_df[valid_lineages]

    # Calculate lineages per cell with improved logic
    lineages_per_cell = lineage_df.sum(axis=1)
    potential_doublets = lineages_per_cell >= cfg.min_lineages_for_doublet

    if cfg.ignore_coexpression_pairs:
        log.info(
            f"Ignoring {len(cfg.ignore_coexpression_pairs)} specific co-expression pairs."
        )
        ignored_set = {tuple(sorted(pair)) for pair in cfg.ignore_coexpression_pairs}
        doublet_indices = potential_doublets[potential_doublets].index
        for cell_idx in doublet_indices:
            expressed_lineages = tuple(
                sorted(lineage_df.columns[lineage_df.loc[cell_idx]])
            )
            if expressed_lineages in ignored_set:
                potential_doublets.loc[cell_idx] = False

    return potential_doublets, lineage_df


def _export_doublet_stats(
    adata: AnnData,
    sample_key: str = "sampleID",
    save_dir: Optional[Union[str, Path]] = None,
    export_csv: bool = True,
    export_xlsx: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Export comprehensive doublet statistics per sample and globally.

    This function generates detailed statistical summaries of doublet detection
    results, including counts, percentages, and score distributions.

    Args:
        adata: AnnData object with doublet predictions
        sample_key: Key for sample identification
        save_dir: Directory to save statistics files
        export_csv: Whether to export as CSV files
        export_xlsx: Whether to export as Excel file

    Returns:
        Dictionary containing sample-wise and global statistics DataFrames
    """
    # Identify all doublet-related columns
    doublet_cols = [
        col
        for col in adata.obs.columns
        if any(
            keyword in col.lower() for keyword in ["doublet", "scrublet", "heuristic"]
        )
    ]

    if not doublet_cols:
        log.warning("No doublet-related columns found in adata.obs")
        return {}

    log.info(f"Found doublet columns: {doublet_cols}")

    # Calculate per-sample statistics
    sample_stats = []
    unique_samples = adata.obs[sample_key].unique()
    if not pd.api.types.is_categorical_dtype(adata.obs[sample_key]):
        unique_samples = sorted(unique_samples)

    for sample in unique_samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs.loc[sample_mask]

        stats = {
            "sample": sample,
            "total_cells": len(sample_data),
        }

        for col in doublet_cols:
            if col in sample_data.columns:
                col_data = sample_data[col].dropna()
                if (
                    pd.api.types.is_numeric_dtype(col_data)
                    and not pd.api.types.is_bool_dtype(col_data)
                    and col_data.nunique() > 2
                ):
                    # Continuous column (scores)
                    stats[f"{col}_mean"] = col_data.mean()
                    stats[f"{col}_median"] = col_data.median()
                    stats[f"{col}_std"] = col_data.std()
                elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                    # Boolean/binary column (predictions)
                    positive_count = col_data.astype(bool).sum()
                    stats[f"{col}_count"] = positive_count
                    stats[f"{col}_percentage"] = (
                        (positive_count / len(sample_data) * 100)
                        if len(sample_data) > 0
                        else 0
                    )

        sample_stats.append(stats)

    sample_df = pd.DataFrame(sample_stats).set_index("sample")

    global_stats = {"metric": "global", "total_cells": adata.n_obs}
    for col in doublet_cols:
        if col in adata.obs.columns:
            col_data = adata.obs[col].dropna()
            if (
                pd.api.types.is_numeric_dtype(col_data)
                and not pd.api.types.is_bool_dtype(col_data)
                and col_data.nunique() > 2
            ):
                global_stats[f"{col}_mean"] = col_data.mean()
                global_stats[f"{col}_median"] = col_data.median()
                global_stats[f"{col}_std"] = col_data.std()
            elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                positive_count = col_data.astype(bool).sum()
                global_stats[f"{col}_count"] = positive_count
                global_stats[f"{col}_percentage"] = (
                    (positive_count / adata.n_obs * 100) if adata.n_obs > 0 else 0
                )

    global_df = pd.DataFrame([global_stats])

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        if export_csv:
            sample_df.to_csv(save_dir / "doublet_stats_per_sample.csv")
            global_df.to_csv(save_dir / "doublet_stats_global.csv", index=False)
            log.info(f"Exported CSV files to {save_dir}")
        if export_xlsx:
            with pd.ExcelWriter(save_dir / "doublet_stats.xlsx") as writer:
                sample_df.to_excel(writer, sheet_name="per_sample")
                global_df.to_excel(writer, sheet_name="global", index=False)
            log.info(f"Exported Excel file to {save_dir / 'doublet_stats.xlsx'}")

    return {"sample": sample_df, "global": global_df}


def _plot_doublet_summary(
    adata: AnnData,
    sample_key: str = "sampleID",
    save_dir: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> None:
    """
    Generates a comprehensive, UMAP-independent summary plot for doublet detection results.

    This function creates a multi-panel figure showing:
    1. A stacked bar plot of doublet counts and percentages per sample, broken down by prediction source.
    2. A scatter plot of doublet scores vs. gene counts, a key diagnostic for algorithmic methods.
    3. An UpSet plot visualizing the co-expression patterns from the heuristic analysis.

    Args:
        adata: AnnData object after running `predict_doublets`.
        sample_key: The key in adata.obs for sample identification.
        save_dir: Directory to save the plot.
        show: Whether to display the plot.
    """
    log.info("Generating UMAP-independent doublet detection summary plot...")

    # --- 0. Check for required columns ---
    required_cols = [FINAL_PRED_COL]
    algo_pred_col_list = [
        c
        for c in adata.obs.columns
        if c.endswith("_predicted") and "heuristic" not in c
    ]
    if not algo_pred_col_list:
        log.warning("No algorithm prediction column found. Skipping summary plot.")
        return
    algo_pred_col = algo_pred_col_list[0]
    algo_score_col = algo_pred_col.replace("_predicted", "_score")
    required_cols.append(algo_pred_col)

    # Heuristic column is optional if use_heuristics is False
    if HEURISTIC_PRED_COL not in adata.obs.columns:
        adata.obs[HEURISTIC_PRED_COL] = False

    if not all(col in adata.obs.columns for col in required_cols):
        raise ValueError("Required columns not found. Run `predict_doublets` first.")

    # Prepare doublet source column for plotting
    source_categories = [
        "Heuristic Only",
        "Algorithm Only",
        "Both Methods",
        "Singleton",
    ]
    conditions = [
        adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL],
        adata.obs[algo_pred_col] & ~adata.obs[HEURISTIC_PRED_COL],
        ~adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL],
    ]
    choices = ["Both Methods", "Algorithm Only", "Heuristic Only"]
    adata.obs["doublet_source"] = np.select(conditions, choices, default="Singleton")
    adata.obs["doublet_source"] = pd.Categorical(
        adata.obs["doublet_source"], categories=source_categories, ordered=True
    )

    # --- Figure 1: Overview (Bar Plot + Scatter Plot) ---
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5), facecolor="white")
    fig1.suptitle("Doublet Detection Summary", fontsize=18, fontweight="bold")

    # Panel 1: Doublet Breakdown Bar Plot
    summary_df = (
        adata.obs.groupby(sample_key)["doublet_source"]
        .value_counts(normalize=True, sort=False)
        .unstack(fill_value=0)
        * 100
    )
    summary_df[source_categories].plot(
        kind="bar",
        stacked=True,
        ax=ax1,
        color={
            "Singleton": "lightgray",
            "Heuristic Only": "coral",
            "Algorithm Only": "skyblue",
            "Both Methods": "darkred",
        },
        width=0.8,
    )
    ax1.set_title("Doublet Breakdown per Sample", fontsize=14)
    ax1.set_ylabel("Percentage of Cells (%)")
    ax1.set_xlabel("Sample")
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right")
    ax1.legend(title="Source", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax1.grid(axis="y", linestyle="--", alpha=0.7)

    # Panel 2: Doublet Score vs. Gene Count Scatter Plot
    if algo_score_col in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
        sns.scatterplot(
            data=adata.obs,
            x="n_genes_by_counts",
            y=algo_score_col,
            hue=FINAL_PRED_COL,
            s=5,
            alpha=0.6,
            palette={True: "red", False: "gray"},
            ax=ax2,
            rasterized=True,
        )
        ax2.set_title(
            f"{algo_score_col.replace('_', ' ').title()} vs. Gene Count", fontsize=14
        )
        ax2.set_xlabel("Number of Genes")
        ax2.set_ylabel("Doublet Score")
        ax2.legend(title="Final Prediction")
    else:
        ax2.text(
            0.5,
            0.5,
            "Score or gene count data not available.",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )

    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        overview_save_path = save_path / "doublet_summary_overview.png"
        fig1.savefig(overview_save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved doublet overview plot to {overview_save_path}")

    # --- Figure 2: Standalone UpSet Plot for Heuristics ---
    if LINEAGE_PRED_KEY in adata.obsm and not adata.obsm[LINEAGE_PRED_KEY].empty:
        lineage_df = adata.obsm[LINEAGE_PRED_KEY]
        coexpressing_cells = lineage_df[lineage_df.sum(axis=1) >= 2]

        if not coexpressing_cells.empty:
            lineage_combinations = coexpressing_cells.groupby(
                list(coexpressing_cells.columns)
            ).size()
            lineage_combinations = lineage_combinations[lineage_combinations > 0]

            if not lineage_combinations.empty:
                fig2 = plt.figure(figsize=(12, 7), facecolor="white")
                upset_plot(lineage_combinations, fig=fig2, element_size=32, show_counts=True)
                fig2.suptitle("Heuristic: Lineage Co-expression Summary", fontsize=16)

                if save_dir:
                    upset_save_path = Path(save_dir) / "doublet_summary_upsetplot.png"
                    fig2.savefig(upset_save_path, dpi=300, bbox_inches="tight")
                    log.info(f"Saved doublet UpSet plot to {upset_save_path}")

                if show:
                    plt.show()  # Show after saving this figure
                plt.close(fig2)

            else:
                log.info(
                    "No co-expressing cells found by heuristics to generate an UpSet plot."
                )
        else:
            log.info(
                "No co-expressing cells found by heuristics to generate an UpSet plot."
            )
    else:
        log.info("Heuristic analysis not performed, skipping UpSet plot.")

    if show:
        plt.show()

    plt.close(fig1)


# --- Main Functions ---
def generate_doublet_rates(
    adata: AnnData,
    sample_key: str = "sampleID",
    rate_per_1000_cells: float = 0.008,
    max_rate: float = 0.20,
    min_rate: float = 0.001,
    chemistry: str = "v3",
) -> Dict[str, float]:
    """
    Automatically generate expected doublet rates based on cell count per sample or platform.

    For 10x Genomics data ('v2', 'v3', 'HT'), this function calculates rates based on the
    guideline that multiplet rate increases with cell load.

    For BD Rhapsody data ('BD'), it applies a fixed, empirically-derived doublet rate, as the
    microwell technology's multiplet rate is largely independent of recovered cell count.

    Args:
        adata: AnnData object containing cell count information.
        sample_key: Column name in adata.obs used to distinguish samples.
        rate_per_1000_cells: (For 10x) Expected doublet rate per 1000 cells.
        max_rate: (For 10x) Maximum doublet rate cap.
        min_rate: (For 10x) Minimum doublet rate floor.
        chemistry: The technology platform ('v2', 'v3', 'HT', 'BD', or 'custom').

    Returns:
        Dictionary mapping sample IDs to calculated doublet rates.
    """
    log.info(
        "Automatically generating doublet rates based on sample chemistry and cell counts..."
    )

    # Define platform-specific rates
    # The value for 'BD' is a fixed rate, not a scaling factor.
    chemistry_rates = {
        "v2": 0.007,  # 10x v2 chemistry (rate per 1000 cells)
        "v3": 0.008,  # 10x v3 chemistry (rate per 1000 cells)
        "HT": 0.016,  # 10x High-throughput (rate per 1000 cells)
        "BD": 0.025,  # BD Rhapsody (fixed rate of 2.5%)
        "custom": rate_per_1000_cells,
    }

    if chemistry not in chemistry_rates:
        log.warning(f"Unknown chemistry '{chemistry}'. Using default 'v3' rate model.")
        chemistry = "v3"

    actual_rate = chemistry_rates[chemistry]
    cell_counts = adata.obs[sample_key].value_counts()
    doublet_rates = {}

    if chemistry == "BD":
        log.info(
            f"Using fixed doublet rate of {actual_rate:.4f} for BD Rhapsody platform."
        )
        for sample, n_cells in cell_counts.items():
            doublet_rates[sample] = actual_rate
            log.info(
                f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {actual_rate:.4f}"
            )
    else:
        # Handle 10x-like linear scaling models
        log.info(
            f"Using {chemistry} chemistry with base rate of {actual_rate} per 1000 cells."
        )
        for sample, n_cells in cell_counts.items():
            if n_cells > 10000 and chemistry in ["v2", "v3"]:
                # Non-linear scaling for very high cell counts in standard 10x kits
                rate = (0.8 * (n_cells / 1000) * actual_rate) + (
                    0.2 * actual_rate * (n_cells / 1000) ** 1.5
                )
            else:
                # Standard linear scaling
                rate = (n_cells / 1000) * actual_rate

            # Apply rate constraints
            rate = max(min_rate, min(rate, max_rate))
            doublet_rates[sample] = rate
            log.info(
                f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate:.4f}"
            )

    return doublet_rates


def create_custom_marker_dict(
    lineage_definitions: Dict[str, Dict], save_path: Optional[Union[str, Path]] = None
) -> Dict[str, MarkerConfig]:
    """
    Create custom marker dictionary from user-defined lineage specifications.

    This function allows users to define their own marker sets with custom
    parameters for specialized doublet detection scenarios.

    Args:
        lineage_definitions: Dictionary defining lineages and their parameters
        save_path: Optional path to save the configuration for future use

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects

    Example:
        lineage_defs = {
            "T_cells": {
                "genes": ["CD3D", "CD3E", "CD8A"],
                "expression_threshold": 0.5,
                "min_genes_required": 1
            },
            "Epithelial": {
                "genes": r"^KRT[0-9]+",  # Regex pattern
                "expression_threshold": 1.0,
                "min_genes_required": 2
            }
        }
        marker_configs = create_custom_marker_dict(lineage_defs)
    """
    config_dict = {}

    for lineage, definition in lineage_definitions.items():
        # Validate required 'genes' field
        if "genes" not in definition:
            raise ValueError(f"Missing 'genes' field for lineage '{lineage}'")

        config_dict[lineage] = MarkerConfig(**definition)

    # Save configuration if requested
    if save_path:
        import json

        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to serializable format
        serializable_dict = {}
        for lineage, config in config_dict.items():
            serializable_dict[lineage] = {
                "genes": config.genes,
                "expression_threshold": config.expression_threshold,
                "min_genes_required": config.min_genes_required,
                "use_raw": config.use_raw,
            }

        with open(save_path, "w") as f:
            json.dump(serializable_dict, f, indent=2)
        log.info(f"Marker configuration saved to {save_path}")

    return config_dict


def run_heuristic_analysis(adata: AnnData, config: DoubletConfig) -> AnnData:
    """
    (Step 1: Calculate) Runs the heuristic analysis and stores results in the AnnData object.

    This function's sole purpose is to perform the marker-based co-expression calculation.
    Results are stored in `adata.obsm['lineage_predictions']` and `adata.obs['heuristic_predicted']`.

    Args:
        adata: The AnnData object.
        config: A DoubletConfig object with heuristic parameters configured.

    Returns:
        The AnnData object, modified with heuristic analysis results.
    """
    log.info("--- Running Heuristic Analysis (Calculation Step) ---")
    if not config.use_heuristics:
        log.warning("`use_heuristics` is False in config. Skipping analysis.")
        adata.obs[HEURISTIC_PRED_COL] = False
        adata.obsm[LINEAGE_PRED_KEY] = pd.DataFrame(index=adata.obs_names)
        return adata

    heuristic_pred, lineage_df = _run_heuristic(adata, config)
    adata.obs[HEURISTIC_PRED_COL] = heuristic_pred
    adata.obsm[LINEAGE_PRED_KEY] = lineage_df

    log.info(
        f"Heuristic analysis complete. Found {heuristic_pred.sum()} potential doublets."
    )
    return adata


def analyze_lineage_coexpression(
    adata: AnnData, save_dir: Optional[str] = None, show: bool = True
):
    """
    (Step 2: Visualize) Visualizes pre-computed lineage co-expression patterns.

    This function's sole purpose is to generate an UpSet plot from results
    already stored in `adata.obsm['lineage_predictions']`. It does NOT perform any calculations.
    You must run `run_heuristic_analysis` first.

    Args:
        adata: AnnData object containing results from `run_heuristic_analysis`.
        save_dir: Directory to save the plot.
        show: Whether to display the plot interactively.
    """
    log.info("--- Visualizing Lineage Co-expression (Visualization Step) ---")
    if LINEAGE_PRED_KEY not in adata.obsm or adata.obsm[LINEAGE_PRED_KEY].empty:
        log.error(
            "No lineage predictions found in `adata.obsm`. Please run `run_heuristic_analysis` first."
        )
        return

    lineage_df = adata.obsm[LINEAGE_PRED_KEY]

    # 1. UpSet Plot for intuitive visualization of intersections
    try:
        # Prepare data for UpSet plot by counting cells in each combination
        lineage_combinations = lineage_df.groupby(list(lineage_df.columns)).size()

        fig = plt.figure(figsize=(15, 7), facecolor="white")
        upset_plot(
            lineage_combinations,
            fig=fig,
            # min_subset_size=min_subset_size,
            show_counts=True,
        )
        plt.suptitle("Co-expression of Cell Lineages", fontsize=16)

        if save_dir:
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path / "lineage_coexpression_upsetplot.png", dpi=300)
        if show:
            plt.show()
        plt.close(fig)

    except ImportError:
        log.warning(
            "`upsetplot` library not installed. Skipping UpSet plot. Please run: pip install upsetplot"
        )
    except Exception as e:
        log.error(f"Failed to generate UpSet plot: {e}")

    # 2. Pairwise Overlap Counts
    lineage_pairs = []
    for i, lineage1 in enumerate(lineage_df.columns):
        for lineage2 in lineage_df.columns[i + 1 :]:
            overlap_count = (lineage_df[lineage1] & lineage_df[lineage2]).sum()
            if overlap_count > 0:
                lineage_pairs.append(
                    {
                        "lineage1": lineage1,
                        "lineage2": lineage2,
                        "overlap_count": overlap_count,
                        "overlap_percent": overlap_count / len(lineage_df) * 100,
                    }
                )

    pairs_df = pd.DataFrame(lineage_pairs).sort_values("overlap_count", ascending=False)
    log.info("Top 10 lineage co-expression pairs:")
    log.info("\n" + pairs_df.head(10).to_string())

    if save_dir:
        pairs_df.to_csv(Path(save_dir) / "lineage_coexpression_pairs.csv", index=False)

    return pairs_df


def predict_doublets(
    adata: AnnData, config: DoubletConfig, sample_key: str = "sampleID", **kwargs
) -> AnnData:
    """
    Enhanced doublet prediction with a clear, config-driven workflow.

    This function serves as the main entry point for doublet detection. For advanced control
    and reproducibility, it is highly recommended to create and pass a `DoubletConfig` object.

    Args:
        adata: AnnData object containing single-cell expression data.
        config: A `DoubletConfig` object that controls the entire workflow.
        sample_key: Key for sample identification in adata.obs.

    Returns:
        AnnData object with doublet predictions added to .obs.
    """
    # === 1. CONFIGURATION SETUP ===
    base_config = DoubletConfig()

    if config is not None:
        config_dict = asdict(config)
        for key, value in config_dict.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)

    if kwargs:
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")

    cfg = base_config 
    cfg.validate()
    log.info("--- Running Final Doublet Prediction Workflow ---")

    # Validate input data
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not found in adata.obs")

    samples = adata.obs[sample_key].unique()
    if len(samples) == 0:
        raise ValueError(f"No samples found for key '{sample_key}'")

    log.info(
        f"Starting doublet prediction for {adata.n_obs} cells across {len(samples)} samples"
    )
    log.info(
        f"Configuration: method={cfg.method}, merge_strategy={cfg.merge_strategy}, "  # --- FIX: Use cfg consistently
        f"use_heuristics={cfg.use_heuristics}"
    )

    # Initialize result columns
    algo_score_col = f"{cfg.method}_score"
    algo_pred_col = f"{cfg.method}_predicted"
    adata.obs[algo_score_col] = np.nan
    adata.obs[algo_pred_col] = False
    adata.obs[HEURISTIC_PRED_COL] = False

    # Use a dispatcher for multi-algorithm support ---
    ALGORITHM_DISPATCHER = {
        "scrublet": _run_scrublet
        # "doubletfinder": _run_doubletfinder # Future-ready
    }
    if cfg.method not in ALGORITHM_DISPATCHER:
        raise ValueError(
            f"Method '{cfg.method}' is not supported. Available: {list(ALGORITHM_DISPATCHER.keys())}"
        )

    # === 2. ALGORITHMIC DETECTION (Per-Sample) ===
    log.info(f"Running {cfg.method} doublet detection...")

    for sample in samples:
        log.info(f"Processing sample '{sample}' with {cfg.method}...")
        sample_mask = adata.obs[sample_key] == sample
        data_view = adata[sample_mask]

        if data_view.n_obs < 50:
            log.warning(
                f"Skipping {sample}: fewer than 50 cells (insufficient for reliable doublet detection)."
            )
            continue

        scores, predicted = ALGORITHM_DISPATCHER[cfg.method](data_view, sample, cfg)

        if scores is not None and predicted is not None:
            adata.obs.loc[sample_mask, algo_score_col] = scores
            adata.obs.loc[sample_mask, algo_pred_col] = predicted

    # === 3. HEURISTIC DETECTION (Global) ===
    if cfg.use_heuristics:
        # Check if heuristic results already exist
        if "lineage_predictions" in adata.obsm:
            log.info(
                "Found existing lineage predictions in adata.obsm. Re-evaluating with current config."
            )
            lineage_df = adata.obsm["lineage_predictions"]
            lineages_per_cell = lineage_df.sum(axis=1)
            heuristic_pred = lineages_per_cell >= cfg.min_lineages_for_doublet

            # --- Re-apply the ignore logic here as well ---
            if cfg.ignore_coexpression_pairs:
                log.info("Applying ignore rules to existing lineage predictions.")
                ignored_set = {
                    tuple(sorted(pair)) for pair in cfg.ignore_coexpression_pairs
                }
                doublet_indices = heuristic_pred[heuristic_pred].index
                for cell_idx in doublet_indices:
                    expressed_lineages = tuple(
                        sorted(lineage_df.columns[lineage_df.loc[cell_idx]])
                    )
                    if expressed_lineages in ignored_set:
                        heuristic_pred.loc[cell_idx] = False
        else:
            # If no results exist, run the heuristic from scratch
            heuristic_pred, lineage_df = _run_heuristic(adata, cfg)
            adata.obsm["lineage_predictions"] = lineage_df

        adata.obs[HEURISTIC_PRED_COL] = heuristic_pred

    # === 4. MERGE RESULTS ===
    log.info("Merging algorithmic and heuristic predictions...")
    adata.obs[FINAL_PRED_COL] = _merge_doublet_predictions(
        adata,
        algorithm_col=algo_pred_col,
        heuristic_col=HEURISTIC_PRED_COL,
        strategy=cfg.merge_strategy,
    )

    # STORE PARAMS
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["doublet_params"] = (
        cfg.__dict__
    )

    # === 5. SUMMARY STATISTICS ===
    log.info("\n" + "=" * 50)
    log.info("DOUBLET DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Algorithm results
    algo_count = adata.obs[algo_pred_col].sum()
    log.info(
        f"Algorithm ({cfg.method}): {algo_count} doublets ({algo_count / total_cells:.2%})"
    )

    # Heuristic results
    if cfg.use_heuristics:
        heur_count = adata.obs[HEURISTIC_PRED_COL].sum()
        log.info(f"Heuristic: {heur_count} doublets ({heur_count / total_cells:.2%})")

        # Overlap analysis
        overlap_count = (adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL]).sum()
        log.info(
            f"Overlap: {overlap_count} doublets ({overlap_count / total_cells:.2%})"
        )

    # Final merged results
    final_count = adata.obs[FINAL_PRED_COL].sum()
    log.info(f"Final merged: {final_count} doublets ({final_count / total_cells:.2%})")

    # Per-sample breakdown
    log.info("\nPer-sample statistics:")
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_total = sample_mask.sum()
        sample_doublets = adata.obs[FINAL_PRED_COL][sample_mask].sum()
        sample_rate = sample_doublets / sample_total
        log.info(
            f"  {sample}: {sample_doublets}/{sample_total} doublets ({sample_rate:.2%})"
        )

    log.info("=" * 50)

    # === 6. Reporting & Visualization ===
    if cfg.plot_summary:
        save_path = Path(cfg.save_dir) if cfg.save_dir else None
        _plot_doublet_summary(adata, sample_key, save_path, cfg.show_plots)

    if cfg.export_stats and cfg.save_dir:
        _export_doublet_stats(adata, sample_key, Path(cfg.save_dir))

    log.info("Doublet prediction workflow completed.")

    return adata
