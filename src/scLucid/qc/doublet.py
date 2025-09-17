"""
Enhanced doublet detection for single-cell RNA-seq data.

This module provides comprehensive functions for identifying potential doublet cells
using an algorithmic method (Scrublet) and a flexible heuristic approach based on
mutually exclusive lineage marker co-expression.
"""

import gc
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns
from anndata import AnnData
from upsetplot import plot as upset_plot

from ..markers import get_marker_manager
from .config import DoubletConfig, MarkerConfig

log = logging.getLogger(__name__)

# --- Use constants for column names for easier maintenance ---
LINEAGE_SCORES_KEY = "lineage_module_scores"
HEURISTIC_SCORE_COL = "heuristic_confidence_score"
HEURISTIC_PRED_COL = "heuristic_predicted"
FINAL_PRED_COL = "predicted_doublet"

__all__ = [
    "generate_doublet_rates",
    "create_custom_marker_dict",
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
        # --- ❗ Select the correct adata object upfront ❗ ---
        # Ensure the intersection is performed on the same data that will be used for scoring.
        adata_for_intersection = (
            adata.raw.to_adata() if cfg.default_use_raw and adata.raw else adata
        )
        log.info(f"Performing marker intersection on {'adata.raw' if cfg.default_use_raw and adata.raw else 'adata'}.")
        
        manager = get_marker_manager(
            species=cfg.marker_species,
            tissue=cfg.marker_tissue,
            case_sensitive=False, # Gene symbols are typically case-insensitive in this context
        )
        # Intersect with the correctly chosen data object
        manager.intersect_with(adata_for_intersection)
        markers_dict = manager.get_doublet_lineage_markers()

        if not markers_dict:
            log.warning(
                "No lineage markers found. Ensure `doublet_lineage = true` is set in your TOML file for desired cell types."
            )
            return {}

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


def _merge_doublet_predictions(
    adata: AnnData,
    algorithm_score_col: str,
    heuristic_score_col: str,
    strategy: str = "weighted_average",
    algo_weight: float = 0.6,
    expected_rate: Optional[Union[float, Dict[str, float]]] = 0.1,
) -> pd.Series:
    """
    Merge algorithmic and heuristic doublet scores for a final, more robust prediction.
    This function combines two continuous score series instead of binary predictions.

    Args:
        adata: AnnData object containing the scores.
        algorithm_score_col: Column name in adata.obs for the algorithm's score (e.g., 'scrublet_score').
        heuristic_score_col: Column name in adata.obs for the heuristic confidence score.
        strategy: The merge strategy ('weighted_average', 'max_score', 'heuristic_boost').
        algo_weight: The weight for the algorithm's score in 'weighted_average' strategy.

    Returns:
        A boolean pandas Series with the final merged doublet predictions.
    """
    algo_scores = adata.obs[algorithm_score_col].fillna(0)
    heur_scores = adata.obs[heuristic_score_col].fillna(0)

    final_score = pd.Series(0.0, index=adata.obs_names)

    if strategy == "weighted_average":
        # A simple weighted average. algo_weight determines the trust in the algorithm.
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)
    elif strategy == "max_score":
        # Takes the highest score from either method, useful if either method is considered reliable on its own.
        final_score = pd.DataFrame({"algo": algo_scores, "heur": heur_scores}).max(
            axis=1
        )
    elif strategy == "heuristic_boost":
        # Uses the algorithm score as a base and the heuristic score as a "booster".
        # This is useful for finding doublets missed by the algorithm but strongly suggested by heuristics.
        final_score = algo_scores + (heur_scores * 0.5)  # Boost factor can be tuned
    else:
        log.warning(
            f"Unknown enhanced merge strategy '{strategy}', falling back to 'weighted_average'."
        )
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)

    # Normalize the final combined score to a [0, 1] range for consistent thresholding.
    if final_score.max() > 0:
        final_score /= final_score.max()

    if expected_rate is None:
        log.warning("expected_doublet_rate is None, using a default of 0.1 for thresholding.")
        expected_rate = 0.1
        
    if isinstance(expected_rate, dict):  # Handle per-sample rates by taking the mean
        expected_rate = np.mean(list(expected_rate.values()))

    threshold = final_score.quantile(1 - expected_rate)
    log.info(
        f"Using a final score threshold of {threshold:.3f} based on expected doublet rate for merged predictions."
    )

    return final_score > threshold


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
) -> Tuple[pd.Series, pd.DataFrame, pd.Series]:
    """
    Runs the full heuristic doublet detection workflow, calculating scores instead of binary predictions.
    This enhanced version uses Scanpy's module scoring to quantify lineage expression and computes
    a confidence score based on co-expression strength and cell complexity (n_genes_by_counts).

    Args:
        adata: AnnData object.
        cfg: The DoubletConfig object.

    Returns:
        A tuple containing:
        - A boolean pd.Series of heuristic doublet predictions based on a score threshold.
        - A pd.DataFrame (`lineage_scores_df`) with module scores for each lineage.
        - A pd.Series (`heuristic_confidence_score`) with the final calculated confidence score, normalized to [0, 1].
    """
    # --- 1. Load Marker Configurations ---
    marker_configs = cfg.marker_configs
    if marker_configs is None:
        marker_configs = _create_doublet_marker_config_from_manager(adata, cfg)

    if not marker_configs:
        log.warning(
            "Heuristics enabled, but no marker configurations were found. Skipping."
        )
        empty_series = pd.Series(0.0, index=adata.obs_names)
        empty_df = pd.DataFrame(index=adata.obs_names)
        return empty_series.astype(bool), empty_df, empty_series

    # --- ❗ In-function Normalization for Accurate Scoring ❗ ---
    # To avoid altering the original adata object at the QC stage, we create a temporary,
    # normalized copy for the sole purpose of running sc.tl.score_genes.
    log.info("Creating temporary normalized data for accurate gene scoring...")
    # Use .raw if available and specified, otherwise use the main adata.X
    source_adata = adata.raw.to_adata() if cfg.default_use_raw and adata.raw else adata.copy()
    
    # Perform standard normalization and log-transformation on the temporary object.
    sc.pp.normalize_total(source_adata, target_sum=1e4)
    sc.pp.log1p(source_adata)
    
    # --- 2. Calculate Module Score for Each Lineage ---
    lineage_scores = {}
    for name, marker_config in marker_configs.items():
        # This simplified version assumes genes are in a list. A real implementation
        # would need to handle the regex case from the original MarkerConfig.
        if marker_config.is_regex:
            log.warning(
                f"Regex markers for lineage '{name}' are not supported in this scoring function. Skipping."
            )
            continue

        valid_genes = [g for g in marker_config.genes if g in source_adata.var_names]

        # sc.tl.score_genes requires at least 2 genes in the list.
        if len(valid_genes) < 2:
            log.warning(
                f"Skipping lineage '{name}': not enough valid genes ({len(valid_genes)}) found for scoring."
            )
            continue

        score_name = f"{name}_score"
        # The control size should not exceed the number of available genes.
        ctrl_size = min(50, source_adata.n_vars - len(valid_genes) - 1)

        sc.tl.score_genes(
            source_adata,
            gene_list=valid_genes,
            score_name=score_name,
            ctrl_size=ctrl_size,
        )
        lineage_scores[name] = source_adata.obs[score_name]

    lineage_scores_df = pd.DataFrame(lineage_scores, index=adata.obs_names).fillna(0)

    if lineage_scores_df.empty:
        log.warning(
            "Heuristic analysis could not proceed because no lineages had sufficient marker genes "
            "present in the data. Please check gene symbols (e.g., case sensitivity) in your marker file "
            "against adata.var_names."
        )
        empty_series = pd.Series(0.0, index=adata.obs_names)
        empty_df = pd.DataFrame(index=adata.obs_names)
        return empty_series.astype(bool), empty_df, empty_series
    
    # --- 3. Compute Heuristic Confidence Score ---
    # Get the top two lineage scores for each cell.
    top_two_scores = lineage_scores_df.apply(
        lambda s: s.nlargest(2).values, axis=1, result_type="expand"
    )
    # Handle cases where a cell has scores for less than 2 lineages
    if top_two_scores.shape[1] == 1:
        top_two_scores[1] = 0.0  # Add a second column of zeros
    top_two_scores.columns = ["score1", "score2"]

    # A cell is a co-expression candidate if its top two scores are both significant.
    # The threshold 0.1 is an empirical value, meaning scores are clearly positive.
    significant_coexpression = (top_two_scores["score1"] > 0.1) & (
        top_two_scores["score2"] > 0.1
    )

    # Use log-transformed gene counts as a weight for cell complexity.
    # Doublets are expected to have more detected genes.
    gene_count_log = np.log1p(adata.obs["n_genes_by_counts"])

    # Calculate the final score: sum of top two scores, weighted by gene counts,
    # and only applied to cells showing significant co-expression.
    heuristic_confidence_score = (
        (top_two_scores["score1"] + top_two_scores["score2"])
        * gene_count_log
        * significant_coexpression
    )

    # Normalize the score to a [0, 1] range to make it comparable to Scrublet's score.
    max_score = heuristic_confidence_score.max()
    if max_score > 0:
        heuristic_confidence_score /= max_score

    # --- 4. Generate a Binary Prediction based on the Score ---
    # This provides a simple True/False output for summary statistics.
    # A high quantile (e.g., 0.95) means only the top 5% of scored cells are flagged.
    # This is an ad-hoc threshold; the continuous score is more informative.
    score_threshold = heuristic_confidence_score[
        heuristic_confidence_score > 0
    ].quantile(0.90)
    if pd.isna(score_threshold):
        score_threshold = 0

    potential_doublets = heuristic_confidence_score > score_threshold

    # The ignore_coexpression_pairs logic can still be applied here if needed,
    # for example, to set scores of certain co-expressing cells to 0.
    
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    adata.uns["sclucid"]["qc"]["doublet_params"]["heuristic_temp_norm"] = {
        "used": True,
        "use_raw": cfg.default_use_raw,
    }
    return potential_doublets, lineage_scores_df, heuristic_confidence_score.fillna(0)


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
    upset_score_threshold: float = 0.1,
    save_dir: Optional[Union[str, Path]] = None,
    show: bool = True,
    plot_bar: bool = True,
    plot_scatter: bool = True,
    plot_upset: bool = True,
) -> None:
    """
    Generates a comprehensive, UMAP-independent summary plot for doublet detection results.

    This function is enhanced to work with the quantitative heuristic scoring workflow.
    It creates a multi-panel figure showing:
    1. A stacked bar plot of doublet counts per sample, broken down by prediction source.
    2. A dual-panel scatter plot comparing algorithmic and heuristic scores against gene counts.
    3. An UpSet plot visualizing co-expression patterns from heuristic analysis, based on binarizing
       the continuous lineage module scores.

    Args:
        adata: AnnData object after running `predict_doublets`.
        sample_key: The key in adata.obs for sample identification.
        upset_score_threshold: The minimum module score required to consider a lineage "positive" for the UpSet plot.
        save_dir: Directory to save the plot.
        show: Whether to display the plot.
    """
    # --- Check if there's anything to plot ---
    if not any([plot_bar, plot_scatter, plot_upset]):
        log.info("All summary plot panels are disabled. Skipping plot generation.")
        return

    log.info("Generating doublet detection summary plot...")

    # --- 1. Data Validation and Preparation ---
    # Find the algorithm's prediction and score columns dynamically.
    algo_pred_col_list = [
        c for c in adata.obs.columns if c.endswith("_predicted") and "heuristic" not in c
    ]
    if not algo_pred_col_list or FINAL_PRED_COL not in adata.obs:
        log.warning(
            f"Required prediction columns not found ('{FINAL_PRED_COL}' or algorithm-specific). "
            "Run `predict_doublets` first. Skipping summary plot."
        )
        return
    algo_pred_col = algo_pred_col_list[0]
    algo_score_col = algo_pred_col.replace("_predicted", "_score")

    # Ensure heuristic column exists, even if analysis was skipped.
    if HEURISTIC_PRED_COL not in adata.obs:
        adata.obs[HEURISTIC_PRED_COL] = False

    # Create a 'doublet_source' column for breakdown visualization.
    source_categories = ["Singleton", "Heuristic Only", "Algorithm Only", "Both Methods"]
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

    # --- 2. Dynamically Create the Figure Layout ---
    # The layout now depends on which plots are enabled.
    n_rows = 0
    if plot_bar or plot_scatter:
        n_rows += 1
    if plot_upset:
        n_rows += 1
        
    if n_rows == 0: # Should be caught by the check above, but as a safeguard.
        return
        
    fig = plt.figure(figsize=(18, 6 * n_rows), facecolor="white")
    
    if n_rows == 2: # All plots enabled
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])
        ax_bar = fig.add_subplot(gs[0, 0]) if plot_bar else None
        gs_scatter_spec = gs[0, 1] if plot_scatter else None
        ax_upset = fig.add_subplot(gs[1, :]) if plot_upset else None
    elif plot_upset: # Only Upset plot
        gs = fig.add_gridspec(1, 1)
        ax_upset = fig.add_subplot(gs[0, 0])
        ax_bar, gs_scatter_spec = None, None
    else: # Only top row (bar and/or scatter)
        gs = fig.add_gridspec(1, 2)
        ax_bar = fig.add_subplot(gs[0, 0]) if plot_bar else None
        gs_scatter_spec = gs[0, 1] if plot_scatter else None
        ax_upset = None

    fig.suptitle("Doublet Detection Summary", fontsize=20, fontweight="bold")
    
    # --- 3. Panel 1: Doublet Breakdown Bar Plot ---
    if plot_bar and ax_bar:
        log.info("Plotting doublet breakdown per sample...")
        summary_df = (
            adata.obs.groupby(sample_key)["doublet_source"]
            .value_counts(normalize=True, sort=False)
            .unstack(fill_value=0)
            * 100
        )
        # Ensure all categories are present for consistent coloring
        for cat in source_categories:
            if cat not in summary_df.columns:
                summary_df[cat] = 0.0
                
        summary_df[source_categories].plot(
            kind="bar",
            stacked=True,
            ax=ax_bar,
            color={
                "Singleton": "lightgray",
                "Heuristic Only": "coral",
                "Algorithm Only": "skyblue",
                "Both Methods": "darkred",
            },
            width=0.8,
        )
        ax_bar.set_title("Doublet Breakdown per Sample", fontsize=14)
        ax_bar.set_ylabel("Percentage of Cells (%)")
        ax_bar.set_xlabel("Sample")
        ax_bar.tick_params(axis="x", labelrotation=45)
        ax_bar.legend(title="Source", bbox_to_anchor=(1.05, 1), loc="upper left")
        ax_bar.grid(axis="y", linestyle="--", alpha=0.7)

    # --- 4. Panel 2: Enhanced Dual Scatter Plots ---
    if plot_scatter and gs_scatter_spec:
        log.info("Plotting diagnostic scatter plots...")
        # This requires a nested GridSpec for the two scatter plots
        gs_scatter = gs[0, 1].subgridspec(1, 2, wspace=0.4)
        ax_scatter_algo = fig.add_subplot(gs_scatter[0])
        ax_scatter_heur = fig.add_subplot(gs_scatter[1])

        # Plot 2a: Algorithm Score vs. Gene Count
        if algo_score_col in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
            sns.scatterplot(
                data=adata.obs, x="n_genes_by_counts", y=algo_score_col,
                hue=FINAL_PRED_COL, s=5, alpha=0.6,
                palette={True: "red", False: "gray"}, ax=ax_scatter_algo, rasterized=True,
            )
            ax_scatter_algo.set_title("Algorithm Score vs. Gene Count", fontsize=12)
            ax_scatter_algo.set_xlabel("Number of Genes")
            ax_scatter_algo.set_ylabel("Algorithm Score")
            ax_scatter_algo.get_legend().remove()
        else:
            ax_scatter_algo.text(0.5, 0.5, "Algorithm score data unavailable.", ha="center", va="center")

        # Plot 2b: Heuristic Score vs. Gene Count
        if HEURISTIC_SCORE_COL in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
            sns.scatterplot(
                data=adata.obs, x="n_genes_by_counts", y=HEURISTIC_SCORE_COL,
                hue=FINAL_PRED_COL, s=5, alpha=0.6,
                palette={True: "red", False: "gray"}, ax=ax_scatter_heur, rasterized=True,
            )
            ax_scatter_heur.set_title("Heuristic Score vs. Gene Count", fontsize=12)
            ax_scatter_heur.set_xlabel("Number of Genes")
            ax_scatter_heur.set_ylabel("Heuristic Score")
            ax_scatter_heur.legend(title="Final Call", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax_scatter_heur.text(0.5, 0.5, "Heuristic score data unavailable.", ha="center", va="center")

    # --- 5. Panel 3: UpSet Plot from Quantitative Scores (THE FIX) ---
    if plot_upset and ax_upset:
        log.info("Plotting heuristic co-expression UpSet plot...")
        if LINEAGE_SCORES_KEY in adata.obsm and not adata.obsm[LINEAGE_SCORES_KEY].empty:
            lineage_scores_df = adata.obsm[LINEAGE_SCORES_KEY]

            # Binarize the continuous scores to create sets for the UpSet plot.
            # This is the key step to adapt the plot to the new quantitative workflow.
            lineage_bool_df = lineage_scores_df > upset_score_threshold
            
            # Filter for cells that are positive for at least two lineages.
            coexpressing_cells = lineage_bool_df[lineage_bool_df.sum(axis=1) >= 2]

            if not coexpressing_cells.empty:
                # Count the number of cells in each co-expression combination.
                lineage_combinations = coexpressing_cells.groupby(
                    list(coexpressing_cells.columns)
                ).size()
                lineage_combinations = lineage_combinations[lineage_combinations > 0]

                if not lineage_combinations.empty:
                    try:
                        upset_plot(lineage_combinations, fig=fig, element_size=32, show_counts=True)
                        ax_upset.set_title(
                            f"Heuristic: Lineage Co-expression (Score > {upset_score_threshold})",
                            fontsize=14, y=0.9 # Adjust title position if needed
                        )
                    except Exception as e:
                        log.error(f"Failed to generate UpSet plot: {e}")
                        ax_upset.text(0.5, 0.5, "Error generating UpSet plot.", ha="center", va="center")
                else:
                    log.info("No cells found with co-expression above the threshold.")
                    ax_upset.text(0.5, 0.5, "No significant co-expression found.", ha="center", va="center")
            else:
                log.info("No cells found with co-expression above the threshold.")
                ax_upset.text(0.5, 0.5, "No significant co-expression found.", ha="center", va="center")
            
            # Remove original subplot axis from the grid to let upsetplot take over
            ax_upset.axis('off')

        else:
            log.info("Heuristic lineage scores not found, skipping UpSet plot.")
            ax_upset.text(0.5, 0.5, "Heuristic scores data unavailable.", ha="center", va="center")
            ax_upset.axis('off')

    # --- 6. Finalization and Saving ---
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        overview_save_path = save_path / "doublet_summary_plot.png"
        fig.savefig(overview_save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved doublet summary plot to {overview_save_path}")

    if show:
        plt.show()

    plt.close(fig)


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


def predict_doublets(
    adata: AnnData, config: DoubletConfig, sample_key: str = "sampleID", **kwargs
) -> AnnData:
    """
    Enhanced doublet prediction with a clear, config-driven workflow.
    This version integrates a quantitative heuristic score with the algorithmic score for improved accuracy.

    Args:
        adata: AnnData object containing single-cell expression data.
        config: A `DoubletConfig` object that controls the entire workflow.
        sample_key: Key for sample identification in adata.obs.

    Returns:
        AnnData object with doublet predictions added to .obs and .obsm.
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
        f"Configuration: method={cfg.method}, merge_strategy={cfg.merge_strategy}, "
    )

    # Initialize result columns
    algo_score_col = f"{cfg.method}_score"
    algo_pred_col = f"{cfg.method}_predicted"
    adata.obs[algo_score_col] = np.nan
    adata.obs[algo_pred_col] = False

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
    if cfg.run_algorithm:
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
    else:
        log.info("Skipping algorithmic detection as per configuration (run_algorithm=False).")

    # === 3. HEURISTIC DETECTION (Global) ===
    adata.obs[HEURISTIC_PRED_COL] = False
    adata.obs[HEURISTIC_SCORE_COL] = 0.0
    if cfg.use_heuristics:
        log.info("Running quantitative heuristic analysis...")
        # Call the new heuristic function and receive its multiple outputs
        heuristic_pred, lineage_scores_df, heuristic_scores = _run_heuristic(adata, cfg)

        # Store all the new results in the AnnData object
        adata.obsm["lineage_module_scores"] = (
            lineage_scores_df  # Store detailed scores in .obsm
        )
        adata.obs[HEURISTIC_PRED_COL] = (
            heuristic_pred  # Store the binary call for simple stats
        )
        adata.obs[HEURISTIC_SCORE_COL] = (
            heuristic_scores  # Store the informative continuous score
        )
        log.info(
            f"Heuristic analysis complete. Found {heuristic_pred.sum()} potential doublets based on score threshold."
        )

    # === 4. MERGE RESULTS ===
    log.info("Merging algorithmic and heuristic scores for final prediction...")
    merged_pred = _merge_doublet_predictions(
        adata,
        algorithm_score_col=algo_score_col,
        heuristic_score_col=HEURISTIC_SCORE_COL,
        strategy=cfg.merge_strategy,
        expected_rate=cfg.expected_doublet_rate,
    )
    adata.obs[FINAL_PRED_COL] = merged_pred

    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    adata.uns["sclucid"]["qc"]["doublet_params"].update({
        "merge_strategy": cfg.merge_strategy,
        "algorithm_weight": cfg.algorithm_weight,
        "expected_doublet_rate": cfg.expected_doublet_rate,
        "method": cfg.method,
    })

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
        _plot_doublet_summary(
            adata=adata,
            sample_key=sample_key,
            save_dir=save_path,
            show=cfg.show_plots,
            plot_bar=cfg.plot_bar,
            plot_scatter=cfg.plot_scatter,
            plot_upset=cfg.plot_upset,
        )

    if cfg.export_stats and cfg.save_dir:
        _export_doublet_stats(adata, sample_key, Path(cfg.save_dir))

    log.info("Doublet prediction workflow completed.")

    return adata
