"""Heuristic doublet detection based on lineage marker co-expression.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from ...utils import get_marker_manager
from ..config import DoubletConfig, MarkerConfig
from .core import (
    FINAL_PRED_COL,
    HEURISTIC_PRED_COL,
    HEURISTIC_SCORE_COL,
    LINEAGE_SCORES_KEY,
    _create_doublet_marker_config_from_manager,
)

log = logging.getLogger(__name__)

def _run_heuristic(adata: AnnData, cfg: DoubletConfig) -> Tuple[pd.Series, pd.DataFrame, pd.Series]:
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
        log.warning("Heuristics enabled, but no marker configurations were found. Skipping.")
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
    significant_coexpression = (top_two_scores["score1"] > 0.1) & (top_two_scores["score2"] > 0.1)

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
    score_threshold = heuristic_confidence_score[heuristic_confidence_score > 0].quantile(0.90)
    if pd.isna(score_threshold):
        score_threshold = 0

    potential_doublets = heuristic_confidence_score > score_threshold

    # The ignore_coexpression_pairs logic can still be applied here if needed,
    # for example, to set scores of certain co-expressing cells to 0.
    if cfg.ignore_coexpression_pairs:
        # top_two_scores 包含了 lineage names，假设你修改了逻辑让它返回 Name 而不是 Score
        # 或者在这里遍历 whitelist
        for lin1, lin2 in cfg.ignore_coexpression_pairs:
            # 检查谱系是否存在于lineage_scores_df中
            if lin1 in lineage_scores_df.columns and lin2 in lineage_scores_df.columns:
                # 找到同时高表达 lin1 和 lin2 的细胞，将其 heuristic_score 强制置 0
                mask = (lineage_scores_df[lin1] > 0.1) & (lineage_scores_df[lin2] > 0.1)
                heuristic_confidence_score[mask] = 0.0
                log.info(
                    f"Ignoring co-expression of {lin1} + {lin2} in {mask.sum()} cells (Allowlist)."
                )
            else:
                missing = []
                if lin1 not in lineage_scores_df.columns:
                    missing.append(lin1)
                if lin2 not in lineage_scores_df.columns:
                    missing.append(lin2)
                log.debug(
                    f"Skipping co-expression pair ({lin1}, {lin2}) - lineages not found: {missing}"
                )

    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    adata.uns["sclucid"]["qc"]["doublet_params"]["heuristic_temp_norm"] = {
        "used": True,
        "use_raw": cfg.default_use_raw,
    }
    return potential_doublets, lineage_scores_df, heuristic_confidence_score.fillna(0)



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
    Generates comprehensive, separate summary plots for doublet detection results.

    This refactored version saves each plot type (bar, scatter, upset) as a
    separate file for improved clarity and usability.

    Args:
        adata: AnnData object after running `predict_doublets`.
        sample_key: The key in adata.obs for sample identification.
        upset_score_threshold: The minimum module score to consider a lineage "positive" for the UpSet plot.
        save_dir: Directory to save the plots.
        show: Whether to display the plots.
        plot_bar, plot_scatter, plot_upset: Booleans to control which plots are generated.
    """
    if not any([plot_bar, plot_scatter, plot_upset]):
        log.info("All summary plot panels are disabled. Skipping plot generation.")
        return

    log.info("Generating doublet detection summary plots...")
    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        log.info(f"Plots will be saved to: {save_path.resolve()}")

    # --- 1. Data Validation and Preparation ---
    algo_pred_col_list = [
        c for c in adata.obs.columns if c.endswith("_predicted") and "heuristic" not in c
    ]
    if not algo_pred_col_list or FINAL_PRED_COL not in adata.obs:
        log.warning(
            f"Required prediction columns not found ('{FINAL_PRED_COL}' or algorithm-specific). "
            "Run `predict_doublets` first. Skipping summary plots."
        )
        return
    algo_pred_col = algo_pred_col_list[0]
    algo_score_col = algo_pred_col.replace("_predicted", "_score")

    if HEURISTIC_PRED_COL not in adata.obs:
        adata.obs[HEURISTIC_PRED_COL] = False

    # --- 2. Panel 1: Doublet Breakdown Bar Plot ---
    if plot_bar:
        log.info("Generating doublet breakdown bar plot...")
        fig_bar, ax_bar = plt.subplots(
            figsize=(max(6, len(adata.obs[sample_key].unique()) * 0.5), 6),
            facecolor="white",
        )

        source_categories = [
            "Singleton",
            "Heuristic Only",
            "Algorithm Only",
            "Both Methods",
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

        summary_df = (
            adata.obs.groupby(sample_key)["doublet_source"]
            .value_counts(normalize=True, sort=False)
            .unstack(fill_value=0)
            * 100
        )
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
        ax_bar.set_title("Doublet Breakdown per Sample", fontsize=16, fontweight="bold")
        ax_bar.set_ylabel("Percentage of Cells (%)")
        ax_bar.set_xlabel("Sample")
        ax_bar.tick_params(axis="x", labelrotation=45)
        plt.setp(ax_bar.get_xticklabels(), ha="right")
        ax_bar.legend(title="Source", bbox_to_anchor=(1.05, 1), loc="upper left")
        ax_bar.grid(axis="y", linestyle="--", alpha=0.7)
        fig_bar.tight_layout()

        if save_dir:
            bar_save_path = save_path / "doublet_summary_bar_plot.png"
            fig_bar.savefig(bar_save_path, dpi=300, bbox_inches="tight")
            log.info(f"Saved bar plot to {bar_save_path}")
        if show:
            plt.show()
        plt.close(fig_bar)

    # --- 3. Panel 2: Enhanced Dual Scatter Plots ---
    if plot_scatter:
        log.info("Generating diagnostic scatter plots...")
        fig_scatter, (ax_scatter_algo, ax_scatter_heur) = plt.subplots(
            1, 2, figsize=(12, 5), facecolor="white"
        )
        fig_scatter.suptitle("Diagnostic Scatter Plots", fontsize=16, fontweight="bold")

        if algo_score_col in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
            sns.scatterplot(
                data=adata.obs,
                x="n_genes_by_counts",
                y=algo_score_col,
                hue=FINAL_PRED_COL,
                s=5,
                alpha=0.6,
                palette={True: "red", False: "gray"},
                ax=ax_scatter_algo,
                rasterized=True,
            )
            ax_scatter_algo.set_title("Algorithm Score vs. Gene Count", fontsize=12)
            ax_scatter_algo.set_xlabel("Number of Genes")
            ax_scatter_algo.set_ylabel("Algorithm Score")
            ax_scatter_algo.legend().remove()
        else:
            ax_scatter_algo.text(
                0.5, 0.5, "Algorithm score data unavailable.", ha="center", va="center"
            )

        if HEURISTIC_SCORE_COL in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns:
            sns.scatterplot(
                data=adata.obs,
                x="n_genes_by_counts",
                y=HEURISTIC_SCORE_COL,
                hue=FINAL_PRED_COL,
                s=5,
                alpha=0.6,
                palette={True: "red", False: "gray"},
                ax=ax_scatter_heur,
                rasterized=True,
            )
            ax_scatter_heur.set_title("Heuristic Score vs. Gene Count", fontsize=12)
            ax_scatter_heur.set_xlabel("Number of Genes")
            ax_scatter_heur.set_ylabel("Heuristic Score")
            ax_scatter_heur.legend(title="Final Call", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax_scatter_heur.text(
                0.5, 0.5, "Heuristic score data unavailable.", ha="center", va="center"
            )

        fig_scatter.tight_layout(rect=[0, 0, 1, 0.96])

        if save_dir:
            scatter_save_path = save_path / "doublet_summary_scatter_plot.png"
            fig_scatter.savefig(scatter_save_path, dpi=300, bbox_inches="tight")
            log.info(f"Saved scatter plot to {scatter_save_path}")
        if show:
            plt.show()
        plt.close(fig_scatter)

    # --- 4. Panel 3: UpSet Plot ---
    if plot_upset:
        log.info("Generating heuristic co-expression UpSet plot...")
        if LINEAGE_SCORES_KEY in adata.obsm and not adata.obsm[LINEAGE_SCORES_KEY].empty:
            lineage_scores_df = adata.obsm[LINEAGE_SCORES_KEY]
            lineage_bool_df = lineage_scores_df > upset_score_threshold
            coexpressing_cells = lineage_bool_df[lineage_bool_df.sum(axis=1) >= 2]

            if not coexpressing_cells.empty:
                lineage_combinations = coexpressing_cells.groupby(
                    list(coexpressing_cells.columns)
                ).size()
                lineage_combinations = lineage_combinations[lineage_combinations > 0]

                if not lineage_combinations.empty:
                    try:
                        # Upsetplot creates its own figure
                        fig_upset = plt.figure(figsize=(12, 7), facecolor="white")
                        upset_plot(
                            lineage_combinations,
                            fig=fig_upset,
                            element_size=32,
                            show_counts=True,
                        )
                        fig_upset.suptitle(
                            f"Heuristic: Lineage Co-expression (Score > {upset_score_threshold})",
                            fontsize=16,
                            fontweight="bold",
                        )

                        if save_dir:
                            upset_save_path = save_path / "doublet_summary_upset_plot.png"
                            fig_upset.savefig(upset_save_path, dpi=300, bbox_inches="tight")
                            log.info(f"Saved UpSet plot to {upset_save_path}")
                        if show:
                            plt.show()
                        plt.close(fig_upset)

                    except Exception as e:
                        log.error(f"Failed to generate UpSet plot: {e}")
                else:
                    log.info(
                        "No cells found with co-expression above the threshold. Skipping UpSet plot."
                    )
            else:
                log.info(
                    "No cells found with co-expression above the threshold. Skipping UpSet plot."
                )
        else:
            log.info("Heuristic lineage scores not found, skipping UpSet plot.")


# --- Main Functions ---
