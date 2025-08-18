"""
Cell filtering utilities for single-cell RNA-seq data.

This module provides functions for identifying and filtering low-quality
cells based on various quality metrics.
"""

import logging
import os
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from anndata import AnnData

from ..utils.utils import identify_outliers

log = logging.getLogger(__name__)

__all__ = [
    "mark_low_quality_cell",
    "filter_cells",
]


def mark_low_quality_cell(
    adata: AnnData,
    sample_key: str = "sampleID",
    min_genes: int = 200,
    nmads: float = 5.0,
    pc_mt: float = 20.0,
    pc_hb: float = 20.0,
    pc_top20_genes: Optional[float] = None,
    use_fixed_top20_threshold: bool = False,
    plot_outliers: bool = False,
    save_dir: Optional[str] = None,
    show: bool = True,
    qc_metrics: Optional[List[Tuple[str, str]]] = None,
    cols_to_plot: Optional[List[str]] = None,
) -> AnnData:
    """
    Identify and mark low-quality cells based on various QC criteria.
    This function focuses on technical artifacts and low-quality cells,
    not doublet detection (which is handled separately in the doublet module).

    Args:
        adata: AnnData object with QC metrics calculated.
        sample_key: Key in adata.obs for sample identification.
        min_genes: Minimum number of genes to keep a cell.
        nmads: Number of MADs for outlier detection.
        pc_mt: Maximum percentage of mitochondrial counts.
        pc_hb: Maximum percentage of hemoglobin counts.
        pc_top20_genes: Fixed threshold for pct_counts_in_top_20_genes if use_fixed_top20_threshold=True.
        use_fixed_top20_threshold: Whether to use fixed threshold for pct_counts_in_top_20_genes instead of MAD.
        plot_outliers: Whether to generate scatter plots highlighting outliers.
        save_dir: Directory to save plots.
        show: Whether to display plots.
        qc_metrics: List of (metric, direction) tuples for QC-based outlier detection.
        cols_to_plot: List of .obs columns to plot as scatter plots if plot_outliers is True.

    Returns:
        AnnData object with boolean columns in .obs marking low-quality cells.
    """
    required = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    if any(col not in adata.obs for col in required):
        raise ValueError(
            "Missing required QC columns. Run calculate_qc_metric() first."
        )

    if qc_metrics is None:
        qc_metrics = [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
        ]
        if use_fixed_top20_threshold and pc_top20_genes is not None:
            log.info(
                f"Using fixed threshold {pc_top20_genes}% for pct_counts_in_top_20_genes in QC"
            )
        elif not use_fixed_top20_threshold:
            qc_metrics.append(("pct_counts_in_top_20_genes", "upper"))
            log.info("Using MAD-based threshold for pct_counts_in_top_20_genes in QC")

    formatted_metrics = [(metric, direction, None) for metric, direction in qc_metrics]

    # Precompute sample indices to improve performance
    sample_indices = {
        sample: adata.obs[sample_key] == sample
        for sample in adata.obs[sample_key].unique()
    }

    # --- Mark cells based on individual criteria ---
    adata.obs["outlier_low_genes"] = adata.obs["n_genes_by_counts"] < min_genes
    log.info(
        f"Cells with low gene counts (< {min_genes}): {adata.obs['outlier_low_genes'].sum()} "
        f"({adata.obs['outlier_low_genes'].sum() / adata.n_obs:.2%})"
    )

    adata.obs["outlier_mt"] = adata.obs["pct_counts_mt"] > pc_mt
    log.info(
        f"Cells with high mitochondrial percentage (> {pc_mt}%): {adata.obs['outlier_mt'].sum()} "
        f"({adata.obs['outlier_mt'].sum() / adata.n_obs:.2%})"
    )

    if "pct_counts_hb" in adata.obs:
        adata.obs["outlier_hb"] = adata.obs["pct_counts_hb"] > pc_hb
        log.info(
            f"Cells with high hemoglobin percentage (> {pc_hb}%): {adata.obs['outlier_hb'].sum()} "
            f"({adata.obs['outlier_hb'].sum() / adata.n_obs:.2%})"
        )

    # --- Mark cells based on MAD outliers (per sample) ---
    adata.obs["outlier_qc_metrics"] = identify_outliers(
        adata, metrics=formatted_metrics, sample_key=sample_key, nmads=nmads
    )

    if use_fixed_top20_threshold and pc_top20_genes is not None:
        adata.obs["outlier_top20_genes"] = (
            adata.obs["pct_counts_in_top_20_genes"] > pc_top20_genes
        )
        log.info(
            f"Cells with high top 20 genes percentage (> {pc_top20_genes}%): {adata.obs['outlier_top20_genes'].sum()} "
            f"({adata.obs['outlier_top20_genes'].sum() / adata.n_obs:.2%})"
        )
        adata.obs["outlier_qc_metrics"] = (
            adata.obs["outlier_qc_metrics"] | adata.obs["outlier_top20_genes"]
        )

    qc_count = adata.obs["outlier_qc_metrics"].sum()
    log.info(
        f"Cells marked as QC matric outliers: {qc_count} ({qc_count / adata.n_obs:.2%})"
    )

    # --- Final Statistics ---
    log.info("\n--- Overall Low-Quality Cell Statistics ---")
    total_cells = adata.n_obs

    # Count cells with multiple issues
    outlier_cols = [col for col in adata.obs.columns if col.startswith("outlier_")]
    adata.obs["outlier_count"] = adata.obs[outlier_cols].sum(axis=1)

    # Report counts per outlier type
    for col in outlier_cols:
        count = adata.obs[col].sum()
        log.info(f"Cells marked as '{col}': {count} ({count / total_cells:.2%})")

    # Report cells with multiple issues
    for n_outliers in range(1, len(outlier_cols) + 1):
        count = (adata.obs["outlier_count"] == n_outliers).sum()
        if count > 0:
            log.info(
                f"Cells with exactly {n_outliers} types of issues: {count} ({count / total_cells:.2%})"
            )

    # Include doublet statistics if available
    if "predicted_doublet_scrublet" in adata.obs:
        count = adata.obs["predicted_doublet_scrublet"].sum()
        log.info(
            f"Cells marked as 'predicted_doublet_scrublet': {count} ({count / total_cells:.2%})"
        )

    if "doublet_expression_pattern" in adata.obs:
        count = adata.obs["doublet_expression_pattern"].sum()
        log.info(
            f"Cells marked as 'doublet_expression_pattern': {count} ({count / total_cells:.2%})"
        )

    # --- Plotting ---
    if plot_outliers:
        if cols_to_plot is None:
            # Define default columns to plot, checking for existence
            default_cols = [
                "outlier_low_genes",
                "outlier_mt",
                "outlier_hb",
                "outlier_qc_metrics",
                "predicted_doublet_scrublet",
                "doublet_expression_pattern",
            ]
            cols_to_plot = [col for col in default_cols if col in adata.obs.columns]

        for sample, sample_mask in sample_indices.items():
            log.info(f"Plotting outlier QC for sample: {sample}")
            data_view = adata[sample_mask]

            if data_view.n_obs == 0:
                log.warning(f"No cells found for sample {sample}, skipping plot")
                continue

            fig, axs = plt.subplots(2, 3, figsize=(11, 7), facecolor="white")
            axs = axs.flatten()
            fig.suptitle(f"Outlier Plots for Sample: {sample}")

            for i, col in enumerate(cols_to_plot):
                if i >= len(axs):
                    break  # Avoid index error if more cols than subplots
                ax = axs[i]
                if col in data_view.obs:
                    # Determine coloring based on data type
                    if data_view.obs[col].dtype == "bool" or set(
                        data_view.obs[col].unique()
                    ).issubset({0, 1, True, False}):
                        colors = data_view.obs[col].map(
                            {
                                False: "#655a5a",
                                True: "#e53639",
                                0: "#655a5a",
                                1: "#e53639",
                            }
                        )
                    else:
                        # Default color for non-boolean data
                        colors = "#655a5a"

                    ax.scatter(
                        data_view.obs["total_counts"],
                        data_view.obs["n_genes_by_counts"],
                        c=colors,
                        s=5,
                        edgecolor="none",
                        rasterized=True,
                    )
                    ax.set_title(col.replace("_", " ").title())
                    ax.set_xlabel("Total Counts")
                    ax.set_ylabel("n_genes_by_counts")

                    # Add count and percentage to title
                    if data_view.obs[col].dtype == "bool" or set(
                        data_view.obs[col].unique()
                    ).issubset({0, 1, True, False}):
                        count = data_view.obs[col].sum()
                        ax.set_title(
                            f"{col.replace('_', ' ').title()}\n{count} cells ({count / data_view.n_obs:.1%})"
                        )
                else:
                    ax.set_visible(False)

            # Hide unused subplots
            for j in range(len(cols_to_plot), len(axs)):
                axs[j].set_visible(False)

            plt.tight_layout(rect=[0, 0, 1, 0.96])
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                plt.savefig(
                    os.path.join(save_dir, f"{sample}_qc_outliers.png"),
                    dpi=300,
                    facecolor="white",
                )
            if show:
                plt.show()
            plt.close(fig)

    return adata


def filter_cells(
    adata: AnnData,
    criteria: Optional[List[str]] = None,
    copy: bool = False,
    filter_by_outlier_low_genes: bool = True,
    filter_by_outlier_mt: bool = True,
    filter_by_outlier_hb: bool = True,
    filter_by_outlier_qc_metrics: bool = True,
    filter_by_predicted_doublet_scrublet: bool = True,
    filter_by_doublet_expression_pattern: bool = True,
) -> Optional[AnnData]:
    """
    Filter cells based on previously calculated QC and doublet boolean flags.

    Args:
        adata: AnnData object with QC metrics calculated by functions like
               `is_low_quality_cell` and `is_doublet`.
        criteria: A list of boolean columns in adata.obs to use for filtering.
                  Cells where ANY of these columns is True will be removed.
                  If None, a default set of criteria is used based on filter_by_* parameters.
        copy: Whether to return a new filtered AnnData object.
        filter_by_*: Boolean flags to enable/disable specific filtering criteria.
                     Only used when criteria=None.

    Returns:
        Filtered AnnData object if copy=True, otherwise filters in place and returns None.
    """
    # If criteria is None, build it from filter_by_* parameters
    if criteria is None:
        criteria = []

        if filter_by_outlier_low_genes and "outlier_low_genes" in adata.obs.columns:
            criteria.append("outlier_low_genes")

        if filter_by_outlier_mt and "outlier_mt" in adata.obs.columns:
            criteria.append("outlier_mt")

        if filter_by_outlier_hb and "outlier_hb" in adata.obs.columns:
            criteria.append("outlier_hb")

        if filter_by_outlier_qc_metrics and "outlier_qc_metrics" in adata.obs.columns:
            criteria.append("outlier_qc_metrics")

        if (
            filter_by_predicted_doublet_scrublet
            and "predicted_doublet_scrublet" in adata.obs.columns
        ):
            criteria.append("predicted_doublet_scrublet")

        if (
            filter_by_doublet_expression_pattern
            and "doublet_expression_pattern" in adata.obs.columns
        ):
            criteria.append("doublet_expression_pattern")

    # Filter out criteria that don't exist in the object
    valid_criteria = [c for c in criteria if c in adata.obs.columns]
    missing_criteria = set(criteria) - set(valid_criteria)
    if missing_criteria:
        log.warning(
            f"Criteria not found in adata.obs and will be ignored: {missing_criteria}"
        )

    if not valid_criteria:
        log.warning("No valid filtering criteria selected. Returning original object.")
        return adata.copy() if copy else None

    initial_cells = adata.n_obs

    # Create a combined mask. A cell is kept if it's False for all criteria.
    combined_mask = pd.Series(False, index=adata.obs_names)

    # Track individual contributions
    criteria_counts = {}
    for col in valid_criteria:
        col_mask = adata.obs[col]
        criteria_counts[col] = col_mask.sum()
        combined_mask |= col_mask

    # The mask for keeping cells is the inverse of the removal mask
    keep_mask = ~combined_mask

    log.info(f"Filtering cells based on: {', '.join(valid_criteria)}")

    # Report individual criteria counts
    for col, count in criteria_counts.items():
        log.info(f"  - {col}: {count} cells ({count / initial_cells:.2%})")

    adata_to_return = adata[keep_mask, :].copy() if copy else None

    if not copy:
        # sc.AnnData._inplace_subset_obs
        adata._inplace_subset_obs(keep_mask)

    final_cells = adata.n_obs if not copy else adata_to_return.n_obs

    log.info(f"  Initial cell count: {initial_cells}")
    log.info(f"  Final cell count:   {final_cells}")
    log.info(
        f"  Cells removed:      {initial_cells - final_cells} ({(initial_cells - final_cells) / initial_cells:.2%})"
    )

    # Calculate cells removed by multiple criteria
    if len(valid_criteria) > 1:
        outlier_counts = adata.obs[valid_criteria].sum(axis=1)
        for i in range(1, len(valid_criteria) + 1):
            count = (outlier_counts == i).sum()
            if count > 0:
                log.info(
                    f"  Cells with exactly {i} types of issues: {count} ({count / initial_cells:.2%})"
                )

    return adata_to_return
