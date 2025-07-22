"""
Cell filtering utilities for single-cell RNA-seq data.

This module provides functions for identifying and filtering low-quality
cells based on various quality metrics.
"""

import os
from typing import List, Literal, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

__all__ = [
    "is_low_quality_cell",
    "filter_cells",
    "_identify_outliers",
]


# --- Helper functions ---#
def _identify_outliers(
    adata: sc.AnnData,
    metric: str,
    nmads: int,
    direction: Literal["both", "upper", "lower"] = "both",
) -> pd.Series:
    """
    Identify outliers based on the given metric and number of median absolute deviations.

    Args:
        adata (AnnData): AnnData object to check for outliers.
        metric (str): The metric to use for outlier detection. Must be a valid column in adata.obs.
        nmads (int): Number of median absolute deviations for outlier detection.
        direction (str): Direction for outlier detection. Options are 'both', 'upper', or 'lower'.
                        'both': Detects values that are too high or too low
                        'upper': Only detects values that are too high
                        'lower': Only detects values that are too low

    Returns:
        outliers (pandas.Series): Boolean mask indicating if a cell is an outlier or not.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Invalid metric '{metric}'. Must be a column in adata.obs.")

    if nmads <= 0:
        raise ValueError(f"nmads must be positive, got {nmads}")

    if direction not in ["both", "upper", "lower"]:
        raise ValueError(
            f"direction must be one of 'both', 'upper', or 'lower', got {direction}"
        )

    values = adata.obs[metric].copy()

    if values.isna().any():
        print(
            f"Warning: {values.isna().sum()} NaN values in {metric}, will be excluded from outlier detection"
        )
        values = values.dropna()

    # Add informative message about distribution
    median = np.median(values)
    mad = median_abs_deviation(values)
    print(f"Distribution info for {metric}: median={median:.2f}, MAD={mad:.2f}")

    if mad == 0:
        print(f"Warning: MAD=0 for {metric}, no outliers will be detected")
        return pd.Series(False, index=adata.obs_names)

    outliers = pd.Series(False, index=adata.obs_names)

    if direction == "both":
        outliers.loc[values.index] = [
            abs(value - median) > nmads * mad for value in values
        ]
    elif direction == "upper":
        outliers.loc[values.index] = [value > median + nmads * mad for value in values]
    elif direction == "lower":
        outliers.loc[values.index] = [value < median - nmads * mad for value in values]

    print(
        f"Identified {outliers.sum()} outliers in {metric} with nmads={nmads} in direction '{direction}'"
    )
    return outliers


# --- Main functions ---#
def is_low_quality_cell(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    min_genes: int = 200,
    nmad: int = 5,
    pc_mt: int = 20,
    pc_hb: int = 20,
    plot_outliers: bool = False,
    save_dir: str = None,
    show: bool = True,
    outlier_metrics: List[Tuple[str, str]] = None,
) -> sc.AnnData:
    """
    Identify low-quality cells based on various quality metrics.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        min_genes (int, optional): Minimum number of genes expressed in a cell. Defaults to 200.
        nmad (int, optional): Number of median absolute deviations for outlier detection. Defaults to 5.
        pc_mt (int, optional): Maximum percentage of mitochondrial counts allowed. Defaults to 20.
        pc_hb (int, optional): Maximum percentage of hemoglobin counts allowed. Defaults to 20.
        plot_outliers (bool, optional): Whether to plot outlier cells. Defaults to False.
        save_dir (str, optional): Directory to save plots. Defaults to None.
        show (bool, optional): Whether to show plots. Defaults to True.
        outlier_metrics (List[Tuple[str, str]], optional): List of (metric, direction) tuples for outlier detection.
                                                         If None, default metrics will be used. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with low-quality cells identified.
    """
    # Check if required columns exist
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns in adata.obs: {missing_cols}. "
            f"Run calculate_qc_metric() first."
        )

    # Define default outlier metrics if not provided
    if outlier_metrics is None:
        outlier_metrics = [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
            ("pct_counts_in_top_20_genes", "upper"),
        ]

    # Precompute sample indexes to avoid repeated filtering
    sample_indices = {}
    for sample in adata.obs[sample_key].unique():
        sample_indices[sample] = adata.obs[sample_key] == sample

    # Mark cells with low gene counts instead of filtering
    adata.obs["low_genes_outlier"] = adata.obs["n_genes_by_counts"] < min_genes
    print(
        f"Cells with low gene counts (< {min_genes}): {adata.obs['low_genes_outlier'].sum()}"
    )

    # Initialize the result column
    adata.obs["outlier"] = False
    adata.obs["mt_outlier"] = False
    adata.obs["hb_outlier"] = False

    for sample, indices in sample_indices.items():
        # Using precomputed indexes
        if sum(indices) == 0:
            continue  # Skip if no cells for this sample

        print(f"QC of low quality cells for sample: {sample}")
        data = adata[indices, :]

        # Identify outlier cells
        outlier_mask = pd.Series(False, index=data.obs_names)
        for metric, direction in outlier_metrics:
            metric_outliers = _identify_outliers(
                data, metric=metric, nmads=nmad, direction=direction
            )
            outlier_mask = outlier_mask | metric_outliers
            print(f"  {metric} outliers ({direction}): {metric_outliers.sum()}")

        data.obs["outlier"] = outlier_mask
        print(f"Outlier cells (combined): {data.obs.outlier.sum()}")

        # Identify cells with high mitochondrial percentage
        data.obs["mt_outlier"] = data.obs["pct_counts_mt"] > pc_mt
        print(
            f"Cells with high mitochondrial percentage (> {pc_mt}%): {data.obs.mt_outlier.sum()}"
        )

        # Identify cells with high hemoglobin percentage
        data.obs["hb_outlier"] = data.obs["pct_counts_hb"] > pc_hb
        print(
            f"Cells with high hemoglobin percentage (> {pc_hb}%): {data.obs.hb_outlier.sum()}"
        )

        adata.obs.loc[data.obs.index, "outlier"] = data.obs["outlier"]
        adata.obs.loc[data.obs.index, "mt_outlier"] = data.obs["mt_outlier"]
        adata.obs.loc[data.obs.index, "hb_outlier"] = data.obs["hb_outlier"]

    # Print the overall statistics for the entire adata object
    print("\nOverall statistics for the entire adata object:")
    total_cells = adata.n_obs
    print(f"Total number of cells: {total_cells}")

    outlier_cells = adata.obs.outlier.sum()
    print(
        f"Outlier cells (based on various metrics): {outlier_cells} ({outlier_cells / total_cells * 100:.2f}%)"
    )

    low_genes_cells = adata.obs.low_genes_outlier.sum()
    print(
        f"Cells with low gene counts (< {min_genes}): {low_genes_cells} ({low_genes_cells / total_cells * 100:.2f}%)"
    )

    mt_outlier_cells = adata.obs.mt_outlier.sum()
    print(
        f"Cells with high mitochondrial percentage (> {pc_mt}%): {mt_outlier_cells} ({mt_outlier_cells / total_cells * 100:.2f}%)"
    )

    hb_outlier_cells = adata.obs.hb_outlier.sum()
    print(
        f"Cells with high hemoglobin percentage (> {pc_hb}%): {hb_outlier_cells} ({hb_outlier_cells / total_cells * 100:.2f}%)"
    )

    doublets_outlier_cells = adata.obs.predicted_doublets_final.sum()
    print(
        f"Predicted doublets: {doublets_outlier_cells} ({doublets_outlier_cells / total_cells * 100:.2f}%)"
    )

    overexpression_outlier_cells = adata.obs.overexpressed_doublets.sum()
    print(
        f"Potential doublets with high overexpression: {overexpression_outlier_cells} ({overexpression_outlier_cells / total_cells * 100:.2f}%)"
    )

    combined_outliers = adata.obs.filter(regex="outlier").sum(axis=1).value_counts()
    for n_outliers, count in sorted(combined_outliers.items()):
        print(
            f"{n_outliers} types of outliers: {count} ({count / total_cells * 100:.2f}%)"
        )

    if plot_outliers:
        for sample in adata.obs[sample_key].unique():
            data = adata[adata.obs[sample_key] == sample, :]

            # Plotting outliers
            print(f"Plotting outliers for sample: {sample}")
            fig, axs = plt.subplots(2, 3, figsize=(11, 7), facecolor="white")
            axs = axs.flatten()
            fig.suptitle(f"Outlier Plots for Sample: {sample}")
            for i, col in enumerate(
                [
                    "outlier",
                    "mt_outlier",
                    "hb_outlier",
                    "low_genes_outlier",
                    "predicted_doublets_final",
                    "overexpressed_doublets",
                ]
            ):
                if col in data.obs:
                    sc.pl.scatter(
                        data,
                        x="total_counts",
                        y="n_genes_by_counts",
                        color=col,
                        ax=axs[i],
                        show=False,
                        title=col.replace("_", " ").title(),
                        # palette=["#bbbbbb", "#e41a1dbf"],
                        color_map=None,
                    )
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                fig.savefig(
                    os.path.join(save_dir, f"{sample}_qc_outliers.png"),
                    dpi=300,
                    facecolor="white",
                )
            if show:
                plt.show()
            plt.close(fig)
    return adata


def filter_cells(
    adata: sc.AnnData,
    filter_by_outliers: bool = True,
    filter_by_low_genes: bool = True,
    filter_by_mt: bool = True,
    filter_by_hb: bool = False,  # Often not default
    filter_by_doublets: bool = True,
    filter_by_overexpression: bool = True,
    copy: bool = False,
) -> sc.AnnData:
    """
    Filter cells based on previously calculated QC and doublet metrics.

    This function acts as a unified wrapper to remove cells marked by
    `is_low_quality_cell` and `is_doublet`.

    Args:
        adata: AnnData object with QC metrics calculated.
        filter_by_*: Flags to determine which criteria to use for filtering.
        copy: Whether to return a copy or filter in place.

    Returns:
        Filtered AnnData object.
    """

    initial_cells = adata.n_obs
    cell_mask = pd.Series(True, index=adata.obs.index)

    reasons = []

    if filter_by_outliers and "outlier" in adata.obs.columns:
        cell_mask &= ~adata.obs["outlier"]
        reasons.append("outliers")

    if filter_by_low_genes and "low_genes_outlier" in adata.obs.columns:
        cell_mask &= ~adata.obs["low_genes_outlier"]
        reasons.append("low gene counts")

    if filter_by_mt and "mt_outlier" in adata.obs.columns:
        cell_mask &= ~adata.obs["mt_outlier"]
        reasons.append("high mitochondrial %")

    if filter_by_hb and "hb_outlier" in adata.obs.columns:
        cell_mask &= ~adata.obs["hb_outlier"]
        reasons.append("high hemoglobin %")

    if filter_by_doublets and "predicted_doublets_final" in adata.obs.columns:
        cell_mask &= ~adata.obs["predicted_doublets_final"]
        reasons.append("predicted doublets")

    if filter_by_overexpression and "overexpressed_doublets" in adata.obs.columns:
        cell_mask &= ~adata.obs["overexpressed_doublets"]
        reasons.append("overexpressed_doublets")

    if not reasons:
        print("No filtering criteria selected. Returning original object.")
        return adata.copy() if copy else adata

    print(f"Filtering cells based on: {', '.join(reasons)}")

    adata_filtered = adata[cell_mask, :].copy() if copy else adata[cell_mask, :]

    print(f"  Initial cell count: {initial_cells}")
    print(f"  Final cell count:   {adata_filtered.n_obs}")
    print(
        f"  Cells removed:      {initial_cells - adata_filtered.n_obs} ({(initial_cells - adata_filtered.n_obs) / initial_cells:.2%})"
    )

    return adata_filtered
