"""
Cell filtering utilities for single-cell RNA-seq data.

This module provides functions for identifying and filtering low-quality
cells based on various quality metrics.
"""

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from typing import List, Tuple

__all__ = [
    "is_low_quality_cell", 
    "filter_low_quality_cells"
]
def is_low_quality_cell(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    min_genes: int = 200,
    nmad: int = 5,
    pc_mt: int = 20,
    pc_hb: int = 20,
    plot_outliers: bool = False,
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
        outlier_metrics (List[Tuple[str, str]], optional): List of (metric, direction) tuples for outlier detection.
                                                         If None, default metrics will be used. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with low-quality cells identified.
    """
    # Check if required columns exist
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in adata.obs: {missing_cols}. "
                         f"Run calculate_qc_metric() first.")
        
    # Define default outlier metrics if not provided
    if outlier_metrics is None:
        outlier_metrics = [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
            ("pct_counts_in_top_20_genes", "upper")
        ]
    
    # Precompute sample indexes to avoid repeated filtering
    sample_indices = {}
    for sample in adata.obs[sample_key].unique():
        sample_indices[sample] = adata.obs[sample_key] == sample
    
    # Mark cells with low gene counts instead of filtering
    adata.obs["low_genes_outlier"] = adata.obs["n_genes_by_counts"] < min_genes
    print(f"Cells with low gene counts (< {min_genes}): {adata.obs['low_genes_outlier'].sum()}")
    
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
            metric_outliers = identify_outliers(data, metric=metric, nmads=nmad, direction=direction)
            outlier_mask = outlier_mask | metric_outliers
            print(f"  {metric} outliers ({direction}): {metric_outliers.sum()}")
        
        data.obs["outlier"] = outlier_mask
        print(
            f"Outlier cells (combined): {data.obs.outlier.sum()}"
        )

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
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.scatter(
                data, 
                x="total_counts", 
                y="n_genes_by_counts",
                color=["outlier", "mt_outlier", "hb_outlier", "low_genes_outlier"],
                title=f"Outliers in sample: {sample}",
                show=False,
                ax=ax
            )
            plt.tight_layout()
            plt.show()
            
    return adata
def filter_low_quality_cells(
    adata: sc.AnnData,
    filter_low_genes: bool = True,
    filter_outliers: bool = True,
    filter_mt: bool = True,
    filter_hb: bool = True,
    filter_doublets: bool = False,  # Default: do not filter doublets
    only_predicted_final: bool = True,  # Only filter final predicted doublets
) -> sc.AnnData:
    """
    Filter low-quality cells based on quality control results.
    
    This function completely removes cells that meet the specified filtering criteria,
    returning a new AnnData object with only high-quality cells.
    
    Args:
        adata: AnnData object with QC metrics calculated (use is_low_quality_cell first)
        filter_low_genes: Whether to filter cells with low gene counts
        filter_outliers: Whether to filter outlier cells
        filter_mt: Whether to filter cells with high mitochondrial content
        filter_hb: Whether to filter cells with high hemoglobin content
        filter_doublets: Whether to filter doublet cells
        only_predicted_final: If filtering doublets, whether to only filter final predicted doublets
        
    Returns:
        adata_filtered: A new AnnData object with only high-quality cells
        
    Example:
        >>> adata = scRNA.qc.calculate_qc_metric(adata, sample_key="batch")
        >>> adata = scRNA.qc.is_low_quality_cell(adata, sample_key="batch")
        >>> adata_filtered = scRNA.qc.filter_low_quality_cells(
        ...     adata, filter_outliers=True, filter_mt=True, filter_doublets=True
        ... )
    """
    mask = np.ones(adata.n_obs, dtype=bool)
    
    if filter_low_genes and "low_genes_outlier" in adata.obs.columns:
        mask &= ~adata.obs["low_genes_outlier"].values
    
    if filter_outliers:
        mask &= ~adata.obs["outlier"].values
    
    if filter_mt:
        mask &= ~adata.obs["mt_outlier"].values
    
    if filter_hb:
        mask &= ~adata.obs["hb_outlier"].values
    
    if filter_doublets:
        if "predicted_doublets_final" not in adata.obs.columns:
            print("Warning: Doublet detection has not been run. Skipping doublet filtering.")
        else:
            if only_predicted_final:
                mask &= ~adata.obs["predicted_doublets_final"].values
            else:
                mask &= ~(adata.obs["predicted_doublets"].values | 
                         adata.obs["predicted_doublets_final"].values | 
                         adata.obs["overexpressed_doublets"].values)
    
    # Calculate the number of filtered cells
    cells_before = adata.n_obs
    adata_filtered = adata[mask].copy()
    cells_after = adata_filtered.n_obs
    cells_filtered = cells_before - cells_after
    
    print(f"Cells before filtering: {cells_before}")
    print(f"Cells after filtering: {cells_after}")
    print(f"Filtered cells: {cells_filtered} ({cells_filtered/cells_before:.2%})")
    
    # Detailed report of filtered cell types
    if filter_low_genes and "low_genes_outlier" in adata.obs.columns:
        n_low_genes = adata.obs["low_genes_outlier"].sum()
        print(f"  Low gene count cells: {n_low_genes} ({n_low_genes/cells_before:.2%})")
    
    if filter_outliers:
        n_outliers = adata.obs["outlier"].sum()
        print(f"  Outlier cells: {n_outliers} ({n_outliers/cells_before:.2%})")
    
    if filter_mt:
        n_mt = adata.obs["mt_outlier"].sum()
        print(f"  High mitochondrial cells: {n_mt} ({n_mt/cells_before:.2%})")
    
    if filter_hb:
        n_hb = adata.obs["hb_outlier"].sum()
        print(f"  High hemoglobin cells: {n_hb} ({n_hb/cells_before:.2%})")
    
    if filter_doublets and "predicted_doublets_final" in adata.obs.columns:
        if only_predicted_final:
            n_doublets = adata.obs["predicted_doublets_final"].sum()
            print(f"  Final predicted doublets: {n_doublets} ({n_doublets/cells_before:.2%})")
        else:
            n_doublets = (adata.obs["predicted_doublets"].values | 
                         adata.obs["predicted_doublets_final"].values | 
                         adata.obs["overexpressed_doublets"].values).sum()
            print(f"  All types of doublets: {n_doublets} ({n_doublets/cells_before:.2%})")
    
    # Add simple visualization of filtering results
    if cells_filtered > 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        categories = ["Low genes", "Outliers", "High MT", "High HB", "Doublets"]
        counts = [
            adata.obs["low_genes_outlier"].sum() if filter_low_genes and "low_genes_outlier" in adata.obs.columns else 0,
            adata.obs["outlier"].sum() if filter_outliers else 0,
            adata.obs["mt_outlier"].sum() if filter_mt else 0,
            adata.obs["hb_outlier"].sum() if filter_hb else 0,
            adata.obs["predicted_doublets_final"].sum() if filter_doublets and "predicted_doublets_final" in adata.obs.columns else 0
        ]
        ax.bar(categories, counts)
        ax.set_ylabel("Number of cells")
        ax.set_title("Filtered cells by category")
        plt.tight_layout()
        plt.show()
    
    return adata_filtered