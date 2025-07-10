"""
Quality control module for single-cell RNA-seq data.
"""

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.stats import median_abs_deviation
from typing import Dict, Literal

__all__ = [
    "calculate_qc_metric", 
    "identify_outliers"
]

def calculate_qc_metric(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    keys: list = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ],
    plot_violin: bool = True,
    plot_scatter: bool = True,
    save_dir: str = None,
    gene_patterns: Dict[str, str] = None,
) -> sc.AnnData:
    """
    Calculate and plot QC metrics for each sample in the AnnData object.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        keys (list, optional): List of QC metrics to plot.
        plot_violin (bool, optional): Whether to plot violin plots for QC metrics. Defaults to True.
        plot_scatter (bool, optional): Whether to plot scatter plot for total_counts vs n_genes_by_counts. Defaults to True.
        save_dir (str, optional): Directory to save plots. If None, plots are not saved to disk. Defaults to None.
        gene_patterns (Dict[str, str], optional): Dictionary of regex patterns for special gene groups. 
                                                If None, default patterns will be used. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with QC metrics added.
    """
    # Check if sample key exists
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample key '{sample_key}' is not in adata.obs columns")
    # Ensure there is at least one sample
    if len(adata.obs[sample_key].unique()) == 0:
        raise ValueError(f"No samples found under sample key '{sample_key}'")
    
    # Default gene patterns if not provided
    if gene_patterns is None:
        gene_patterns = {
            "mt": r"^(MT-|mt-|Mt-)",  # mitochondrial genes
            "ribo": r"^(RPS|RPL|Rps|Rpl|Gm\d+)",  # ribosomal genes
            "hb": r"^(HB|Hb)[^(P|p)]"  # hemoglobin genes
        }
    
    adata.obs["total_counts"] = 0
    adata.obs["log1p_total_counts"] = 0
    adata.obs["n_genes_by_counts"] = 0
    adata.obs["log1p_n_genes_by_counts"] = 0
    adata.obs["pct_counts_in_top_20_genes"] = 0
    adata.obs["pct_counts_mt"] = 0.0
    adata.obs["pct_counts_ribo"] = 0.0
    adata.obs["pct_counts_hb"] = 0.0

    for sample in adata.obs[sample_key].unique():
        print(f"Begin of QC metric calculation and QC plot for sample: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        # Calculate the QC covariates or metric using provided patterns
        for gene_type, pattern in gene_patterns.items():
            data.var[gene_type] = data.var_names.str.contains(
                pattern, regex=True, na=False
            )
        
        sc.pp.calculate_qc_metrics(
            data,
            qc_vars=list(gene_patterns.keys()),
            inplace=True,
            percent_top=[20],
            log1p=True,
        )

        # Plot the QC covariates per sample
        if plot_violin:
            fig, axes = plt.subplots(nrows=1, ncols=len(keys), figsize=(12, 3))
            if len(keys) == 1:  # Handle single key case
                axes = [axes]
            for i, ax in enumerate(axes):
                sc.pl.violin(data, keys[i], ax=ax, jitter=0.4, show=False)
                ax.set_title(f"Sample: {sample}", fontsize=10)
            plt.tight_layout()
            plt.show()
            if save_dir:
                import os
                os.makedirs(save_dir, exist_ok=True)
                plt.savefig(os.path.join(save_dir, f"{sample}_qc_violin.png"), dpi=300)
                plt.close()

        if plot_scatter:
            # Standard scatter plot
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.scatter(
                data,
                x="total_counts",
                y="n_genes_by_counts",
                color="pct_counts_mt",
                title=f"Sample: {sample} - Basic QC",
                show=False,
                ax=ax,
                legend_loc="right margin",
            )
            plt.tight_layout()
            plt.show()
            if save_dir:
                import os
                os.makedirs(save_dir, exist_ok=True)
                plt.savefig(os.path.join(save_dir, f"{sample}_qc_scatter.png"), dpi=300)
                plt.close()
            
            # Additional mt vs ribo scatter plot for deeper analysis
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.scatter(
                data,
                x="pct_counts_mt",
                y="pct_counts_ribo",
                color="total_counts",
                title=f"Sample: {sample} - MT vs Ribo",
                show=False,
                ax=ax,
                legend_loc="right margin",
            )
            plt.tight_layout()
            plt.show()
            if save_dir:
                import os
                os.makedirs(save_dir, exist_ok=True)
                plt.savefig(os.path.join(save_dir, f"{sample}_mt_ribo_scatter.png"), dpi=300)
                plt.close()
        
        metric_cols = [
            "total_counts",
            "log1p_total_counts",
            "n_genes_by_counts",
            "log1p_n_genes_by_counts",
            "pct_counts_in_top_20_genes",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
        ]
        adata.obs.loc[data.obs.index, metric_cols] = data.obs[metric_cols]
        print("Done.")
                
    return adata


def identify_outliers(
    adata: sc.AnnData, 
    metric: str, 
    nmads: int,
    direction: Literal["both", "upper", "lower"] = "both"
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
        raise ValueError(f"direction must be one of 'both', 'upper', or 'lower', got {direction}")
    
    values = adata.obs[metric].copy()
    
    if values.isna().any():
        print(f"Warning: {values.isna().sum()} NaN values in {metric}, will be excluded from outlier detection")
        values = values.dropna()
    
    median = np.median(values)
    mad = median_abs_deviation(values)
    
    if mad == 0:
        print(f"Warning: MAD=0 for {metric}, no outliers will be detected")
        return pd.Series(False, index=adata.obs_names)
    
    outliers = pd.Series(False, index=adata.obs_names)
    
    if direction == "both":
        outliers.loc[values.index] = [abs(value - median) > nmads * mad for value in values]
    elif direction == "upper":
        outliers.loc[values.index] = [value > median + nmads * mad for value in values]
    elif direction == "lower":
        outliers.loc[values.index] = [value < median - nmads * mad for value in values]
    
    print(f"Identified {outliers.sum()} outliers in {metric} with nmads={nmads} in direction '{direction}'")
    return outliers