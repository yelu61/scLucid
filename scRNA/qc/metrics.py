"""
Quality control metrics calculation for single-cell RNA-seq data.

This module provides functions for calculating standard QC metrics and 
identifying outlier cells based on statistical methods.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import re

from typing import Dict, Literal

__all__ = [
    "calculate_qc_metric"
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

    total_samples = len(adata.obs[sample_key].unique())
    for i, sample in enumerate(adata.obs[sample_key].unique()):
        print(f"Processing sample {i+1}/{total_samples}: {sample}")
        data = adata[adata.obs[sample_key] == sample, :].copy()

        # Calculate the QC covariates or metric using provided patterns
        for gene_type, pattern in gene_patterns.items():
            data.var[gene_type] = data.var_names.str.contains(
                pattern, regex=True, na=False, flags=re.IGNORECASE
            )
        
        sc.pp.calculate_qc_metrics(
            data,
            qc_vars=list(gene_patterns.keys()),
            inplace=True,
            percent_top=[20],
            log1p=True,
        )

        # Plot the QC covariates per sample
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)
        plt.style.use('default')
        
        if plot_violin:
            fig, axes = plt.subplots(nrows=1, ncols=len(keys), figsize=(12, 3), facecolor='white')
            if len(keys) == 1:  # Handle single key case
                axes = [axes]
            for i, ax in enumerate(axes):
                ax.set_facecolor('white') 
                sc.pl.violin(data, keys[i], ax=ax, jitter=0.4, show=False)
                ax.set_title(f"Sample: {sample}", fontsize=10)
            plt.tight_layout()
            plt.show()
            if save_dir:
                import os
                os.makedirs(save_dir, exist_ok=True)
                fig.savefig(os.path.join(save_dir, f"{sample}_qc_violin.png"), dpi=300, facecolor='white')
                plt.close(fig)

        if plot_scatter:
            # Standard scatter plot
            fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')
            ax.set_facecolor('white')
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
                fig.savefig(os.path.join(save_dir, f"{sample}_qc_scatter.png"), dpi=300, facecolor='white')
                plt.close(fig)
            
            # Additional mt vs ribo scatter plot for deeper analysis
            fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')
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
                fig.savefig(os.path.join(save_dir, f"{sample}_mt_ribo_scatter.png"), dpi=300, facecolor='white')
                plt.close(fig)

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