"""
Quality control metrics calculation for single-cell RNA-seq data.

This module provides functions for calculating standard QC metrics and 
identifying outlier cells based on statistical methods.
"""

import matplotlib.pyplot as plt
import scanpy as sc

from typing import Optional, Dict

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
    gene_patterns: Optional[Dict[str, str]] = None,
    plot_violin: bool = True,
    plot_scatter: bool = True,
    save_dir: str = None,
    show: bool = True,
) -> sc.AnnData:
    """
    Calculate and plot QC metrics for each sample in the AnnData object.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        keys (list, optional): List of QC metrics to plot.
        gene_patterns (dict, optional): Dictionary of gene patterns to identify mitochondrial, ribosomal, and hemoglobin genes.
            Defaults to None, which uses standard patterns.
            If None, the following patterns are used:
                - "mt": r"^(MT|Mt|mt)-"
                - "ribo": r"^(RP[SL]|Rp[sl])"
                - "hb": r"^(HB|hb)[^(P|p)]"
        plot_violin (bool, optional): Whether to plot violin plots for QC metrics. Defaults to True.
        plot_scatter (bool, optional): Whether to plot scatter plot for total_counts vs n_genes_by_counts. Defaults to True.
        save_dir (str, optional): Directory to save plots. If None, plots are not saved to disk. Defaults to None.
        show (bool, optional): Whether to display the plots. Defaults to True.
            If False, plots are saved to disk if save_dir is provided.
            If True, plots are displayed in the notebook.
    Returns:
        adata (AnnData): AnnData object with QC metrics added.
    """
    # Check if sample key exists
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample key '{sample_key}' is not in adata.obs columns")
    # Ensure there is at least one sample
    if len(adata.obs[sample_key].unique()) == 0:
        raise ValueError(f"No samples found under sample key '{sample_key}'")

    # --- Main Calculation (Optimized) ---
    if gene_patterns is None:
        # Default gene patterns if none provided
        gene_patterns = {
            "mt": r"^(MT|Mt|mt)-",
            "ribo": r"^(RP[SL]|Rp[sl])",
            "hb": r"^(HB|hb)[^(P|p)]"
        }

    print("Calculating QC metrics for the entire dataset...")
    # Calculate metrics for all cells at once
    for gene_type, pattern in gene_patterns.items():
        adata.var[gene_type] = adata.var_names.str.contains(pattern, regex=True, na=False)
    
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=list(gene_patterns.keys()),
        inplace=True,
        percent_top=[20],
        log1p=True
    )
    print("QC metrics calculation complete.")

    # --- Plotting (Remains per-sample) ---
    if plot_violin or plot_scatter:
        if save_dir:
            import os
            os.makedirs(save_dir, exist_ok=True) # Create directory once
            
        samples = adata.obs[sample_key].unique()
        for sample in samples:
            print(f"Plotting QC for sample: {sample}")
            data_view = adata[adata.obs[sample_key] == sample]
            
            if plot_violin:
                fig, axes = plt.subplots(nrows=1, ncols=len(keys), figsize=(12, 3), facecolor='white')
                if len(keys) == 1:  # Handle single key case
                    axes = [axes]
                for i, ax in enumerate(axes):
                    ax.set_facecolor('white') 
                    sc.pl.violin(data_view, keys[i], ax=ax, jitter=0.4, show=False)
                    ax.set_title(f"Sample: {sample}", fontsize=10)
                plt.tight_layout()
                if save_dir:
                    fig.savefig(os.path.join(save_dir, f"{sample}_qc_violin.png"), dpi=300, facecolor='white')
                if show:
                    plt.show()
                plt.close(fig)
            if plot_scatter:
                fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')
                ax.set_facecolor('white')
                sc.pl.scatter(
                    data_view,
                    x="total_counts",
                    y="n_genes_by_counts",
                    color="pct_counts_mt",
                    title=f"Sample: {sample} - Basic QC",
                    show=False,
                    ax=ax,
                    legend_loc="right margin",
                )
                plt.tight_layout()
                if save_dir:
                    fig.savefig(os.path.join(save_dir, f"{sample}_qc_scatter.png"), dpi=300, facecolor='white')
                if show:
                    plt.show()
                plt.close(fig)
                           
    return adata