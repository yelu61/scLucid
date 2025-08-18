"""
Quality control metrics calculation for single-cell RNA-seq data.

This module provides functions for calculating standard QC metrics and
identifying outlier cells based on statistical methods.
"""

import logging
import os
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
log = logging.getLogger(__name__)

__all__ = [
    "calculate_qc_metric",
    "_plot_top20_genes_distribution",
    "_plot_qc_violin",
    "_plot_qc_scatter",
]


def calculate_qc_metric(
    adata: AnnData,
    sample_key: str = "sampleID",
    keys: List[str] = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ],
    gene_patterns: Optional[Dict[str, str]] = None,
    plot_violin: bool = True,
    plot_scatter: bool = True,
    plot_top20: bool = True,  # New parameter to control top20 genes plotting
    save_dir: Optional[str] = None,
    show: bool = True,
) -> AnnData:
    """
    Calculate and plot QC metrics for each sample in the AnnData object.

    Args:
        adata: AnnData object containing single-cell data.
        sample_key: The key in adata.obs to identify different samples.
        keys: List of QC metrics to plot.
        gene_patterns: Dictionary of gene patterns to identify specific gene sets.
            If None, defaults are used for mitochondrial, ribosomal, and hemoglobin genes.
        plot_violin: Whether to plot violin plots for QC metrics.
        plot_scatter: Whether to plot scatter plot for total_counts vs n_genes_by_counts.
        plot_top20: Whether to plot distribution of pct_counts_in_top_20_genes.
        save_dir: Directory to save plots. If None, plots are not saved.
        show: Whether to display the plots.

    Returns:
        AnnData object with QC metrics added to .obs and .var.
    """
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not in adata.obs columns.")
    if not adata.obs[sample_key].nunique():
        raise ValueError(f"No samples found under sample key '{sample_key}'.")

    # --- Main Calculation ---
    if gene_patterns is None:
        gene_patterns = {
            "mt": r"^(MT|Mt|mt)-",
            "ribo": r"^(RP[SL]|Rp[sl])",
            "hb": r"^(HB|hb)[^(P|p)]",
        }

    log.info("Calculating QC metrics for the entire dataset...")
    for gene_type, pattern in gene_patterns.items():
        adata.var[gene_type] = adata.var_names.str.contains(
            pattern, regex=True, na=False
        )

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=list(gene_patterns.keys()),
        inplace=True,
        percent_top=[20],
        log1p=True,
    )
    log.info("QC metrics calculation complete.")

    # --- Plotting (Remains per-sample) ---
    if plot_violin or plot_scatter or plot_top20:
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

        samples = adata.obs[sample_key].unique()

        # Plot top 20 genes distribution if requested
        if plot_top20:
            _plot_top20_genes_distribution(
                adata, sample_key=sample_key, save_dir=save_dir, show=show
            )

        for sample in samples:
            log.info(f"Plotting QC for sample: {sample}")
            data_view = adata[adata.obs[sample_key] == sample]

            if plot_violin:
                _plot_qc_violin(data_view, keys, sample, save_dir=save_dir, show=show)

            if plot_scatter:
                _plot_qc_scatter(data_view, sample, save_dir=save_dir, show=show)

    return adata


def _plot_top20_genes_distribution(
    adata: AnnData,
    sample_key: str = "sampleID",
    thresholds: List[float] = [50, 60, 65, 70],
    save_dir: Optional[str] = None,
    show: bool = True,
) -> Optional[Tuple[plt.Figure, plt.Figure, plt.Figure]]:
    """
    Plot detailed distribution of pct_counts_in_top_20_genes to help determine appropriate thresholds.

    Args:
        adata: AnnData object with QC metrics calculated.
        sample_key: The key in adata.obs to identify different samples.
        thresholds: List of thresholds to mark in the plots.
        save_dir: Directory to save plots. If None, plots are not saved.
        show: Whether to display the plots.

    Returns:
        Tuple of Figure objects for the three plots (histogram, boxplot, scatter) if show=True
    """
    if "pct_counts_in_top_20_genes" not in adata.obs.columns:
        log.warning(
            "pct_counts_in_top_20_genes not found in adata.obs. "
            "Please run calculate_qc_metric first."
        )
        return None

    log.info("Analyzing distribution of pct_counts_in_top_20_genes...")

    # Calculate summary statistics
    stats = adata.obs.groupby(sample_key)["pct_counts_in_top_20_genes"].describe()
    log.info(f"Summary statistics for pct_counts_in_top_20_genes by sample:\n{stats}")

    # For each threshold, calculate percentage of cells above it
    for threshold in thresholds:
        counts = adata.obs.groupby(sample_key)["pct_counts_in_top_20_genes"].apply(
            lambda x: (x > threshold).sum()
        )
        percentages = counts / adata.obs.groupby(sample_key).size() * 100
        threshold_stats = pd.DataFrame({"counts": counts, "percentage": percentages})
        log.info(f"Cells above {threshold}% threshold:\n{threshold_stats}")

    # 1. Histogram with density curves by sample
    fig1, ax1 = plt.subplots(figsize=(12, 8), facecolor="white")

    samples = adata.obs[sample_key].unique()
    for sample in samples:
        sample_data = adata[adata.obs[sample_key] == sample].obs
        sns.kdeplot(
            data=sample_data,
            x="pct_counts_in_top_20_genes",
            label=f"Sample {sample}",
            ax=ax1,
        )

    # Add vertical lines for thresholds
    for threshold in thresholds:
        ax1.axvline(
            x=threshold,
            linestyle="--",
            color="red",
            alpha=0.7,
            label=f"{threshold}% threshold",
        )

    # Add percentile markers (90th, 95th, 99th)
    percentiles = [90, 95, 99]
    for p in percentiles:
        val = np.percentile(adata.obs["pct_counts_in_top_20_genes"], p)
        ax1.axvline(
            x=val,
            linestyle=":",
            color="green",
            alpha=0.7,
            label=f"{p}th percentile: {val:.1f}%",
        )

    ax1.set_title(
        "Distribution of Percentage Counts in Top 20 Genes by Sample", fontsize=14
    )
    ax1.set_xlabel("Percentage of counts in top 20 genes (%)", fontsize=12)
    ax1.set_ylabel("Density", fontsize=12)
    ax1.legend(title="Sample / Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            os.path.join(save_dir, "pct_top20_genes_distribution.png"),
            dpi=300,
            bbox_inches="tight",
        )

    # 2. Boxplot by sample
    fig2, ax2 = plt.subplots(figsize=(12, 8), facecolor="white")

    sns.boxplot(data=adata.obs, x=sample_key, y="pct_counts_in_top_20_genes", ax=ax2)

    # Add horizontal lines for thresholds
    for threshold in thresholds:
        ax2.axhline(
            y=threshold,
            linestyle="--",
            color="red",
            alpha=0.7,
            label=f"{threshold}% threshold",
        )

    ax2.set_title("Boxplot of Percentage Counts in Top 20 Genes by Sample", fontsize=14)
    ax2.set_xlabel("Sample", fontsize=12)
    ax2.set_ylabel("Percentage of counts in top 20 genes (%)", fontsize=12)
    ax2.legend(title="Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45)
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            os.path.join(save_dir, "pct_top20_genes_boxplot.png"),
            dpi=300,
            bbox_inches="tight",
        )

    # 3. Scatter plot: pct_counts_in_top_20_genes vs n_genes_by_counts
    fig3, ax3 = plt.subplots(figsize=(10, 8), facecolor="white")

    sc.pl.scatter(
        adata,
        x="n_genes_by_counts",
        y="pct_counts_in_top_20_genes",
        color=sample_key,
        ax=ax3,
        show=False,
    )

    # Add horizontal lines for thresholds
    for threshold in thresholds:
        ax3.axhline(
            y=threshold,
            linestyle="--",
            color="red",
            alpha=0.7,
            label=f"{threshold}% threshold",
        )

    ax3.set_title(
        "Relationship between Gene Counts and Top 20 Genes Percentage", fontsize=14
    )
    ax3.set_xlabel("Number of genes detected", fontsize=12)
    ax3.set_ylabel("Percentage of counts in top 20 genes (%)", fontsize=12)
    plt.legend(title="Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            os.path.join(save_dir, "pct_top20_genes_scatter.png"),
            dpi=300,
            bbox_inches="tight",
        )

    log.info("Analysis of pct_counts_in_top_20_genes complete.")

    if show:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)

    return (fig1, fig2, fig3) if show else None


def _plot_qc_violin(data, keys, sample, save_dir=None, show=False):
    fig, axs = plt.subplots(1, len(keys), figsize=(15, 4), facecolor="white")
    if len(keys) == 1:  # Ensure axs is always iterable
        axs = [axs]

    for ax, key in zip(axs, keys):
        sc.pl.violin(data, key, ax=ax, show=False)
        ax.set_title(key.replace("_", " ").title())
        ax.set_ylabel(key)
        plt.setp(ax.get_xticklabels(), visible=True)

    fig.suptitle(f"QC Metrics for Sample: {sample}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_dir:
        plt.savefig(
            os.path.join(save_dir, f"{sample}_qc_violin.png"),
            dpi=300,
            bbox_inches="tight",
        )
    if show:
        plt.show()
    plt.close(fig)


def _plot_qc_scatter(data, sample, save_dir=None, show=False):
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    sc.pl.scatter(
        data,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        ax=ax,
        show=False,
    )

    ax.set_title(f"Sample: {sample} - Basic QC")
    ax.set_xlabel("Total Counts")
    ax.set_ylabel("Number of Genes")

    for im in ax.get_images():
        if im.get_cmap():
            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label("% Mitochondrial")
            break

    plt.tight_layout()

    if save_dir:
        plt.savefig(
            os.path.join(save_dir, f"{sample}_qc_scatter.png"),
            dpi=300,
            bbox_inches="tight",
        )
    if show:
        plt.show()
    plt.close(fig)
