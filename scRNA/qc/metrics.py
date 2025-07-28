"""
Quality control metrics calculation for single-cell RNA-seq data.

This module provides functions for calculating standard QC metrics and
identifying outlier cells based on statistical methods.
"""

import logging
import os
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
log = logging.getLogger(__name__)

__all__ = ["calculate_qc_metric"]


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
        save_dir: Directory to save plots. If None, plots are not saved.
        show: Whether to display the plots.

    Returns:
        AnnData object with QC metrics added to .obs and .var.
    """
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not in adata.obs columns.")
    if not adata.obs[sample_key].nunique():
        raise ValueError(f"No samples found under sample key '{sample_key}'.")

    # --- Main Calculation (Optimized) ---
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
    if plot_violin or plot_scatter:
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

        samples = adata.obs[sample_key].unique()
        for sample in samples:
            log.info(f"Plotting QC for sample: {sample}")
            data_view = adata[adata.obs[sample_key] == sample]

            if plot_violin:
                fig, axs = plt.subplots(
                    1, len(keys), figsize=(15, 4), facecolor="white"
                )
                if len(keys) == 1:  # Ensure axs is always iterable
                    axs = [axs]

                for ax, key in zip(axs, keys):
                    sc.pl.violin(data_view, key, ax=ax, show=False)
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

            if plot_scatter:
                fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
                sc.pl.scatter(
                    data_view,
                    x="total_counts",
                    y="n_genes_by_counts",
                    color="pct_counts_mt",
                    ax=ax,
                    show=False,
                )

                ax.set_title(f"Sample: {sample} - Basic QC")
                ax.set_xlabel("Total Counts")
                ax.set_ylabel("Number of Genes")

                # The colorbar is now part of the figure, find it and label it
                # This is a bit more robust than creating a new one
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

    return adata
