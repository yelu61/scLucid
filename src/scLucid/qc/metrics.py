"""
Quality control metrics calculation for single-cell RNA-seq data.

This module provides functions for calculating and visualizing quality control
metrics for single-cell datasets, with robust parameter validation,
auto-detection of sample keys, and user-friendly logging.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from .config import MetricsReportingConfig

log = logging.getLogger(__name__)

__all__ = ["calculate_qc_metric"]


# --- Helper Functions ---
def _find_sample_key(adata: AnnData, sample_key: Optional[str] = None) -> str:
    """
    Automatically detect the sample key column in adata.obs.

    Returns the first matching column from a list of common names, or raises an error.
    """
    if sample_key and sample_key in adata.obs.columns:
        return sample_key
    for key in ["sampleID", "sample", "Sample", "batch", "Batch"]:
        if key in adata.obs.columns:
            log.info(f"Auto-detected sample/group column: {key}")
            return key
    raise ValueError(
        "No valid sample key found in adata.obs. Please provide a sample_key argument."
    )


def _export_qc_stats(
    adata: AnnData,
    sample_key: str,
    percent_top_cols: List[str],
    outdir: Optional[str] = None,
    export_csv: bool = True,
    export_xlsx: bool = False,
    print_summary: bool = True,
) -> Dict[str, pd.DataFrame]:
    """
    Export per-sample and global QC statistics as CSV/XLSX tables.

    Args:
        adata: AnnData object with calculated QC metrics.
        sample_key: obs key used for grouping samples.
        percent_top_cols: obs keys for percent of top genes.
        outdir: Directory for saving summary tables.
        export_csv: Whether to save CSV files.
        export_xlsx: Whether to save XLSX file.
        print_summary: Whether to print summaries to the log.

    Returns:
        Dictionary with per-sample and global summary DataFrames.
    """
    # Define standard metrics to look for
    standard_metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ]

    # Filter metrics to only include columns that actually exist in adata.obs
    metrics = [metric for metric in standard_metrics if metric in adata.obs.columns]
    if percent_top_cols:
        metrics.extend([col for col in percent_top_cols if col in adata.obs.columns])

    if not metrics:
        log.warning("No QC metrics found to export.")
        return {}

    # Calculate per-sample and global statistics
    sample_stats = adata.obs.groupby(sample_key)[metrics].describe()
    global_stats = adata.obs[metrics].describe().T

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        if export_csv:
            sample_stats.to_csv(Path(outdir) / "qc_stats_per_sample.csv")
            global_stats.to_csv(Path(outdir) / "qc_stats_global.csv")
        if export_xlsx:
            with pd.ExcelWriter(Path(outdir) / "qc_stats.xlsx") as writer:
                sample_stats.to_excel(writer, sheet_name="per_sample")
                global_stats.to_excel(writer, sheet_name="global")

    if print_summary:
        log.info("==== Global QC Metrics Summary ====")
        log.info(global_stats)
        log.info("==== Per Sample QC Metrics (head) ====")
        log.info(sample_stats.head(10))

    return {"sample": sample_stats, "global": global_stats}


def _detect_qc_outliers(
    adata: AnnData,
    sample_key: str,
    percent_top_cols: List[str],
    outdir: Optional[str] = None,
) -> List[str]:
    """
    Detect and log potential outlier samples/metrics.

    Args:
        adata: AnnData object with QC metrics.
        sample_key: obs key for grouping.
        percent_top_cols: obs keys for percent of top genes.
        outdir: Directory to save tips/warnings.

    Returns:
        List of textual warnings or suggestions.
    """
    tips = []
    warn = log.warning

    # Mitochondrial percentage warning
    if "pct_counts_mt" in adata.obs.columns:
        mt_stats = adata.obs.groupby(sample_key)["pct_counts_mt"].mean()
        for s, v in mt_stats.items():
            if v > 15:
                warn(
                    f"[QC Alert] Sample {s} has a high mean mitochondrial content: {v:.1f}%. Recommended is typically ≤ 15%."
                )
                tips.append(f"{s} mt_mean={v:.1f}%")

    # Loop through all top gene columns for outlier detection
    for col in percent_top_cols:
        if col in adata.obs.columns:
            m = re.search(r"pct_counts_in_top_(\d+)_genes", col)
            top_n = int(m.group(1)) if m else "N"
            top_gene_stats = adata.obs.groupby(sample_key)[col].mean()
            for s, v in top_gene_stats.items():
                if v > 60:
                    warn(
                        f"[QC Alert] Sample {s} has a high mean Top-{top_n} gene fraction: {v:.1f}%. "
                        "This could indicate low-complexity libraries or technical artifacts."
                    )
                    tips.append(f"{s} top{top_n}_mean={v:.1f}%")

    # Detected gene number warning
    if "n_genes_by_counts" in adata.obs.columns:
        gene_stats = adata.obs.groupby(sample_key)["n_genes_by_counts"].median()
        for s, v in gene_stats.items():
            if v < 1000:
                warn(
                    f"[QC Alert] Sample {s} has a low median of genes detected: {int(v)}. Possible low quality."
                )
                tips.append(f"{s} median_n_genes={int(v)}")

    if outdir and tips:
        with open(Path(outdir) / "qc_outlier_tips.txt", "w") as f:
            f.write("\n".join(tips))
    return tips


def _sample_for_plotting(
    adata: AnnData,
    max_cells: int = 50000,
    sample_key: str = "sampleID",
    random_state: int = 61,
) -> AnnData:
    """
    Sample cells for plotting when dataset is too large.

    Args:
        adata: Input AnnData object
        max_cells: Maximum number of cells to sample for plotting
        sample_key: Key to stratify sampling by samples
        random_state: Random state for reproducibility

    Returns:
        Sampled AnnData object
    """
    if adata.n_obs <= max_cells:
        return adata

    log.info(
        f"Large dataset detected ({adata.n_obs} cells), sampling {max_cells} cells for visualization..."
    )

    # Stratified sampling by sample to maintain representation
    if sample_key in adata.obs.columns:
        sampled_indices = []
        samples = adata.obs[sample_key].unique()
        cells_per_sample = max_cells // len(samples)

        for sample in samples:
            sample_mask = adata.obs[sample_key] == sample
            sample_indices = np.where(sample_mask)[0]

            if len(sample_indices) <= cells_per_sample:
                sampled_indices.extend(sample_indices)
            else:
                np.random.seed(random_state)
                selected = np.random.choice(
                    sample_indices, size=cells_per_sample, replace=False
                )
                sampled_indices.extend(selected)

        # If we still have room, randomly sample more cells
        remaining_slots = max_cells - len(sampled_indices)
        if remaining_slots > 0:
            all_indices = set(range(adata.n_obs))
            available_indices = list(all_indices - set(sampled_indices))
            if available_indices:
                np.random.seed(random_state)
                additional = np.random.choice(
                    available_indices,
                    size=min(remaining_slots, len(available_indices)),
                    replace=False,
                )
                sampled_indices.extend(additional)

    else:
        # Simple random sampling if no sample key
        np.random.seed(random_state)
        sampled_indices = np.random.choice(adata.n_obs, size=max_cells, replace=False)

    return adata[sampled_indices].copy()


def _plot_top_genes_distribution(
    adata: AnnData,
    sample_key: str = "sampleID",
    percent_top_col: str = "pct_counts_in_top_20_genes",
    thresholds: List[float] = [50, 60, 65, 70],
    save_dir: Optional[str] = None,
    show: bool = True,
    max_cells_for_plotting: int = 50000,
    random_state: int = 61,
) -> Optional[Tuple[plt.Figure, plt.Figure, plt.Figure]]:
    """
    Plot detailed distribution of pct_counts_in_top_X_genes to help determine appropriate thresholds.
    Optimized for large datasets with intelligent sampling.

    Args:
        adata: AnnData object with QC metrics calculated.
        sample_key: The key in adata.obs to identify different samples.
        percent_top_col: Column name in adata.obs containing the percentage of counts in top X genes.
        thresholds: List of thresholds to mark in the plots.
        save_dir: Directory to save plots. If None, plots are not saved.
        show: Whether to display the plots.
        max_cells_for_plotting: Maximum number of cells to use for plotting (optimization for large datasets).
        random_state: Random state for reproducible sampling.
        percent_top_col: Column name in adata.obs containing the percentage of counts in top X genes.

    Returns:
        Tuple of Figure objects for the three plots (histogram, boxplot, scatter) if show=True
    """
    if percent_top_col not in adata.obs.columns:
        log.warning(
            f"{percent_top_col} not found in adata.obs. Skipping its distribution plot."
        )
        return

    m = re.search(r"pct_counts_in_top_(\d+)_genes", percent_top_col)
    main_top_n = int(m.group(1)) if m else "N"

    log.info(
        f"Plotting detailed distribution for 'pct_counts_in_top_{main_top_n}_genes'..."
    )

    adata_plot = _sample_for_plotting(
        adata,
        max_cells=max_cells_for_plotting,
        sample_key=sample_key,
        random_state=random_state,
    )

    if adata_plot.n_obs < adata.n_obs:
        log.info(
            f"Using {adata_plot.n_obs} sampled cells for visualization, "
            f"but statistics computed on all {adata.n_obs} cells"
        )

    # Calculate summary statistics on full dataset
    if sample_key in adata.obs.columns:
        stats = adata.obs.groupby(sample_key)[percent_top_col].describe()
        log.info(f"Summary statistics for {percent_top_col} by sample:\n{stats}")

        # For each threshold, calculate percentage of cells above it
        for threshold in thresholds:
            counts = adata.obs.groupby(sample_key)[percent_top_col].apply(
                lambda x: (x > threshold).sum()
            )
            percentages = counts / adata.obs.groupby(sample_key).size() * 100
            threshold_stats = pd.DataFrame(
                {"counts": counts, "percentage": percentages}
            )
            log.info(f"Cells above {threshold}% threshold:\n{threshold_stats}")
    else:
        stats = adata.obs[percent_top_col].describe()
        log.info(f"Summary statistics for {percent_top_col}:\n{stats}")

    # 1. Histogram with density curves by sample
    fig1, ax1 = plt.subplots(figsize=(12, 8), facecolor="white")

    if sample_key in adata_plot.obs.columns:
        samples = adata_plot.obs[sample_key].unique()
        for sample in samples:
            sample_data = adata_plot[adata_plot.obs[sample_key] == sample].obs
            if len(sample_data) > 0:  # Only plot if sample has data
                sns.kdeplot(
                    data=sample_data,
                    x=percent_top_col,
                    label=f"Sample {sample}",
                    ax=ax1,
                )
    else:
        sns.kdeplot(
            data=adata_plot.obs,
            x=percent_top_col,
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

    # Add percentile markers (90th, 95th, 99th) - use full dataset for percentiles
    percentiles = [90, 95, 99]
    for p in percentiles:
        val = np.percentile(adata.obs[percent_top_col], p)
        ax1.axvline(
            x=val,
            linestyle=":",
            color="green",
            alpha=0.7,
            label=f"{p}th percentile: {val:.1f}%",
        )

    ax1.set_title(
        "Distribution of Percentage Counts in Top Genes by Sample", fontsize=14
    )
    ax1.set_xlabel(f"Percentage of counts in top {main_top_n} genes (%)", fontsize=12)
    ax1.set_ylabel("Density", fontsize=12)
    ax1.legend(title="Sample / Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            Path(save_dir) / f"pct_top_{main_top_n}_genes_distribution.png",
            dpi=300,
            bbox_inches="tight",
        )

    # 2. Boxplot by sample
    fig2, ax2 = plt.subplots(figsize=(12, 8), facecolor="white")

    if sample_key in adata_plot.obs.columns:
        sns.boxplot(data=adata_plot.obs, x=sample_key, y=percent_top_col, ax=ax2)
    else:
        sns.boxplot(data=adata_plot.obs, y=percent_top_col, ax=ax2)

    # Add horizontal lines for thresholds
    for threshold in thresholds:
        ax2.axhline(
            y=threshold,
            linestyle="--",
            color="red",
            alpha=0.7,
            label=f"{threshold}% threshold",
        )

    ax2.set_title("Boxplot of Percentage Counts in Top Genes by Sample", fontsize=14)
    if sample_key in adata_plot.obs.columns:
        ax2.set_xlabel("Sample", fontsize=12)
        plt.xticks(rotation=45)
    ax2.set_ylabel("Percentage of counts in top genes (%)", fontsize=12)
    ax2.legend(title="Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            Path(save_dir) / f"pct_top_{main_top_n}_genes_boxplot.png",
            dpi=300,
            bbox_inches="tight",
        )

    # 3. Scatter plot: pct_counts_in_top_X_genes vs n_genes_by_counts
    fig3, ax3 = plt.subplots(figsize=(10, 8), facecolor="white")

    if sample_key in adata_plot.obs.columns:
        sc.pl.scatter(
            adata_plot,
            x="n_genes_by_counts",
            y=percent_top_col,
            color=sample_key,
            ax=ax3,
            show=False,
        )
    else:
        ax3.scatter(
            adata_plot.obs["n_genes_by_counts"],
            adata_plot.obs[percent_top_col],
            alpha=0.6,
            s=1,
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
        "Relationship between Gene Counts and Top Genes Percentage", fontsize=14
    )
    ax3.set_xlabel("Number of genes detected", fontsize=12)
    ax3.set_ylabel("Percentage of counts in top genes (%)", fontsize=12)
    ax3.legend(title="Threshold", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            Path(save_dir) / f"pct_top_{main_top_n}_genes_scatter.png",
            dpi=300,
            bbox_inches="tight",
        )

    log.info(f"Analysis of {percent_top_col} complete.")

    if show:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)

    # Make sure to close figures if not showing
    if not show:
        plt.close("all")

    return (fig1, fig2, fig3) if show else None


def _plot_qc_violin(
    adata_view: AnnData,
    keys: List[str],
    sample: str,
    save_dir: Optional[str] = None,
    show: bool = False,
    max_cells_for_plotting: int = 10000,
) -> None:
    """
    Internal helper to plot QC violin plots for a sample.
    Optimized for large datasets.
    """
    # Sample data if too large
    if adata_view.n_obs > max_cells_for_plotting:
        log.info(
            f"Sampling {max_cells_for_plotting} cells from {adata_view.n_obs} for violin plot"
        )
        indices = np.random.choice(
            adata_view.n_obs, max_cells_for_plotting, replace=False
        )
        adata_plot = adata_view[indices]
    else:
        adata_plot = adata_view

    fig, axs = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4), facecolor="white")
    axs = [axs] if len(keys) == 1 else axs

    for ax, key in zip(axs, keys):
        if key in adata_plot.obs.columns:
            sc.pl.violin(adata_plot, key, ax=ax, show=False)
            ax.set_title(key.replace("_", " ").title())
            ax.set_ylabel(key)
            plt.setp(ax.get_xticklabels(), visible=True)
        else:
            ax.text(
                0.5,
                0.5,
                f"'{key}' not found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(f"{key} (not available)")

    fig.suptitle(f"QC Metrics for Sample: {sample}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_dir:
        plt.savefig(
            Path(save_dir) / f"{sample}_qc_violin.png",
            dpi=300,
            bbox_inches="tight",
        )
    if show:
        plt.show()
    plt.close(fig)


def _plot_qc_scatter(
    adata_view: AnnData,
    sample: str,
    save_dir: Optional[str] = None,
    show: bool = False,
    max_cells_for_plotting: int = 10000,
) -> None:
    """
    Plot total_counts vs n_genes_by_counts scatter for a sample.
    """
    # Sample data if too large
    if adata_view.n_obs > max_cells_for_plotting:
        log.info(
            f"Sampling {max_cells_for_plotting} cells from {adata_view.n_obs} for scatter plot"
        )
        indices = np.random.choice(
            adata_view.n_obs, max_cells_for_plotting, replace=False
        )
        adata_plot = adata_view[indices]
    else:
        adata_plot = adata_view

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")

    if "pct_counts_mt" in adata_plot.obs.columns:
        sc.pl.scatter(
            adata_plot,
            x="total_counts",
            y="n_genes_by_counts",
            color="pct_counts_mt",
            ax=ax,
            show=False,
        )

        # Add colorbar
        for im in ax.get_images():
            if im.get_cmap():
                cbar = fig.colorbar(im, ax=ax)
                cbar.set_label("% Mitochondrial")
                break
    else:
        ax.scatter(
            adata_plot.obs["total_counts"],
            adata_plot.obs["n_genes_by_counts"],
            alpha=0.6,
            s=1,
        )

    ax.set_title(f"Sample: {sample} - Basic QC")
    ax.set_xlabel("Total Counts")
    ax.set_ylabel("Number of Genes")
    plt.tight_layout()

    if save_dir:
        plt.savefig(
            Path(save_dir) / f"{sample}_qc_scatter.png",
            dpi=300,
            bbox_inches="tight",
        )
    if show:
        plt.show()
    plt.close(fig)


def _get_default_gene_patterns() -> Dict[str, str]:
    """Get default gene patterns for standard gene sets."""
    return {
        "mt": r"^(MT|Mt|mt)-",
        "ribo": r"^(RP[SL]|Rp[sl])",
        "hb": r"^(HB|hb)[^(P|p)]",
    }


def _validate_gene_patterns(
    gene_patterns: Dict[str, str], var_names: pd.Index
) -> Dict[str, str]:
    """
    Validate gene patterns and warn if no genes match.

    Args:
        gene_patterns: Dictionary of gene patterns
        var_names: Gene names from AnnData

    Returns:
        Validated gene patterns dictionary
    """
    validated_patterns = {}

    for gene_type, pattern in gene_patterns.items():
        try:
            matches = var_names.str.contains(pattern, regex=True, na=False)
            n_matches = matches.sum()

            if n_matches == 0:
                log.warning(
                    f"No genes found matching pattern '{pattern}' for gene type '{gene_type}'"
                )
            else:
                log.info(
                    f"Found {n_matches} genes matching pattern '{pattern}' for gene type '{gene_type}'"
                )
                validated_patterns[gene_type] = pattern

        except Exception as e:
            log.error(
                f"Invalid regex pattern '{pattern}' for gene type '{gene_type}': {e}"
            )

    return validated_patterns


# --- Main Functions ---
def calculate_qc_metric(
    adata: AnnData,
    sample_key: Optional[str] = None,
    extra_gene_sets: Optional[Dict[str, Union[str, List[str]]]] = None,
    max_cells_for_plotting: int = 50000,
    random_state: int = 61,
    percent_top: Optional[Union[int, List[int]]] = None,
    calculate_cell_cycle: bool = False,
    cell_cycle_species: str = "human",
    cell_cycle_kwargs: Optional[dict] = None,
    reporting_config: Optional[MetricsReportingConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Calculate and plot QC metrics for each sample in the AnnData object.

    Args:
        adata: AnnData object containing single-cell data.
        sample_key: The key in adata.obs to identify different samples.
        extra_gene_sets: Dictionary of additional gene sets to calculate metrics for.
            Can contain either regex patterns (str) or gene lists (List[str]).
            Example: {
                'stress': r'^(HSP|HSPA|HSPB)',
                'cell_cycle': ['MKI67', 'TOP2A', 'PCNA'],
                'custom_set': r'^CUSTOM'
            }
        show: Whether to display the plots.
        max_cells_for_plotting: Maximum number of cells to use for plotting (large dataset optimization).
        random_state: Random state for reproducible sampling.
        percent_top: List of top gene percentages to calculate. Default is [20].
        calculate_cell_cycle: Also compute cell cycle scores.
        cell_cycle_species: Species for cell cycle scoring.
        cell_cycle_kwargs: Extra kwargs for cell cycle scoring function.

    Returns:
        AnnData object with QC metrics added to .obs and .var.
    """
    report_config = MetricsReportingConfig()

    if reporting_config is not None:
        config_dict = reporting_config.to_dict()  # Pydantic's built-in serialization
        for key, value in config_dict.items():
            if hasattr(report_config, key):
                setattr(report_config, key, value)

    if kwargs:
        for key, value in kwargs.items():
            if hasattr(report_config, key):
                setattr(report_config, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")

    cfg = report_config

    # --- Basic sanity checks ---
    sample_key = _find_sample_key(adata, sample_key)

    if adata.n_obs == 0:
        raise ValueError("Input AnnData object is empty")

    if adata.X is None:
        raise ValueError("Input AnnData object has no expression matrix")

    log.info(
        f"Calculating QC metrics for {adata.n_obs} cells and {adata.n_vars} genes..."
    )

    # --- Prepare gene sets for QC ---
    gene_patterns: Dict[str, str] = {}

    # Add standard gene patterns if requested
    if cfg.include_standard_qc:
        gene_patterns.update(_get_default_gene_patterns())
        log.info(
            "Including standard QC gene sets: mitochondrial, ribosomal, hemoglobin"
        )

    # Add extra gene sets
    if extra_gene_sets:
        for gene_set_name, gene_set_def in extra_gene_sets.items():
            if isinstance(gene_set_def, str):
                # It's a regex pattern
                gene_patterns[gene_set_name] = gene_set_def
            elif isinstance(gene_set_def, list):
                # It's a gene list - create boolean mask directly
                gene_mask = adata.var_names.isin(gene_set_def)
                n_found = gene_mask.sum()
                if n_found == 0:
                    log.warning(
                        f"No genes found for gene set '{gene_set_name}' from provided list"
                    )
                else:
                    log.info(
                        f"Found {n_found}/{len(gene_set_def)} genes for gene set '{gene_set_name}'"
                    )
                    adata.var[gene_set_name] = gene_mask
            else:
                log.warning(
                    f"Invalid gene set definition for '{gene_set_name}'. Must be str (regex) or List[str] (gene list)"
                )

    # Validate and apply regex patterns
    if gene_patterns:
        validated_patterns = _validate_gene_patterns(gene_patterns, adata.var_names)

        for gene_type, pattern in validated_patterns.items():
            adata.var[gene_type] = adata.var_names.str.contains(
                pattern, regex=True, na=False
            )

    # --- Dynamically determine the primary top-gene column name ---
    if percent_top is None:
        percent_top_list = [20]
    elif isinstance(percent_top, int):
        percent_top_list = [percent_top]
    else:
        percent_top_list = sorted(list(set(percent_top)))  # Ensure unique and sorted

    percent_top_cols = [f"pct_counts_in_top_{n}_genes" for n in percent_top_list]
    log.info(
        f"Will calculate and analyze top gene percentages for: Top {percent_top_list}"
    )

    # --- Calculate QC metrics with scanpy ---
    qc_vars = [
        col
        for col in adata.var.columns
        if col
        in list(gene_patterns.keys())
        + [
            name
            for name in (extra_gene_sets or {}).keys()
            if isinstance(extra_gene_sets[name], list)
        ]
    ]

    if qc_vars:
        log.info(f"Calculating QC metrics for gene sets: {qc_vars}")
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            inplace=True,
            percent_top=percent_top_list,
            log1p=True,
        )
    else:
        log.warning("No valid gene sets found. Calculating basic QC metrics only.")
        sc.pp.calculate_qc_metrics(
            adata,
            inplace=True,
            percent_top=percent_top_list,
            log1p=True,
        )

    log.info("QC metrics calculation complete.")

    # --- Centralized .uns storage ---
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("metrics", {})
    adata.uns["sclucid"]["qc"]["metrics"]["params"] = {
        "sample_key": sample_key,
        "include_standard_qc": cfg.include_standard_qc,
        "extra_gene_sets_provided": bool(extra_gene_sets),
        "percent_top_calculated": percent_top_list,
        "scanpy_version": getattr(sc, "__version__", "unknown"),
    }

    # --- Optional: Cell cycle scoring ---
    if calculate_cell_cycle:
        if cell_cycle_kwargs is None:
            cell_cycle_kwargs = {}
        try:
            from .cycle import score_cell_cycle
        except ImportError:
            raise ImportError("cycle.py not found or not importable in the current module.")
        log.info("Automatic cell cycle scoring enabled.")
        adata = score_cell_cycle(
            adata,
            species=cell_cycle_species,
            plot=cfg.show_plots,
            save_dir=cfg.save_dir,
            **cell_cycle_kwargs,
        )
        log.info("Cell cycle scores and phase added to .obs.")

    # --- Export parameters and QC statistics ---
    if cfg.export_stats and cfg.save_dir:
        _export_qc_stats(
            adata,
            sample_key=sample_key,
            percent_top_cols=percent_top_cols,
            outdir=cfg.save_dir,
            export_csv=True,
            export_xlsx=cfg.export_xlsx,
            print_summary=cfg.print_stats,
        )

    # --- Detect and log outlier warnings ---
    _detect_qc_outliers(
        adata,
        sample_key=sample_key,
        percent_top_cols=percent_top_cols,
        outdir=cfg.save_dir if cfg.save_dir else None,
    )

    # --- Plotting ---
    if cfg.save_dir:
        Path(cfg.save_dir).mkdir(parents=True, exist_ok=True)

    if cfg.plot_top_genes:
        for col_name in percent_top_cols:
            _plot_top_genes_distribution(
                adata,
                sample_key=sample_key,
                percent_top_col=col_name,
                save_dir=cfg.save_dir,
                show=cfg.show_plots,
                max_cells_for_plotting=max_cells_for_plotting,
                random_state=random_state,
            )

    if cfg.plot_violin or cfg.plot_scatter:
        # Determine which metrics to plot
        keys_to_plot = ["total_counts", "n_genes_by_counts"]
        # Gather all potential pct_counts_* columns that exist in adata.obs
        potential_metrics = [
            c for c in adata.obs.columns if c.startswith("pct_counts_")
        ]
        keys_to_plot.extend(
            [metric for metric in potential_metrics if metric in adata.obs.columns]
        )
        # We don't want to plot all 50 top gene %s, so we exclude them here.
        keys_to_plot = [k for k in keys_to_plot if "in_top" not in k]

        # Limit to reasonable number for plotting
        if len(keys_to_plot) > 6:
            log.info(
                f"Too many metrics ({len(keys_to_plot)}), selecting first 6 for plotting"
            )
            keys_to_plot = keys_to_plot[:6]

        for sample in adata.obs[sample_key].unique():
            log.info(f"Plotting QC for sample: {sample}")
            adata_view = adata[adata.obs[sample_key] == sample]

            if adata_view.n_obs == 0:
                log.warning(f"No cells found for sample {sample}, skipping plots")
                continue

            if cfg.plot_violin:
                _plot_qc_violin(
                    adata_view,
                    keys_to_plot,
                    sample,
                    save_dir=cfg.save_dir,
                    show=cfg.show_plots,
                    max_cells_for_plotting=max_cells_for_plotting
                    // max(1, len(adata.obs[sample_key].unique())),
                )
            if cfg.plot_scatter:
                _plot_qc_scatter(
                    adata_view,
                    sample,
                    save_dir=cfg.save_dir,
                    show=cfg.show_plots,
                    max_cells_for_plotting=max_cells_for_plotting
                    // max(1, len(adata.obs[sample_key].unique())),
                )

    return adata
