"""
Cell filtering and quality control logic for single-cell RNA-seq data.

This module provides marking, filtering, and reporting functions for low-quality
cells, doublets, and custom outliers, with flexible logic and clear outputs.
"""

import logging
import re
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from ..utils.utils import identify_outliers
from .config import FilterConfig, MarkingConfig, QCThresholds

log = logging.getLogger(__name__)

__all__ = [
    "suggest_qc_thresholds",
    "mark_low_quality_cell",
    "filter_cells",
    "generate_qc_report",
]


# --- Helper Functions ---
def _safe_threshold_check(
    data: pd.Series, threshold: Optional[float], operator: str, name: str
) -> pd.Series:
    """
    Safely apply threshold checks with None handling.

    Args:
        data: Data series to check
        threshold: Threshold value (can be None)
        operator: Comparison operator ('>', '<', '>=', '<=')
        name: Threshold name for logging

    Returns:
        Boolean series indicating threshold violations
    """
    if threshold is None:
        log.debug(f"Skipping {name} threshold check (threshold is None)")
        return pd.Series(False, index=data.index)

    if operator == ">":
        result = data > threshold
    elif operator == "<":
        result = data < threshold
    elif operator == ">=":
        result = data >= threshold
    elif operator == "<=":
        result = data <= threshold
    else:
        raise ValueError(f"Unsupported operator: {operator}")

    count = result.sum()
    percentage = count / len(data) * 100
    log.info(
        f"Cells failing {name} ({operator} {threshold}): {count} ({percentage:.2f}%)"
    )

    return result


def _plot_qc_outliers(
    adata: AnnData,
    sample_indices: Dict[str, pd.Series],
    cols_to_plot: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
    show: bool = True,
):
    """
    Generate QC outlier visualization plots.

    Args:
        adata: AnnData object with outlier annotations
        sample_indices: Dictionary mapping sample names to boolean masks
        cols_to_plot: List of columns to plot
        save_dir: Directory to save plots
        show: Whether to display plots
    """
    if cols_to_plot is None:
        # Define default columns to plot, checking for existence
        default_cols = [
            "outlier_min_genes",
            "outlier_mt",
            "outlier_hb",
            "outlier_qc_metrics",
        ]

        # Add doublet columns if available
        doublet_cols = [
            "scrublet_predicted",
            "heuristic_predicted",
            "predicted_doublet",
        ]
        for col in doublet_cols:
            if col in adata.obs.columns:
                default_cols.append(col)

        # Add custom outlier columns
        custom_cols = [
            col for col in adata.obs.columns if col.startswith("outlier_custom_")
        ]
        default_cols.extend(custom_cols)

        cols_to_plot = [col for col in default_cols if col in adata.obs.columns]

    if not cols_to_plot:
        log.warning("No valid columns found for plotting")
        return

    for sample, sample_mask in sample_indices.items():
        log.info(f"Plotting QC outliers for sample: {sample}")
        data_view = adata[sample_mask]

        if data_view.n_obs == 0:
            log.warning(f"No cells found for sample {sample}, skipping plot")
            continue

        # Calculate subplot layout
        n_plots = len(cols_to_plot)
        n_cols = min(3, n_plots)
        n_rows = (n_plots + n_cols - 1) // n_cols

        fig, axs = plt.subplots(
            n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows), facecolor="white"
        )
        if n_plots == 1:
            axs = [axs]
        elif n_rows == 1:
            axs = axs if isinstance(axs, np.ndarray) else [axs]
        else:
            axs = axs.flatten()

        fig.suptitle(f"QC Outlier Analysis for Sample: {sample}", fontsize=14, y=0.98)

        for i, col in enumerate(cols_to_plot):
            if i >= len(axs):
                break

            ax = axs[i]

            if col in data_view.obs.columns:
                # Determine coloring based on data type
                col_data = data_view.obs[col]

                if col_data.dtype == "bool" or set(col_data.unique()).issubset(
                    {0, 1, True, False}
                ):
                    # Boolean data - color by outlier status
                    colors = col_data.map(
                        {
                            False: "#637b8a",  # Blue for normal cells
                            True: "#d62728",  # Red for outliers
                            0: "#637b8a",
                            1: "#d62728",
                        }
                    )

                    # Count and percentage for title
                    outlier_count = col_data.sum()
                    outlier_pct = outlier_count / len(col_data) * 100
                    title = f"{col.replace('_', ' ').title()}\n{outlier_count} cells ({outlier_pct:.1f}%)"

                else:
                    # Continuous data - use value-based coloring
                    colors = col_data
                    title = col.replace("_", " ").title()

                # Create scatter plot
                scatter = ax.scatter(
                    data_view.obs["total_counts"],
                    data_view.obs["n_genes_by_counts"],
                    c=colors,
                    s=8,
                    alpha=0.7,
                    edgecolors="none",
                    rasterized=True,
                    cmap="viridis"
                    if not isinstance(colors, pd.Series) or colors.dtype != "object"
                    else None,
                )

                ax.set_title(title, fontsize=10)
                ax.set_xlabel("Total Counts", fontsize=9)
                ax.set_ylabel("Number of Genes", fontsize=9)
                ax.tick_params(labelsize=8)

                # Add colorbar for continuous data
                if not (
                    col_data.dtype == "bool"
                    or set(col_data.unique()).issubset({0, 1, True, False})
                ):
                    plt.colorbar(scatter, ax=ax, shrink=0.8)

            else:
                ax.set_visible(False)

        # Hide unused subplots
        for j in range(len(cols_to_plot), len(axs)):
            axs[j].set_visible(False)

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save_dir:
            Path(save_dir).mkdir(parents=True, exist_ok=True)
            filename = f"{sample}_qc_outliers.png"
            filepath = Path(save_dir) / filename
            plt.savefig(filepath, dpi=300, facecolor="white", bbox_inches="tight")
            log.info(f"Saved QC outlier plot to {filepath}")

        if show:
            plt.show()
        else:
            plt.close(fig)


def _plot_before_after_comparison(
    adata_before: AnnData,
    adata_after: AnnData,
    save_dir: str,
    sample_key: str,
    qc_metrics: List[str],
) -> None:
    """Plot before/after filtering comparison."""

    # Cell count comparison
    before_counts = adata_before.obs[sample_key].value_counts()
    after_counts = adata_after.obs[sample_key].value_counts()

    comparison_df = pd.DataFrame(
        {
            "before": before_counts,
            "after": after_counts.reindex(before_counts.index, fill_value=0),
        }
    ).fillna(0)

    comparison_df["removed"] = comparison_df["before"] - comparison_df["after"]
    comparison_df["retention_rate"] = (
        comparison_df["after"] / comparison_df["before"] * 100
    )

    # Plot cell count comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Bar plot of cell counts
    x = range(len(comparison_df))
    width = 0.35

    ax1.bar(
        [i - width / 2 for i in x],
        comparison_df["before"],
        width,
        label="Before",
        alpha=0.8,
    )
    ax1.bar(
        [i + width / 2 for i in x],
        comparison_df["after"],
        width,
        label="After",
        alpha=0.8,
    )

    ax1.set_xlabel("Sample")
    ax1.set_ylabel("Cell Count")
    ax1.set_title("Cell Counts Before/After Filtering")
    ax1.set_xticks(x)
    ax1.set_xticklabels(comparison_df.index, rotation=45)
    ax1.legend()

    # Retention rate plot
    ax2.bar(x, comparison_df["retention_rate"], alpha=0.8, color="green")
    ax2.set_xlabel("Sample")
    ax2.set_ylabel("Retention Rate (%)")
    ax2.set_title("Cell Retention Rate by Sample")
    ax2.set_xticks(x)
    ax2.set_xticklabels(comparison_df.index, rotation=45)
    ax2.axhline(y=80, color="red", linestyle="--", alpha=0.7, label="80% threshold")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(
        Path(save_dir) / "filtering_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Save comparison statistics
    comparison_df.to_csv(Path(save_dir) / "filtering_comparison_stats.csv")


# --- Main Functions ---
def suggest_qc_thresholds(
    adata: AnnData,
    method: Literal["mad", "iqr", "percentile"] = "mad",
    mad_multipliers: Union[float, List[float]] = [3.0, 4.0, 5.0],
    iqr_multiplier: float = 1.5,
    percentile_range: Tuple[float, float] = (2.5, 97.5),
    plot_distributions: bool = True,
    save_dir: Optional[str] = None,
) -> QCThresholds:
    """
    Automatically suggest QC thresholds based on data distribution and generate informative plots.

    This function analyzes the distribution of QC metrics and suggests reasonable
    thresholds. The generated plots now include the specific threshold values in the
    legend for clarity.

    Args:
        adata: AnnData object with calculated QC metrics.
        method: Method for threshold suggestion ("mad", "iqr", "percentile").
        mad_multipliers: A single multiplier or a list for MAD-based thresholds.
        iqr_multiplier: Multiplier for IQR-based thresholds.
        percentile_range: Percentile range for threshold suggestion.
        plot_distributions: Whether to plot distribution analysis.
        save_dir: Directory to save plots.

    Returns:
        QCThresholds object with suggested values (based on the first MAD multiplier if applicable).
    """
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required QC columns: {missing_cols}")

    log.info(f"Suggesting QC thresholds using '{method}' method...")

    if isinstance(mad_multipliers, (int, float)):
        mad_multipliers = [mad_multipliers]

    suggestions = {}

    # Define which metrics to analyze
    metrics = {
        "n_genes_by_counts": "Gene counts per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial percentage",
    }
    if "pct_counts_hb" in adata.obs.columns:
        metrics["pct_counts_hb"] = "Hemoglobin percentage"

    top_gene_cols = [
        col
        for col in adata.obs.columns
        if re.match(r"pct_counts_in_top_\d+_genes", col)
    ]
    for col in top_gene_cols:
        metrics[col] = (
            col.replace("_", " ")
            .replace("pct counts in ", "")
            .replace(" genes", "")
            .title()
        )

    if plot_distributions:
        n_metrics = len(metrics)
        n_cols = min(2, n_metrics)
        n_rows = (n_metrics + n_cols - 1) // n_cols
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows), constrained_layout=True
        )
        axes = np.array(axes).flatten()

    for i, (metric, title) in enumerate(metrics.items()):
        data = adata.obs[metric].dropna()
        ax = axes[i] if plot_distributions and i < len(axes) else None

        # --- Centralized threshold calculation logic ---
        # This part calculates bounds for all multipliers and stores them for plotting
        plot_lines = []
        is_count_metric = metric in ["n_genes_by_counts", "total_counts"]

        def _get_metric_key(metric):
            """Map metrics to QCThresholds parameter names"""
            if "top" in metric:
                match = re.search(r"pct_counts_in_top_(\d+)_genes", metric)
                if match:
                    # Use the format pc_top_20_genes instead of pc_top20_genes
                    return f"pc_top_{match.group(1)}_genes"
                else:
                    match = re.search(r"(\d+)", metric)
                    return f"pc_top_{match.group(1)}_genes" if match else "pc_top_genes"
            else:
                parts = metric.split("_")
                if len(parts) > 0:
                    return f"pc_{parts[-1]}"
                else:
                    return f"pc_{metric}"

        if method == "mad":
            median_val = data.median()
            mad_val = np.median(np.abs(data - median_val))
            if mad_val == 0:
                log.warning(
                    f"MAD for metric '{metric}' is zero. MAD-based thresholds may be unreliable."
                )
                mad_val = 1e-5  # Small value to avoid division by zero

            for idx, multiplier in enumerate(mad_multipliers):
                upper_bound = median_val + multiplier * mad_val
                if is_count_metric:
                    lower_bound = max(0, median_val - multiplier * mad_val)
                    plot_lines.append(
                        {
                            "val": lower_bound,
                            "label": f"Min (MAD x{multiplier})",
                            "color": "red",
                        }
                    )
                    plot_lines.append(
                        {
                            "val": upper_bound,
                            "label": f"Max (MAD x{multiplier})",
                            "color": "orange",
                        }
                    )
                    if idx == 0:  # Store suggestion from primary multiplier
                        suggestions[
                            f"min_{'genes' if 'genes' in metric else 'counts'}"
                        ] = int(lower_bound)
                        suggestions[
                            f"max_{'genes' if 'genes' in metric else 'counts'}"
                        ] = int(upper_bound)
                else:  # Percentage metric
                    plot_lines.append(
                        {
                            "val": upper_bound,
                            "label": f"Max (MAD x{multiplier})",
                            "color": "red",
                        }
                    )
                    if idx == 0:
                        metric_key = _get_metric_key(metric)
                        if metric_key:
                            suggestions[metric_key] = min(100.0, upper_bound)
                        else:
                            log.warning(
                                f"Could not generate suggestion key for metric '{metric}'"
                            )

        elif method == "iqr":
            q25, q75 = data.quantile([0.25, 0.75])
            iqr = q75 - q25
            upper_bound = q75 + iqr_multiplier * iqr
            if is_count_metric:
                lower_bound = max(0, q25 - iqr_multiplier * iqr)
                suggestions[f"min_{'genes' if 'genes' in metric else 'counts'}"] = int(
                    lower_bound
                )
                suggestions[f"max_{'genes' if 'genes' in metric else 'counts'}"] = int(
                    upper_bound
                )
                plot_lines.append(
                    {"val": lower_bound, "label": "Min (IQR)", "color": "red"}
                )
                plot_lines.append(
                    {"val": upper_bound, "label": "Max (IQR)", "color": "orange"}
                )
            else:
                metric_key = _get_metric_key(metric)
                plot_lines.append(
                    {"val": upper_bound, "label": "Max (IQR)", "color": "red"}
                )

        elif method == "percentile":
            upper_bound = data.quantile(percentile_range[1] / 100)
            if is_count_metric:
                lower_bound = data.quantile(percentile_range[0] / 100)
                suggestions[f"min_{'genes' if 'genes' in metric else 'counts'}"] = int(
                    lower_bound
                )
                suggestions[f"max_{'genes' if 'genes' in metric else 'counts'}"] = int(
                    upper_bound
                )
                plot_lines.append(
                    {
                        "val": lower_bound,
                        "label": f"Min ({percentile_range[0]}th %ile)",
                        "color": "red",
                    }
                )
                plot_lines.append(
                    {
                        "val": upper_bound,
                        "label": f"Max ({percentile_range[1]}th %ile)",
                        "color": "orange",
                    }
                )
            else:
                metric_key = _get_metric_key(metric)
                plot_lines.append(
                    {
                        "val": upper_bound,
                        "label": f"Max ({percentile_range[1]}th %ile)",
                        "color": "red",
                    }
                )

        # --- Plotting logic with dynamic labels ---
        if plot_distributions and ax is not None:
            ax.hist(data, bins=50, alpha=0.75, edgecolor="black")
            ax.set_title(title, fontsize=14, fontweight="bold")
            ax.set_xlabel(metric.replace("_", " ").title(), fontsize=12)
            ax.set_ylabel("Frequency", fontsize=12)

            for line in plot_lines:
                # Format label with the calculated value
                if is_count_metric:
                    formatted_label = f"{line['label']}: {line['val']:.0f}"
                else:  # Percentage
                    formatted_label = f"{line['label']}: {line['val']:.1f}%"

                ax.axvline(
                    x=line["val"],
                    color=line["color"],
                    linestyle="--",
                    alpha=0.8,
                    linewidth=1.5,
                    label=formatted_label,
                )

            # Create a clean legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))  # Removes duplicate labels
            ax.legend(by_label.values(), by_label.keys(), loc="upper right")
            ax.grid(axis="y", linestyle="--", alpha=0.6)

    if plot_distributions:
        for j in range(len(metrics), len(axes)):
            axes[j].set_visible(False)  # Hide unused subplots

        fig.suptitle(
            "Suggested QC Thresholds from Data Distribution",
            fontsize=18,
            fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for suptitle

        if save_dir:
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            plt.savefig(
                save_path / "qc_threshold_suggestions.png",
                dpi=300,
                bbox_inches="tight",
            )
        plt.show()

    # Create QCThresholds object with suggestions
    # Log which parameters we found and which we're passing to QCThresholds
    log.info("Parameters extracted from metrics:")
    pc_top_genes_dict = {}
    for k, v in suggestions.items():
        if k.startswith("pc_top_"):
            pc_top_genes_dict[k] = v

    # Filter out any parameters that the QCThresholds class doesn't accept
    # This requires knowing the valid parameters for QCThresholds
    # For now, just log what we're passing
    final_kwargs = {
        "min_genes": suggestions.get("min_genes"),
        "max_genes": suggestions.get("max_genes"),
        "min_counts": suggestions.get("min_counts"),
        "max_counts": suggestions.get("max_counts"),
        "pc_mt": suggestions.get("pc_mt"),
        "pc_hb": suggestions.get("pc_hb"),
        "pc_top_genes": pc_top_genes_dict,
    }

    final_kwargs = {k: v for k, v in final_kwargs.items() if v is not None}

    suggested_thresholds = QCThresholds(**final_kwargs)

    log.info("Suggested QC thresholds (based on primary multiplier/setting):")
    for param, value in suggested_thresholds.to_dict().items():
        if value is not None and value != {}:
            log.info(f"  {param}: {value}")

    return suggested_thresholds


def mark_low_quality_cell(
    adata: AnnData,
    sample_key: str = "sampleID",
    config: Optional[MarkingConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Identifies and marks low-quality cells using a configuration-driven workflow.

    Args:
        adata: AnnData object with QC metrics calculated.
        sample_key: Key in adata.obs for sample identification.
        config: A MarkingConfig object with all parameters.
        **kwargs: Additional parameters to override defaults in the config.

    Returns:
        AnnData object with boolean columns in .obs marking low-quality cells.
    """
    # === CONFIGURATION SETUP ===
    base_config = MarkingConfig()
    if config is not None:
        config_dict = asdict(config)
        if "thresholds" in config_dict:
            # Update the default thresholds object field by field
            for th_key, th_value in config_dict["thresholds"].items():
                if hasattr(base_config.thresholds, th_key):
                    setattr(base_config.thresholds, th_key, th_value)
            # Remove 'thresholds' so it's not processed again
            del config_dict["thresholds"]

        # Update the rest of the config fields
        for key, value in config_dict.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)

    if kwargs:
        # Override with any specific kwargs
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            # Allow nested threshold overrides like `min_genes=300`
            elif hasattr(base_config.thresholds, key):
                setattr(base_config.thresholds, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")

    cfg = base_config
    thresholds = cfg.thresholds

    # Check required QC columns
    required = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Missing required QC columns: {missing}. Run calculate_qc_metric() first."
        )

    log.info("Marking low-quality cells with the following thresholds:")
    for param, value in thresholds.to_dict().items():
        if value is not None:
            log.info(f"  {param}: {value}")

    # Precompute sample indices for efficiency
    sample_indices = {
        sample: adata.obs[sample_key] == sample
        for sample in adata.obs[sample_key].unique()
    }

    log.info(f"Processing {len(sample_indices)} samples with {adata.n_obs} total cells")

    # === THRESHOLD CHECKS ===

    # Gene count thresholds
    adata.obs["outlier_min_genes"] = _safe_threshold_check(
        adata.obs["n_genes_by_counts"], thresholds.min_genes, "<", "min_genes"
    )

    adata.obs["outlier_max_genes"] = _safe_threshold_check(
        adata.obs["n_genes_by_counts"], thresholds.max_genes, ">", "max_genes"
    )

    # Total count thresholds
    adata.obs["outlier_min_counts"] = _safe_threshold_check(
        adata.obs["total_counts"], thresholds.min_counts, "<", "min_counts"
    )

    adata.obs["outlier_max_counts"] = _safe_threshold_check(
        adata.obs["total_counts"], thresholds.max_counts, ">", "max_counts"
    )

    # Mitochondrial percentage threshold
    adata.obs["outlier_mt"] = _safe_threshold_check(
        adata.obs["pct_counts_mt"], thresholds.pc_mt, ">", "mitochondrial_percentage"
    )

    # Hemoglobin percentage threshold (if available)
    if "pct_counts_hb" in adata.obs.columns:
        adata.obs["outlier_hb"] = _safe_threshold_check(
            adata.obs["pct_counts_hb"], thresholds.pc_hb, ">", "hemoglobin_percentage"
        )
    else:
        adata.obs["outlier_hb"] = False
        log.info("Hemoglobin percentage not available, setting outlier_hb to False")

    # === MAD-BASED OUTLIER DETECTION ===

    # Format metrics for identify_outliers function
    formatted_metrics = [
        (metric, direction, None) for metric, direction in cfg.qc_metrics_mad
    ]

    # Run MAD-based outlier detection
    adata.obs["outlier_qc_metrics"] = identify_outliers(
        adata, metrics=formatted_metrics, sample_key=sample_key, nmads=thresholds.nmads
    )

    # Handle fixed top gene thresholds if specified
    if thresholds.use_fixed_top_gene_threshold:
        for metric_key, threshold_value in thresholds.pc_top_genes.items():
            # Construct the column name from the key, e.g., pc_top_20_genes -> pct_counts_in_top_20_genes
            col_name = f"pct_counts_in_{metric_key.split('pc_')[-1]}"
            outlier_col_name = (
                f"outlier_{metric_key.split('pc_')[-1]}"  # e.g., outlier_top_20_genes
            )

            if col_name in adata.obs.columns:
                adata.obs[outlier_col_name] = _safe_threshold_check(
                    adata.obs[col_name],
                    threshold_value,
                    ">",
                    f"fixed {metric_key}",
                )
                # Combine with other QC metrics for a unified outlier flag
                adata.obs["outlier_qc_metrics"] = (
                    adata.obs["outlier_qc_metrics"] | adata.obs[outlier_col_name]
                )
            else:
                log.warning(
                    f"Fixed threshold provided for '{metric_key}', but column '{col_name}' not found in data."
                )

    qc_count = adata.obs["outlier_qc_metrics"].sum()
    log.info(
        f"Cells marked as QC metric outliers: {qc_count} ({qc_count / adata.n_obs:.2%})"
    )

    # === CUSTOM OUTLIER DETECTION ===

    if cfg.custom_outlier_functions:
        log.info(
            f"Running {len(cfg.custom_outlier_functions)} custom outlier detection functions..."
        )
        for func_name, func in cfg.custom_outlier_functions.items():
            try:
                custom_outliers = func(adata)
                if not isinstance(custom_outliers, pd.Series):
                    raise ValueError(
                        f"Custom function {func_name} must return a pandas Series"
                    )
                if len(custom_outliers) != adata.n_obs:
                    raise ValueError(
                        f"Custom function {func_name} returned wrong length"
                    )

                col_name = f"outlier_custom_{func_name}"
                adata.obs[col_name] = custom_outliers.astype(bool)

                custom_count = custom_outliers.sum()
                log.info(
                    f"Custom outliers ({func_name}): {custom_count} ({custom_count / adata.n_obs:.2%})"
                )

            except Exception as e:
                log.error(f"Error in custom outlier function {func_name}: {e}")

    # === Store parameters in the unified namespace ===
    if "sclucid" not in adata.uns:
        adata.uns["sclucid"] = {}
    if "qc" not in adata.uns["sclucid"]:
        adata.uns["sclucid"]["qc"] = {}

    adata.uns["sclucid"]["qc"]["marking_params"] = {
        "thresholds": thresholds.to_dict(),
    }

    # === SUMMARY STATISTICS ===

    log.info("\n" + "=" * 50)
    log.info("LOW-QUALITY CELL DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Count cells with different types of issues
    outlier_cols = [col for col in adata.obs.columns if col.startswith("outlier_")]
    adata.obs["outlier_count"] = adata.obs[outlier_cols].sum(axis=1)

    # Report counts per outlier type
    for col in outlier_cols:
        count = adata.obs[col].sum()
        percentage = count / total_cells * 100
        log.info(
            f"{col.replace('outlier_', '').replace('_', ' ').title()}: {count} cells ({percentage:.2f}%)"
        )

    # Report cells with multiple issues
    max_issues = adata.obs["outlier_count"].max()
    for n_outliers in range(1, int(max_issues) + 1):
        count = (adata.obs["outlier_count"] == n_outliers).sum()
        if count > 0:
            percentage = count / total_cells * 100
            log.info(
                f"Cells with exactly {n_outliers} types of issues: {count} ({percentage:.2f}%)"
            )

    # Include doublet statistics if available
    doublet_cols = ["predicted_doublet", "scrublet_predicted", "heuristic_predicted"]
    for doublet_col in doublet_cols:
        if doublet_col in adata.obs.columns:
            count = adata.obs[doublet_col].sum()
            percentage = count / total_cells * 100
            log.info(
                f"{doublet_col.replace('_', ' ').title()}: {count} cells ({percentage:.2f}%)"
            )

    log.info("=" * 50)

    # === VISUALIZATION ===

    if cfg.plot_outliers:
        _plot_qc_outliers(
            adata, sample_indices, cfg.cols_to_plot, cfg.save_dir, cfg.show
        )

    return adata


def filter_cells(
    adata: AnnData,
    config: Optional[FilterConfig] = None,
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Enhanced cell filtering with flexible logical combinations and detailed reporting.

    This function filters cells based on previously calculated QC and doublet boolean flags.
    It supports multiple combination strategies and provides comprehensive statistics.

    Args:
        adata: AnnData object with QC metrics calculated
        config: FilterConfig object with filtering parameters
        copy: Whether to return a new filtered AnnData object

    Returns:
        Filtered AnnData object if copy=True, otherwise filters in place and returns None

    Example:
        # Basic usage
        adata_filtered = filter_cells(adata, copy=True)

        # Custom logic - remove cells with both MT and QC issues
        config = FilterConfig(
            combination_logic="custom",
            custom_logic_expr="outlier_mt & outlier_qc_metrics"
        )
        filter_cells(adata, config=config)

        # Threshold-based - remove cells with at least 2 issues
        config = FilterConfig(
            combination_logic="threshold",
            min_criteria_for_removal=2
        )
        filter_cells(adata, config=config)
    """
    # === 1. CONFIGURATION SETUP ===
    base_config = FilterConfig()
    if config is not None:
        base_config.__dict__.update(asdict(config))
    if kwargs:
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")
    cfg = base_config
    cfg.validate()

    # --- Use cfg.criteria_to_filter instead of building a new list ---
    criteria = cfg.criteria_to_filter

    # === 2. FILTERING ===
    # Build criteria list if not provided
    if criteria is None:
        criteria = []

        # Map config attributes to column names
        criteria_mapping = {
            "filter_by_outlier_min_genes": "outlier_min_genes",
            "filter_by_outlier_max_genes": "outlier_max_genes",
            "filter_by_outlier_min_counts": "outlier_min_counts",
            "filter_by_outlier_max_counts": "outlier_max_counts",
            "filter_by_outlier_mt": "outlier_mt",
            "filter_by_outlier_hb": "outlier_hb",
            "filter_by_outlier_qc_metrics": "outlier_qc_metrics",
            "filter_by_scrublet_predicted": "scrublet_predicted",
            "filter_by_heuristic_predicted": "heuristic_predicted",
            "filter_by_predicted_doublet": "predicted_doublet",
        }

        for config_attr, col_name in criteria_mapping.items():
            if getattr(cfg, config_attr, False) and col_name in adata.obs.columns:
                criteria.append(col_name)

    # Filter out criteria that don't exist
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
    log.info(f"Starting cell filtering with {initial_cells} cells")
    log.info(f"Using criteria: {', '.join(valid_criteria)}")
    log.info(f"Combination logic: {cfg.combination_logic}")

    # Apply metadata filters first if specified
    metadata_mask = pd.Series(True, index=adata.obs_names)
    if cfg.metadata_filters:
        for key, value in cfg.metadata_filters.items():
            if key in adata.obs.columns:
                if isinstance(value, list):
                    metadata_mask &= adata.obs[key].isin(value)
                else:
                    metadata_mask &= adata.obs[key] == value
                log.info(f"Applied metadata filter {key}={value}")
            else:
                log.warning(f"Metadata key '{key}' not found in adata.obs")

    # Calculate individual criteria masks
    criteria_masks = {}
    criteria_counts = {}

    for col in valid_criteria:
        col_mask = adata.obs[col].fillna(False).astype(bool)
        criteria_masks[col] = col_mask
        criteria_counts[col] = col_mask.sum()

    # Apply combination logic
    if cfg.combination_logic == "any":
        # Remove if ANY criterion is true (default behavior)
        combined_removal_mask = pd.Series(False, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask |= col_mask

    elif cfg.combination_logic == "all":
        # Remove only if ALL criteria are true
        combined_removal_mask = pd.Series(True, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask &= col_mask

    elif cfg.combination_logic == "custom":
        # Use custom expression
        try:
            # Create a namespace with all criteria for evaluation
            namespace = {col: criteria_masks[col] for col in valid_criteria}
            combined_removal_mask = eval(
                cfg.custom_logic_expr, {"__builtins__": {}}, namespace
            )

            if not isinstance(combined_removal_mask, pd.Series):
                combined_removal_mask = pd.Series(
                    combined_removal_mask, index=adata.obs_names
                )

        except Exception as e:
            log.error(f"Error evaluating custom logic expression: {e}")
            raise ValueError(
                f"Invalid custom logic expression: {cfg.custom_logic_expr}"
            )

    elif cfg.combination_logic == "threshold":
        # Remove if at least min_criteria_for_removal criteria are true
        criteria_sum = pd.Series(0, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            criteria_sum += col_mask.astype(int)
        combined_removal_mask = criteria_sum >= cfg.min_criteria_for_removal

    else:
        raise ValueError(f"Unknown combination logic: {cfg.combination_logic}")

    # Apply metadata filter
    combined_removal_mask = combined_removal_mask & metadata_mask

    # Calculate final keep mask
    keep_mask = ~combined_removal_mask

    # Report statistics
    log.info("\n" + "=" * 40)
    log.info("CELL FILTERING STATISTICS")
    log.info("=" * 40)

    # Individual criteria counts
    total_cells = len(adata.obs_names)
    for col, count in criteria_counts.items():
        percentage = count / total_cells * 100
        log.info(f"{col}: {count} cells ({percentage:.2f}%)")

    # Final filtering results
    removed_count = combined_removal_mask.sum()
    kept_count = keep_mask.sum()

    log.info("\nFiltering results:")
    log.info(f"  Initial cells: {initial_cells}")
    log.info(f"  Cells removed: {removed_count} ({removed_count / initial_cells:.2%})")
    log.info(f"  Cells retained: {kept_count} ({kept_count / initial_cells:.2%})")

    # Analyze overlap between criteria
    if len(valid_criteria) > 1:
        log.info("\nOverlap analysis:")

        # Count cells with multiple issues
        criteria_sum = pd.Series(0, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            criteria_sum += col_mask.astype(int)

        for i in range(1, len(valid_criteria) + 1):
            count = (criteria_sum == i).sum()
            if count > 0:
                percentage = count / total_cells * 100
                log.info(
                    f"  Cells with exactly {i} issues: {count} ({percentage:.2f}%)"
                )

    log.info("=" * 40)

    # === Store filtering results in the unified namespace ===
    if "sclucid" not in adata.uns:
        adata.uns["sclucid"] = {}
    if "qc" not in adata.uns["sclucid"]:
        adata.uns["sclucid"]["qc"] = {}

    # Store stats in a dictionary
    stats = {
        "initial_cells": initial_cells,
        "final_cells": kept_count,
        "removed_cells": removed_count,
        "removed_fraction": removed_count / initial_cells if initial_cells > 0 else 0,
        "criteria_used": valid_criteria,
        "combination_logic": config.combination_logic,
        "criteria_counts": criteria_counts,
        "config": config.to_dict() if config else {"criteria": criteria},
    }

    # Decide where to save the results. If this is a copy, save to the new object.
    target_adata = adata
    if copy:
        # If we return a copy, the stats should be in the new object's .uns
        # The filtering logic already returns a copy, so we just add the .uns dict to it.
        adata_filtered = adata[keep_mask, :].copy()
        target_adata = adata_filtered

    if "sclucid" not in target_adata.uns:
        target_adata.uns["sclucid"] = {}
    if "qc" not in target_adata.uns["sclucid"]:
        target_adata.uns["sclucid"]["qc"] = {}

    target_adata.uns["sclucid"]["qc"]["filtering_results"] = stats

    # Perform filtering
    if copy:
        # Create the filtered copy ONCE.
        adata_filtered = adata[keep_mask, :].copy()
        # Ensure the .uns structure exists and add the results.
        adata_filtered.uns.setdefault("sclucid", {}).setdefault("qc", {})[
            "filtering_results"
        ] = stats
        return adata_filtered
    else:
        # Filter in place.
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
            "filtering_results"
        ] = stats
        adata._inplace_subset_obs(keep_mask)
        return None


def generate_qc_report(
    adata: AnnData,
    save_dir: str,
    sample_key: str = "sampleID",
    include_before_after: bool = True,
    adata_before: Optional[AnnData] = None,
) -> None:
    """
    Generate comprehensive QC report with visualizations.

    This function creates a detailed report of quality control metrics
    and filtering results, including before/after comparisons and
    statistical summaries.

    Args:
        adata: AnnData object (after filtering)
        save_dir: Directory to save report files
        sample_key: Key for sample identification
        include_before_after: Whether to include before/after comparison
        adata_before: AnnData object before filtering (for comparison)
    """
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    log.info(f"Generating QC report in {save_dir}")

    if include_before_after and adata_before is None:
        log.warning(
            "`adata_before` not provided, cannot generate before/after comparison plots."
        )
        include_before_after = False

    # QC metrics to analyze
    qc_metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    if "sclucid" in adata.uns and "metrics" in adata.uns["sclucid"]["qc"]:
        # This makes the report automatically adapt to what was calculated
        params = adata.uns["sclucid"]["qc"]["metrics"]["params"]
        if params.get("extra_gene_sets_provided"):
            # Logic to find the pct_counts_* columns from the params
            pass  # You can add logic here to make it even smarter

    # 1. Summary statistics table
    summary_stats = []

    samples = adata.obs[sample_key].unique()
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs[sample_mask]

        stats = {"sample": sample, "n_cells": len(sample_data)}

        for metric in qc_metrics:
            if metric in sample_data.columns:
                stats[f"{metric}_mean"] = sample_data[metric].mean()
                stats[f"{metric}_median"] = sample_data[metric].median()
                stats[f"{metric}_std"] = sample_data[metric].std()

        summary_stats.append(stats)

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(Path(save_dir) / "qc_summary_statistics.csv", index=False)

    # 2. QC distributions plot
    n_metrics = len(qc_metrics)
    n_samples = len(samples)

    fig, axes = plt.subplots(n_metrics, 1, figsize=(12, 4 * n_metrics))
    if n_metrics == 1:
        axes = [axes]

    for i, metric in enumerate(qc_metrics):
        ax = axes[i]

        # Box plot by sample
        sample_data = []
        sample_labels = []

        for sample in samples:
            sample_mask = adata.obs[sample_key] == sample
            sample_values = adata.obs.loc[sample_mask, metric].dropna()
            sample_data.append(sample_values)
            sample_labels.append(f"{sample}\n(n={len(sample_values)})")

        ax.boxplot(sample_data, labels=sample_labels)
        ax.set_title(f"{metric.replace('_', ' ').title()} Distribution by Sample")
        ax.set_ylabel(metric.replace("_", " ").title())

        if len(samples) > 5:
            ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(Path(save_dir) / "qc_distributions.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 3. Before/after comparison if requested
    if include_before_after and adata_before is not None:
        _plot_before_after_comparison(
            adata_before, adata, save_dir, sample_key, qc_metrics
        )

    # 4. Outlier summary
    outlier_cols = [col for col in adata.obs.columns if col.startswith("outlier_")]
    if outlier_cols:
        outlier_summary = []

        for sample in samples:
            sample_mask = adata.obs[sample_key] == sample
            sample_data = adata.obs[sample_mask]

            stats = {"sample": sample, "n_cells": len(sample_data)}

            for col in outlier_cols:
                if col in sample_data.columns:
                    count = sample_data[col].sum()
                    stats[f"{col}_count"] = count
                    stats[f"{col}_percentage"] = count / len(sample_data) * 100

            outlier_summary.append(stats)

        outlier_df = pd.DataFrame(outlier_summary)
        outlier_df.to_csv(Path(save_dir) / "outlier_summary.csv", index=False)

    log.info("QC report generation completed")
