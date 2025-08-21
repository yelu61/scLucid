"""
Enhanced cell filtering utilities for single-cell RNA-seq data.

This module provides comprehensive functions for identifying and filtering low-quality
cells based on various quality metrics. It supports flexible filtering logic,
automatic parameter suggestion, custom outlier detection functions, and detailed
reporting with visualizations.
"""

import logging
import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from ..utils.utils import identify_outliers

log = logging.getLogger(__name__)

__all__ = [
    "FilterConfig",
    "QCThresholds",
    "mark_low_quality_cell",
    "filter_cells",
    "suggest_qc_thresholds",
    "generate_qc_report",
]


@dataclass
class QCThresholds:
    """
    Configuration class for QC thresholds with automatic validation.

    This class encapsulates all quality control thresholds and provides
    validation to ensure parameters are reasonable and consistent.
    """

    # Gene count thresholds
    min_genes: Optional[int] = 200
    max_genes: Optional[int] = None

    # Total count thresholds
    min_counts: Optional[int] = None
    max_counts: Optional[int] = None

    # Percentage thresholds
    pc_mt: Optional[float] = 20.0
    pc_hb: Optional[float] = 20.0
    pc_top20_genes: Optional[float] = None

    # MAD-based outlier detection
    nmads: float = 5.0

    # Custom threshold validation
    use_fixed_top20_threshold: bool = False

    def __post_init__(self):
        """Validate threshold parameters after initialization."""
        self.validate()

    def validate(self):
        """
        Validate threshold parameters for consistency and reasonableness.

        Raises:
            ValueError: If any threshold values are invalid or inconsistent
        """
        # Validate gene count thresholds
        if self.min_genes is not None and self.min_genes < 0:
            raise ValueError("min_genes must be non-negative")

        if self.max_genes is not None and self.max_genes < 0:
            raise ValueError("max_genes must be non-negative")

        if (
            self.min_genes is not None
            and self.max_genes is not None
            and self.min_genes >= self.max_genes
        ):
            raise ValueError("min_genes must be less than max_genes")

        # Validate count thresholds
        if self.min_counts is not None and self.min_counts < 0:
            raise ValueError("min_counts must be non-negative")

        if self.max_counts is not None and self.max_counts < 0:
            raise ValueError("max_counts must be non-negative")

        if (
            self.min_counts is not None
            and self.max_counts is not None
            and self.min_counts >= self.max_counts
        ):
            raise ValueError("min_counts must be less than max_counts")

        # Validate percentage thresholds
        for param_name, value in [
            ("pc_mt", self.pc_mt),
            ("pc_hb", self.pc_hb),
            ("pc_top20_genes", self.pc_top20_genes),
        ]:
            if value is not None and not (0 <= value <= 100):
                raise ValueError(f"{param_name} must be between 0 and 100")

        # Validate MAD parameter
        if self.nmads <= 0:
            raise ValueError("nmads must be positive")

    def to_dict(self) -> Dict[str, Any]:
        """Convert thresholds to dictionary format."""
        return {
            "min_genes": self.min_genes,
            "max_genes": self.max_genes,
            "min_counts": self.min_counts,
            "max_counts": self.max_counts,
            "pc_mt": self.pc_mt,
            "pc_hb": self.pc_hb,
            "pc_top20_genes": self.pc_top20_genes,
            "nmads": self.nmads,
            "use_fixed_top20_threshold": self.use_fixed_top20_threshold,
        }


@dataclass
class FilterConfig:
    """
    Configuration class for cell filtering with flexible logic combinations.

    This class defines how different QC criteria should be combined and
    which filtering operations to perform.
    """

    # Basic filtering criteria
    filter_by_outlier_min_genes: bool = True
    filter_by_outlier_max_genes: bool = True
    filter_by_outlier_min_counts: bool = True
    filter_by_outlier_max_counts: bool = True
    filter_by_outlier_mt: bool = True
    filter_by_outlier_hb: bool = True
    filter_by_outlier_qc_metrics: bool = True
    filter_by_scrublet_predicted: bool = True
    filter_by_heuristic_predicted: bool = True
    filter_by_predicted_doublet: bool = True

    # Logical combination strategy
    combination_logic: Literal["any", "all", "custom"] = "any"
    custom_logic_expr: Optional[str] = (
        None  # e.g., "(outlier_mt & outlier_qc_metrics) | predicted_doublet"
    )

    # Metadata-based filtering
    metadata_filters: Optional[Dict[str, Any]] = None

    # Minimum criteria for removal (when combination_logic="threshold")
    min_criteria_for_removal: int = 1

    def validate(self):
        """Validate filter configuration."""
        if self.combination_logic == "custom" and not self.custom_logic_expr:
            raise ValueError(
                "custom_logic_expr must be provided when combination_logic='custom'"
            )

        if self.min_criteria_for_removal < 1:
            raise ValueError("min_criteria_for_removal must be at least 1")


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
            os.makedirs(save_dir, exist_ok=True)
            filename = f"{sample}_qc_outliers.png"
            filepath = os.path.join(save_dir, filename)
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
        os.path.join(save_dir, "filtering_comparison.png"), dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Save comparison statistics
    comparison_df.to_csv(os.path.join(save_dir, "filtering_comparison_stats.csv"))


# --- Main Functions ---
def suggest_qc_thresholds(
    adata: AnnData,
    sample_key: str = "sampleID",
    method: Literal["mad", "iqr", "percentile"] = "mad",
    mad_multiplier: float = 3.0,
    iqr_multiplier: float = 1.5,
    percentile_range: Tuple[float, float] = (2.5, 97.5),
    plot_distributions: bool = True,
    save_dir: Optional[str] = None,
) -> QCThresholds:
    """
    Automatically suggest QC thresholds based on data distribution.

    This function analyzes the distribution of QC metrics and suggests
    reasonable thresholds based on statistical outlier detection methods.

    Args:
        adata: AnnData object with calculated QC metrics
        sample_key: Key for sample identification
        method: Method for threshold suggestion ("mad", "iqr", "percentile")
        mad_multiplier: Multiplier for MAD-based thresholds
        iqr_multiplier: Multiplier for IQR-based thresholds
        percentile_range: Percentile range for threshold suggestion
        plot_distributions: Whether to plot distribution analysis
        save_dir: Directory to save plots

    Returns:
        QCThresholds object with suggested values
    """
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required QC columns: {missing_cols}")

    log.info(f"Suggesting QC thresholds using {method} method...")

    # Initialize suggestions dictionary
    suggestions = {}

    # Analyze each metric
    metrics = {
        "n_genes_by_counts": "Gene counts per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial percentage",
    }

    # Add hemoglobin if available
    if "pct_counts_hb" in adata.obs.columns:
        metrics["pct_counts_hb"] = "Hemoglobin percentage"

    # Add top20 genes if available
    if "pct_counts_in_top_20_genes" in adata.obs.columns:
        metrics["pct_counts_in_top_20_genes"] = "Top 20 genes percentage"

    if plot_distributions:
        n_metrics = len(metrics)
        fig, axes = plt.subplots(2, (n_metrics + 1) // 2, figsize=(15, 8))
        if n_metrics == 1:
            axes = [axes]
        elif n_metrics <= 2:
            axes = axes.flatten()
        else:
            axes = axes.flatten()

    for i, (metric, title) in enumerate(metrics.items()):
        data = adata.obs[metric].dropna()

        if method == "mad":
            median_val = data.median()
            mad_val = np.median(np.abs(data - median_val))

            if metric in ["n_genes_by_counts", "total_counts"]:
                # For count data, suggest lower bound only
                lower_bound = max(0, median_val - mad_multiplier * mad_val)
                upper_bound = median_val + mad_multiplier * mad_val
                suggestions[
                    f"min_{metric.split('_')[0] if 'genes' in metric else 'counts'}"
                ] = int(lower_bound)
                if metric == "n_genes_by_counts":
                    suggestions["max_genes"] = int(upper_bound)
                else:
                    suggestions["max_counts"] = int(upper_bound)
            else:
                # For percentage data, suggest upper bound only
                upper_bound = median_val + mad_multiplier * mad_val
                suggestions[f"pc_{metric.split('_')[-1]}"] = min(100, upper_bound)

        elif method == "iqr":
            q25, q75 = data.quantile([0.25, 0.75])
            iqr = q75 - q25

            if metric in ["n_genes_by_counts", "total_counts"]:
                lower_bound = max(0, q25 - iqr_multiplier * iqr)
                upper_bound = q75 + iqr_multiplier * iqr
                suggestions[
                    f"min_{metric.split('_')[0] if 'genes' in metric else 'counts'}"
                ] = int(lower_bound)
                if metric == "n_genes_by_counts":
                    suggestions["max_genes"] = int(upper_bound)
                else:
                    suggestions["max_counts"] = int(upper_bound)
            else:
                upper_bound = q75 + iqr_multiplier * iqr
                suggestions[f"pc_{metric.split('_')[-1]}"] = min(100, upper_bound)

        elif method == "percentile":
            if metric in ["n_genes_by_counts", "total_counts"]:
                lower_bound = data.quantile(percentile_range[0] / 100)
                upper_bound = data.quantile(percentile_range[1] / 100)
                suggestions[
                    f"min_{metric.split('_')[0] if 'genes' in metric else 'counts'}"
                ] = int(lower_bound)
                if metric == "n_genes_by_counts":
                    suggestions["max_genes"] = int(upper_bound)
                else:
                    suggestions["max_counts"] = int(upper_bound)
            else:
                upper_bound = data.quantile(percentile_range[1] / 100)
                suggestions[f"pc_{metric.split('_')[-1]}"] = min(100, upper_bound)

        # Plot distribution with suggested thresholds
        if plot_distributions and i < len(axes):
            ax = axes[i]
            ax.hist(data, bins=50, alpha=0.7, edgecolor="black")
            ax.set_title(title)
            ax.set_xlabel(metric.replace("_", " ").title())
            ax.set_ylabel("Frequency")

            # Add threshold lines
            if metric in ["n_genes_by_counts", "total_counts"]:
                param_name = "min_genes" if "genes" in metric else "min_counts"
                if param_name in suggestions:
                    ax.axvline(
                        suggestions[param_name],
                        color="red",
                        linestyle="--",
                        label=f"Min: {suggestions[param_name]}",
                    )

                param_name = "max_genes" if "genes" in metric else "max_counts"
                if param_name in suggestions:
                    ax.axvline(
                        suggestions[param_name],
                        color="orange",
                        linestyle="--",
                        label=f"Max: {suggestions[param_name]}",
                    )
            else:
                param_name = f"pc_{metric.split('_')[-1]}"
                if param_name in suggestions:
                    ax.axvline(
                        suggestions[param_name],
                        color="red",
                        linestyle="--",
                        label=f"Max: {suggestions[param_name]:.1f}%",
                    )

            ax.legend()

    if plot_distributions:
        # Hide unused subplots
        for j in range(len(metrics), len(axes)):
            axes[j].set_visible(False)

        plt.tight_layout()

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(
                os.path.join(save_dir, "qc_threshold_suggestions.png"),
                dpi=300,
                bbox_inches="tight",
            )
        plt.show()

    # Create QCThresholds object with suggestions
    threshold_kwargs = {
        "min_genes": suggestions.get("min_genes"),
        "max_genes": suggestions.get("max_genes"),
        "min_counts": suggestions.get("min_counts"),
        "max_counts": suggestions.get("max_counts"),
        "pc_mt": suggestions.get("pc_mt"),
        "pc_hb": suggestions.get("pc_hb"),
        "pc_top20_genes": suggestions.get("pc_top20_genes"),
    }

    # Remove None values
    threshold_kwargs = {k: v for k, v in threshold_kwargs.items() if v is not None}

    suggested_thresholds = QCThresholds(**threshold_kwargs)

    log.info("Suggested QC thresholds:")
    for param, value in suggested_thresholds.to_dict().items():
        if value is not None:
            log.info(f"  {param}: {value}")

    return suggested_thresholds


def mark_low_quality_cell(
    adata: AnnData,
    sample_key: str = "sampleID",
    thresholds: Optional[QCThresholds] = None,
    qc_metrics: Optional[List[Tuple[str, str]]] = None,
    custom_outlier_functions: Optional[
        Dict[str, Callable[[AnnData], pd.Series]]
    ] = None,
    plot_outliers: bool = False,
    save_dir: Optional[str] = None,
    show: bool = True,
    cols_to_plot: Optional[List[str]] = None,
    # Backward compatibility parameters
    min_genes: Optional[int] = None,
    max_genes: Optional[int] = None,
    min_counts: Optional[int] = None,
    max_counts: Optional[int] = None,
    nmads: Optional[float] = None,
    pc_mt: Optional[float] = None,
    pc_hb: Optional[float] = None,
    pc_top20_genes: Optional[float] = None,
    use_fixed_top20_threshold: Optional[bool] = None,
) -> AnnData:
    """
    Enhanced function to identify and mark low-quality cells with robust parameter handling.

    This function identifies low-quality cells based on various QC criteria including
    gene counts, total counts, mitochondrial percentages, and custom metrics. It supports
    both fixed thresholds and MAD-based outlier detection, with flexible parameter handling.

    Args:
        adata: AnnData object with QC metrics calculated
        sample_key: Key in adata.obs for sample identification
        thresholds: QCThresholds object with all threshold parameters
        qc_metrics: List of (metric, direction) tuples for MAD-based outlier detection
        custom_outlier_functions: Dictionary of custom outlier detection functions
        plot_outliers: Whether to generate scatter plots highlighting outliers
        save_dir: Directory to save plots
        show: Whether to display plots
        cols_to_plot: List of .obs columns to plot as scatter plots

        # Backward compatibility parameters (will be used if thresholds is None)
        min_genes: Minimum number of genes to keep a cell
        max_genes: Maximum number of genes to keep a cell
        min_counts: Minimum total counts to keep a cell
        max_counts: Maximum total counts to keep a cell
        nmads: Number of MADs for outlier detection
        pc_mt: Maximum percentage of mitochondrial counts
        pc_hb: Maximum percentage of hemoglobin counts
        pc_top20_genes: Fixed threshold for pct_counts_in_top_20_genes
        use_fixed_top20_threshold: Whether to use fixed threshold for top20 genes

    Returns:
        AnnData object with boolean columns in .obs marking low-quality cells

    Raises:
        ValueError: If required QC metrics are missing
    """
    # Check required QC columns
    required = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing = [col for col in required if col not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Missing required QC columns: {missing}. Run calculate_qc_metric() first."
        )

    # Handle thresholds - use provided object or create from individual parameters
    if thresholds is None:
        # Use individual parameters for backward compatibility
        thresholds = QCThresholds(
            min_genes=min_genes if min_genes is not None else 200,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            pc_mt=pc_mt if pc_mt is not None else 20.0,
            pc_hb=pc_hb if pc_hb is not None else 20.0,
            pc_top20_genes=pc_top20_genes,
            nmads=nmads if nmads is not None else 5.0,
            use_fixed_top20_threshold=use_fixed_top20_threshold
            if use_fixed_top20_threshold is not None
            else False,
        )

    log.info("Marking low-quality cells with the following thresholds:")
    for param, value in thresholds.to_dict().items():
        if value is not None:
            log.info(f"  {param}: {value}")

    # Set up default QC metrics for MAD-based detection
    if qc_metrics is None:
        qc_metrics = [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
        ]

        if (
            not thresholds.use_fixed_top20_threshold
            and "pct_counts_in_top_20_genes" in adata.obs.columns
        ):
            qc_metrics.append(("pct_counts_in_top_20_genes", "upper"))
            log.info("Using MAD-based threshold for pct_counts_in_top_20_genes")
        elif (
            thresholds.use_fixed_top20_threshold
            and thresholds.pc_top20_genes is not None
        ):
            log.info(
                f"Using fixed threshold {thresholds.pc_top20_genes}% for pct_counts_in_top_20_genes"
            )

    # Precompute sample indices for efficiency
    sample_indices = {
        sample: adata.obs[sample_key] == sample
        for sample in adata.obs[sample_key].unique()
    }

    log.info(f"Processing {len(sample_indices)} samples with {adata.n_obs} total cells")

    # === FIXED THRESHOLD CHECKS ===

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
    formatted_metrics = [(metric, direction, None) for metric, direction in qc_metrics]

    # Run MAD-based outlier detection
    adata.obs["outlier_qc_metrics"] = identify_outliers(
        adata, metrics=formatted_metrics, sample_key=sample_key, nmads=thresholds.nmads
    )

    # Handle fixed top20 genes threshold if specified
    if thresholds.use_fixed_top20_threshold and thresholds.pc_top20_genes is not None:
        if "pct_counts_in_top_20_genes" in adata.obs.columns:
            adata.obs["outlier_top20_genes"] = _safe_threshold_check(
                adata.obs["pct_counts_in_top_20_genes"],
                thresholds.pc_top20_genes,
                ">",
                "top20_genes_percentage",
            )
            # Combine with other QC metrics
            adata.obs["outlier_qc_metrics"] = (
                adata.obs["outlier_qc_metrics"] | adata.obs["outlier_top20_genes"]
            )
        else:
            log.warning(
                "pct_counts_in_top_20_genes not found, skipping top20 genes threshold"
            )

    qc_count = adata.obs["outlier_qc_metrics"].sum()
    log.info(
        f"Cells marked as QC metric outliers: {qc_count} ({qc_count / adata.n_obs:.2%})"
    )

    # === CUSTOM OUTLIER DETECTION ===

    if custom_outlier_functions:
        log.info(
            f"Running {len(custom_outlier_functions)} custom outlier detection functions..."
        )
        for func_name, func in custom_outlier_functions.items():
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
    if "scrnatk" not in adata.uns:
        adata.uns["scrnatk"] = {}
    if "qc" not in adata.uns["scrnatk"]:
        adata.uns["scrnatk"]["qc"] = {}

    adata.uns["scrnatk"]["qc"]["marking_params"] = {
        "thresholds": thresholds.to_dict(),
        "qc_metrics_for_mad": qc_metrics,
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

    if plot_outliers:
        _plot_qc_outliers(adata, sample_indices, cols_to_plot, save_dir, show)

    return adata


def filter_cells(
    adata: AnnData,
    config: Optional[FilterConfig] = None,
    criteria: Optional[List[str]] = None,
    copy: bool = False,
    # Backward compatibility parameters
    filter_by_outlier_min_genes: Optional[bool] = None,
    filter_by_outlier_mt: Optional[bool] = None,
    filter_by_outlier_hb: Optional[bool] = None,
    filter_by_outlier_qc_metrics: Optional[bool] = None,
    filter_by_scrublet_predicted: Optional[bool] = None,
    filter_by_heuristic_predicted: Optional[bool] = None,
    filter_by_predicted_doublet: Optional[bool] = None,
) -> Optional[AnnData]:
    """
    Enhanced cell filtering with flexible logical combinations and detailed reporting.

    This function filters cells based on previously calculated QC and doublet boolean flags.
    It supports multiple combination strategies and provides comprehensive statistics.

    Args:
        adata: AnnData object with QC metrics calculated
        config: FilterConfig object with filtering parameters
        criteria: Custom list of boolean columns for filtering (overrides config)
        copy: Whether to return a new filtered AnnData object

        # Backward compatibility parameters
        filter_by_*: Individual boolean flags for specific criteria

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
    # Handle configuration - use provided config or create from individual parameters
    if config is None:
        config = FilterConfig(
            filter_by_outlier_min_genes=filter_by_outlier_min_genes
            if filter_by_outlier_min_genes is not None
            else True,
            filter_by_outlier_mt=filter_by_outlier_mt
            if filter_by_outlier_mt is not None
            else True,
            filter_by_outlier_hb=filter_by_outlier_hb
            if filter_by_outlier_hb is not None
            else True,
            filter_by_outlier_qc_metrics=filter_by_outlier_qc_metrics
            if filter_by_outlier_qc_metrics is not None
            else True,
            filter_by_scrublet_predicted=filter_by_scrublet_predicted
            if filter_by_scrublet_predicted is not None
            else True,
            filter_by_heuristic_predicted=filter_by_heuristic_predicted
            if filter_by_heuristic_predicted is not None
            else True,
            filter_by_predicted_doublet=filter_by_predicted_doublet
            if filter_by_predicted_doublet is not None
            else True,
        )

    config.validate()

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
            if getattr(config, config_attr, False) and col_name in adata.obs.columns:
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
    log.info(f"Combination logic: {config.combination_logic}")

    # Apply metadata filters first if specified
    metadata_mask = pd.Series(True, index=adata.obs_names)
    if config.metadata_filters:
        for key, value in config.metadata_filters.items():
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
    if config.combination_logic == "any":
        # Remove if ANY criterion is true (default behavior)
        combined_removal_mask = pd.Series(False, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask |= col_mask

    elif config.combination_logic == "all":
        # Remove only if ALL criteria are true
        combined_removal_mask = pd.Series(True, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask &= col_mask

    elif config.combination_logic == "custom":
        # Use custom expression
        try:
            # Create a namespace with all criteria for evaluation
            namespace = {col: criteria_masks[col] for col in valid_criteria}
            combined_removal_mask = eval(
                config.custom_logic_expr, {"__builtins__": {}}, namespace
            )

            if not isinstance(combined_removal_mask, pd.Series):
                combined_removal_mask = pd.Series(
                    combined_removal_mask, index=adata.obs_names
                )

        except Exception as e:
            log.error(f"Error evaluating custom logic expression: {e}")
            raise ValueError(
                f"Invalid custom logic expression: {config.custom_logic_expr}"
            )

    elif config.combination_logic == "threshold":
        # Remove if at least min_criteria_for_removal criteria are true
        criteria_sum = pd.Series(0, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            criteria_sum += col_mask.astype(int)
        combined_removal_mask = criteria_sum >= config.min_criteria_for_removal

    else:
        raise ValueError(f"Unknown combination logic: {config.combination_logic}")

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
    if "scrnatk" not in adata.uns:
        adata.uns["scrnatk"] = {}
    if "qc" not in adata.uns["scrnatk"]:
        adata.uns["scrnatk"]["qc"] = {}

    # Store stats in a dictionary
    stats = {
        "initial_cells": initial_cells,
        "final_cells": kept_count,
        "removed_cells": removed_count,
        "removed_fraction": removed_count / initial_cells if initial_cells > 0 else 0,
        "criteria_used": valid_criteria,
        "combination_logic": config.combination_logic,
        "criteria_counts": criteria_counts,
    }

    # Decide where to save the results. If this is a copy, save to the new object.
    target_adata = adata
    if copy:
        # If we return a copy, the stats should be in the new object's .uns
        # The filtering logic already returns a copy, so we just add the .uns dict to it.
        adata_filtered = adata[keep_mask, :].copy()
        target_adata = adata_filtered

    if "scrnatk" not in target_adata.uns:
        target_adata.uns["scrnatk"] = {}
    if "qc" not in target_adata.uns["scrnatk"]:
        target_adata.uns["scrnatk"]["qc"] = {}

    target_adata.uns["scrnatk"]["qc"]["filtering_results"] = stats

    # Perform filtering
    if copy:
        adata_filtered = adata[keep_mask, :].copy()
        return adata_filtered
    else:
        # Filter in place
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
    os.makedirs(save_dir, exist_ok=True)
    log.info(f"Generating QC report in {save_dir}")

    if include_before_after and adata_before is None:
        log.warning(
            "`adata_before` not provided, cannot generate before/after comparison plots."
        )
        include_before_after = False

    # QC metrics to analyze
    qc_metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    if "scrnatk" in adata.uns and "metrics" in adata.uns["scrnatk"]["qc"]:
        # This makes the report automatically adapt to what was calculated
        params = adata.uns["scrnatk"]["qc"]["metrics"]["params"]
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
    summary_df.to_csv(os.path.join(save_dir, "qc_summary_statistics.csv"), index=False)

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
    plt.savefig(
        os.path.join(save_dir, "qc_distributions.png"), dpi=300, bbox_inches="tight"
    )
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
        outlier_df.to_csv(os.path.join(save_dir, "outlier_summary.csv"), index=False)

    log.info("QC report generation completed")
