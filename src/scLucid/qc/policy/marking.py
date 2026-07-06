"""QC evidence marking helpers used by the threshold decision chain."""

import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.plotting.plotting_utils import _is_interactive_backend

from .adaptive_threshold import AdaptiveThresholdCalculator, compute_mad_bounds
from ..artifacts import record_mark_evidence
from ..config import MarkingConfig

log = logging.getLogger(__name__)

__all__: list[str] = []


def audit_filtering(
    adata_before: AnnData,
    adata_after: AnnData,
    sample_key: str,
    group_key: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """Audit cell retention after filtering.

    Computes per-sample before/after counts, number removed, and retention rate.
    If ``group_key`` is provided, also computes a per-group summary.

    Parameters
    ----------
    adata_before, adata_after
        AnnData objects before and after filtering. ``adata_after`` should be a
        subset of ``adata_before``.
    sample_key
        Column in ``.obs`` identifying samples.
    group_key
        Optional column in ``.obs`` for a biological group summary.

    Returns:
    -------
    Dict[str, pd.DataFrame]
        ``{"sample": sample_df, "group": group_df}``. ``group`` is only present
        when ``group_key`` is provided.

    Examples:
    --------
    >>> audit = scl.qc.audit_filtering(
    ...     adata_before, adata_after, sample_key="sampleID", group_key="group"
    ... )
    >>> print(audit["sample"])
    """
    if sample_key not in adata_before.obs.columns:
        raise ValueError(f"sample_key '{sample_key}' not found in adata_before.obs")
    if sample_key not in adata_after.obs.columns:
        raise ValueError(f"sample_key '{sample_key}' not found in adata_after.obs")

    before = adata_before.obs.groupby(sample_key, observed=True).size()
    after = adata_after.obs.groupby(sample_key, observed=True).size()
    sample_audit = pd.DataFrame({"before": before, "after": after}).fillna(0).astype(int)
    sample_audit["removed"] = sample_audit["before"] - sample_audit["after"]
    sample_audit["retention_rate"] = (sample_audit["after"] / sample_audit["before"]).round(3)

    print("=== Retention by sample ===")
    print(sample_audit.reset_index().to_string(index=False))

    result: Dict[str, pd.DataFrame] = {"sample": sample_audit.reset_index()}

    if group_key and group_key in adata_before.obs.columns and group_key in adata_after.obs.columns:
        before_g = adata_before.obs.groupby(group_key, observed=True).size()
        after_g = adata_after.obs.groupby(group_key, observed=True).size()
        group_audit = pd.DataFrame({"before": before_g, "after": after_g}).fillna(0).astype(int)
        group_audit["removed"] = group_audit["before"] - group_audit["after"]
        group_audit["retention_rate"] = (group_audit["after"] / group_audit["before"]).round(3)
        print("=== Retention by group ===")
        print(group_audit.reset_index().to_string(index=False))
        result["group"] = group_audit.reset_index()

    return result


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
    log.info(f"Cells failing {name} ({operator} {threshold}): {count} ({percentage:.2f}%)")

    return result


def _identify_outliers_subset(
    obs_subset: pd.DataFrame,
    metrics: List[Tuple[str, str, Optional[float]]],
    nmads: float = 4.0,
    group_name: str = "global",
) -> pd.Series:
    """
    Internal helper function to identify outliers on a subset of data.
    """
    subset_outliers = pd.Series(False, index=obs_subset.index)

    for metric, direction, threshold in metrics:
        if metric not in obs_subset.columns:
            log.warning(f"Metric '{metric}' not found in data for group '{group_name}', skipping.")
            continue

        values = obs_subset[metric]
        metric_outliers = pd.Series(False, index=obs_subset.index)

        if threshold is not None:
            # Use fixed threshold
            if direction == "upper":
                metric_outliers = values > threshold
            elif direction == "lower":
                metric_outliers = values < threshold
            elif direction == "both":
                # For fixed threshold, 'both' is not meaningful.
                # A user should provide two separate tuples for upper and lower bounds.
                log.warning(
                    f"Direction 'both' with a fixed threshold is ambiguous for '{metric}'. "
                    "Please provide separate 'upper' and 'lower' tuples if needed. Skipping."
                )
                continue
            else:
                log.warning(f"Invalid direction '{direction}' for '{metric}', skipping.")
                continue
        else:
            # Calculate threshold using canonical MAD implementation
            lower_bound, upper_bound = compute_mad_bounds(
                values.values, nmads=nmads, direction=direction
            )

            # Detect degenerate case (all values identical → MAD == 0)
            if lower_bound == upper_bound:
                log.warning(
                    f"MAD is zero for '{metric}' in group '{group_name}'. "
                    "Cannot perform outlier detection for this metric."
                )
                continue

            if direction == "upper":
                metric_outliers = values > upper_bound
            elif direction == "lower":
                metric_outliers = values < lower_bound
            elif direction == "both":
                metric_outliers = (values > upper_bound) | (values < lower_bound)
            else:
                log.warning(f"Invalid direction '{direction}' for '{metric}', skipping.")
                continue

        outlier_count = metric_outliers.sum()
        if outlier_count > 0:
            log.info(
                f"  - Group '{group_name}': Identified {outlier_count} outliers "
                f"({outlier_count / len(values):.2%}) for metric '{metric}' (direction: {direction})"
            )

        subset_outliers |= metric_outliers

    return subset_outliers


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
        custom_cols = [col for col in adata.obs.columns if col.startswith("outlier_custom_")]
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

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows), facecolor="white")
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

                if col_data.dtype == "bool" or set(col_data.unique()).issubset({0, 1}):
                    # Boolean data - color by outlier status
                    colors = col_data.map(
                        {
                            False: "#637b8a",  # Blue for normal cells
                            True: "#d62728",  # Red for outliers
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
                    cmap=(
                        "viridis"
                        if not isinstance(colors, pd.Series) or colors.dtype != "object"
                        else None
                    ),
                )

                ax.set_title(title, fontsize=10)
                ax.set_xlabel("Total Counts", fontsize=9)
                ax.set_ylabel("Number of Genes", fontsize=9)
                ax.tick_params(labelsize=8)

                # Add colorbar for continuous data
                if not (
                    col_data.dtype == "bool" or set(col_data.unique()).issubset({0, 1})
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

        if show and _is_interactive_backend():
            plt.show()
        else:
            plt.close(fig)


# --- Main Functions ---
def identify_outliers(
    adata: AnnData,
    metrics: List[Tuple[str, str, Optional[float]]],
    sample_key: Optional[str] = None,
    nmads: float = 4.0,
) -> pd.Series:
    """
    Identify outliers based on metrics using median absolute deviation (MAD) or fixed thresholds.

    This function can process multiple metrics and optionally group by sample for per-group
    outlier detection.

    Args:
        adata: AnnData object to check for outliers.
        metrics: List of tuples for outlier detection. Each tuple is (metric, direction, threshold).
                 - metric (str): Column name in `adata.obs`.
                 - direction (str): 'upper', 'lower', or 'both'.
                 - threshold (float, optional): If provided, this fixed value is used as the threshold.
                   If None, the threshold is calculated dynamically using MAD.
        sample_key: If provided, outliers will be identified separately per sample group.
        nmads: Number of median absolute deviations for dynamic outlier detection.

    Returns:
        Boolean pd.Series indicating if a cell is an outlier for any of the specified metrics.
    """
    if not metrics:
        return pd.Series(False, index=adata.obs_names)

    final_outliers = pd.Series(False, index=adata.obs_names)

    if sample_key and sample_key in adata.obs.columns:
        log.info(f"Identifying outliers per group in '{sample_key}'...")
        for sample_id, group_df in adata.obs.groupby(sample_key, observed=False):
            group_outliers = _identify_outliers_subset(
                group_df, metrics, nmads, group_name=str(sample_id)
            )
            final_outliers[group_outliers.index] = group_outliers
    else:
        log.info("Identifying outliers on the entire dataset...")
        global_outliers = _identify_outliers_subset(adata.obs, metrics, nmads, group_name="global")
        final_outliers = global_outliers

    total_count = final_outliers.sum()
    log.info(
        f"Total unique outliers identified: {total_count} ({total_count / len(final_outliers):.2%})"
    )

    return final_outliers


def mark_low_quality_cell(
    adata: AnnData,
    sample_key: str = "sampleID",
    config: Optional[MarkingConfig] = None,
    sample_thresholds: Optional[Dict[str, Dict[str, Any]]] = None,
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
    cfg = MarkingConfig()
    if config is not None:
        cfg = config.model_copy(deep=True)
        # Allow threshold overrides via kwargs (e.g., min_genes=300)
        threshold_overrides = {
            key: value for key, value in kwargs.items()
            if hasattr(cfg.thresholds, key)
        }
        top_level_overrides = {
            key: value for key, value in kwargs.items()
            if key not in threshold_overrides and hasattr(cfg, key)
        }
        if threshold_overrides:
            cfg = cfg.model_copy(
                update={
                    "thresholds": cfg.thresholds.model_copy(update=threshold_overrides)
                },
                deep=True,
            )
        if top_level_overrides:
            cfg = cfg.model_copy(update=top_level_overrides, deep=True)
    elif kwargs:
        threshold_overrides = {
            key: value for key, value in kwargs.items()
            if hasattr(cfg.thresholds, key)
        }
        top_level_overrides = {
            key: value for key, value in kwargs.items()
            if key not in threshold_overrides and hasattr(cfg, key)
        }
        if threshold_overrides:
            cfg = cfg.model_copy(
                update={
                    "thresholds": cfg.thresholds.model_copy(update=threshold_overrides)
                },
                deep=True,
            )
        if top_level_overrides:
            cfg = cfg.model_copy(update=top_level_overrides, deep=True)

    thresholds = cfg.thresholds

    # Check required QC columns. pct_counts_mt is only required when the user
    # has configured an MT threshold; datasets without mitochondrial genes
    # (e.g., synthetic test fixtures) may legitimately lack this column.
    required = ["total_counts", "n_genes_by_counts"]
    if thresholds.pc_mt is not None:
        required.append("pct_counts_mt")
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
        sample: adata.obs[sample_key] == sample for sample in adata.obs[sample_key].unique()
    }

    log.info(f"Processing {len(sample_indices)} samples with {adata.n_obs} total cells")

    # === THRESHOLD CHECKS ===

    # --- Helper to resolve per-sample thresholds with both lower/upper support ---
    def _resolve_sample_threshold(
        metric_name: str,
        global_lower_val: Optional[float],
        global_upper_val: Optional[float],
    ) -> Tuple[pd.Series, pd.Series]:
        """Return (lower_outliers, upper_outliers) for a metric."""
        lower_out = pd.Series(False, index=adata.obs_names)
        upper_out = pd.Series(False, index=adata.obs_names)
        has_sample_th = (
            sample_thresholds is not None
            and metric_name in next(iter(sample_thresholds.values()), {})
        )
        if has_sample_th:
            for sample, idx in sample_indices.items():
                st = sample_thresholds.get(sample, {}).get(metric_name, {})
                th_low = st.get("lower")
                th_up = st.get("upper")
                if th_low is not None:
                    lower_out.loc[idx] = adata.obs.loc[idx, metric_name] < th_low
                if th_up is not None:
                    upper_out.loc[idx] = adata.obs.loc[idx, metric_name] > th_up
        else:
            if global_lower_val is not None:
                lower_out = _safe_threshold_check(
                    adata.obs[metric_name], global_lower_val, "<", f"{metric_name}_lower"
                )
            if global_upper_val is not None:
                upper_out = _safe_threshold_check(
                    adata.obs[metric_name], global_upper_val, ">", f"{metric_name}_upper"
                )
        return lower_out, upper_out

    # Gene count thresholds
    outlier_min_genes, outlier_max_genes = _resolve_sample_threshold(
        "n_genes_by_counts",
        thresholds.min_genes,
        thresholds.max_genes,
    )
    adata.obs["outlier_min_genes"] = outlier_min_genes
    adata.obs["outlier_max_genes"] = outlier_max_genes

    # Total count thresholds
    outlier_min_counts, outlier_max_counts = _resolve_sample_threshold(
        "total_counts",
        thresholds.min_counts,
        thresholds.max_counts,
    )
    adata.obs["outlier_min_counts"] = outlier_min_counts
    adata.obs["outlier_max_counts"] = outlier_max_counts

    # Mitochondrial percentage threshold
    outlier_mt_lower, outlier_mt = _resolve_sample_threshold(
        "pct_counts_mt",
        None,  # no lower bound for MT%
        thresholds.pc_mt,
    )
    adata.obs["outlier_mt"] = outlier_mt | outlier_mt_lower

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
    formatted_metrics = [(metric, direction, None) for metric, direction in cfg.qc_metrics_mad]

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
    log.info(f"Cells marked as QC metric outliers: {qc_count} ({qc_count / adata.n_obs:.2%})")

    # === CUSTOM OUTLIER DETECTION ===

    if cfg.custom_outlier_functions:
        log.info(
            f"Running {len(cfg.custom_outlier_functions)} custom outlier detection functions..."
        )
        for func_name, func in cfg.custom_outlier_functions.items():
            try:
                custom_outliers = func(adata)
                if not isinstance(custom_outliers, pd.Series):
                    raise ValueError(f"Custom function {func_name} must return a pandas Series")
                if len(custom_outliers) != adata.n_obs:
                    raise ValueError(f"Custom function {func_name} returned wrong length")

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
        "sample_thresholds": sample_thresholds if sample_thresholds else {},
    }

    # === SUMMARY STATISTICS ===

    log.info("\n" + "=" * 50)
    log.info("LOW-QUALITY CELL DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Count cells with different types of issues
    outlier_cols = [col for col in adata.obs.columns if col.startswith("outlier_")]
    adata.obs["outlier_count"] = adata.obs[outlier_cols].sum(axis=1).astype(int)

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
            log.info(f"{doublet_col.replace('_', ' ').title()}: {count} cells ({percentage:.2f}%)")

    log.info("=" * 50)

    # === VISUALIZATION ===

    if cfg.plot_outliers:
        _plot_qc_outliers(adata, sample_indices, cfg.cols_to_plot, cfg.save_dir, cfg.show_plots)

    # === Final Type Casting for Robustness ===
    log.info("Finalizing data types for all 'outlier_' columns to ensure save compatibility.")
    outlier_cols_to_cast = [col for col in adata.obs.columns if col.startswith("outlier_")]

    for col in outlier_cols_to_cast:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].fillna(False).astype(bool)

    record_mark_evidence(
        adata,
        evidence_columns=outlier_cols_to_cast + ["outlier_count"],
        thresholds=thresholds.to_dict(),
        sample_key=sample_key,
        sample_thresholds=sample_thresholds,
    )

    return adata


def mark_low_quality_cells_adaptive(
    adata: AnnData,
    batch_key: str = "sampleID",
    metrics: List[str] = ["n_genes_by_counts", "pct_counts_mt"],
    method: str = "hierarchical",
    **kwargs,
) -> AnnData:
    """
    Enhanced version of mark_low_quality_cell with batch-aware thresholds.

    This is particularly useful for datasets with strong batch effects
    (e.g., multi-center studies, fresh vs. frozen samples).
    """
    warnings.warn(
        "mark_low_quality_cells_adaptive is a legacy helper. Prefer "
        "recommend_qc_thresholds/decide_qc_thresholds followed by "
        "build_qc_decisions or run_qc(...), which records unified "
        "qc_decision/qc_reason evidence.",
        FutureWarning,
        stacklevel=2,
    )
    calculator = AdaptiveThresholdCalculator(adata, batch_key)

    for metric in metrics:
        log.info(f"Calculating adaptive thresholds for {metric}...")

        # Analyze batch effects (side effect: logs batch effect info)
        calculator._calculate_batch_effects(metric)

        # Get adaptive thresholds
        thresholds = calculator._suggest_adaptive_thresholds(metric, method=method)

        # Apply batch-specific thresholds
        outlier_mask = pd.Series(False, index=adata.obs_names)

        for batch, batch_thresholds in thresholds.items():
            batch_mask = adata.obs[batch_key] == batch
            values = adata.obs.loc[batch_mask, metric]

            if metric in ["n_genes_by_counts", "total_counts"]:
                # Lower threshold for count metrics
                batch_outliers = values < batch_thresholds["lower"]
            else:
                # Upper threshold for percentage metrics
                batch_outliers = values > batch_thresholds["upper"]

            outlier_mask.loc[batch_mask] = batch_outliers.astype(bool).to_numpy()

        adata.obs[f"outlier_{metric}_adaptive"] = outlier_mask

        log.info(
            f"Marked {outlier_mask.sum()} cells as outliers for {metric} "
            f"({outlier_mask.sum() / len(outlier_mask):.2%})"
        )

    return adata
