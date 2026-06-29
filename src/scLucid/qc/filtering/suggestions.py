"""QC threshold suggestion and report generation.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
import re
import warnings
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from ..adaptive_threshold import compute_mad_bounds
from ..config import QCThresholds

log = logging.getLogger(__name__)

def suggest_qc_thresholds(
    adata: AnnData,
    method: Literal["mad", "iqr", "percentile"] = "mad",
    mad_multipliers: Union[float, List[float]] = [3.0, 4.0, 5.0],
    iqr_multiplier: float = 1.5,
    percentile_range: Tuple[float, float] = (2.5, 97.5),
    plot_distributions: bool = True,
    save_dir: Optional[str] = None,
) -> Tuple[pd.DataFrame, QCThresholds]:
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
        Tuple containing:
        - pd.DataFrame: A DataFrame with QC metrics as rows and suggestion levels
                        (e.g., 'mad_x3.0') as columns.
        - QCThresholds: A QCThresholds object with suggested values based on the
                        first MAD multiplier or the default setting, for convenience.
    """
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required QC columns: {missing_cols}")

    log.info(f"Suggesting QC thresholds using '{method}' method...")

    if isinstance(mad_multipliers, (int, float)):
        mad_multipliers = [mad_multipliers]

    all_suggestions = {}

    # Define which metrics to analyze
    metrics = {
        "n_genes_by_counts": "Gene counts per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial percentage",
    }
    if "pct_counts_hb" in adata.obs.columns:
        metrics["pct_counts_hb"] = "Hemoglobin percentage"

    top_gene_cols = [
        col for col in adata.obs.columns if re.match(r"pct_counts_in_top_\d+_genes", col)
    ]
    for col in top_gene_cols:
        metrics[col] = (
            col.replace("_", " ").replace("pct counts in ", "").replace(" genes", "").title()
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

        metric_map = {
            "n_genes_by_counts": ("min_genes", "max_genes"),
            "total_counts": ("min_counts", "max_counts"),
            "pct_counts_mt": "pc_mt",
            "pct_counts_hb": "pc_hb",
        }
        # Dynamically add top gene cols to map
        for col in top_gene_cols:
            metric_map[col] = f"pc_{col.split('pct_counts_in_')[-1]}"

        if method == "mad":
            # Use canonical compute_mad_bounds for consistency with the rest of the QC module.
            # This includes the 1.4826 scale factor and proper MAD=0 handling.
            for multiplier in mad_multipliers:
                level_name = f"mad_x{multiplier}"
                all_suggestions.setdefault(level_name, {})

                lower_bound, upper_bound = compute_mad_bounds(
                    data.values, nmads=multiplier, direction="both"
                )

                # When MAD is zero, compute_mad_bounds returns bounds == median.
                # Fall back to percentile-based bounds for more robust thresholds.
                if lower_bound == upper_bound:
                    log.warning(
                        f"MAD for metric '{metric}' is zero at MAD x{multiplier}. "
                        f"Falling back to percentile-based bounds."
                    )
                    lower_pct, upper_pct = data.quantile([0.05, 0.95])
                    lower_bound = float(lower_pct)
                    upper_bound = float(upper_pct)

                if is_count_metric:
                    lower_bound = max(0.0, lower_bound)
                    min_key, max_key = metric_map[metric]
                    all_suggestions[level_name][min_key] = int(lower_bound)
                    all_suggestions[level_name][max_key] = int(upper_bound)
                else:  # Percentage metric
                    key = metric_map.get(metric)
                    if key:
                        all_suggestions[level_name][key] = min(100.0, upper_bound)

                if is_count_metric:
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
                        "color": "orange" if is_count_metric else "red",
                    }
                )

        elif method == "iqr":
            level_name = f"iqr_x{iqr_multiplier}"
            all_suggestions.setdefault(level_name, {})
            q25, q75 = data.quantile([0.25, 0.75])
            iqr = q75 - q25
            upper_bound = q75 + iqr_multiplier * iqr

            if is_count_metric:
                lower_bound = max(0, q25 - iqr_multiplier * iqr)
                min_key, max_key = metric_map[metric]
                all_suggestions[level_name][min_key] = int(lower_bound)
                all_suggestions[level_name][max_key] = int(upper_bound)
                # Add lines for plotting
                plot_lines.append(
                    {
                        "val": lower_bound,
                        "label": f"Min (IQR x{iqr_multiplier})",
                        "color": "red",
                    }
                )
                plot_lines.append(
                    {
                        "val": upper_bound,
                        "label": f"Max (IQR x{iqr_multiplier})",
                        "color": "orange",
                    }
                )
            else:  # Percentage metric
                key = metric_map.get(metric)
                if key:
                    all_suggestions[level_name][key] = min(100.0, upper_bound)
                # Add line for plotting
                plot_lines.append(
                    {
                        "val": upper_bound,
                        "label": f"Max (IQR x{iqr_multiplier})",
                        "color": "red",
                    }
                )

        elif method == "percentile":
            level_name = f"percentile_{percentile_range[0]}-{percentile_range[1]}"
            all_suggestions.setdefault(level_name, {})
            upper_bound = data.quantile(percentile_range[1] / 100)

            if is_count_metric:
                lower_bound = data.quantile(percentile_range[0] / 100)
                min_key, max_key = metric_map[metric]
                all_suggestions[level_name][min_key] = int(lower_bound)
                all_suggestions[level_name][max_key] = int(upper_bound)
                # Add lines for plotting
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
            else:  # Percentage metric
                key = metric_map.get(metric)
                if key:
                    all_suggestions[level_name][key] = min(100.0, upper_bound)
                # Add line for plotting
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
    suggested_thresholds_df = pd.DataFrame.from_dict(all_suggestions, orient="index")

    # reorder columns
    cols_order = [
        "min_genes",
        "max_genes",
        "min_counts",
        "max_counts",
        "pc_mt",
        "pc_hb",
    ]
    top_gene_cols_sorted = sorted(
        [c for c in suggested_thresholds_df.columns if c.startswith("pc_top_")]
    )
    final_cols = [
        c for c in cols_order if c in suggested_thresholds_df.columns
    ] + top_gene_cols_sorted
    suggested_thresholds_df = suggested_thresholds_df[final_cols]

    # Create default thresholds object
    default_thresholds_obj = QCThresholds()
    if not suggested_thresholds_df.empty:
        default_series = suggested_thresholds_df.iloc[0]
        pc_top_genes_dict = {k: v for k, v in default_series.items() if k.startswith("pc_top_")}

        final_kwargs = {
            k: v for k, v in default_series.items() if not k.startswith("pc_top_") and pd.notna(v)
        }
        final_kwargs["pc_top_genes"] = pc_top_genes_dict

        default_thresholds_obj = QCThresholds(**final_kwargs)

    log.info("Comparison of recommended QC thresholds:")
    log.info("\n" + suggested_thresholds_df.to_string())

    return suggested_thresholds_df, default_thresholds_obj



def resolve_qc_thresholds(
    *,
    intelligent: Optional[Dict[str, Any]] = None,
    mad: Optional[Dict[str, Any]] = None,
    manual: Optional[Dict[str, Any]] = None,
    policy: Literal["intelligent_then_mad", "mad_then_intelligent", "manual_override"] = "intelligent_then_mad",
) -> QCThresholds:
    """Merge QC threshold sources into a single ``QCThresholds`` object.

    This is the canonical helper for combining data-driven suggestions (Intelligent
    QC, MAD) with manual project overrides. Manual overrides are treated as
    safety bounds:

    - ``manual.min_genes`` / ``manual.min_counts`` act as floors.
    - ``manual.max_genes`` / ``manual.max_counts`` / ``manual.pc_mt`` act as
      ceilings.

    Parameters
    ----------
    intelligent
        Threshold dict from ``recommend_intelligent_qc`` (keys like ``min_genes``,
        ``min_counts``, ``pc_mt``).
    mad
        Threshold dict from ``suggest_qc_thresholds`` (a ``QCThresholds.to_dict()``).
    manual
        Manual overrides. Missing keys are ignored.
    policy
        Which data-driven source takes priority when both are provided.

    Returns:
    -------
    QCThresholds
        Resolved thresholds.

    Examples:
    --------
    >>> thresholds = resolve_qc_thresholds(
    ...     intelligent={"min_genes": 300, "pc_mt": 8.0},
    ...     mad={"min_genes": 200, "max_genes": 8000, "pc_mt": 10.0},
    ...     manual={"min_genes": 500, "pc_mt": 12.0},
    ... )
    """
    intelligent = intelligent or {}
    mad = mad or {}
    manual = manual or {}

    def _pick_data(key: str) -> Any:
        i_val = intelligent.get(key)
        m_val = mad.get(key)
        if policy == "manual_override":
            return manual.get(key) if manual.get(key) is not None else (i_val if i_val is not None else m_val)
        if policy == "intelligent_then_mad":
            return i_val if i_val is not None else m_val
        # mad_then_intelligent
        return m_val if m_val is not None else i_val

    min_genes = _pick_data("min_genes")
    max_genes = _pick_data("max_genes")
    min_counts = _pick_data("min_counts")
    max_counts = _pick_data("max_counts")
    pc_mt = _pick_data("pc_mt")
    pc_hb = _pick_data("pc_hb")
    pc_top_genes = _pick_data("pc_top_genes") or {}

    # Manual overrides as floors/ceilings
    if manual.get("min_genes") is not None and min_genes is not None:
        min_genes = max(min_genes, manual["min_genes"])
    elif manual.get("min_genes") is not None:
        min_genes = manual["min_genes"]

    if manual.get("min_counts") is not None and min_counts is not None:
        min_counts = max(min_counts, manual["min_counts"])
    elif manual.get("min_counts") is not None:
        min_counts = manual["min_counts"]

    if manual.get("max_genes") is not None and max_genes is not None:
        max_genes = min(max_genes, manual["max_genes"])
    elif manual.get("max_genes") is not None:
        max_genes = manual["max_genes"]

    if manual.get("max_counts") is not None and max_counts is not None:
        max_counts = min(max_counts, manual["max_counts"])
    elif manual.get("max_counts") is not None:
        max_counts = manual["max_counts"]

    if manual.get("pc_mt") is not None and pc_mt is not None:
        pc_mt = min(pc_mt, manual["pc_mt"])
    elif manual.get("pc_mt") is not None:
        pc_mt = manual["pc_mt"]

    if manual.get("pc_hb") is not None and pc_hb is not None:
        pc_hb = min(pc_hb, manual["pc_hb"])
    elif manual.get("pc_hb") is not None:
        pc_hb = manual["pc_hb"]

    kwargs: Dict[str, Any] = {}
    if min_genes is not None:
        kwargs["min_genes"] = min_genes
    if max_genes is not None:
        kwargs["max_genes"] = max_genes
    if min_counts is not None:
        kwargs["min_counts"] = min_counts
    if max_counts is not None:
        kwargs["max_counts"] = max_counts
    if pc_mt is not None:
        kwargs["pc_mt"] = pc_mt
    if pc_hb is not None:
        kwargs["pc_hb"] = pc_hb
    if pc_top_genes:
        kwargs["pc_top_genes"] = pc_top_genes

    return QCThresholds(**kwargs)


def generate_qc_report(
    adata: AnnData,
    save_dir: str,
    sample_key: str = "sampleID",
    include_before_after: bool = True,
    adata_before: Optional[AnnData] = None,
) -> None:
    """Compatibility wrapper for the reporting-layer implementation."""
    warnings.warn(
        "scLucid.qc.filtering.generate_qc_report is deprecated; use "
        "scLucid.qc.generate_qc_report or scLucid.qc.reporting.generate_qc_report.",
        FutureWarning,
        stacklevel=2,
    )
    from ..reporting import generate_qc_report as _generate_qc_report

    return _generate_qc_report(
        adata=adata,
        save_dir=save_dir,
        sample_key=sample_key,
        include_before_after=include_before_after,
        adata_before=adata_before,
    )
