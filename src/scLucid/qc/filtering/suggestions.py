"""QC threshold suggestion and report generation.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from scipy import stats

from ..adaptive_threshold import compute_mad_bounds
from ..config import QCThresholds
from .core import _plot_before_after_comparison

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
            median_val = data.median()
            mad_val = np.median(np.abs(data - median_val))
            if mad_val == 0:
                log.warning(
                    f"MAD for metric '{metric}' is zero. MAD-based thresholds may be unreliable."
                )
                mad_val = 1e-5  # Small value to avoid division by zero

            for multiplier in mad_multipliers:
                level_name = f"mad_x{multiplier}"
                all_suggestions.setdefault(level_name, {})

                upper_bound = median_val + multiplier * mad_val
                if is_count_metric:
                    lower_bound = max(0, median_val - multiplier * mad_val)
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
        log.warning("`adata_before` not provided, cannot generate before/after comparison plots.")
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

        ax.boxplot(sample_data, tick_labels=sample_labels)
        ax.set_title(f"{metric.replace('_', ' ').title()} Distribution by Sample")
        ax.set_ylabel(metric.replace("_", " ").title())

        if len(samples) > 5:
            ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(Path(save_dir) / "qc_distributions.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 3. Before/after comparison if requested
    if include_before_after and adata_before is not None:
        _plot_before_after_comparison(adata_before, adata, save_dir, sample_key, qc_metrics)

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

    qc_trace = adata.uns.get("sclucid", {}).get("qc", {})
    trace_context = qc_trace.get("context", {}).get("data", {})
    trace_recommendation = qc_trace.get("recommendation", {}).get("data", {})
    trace_warnings = qc_trace.get("warnings", {}).get("data", [])
    trace_filtering = qc_trace.get("filtering_summary", {}).get("data", {})
    trace_thresholds = qc_trace.get("sample_thresholds", {}).get("data", {})
    trace_tumor_flags = qc_trace.get("tumor_aware_flags", {}).get("data", {})

    report_summary = {
        "dataset_shape_after": [adata.n_obs, adata.n_vars],
        "dataset_shape_before": (
            [adata_before.n_obs, adata_before.n_vars] if adata_before is not None else None
        ),
        "context": trace_context,
        "recommendation": trace_recommendation,
        "filtering_summary": trace_filtering,
        "tumor_aware_flags": trace_tumor_flags,
        "warnings": trace_warnings,
        "sample_thresholds": trace_thresholds,
    }
    (Path(save_dir) / "qc_summary.json").write_text(
        json.dumps(report_summary, indent=2, default=str)
    )

    md_lines = [
        "# QC Summary",
        "",
        f"- **Cells before**: {adata_before.n_obs if adata_before is not None else 'NA'}",
        f"- **Cells after**: {adata.n_obs}",
        f"- **Genes**: {adata.n_vars}",
        f"- **Threshold mode**: {trace_context.get('threshold_mode', 'NA')}",
        f"- **Strategy**: {trace_recommendation.get('overall_strategy', 'NA')}",
        f"- **Overall confidence**: {trace_recommendation.get('overall_confidence', 'NA')}",
        f"- **Tissue type**: {trace_context.get('tissue_type', 'NA')}",
        "",
        "## Filtering",
        "",
        f"- **Criteria used**: {', '.join(trace_filtering.get('criteria_used', [])) if trace_filtering else 'NA'}",
        f"- **Removed cells**: {trace_filtering.get('removed_cells', 'NA')}",
        f"- **Removed fraction**: {trace_filtering.get('removed_fraction', 'NA')}",
        "",
        "## Concerns",
        "",
    ]

    concerns = trace_recommendation.get("concerns", []) if trace_recommendation else []
    if concerns:
        md_lines.extend([f"- {concern}" for concern in concerns])
    else:
        md_lines.append("- None")

    md_lines.extend(["", "## Warnings", ""])
    if trace_warnings:
        md_lines.extend([f"- {warning}" for warning in trace_warnings])
    else:
        md_lines.append("- None")

    if trace_tumor_flags:
        md_lines.extend(["", "## Tumor-aware Flags", "", "```json"])
        md_lines.append(json.dumps(trace_tumor_flags, indent=2, default=str))
        md_lines.append("```")

    (Path(save_dir) / "qc_summary.md").write_text("\n".join(md_lines))

    try:
        from ..reporting import generate_qc_html_report

        generate_qc_html_report(
            adata,
            output_path=str(Path(save_dir) / "qc_report.html"),
            adata_before=adata_before,
            title="scLucid Quality Control Report",
        )
    except Exception as exc:
        log.warning(f"Enhanced QC HTML report generation skipped: {exc}")

    log.info("QC report generation completed")
