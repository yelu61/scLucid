"""
Enhanced HTML report generation for QC results.

This module provides comprehensive, publication-ready HTML reports
with embedded visualizations and recommendations.
"""

import json
import logging
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from .benchmark import compute_retention_metrics
from .filtering.core import _plot_before_after_comparison

try:
    import plotly.graph_objects as go
except ImportError:
    go = None  # type: ignore[misc,assignment]

log = logging.getLogger(__name__)


def _infer_report_qc_metrics(adata: AnnData) -> List[str]:
    """Return QC metrics that are actually present and useful for reporting."""
    preferred = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"]
    metrics = [metric for metric in preferred if metric in adata.obs.columns]

    extra_pct_metrics = sorted(
        col
        for col in adata.obs.columns
        if col.startswith("pct_counts_") and col not in metrics
    )
    metrics.extend(extra_pct_metrics)
    return metrics


def generate_qc_report(
    adata: AnnData,
    save_dir: str,
    sample_key: str = "sampleID",
    include_before_after: bool = True,
    adata_before: Optional[AnnData] = None,
) -> None:
    """Generate file-based QC report artifacts.

    This is the canonical implementation for summary CSVs and static plots.
    Filtering modules should mark, filter, and audit cells; report rendering
    belongs in this reporting layer.
    """
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    log.info(f"Generating QC report in {save_dir}")

    if sample_key not in adata.obs.columns:
        raise KeyError(f"sample_key {sample_key!r} not found in adata.obs.")

    if include_before_after and adata_before is None:
        log.warning("`adata_before` not provided, cannot generate before/after comparison plots.")
        include_before_after = False

    qc_metrics = _infer_report_qc_metrics(adata)
    if not qc_metrics:
        raise ValueError("No QC metric columns found in adata.obs for report generation.")

    summary_stats = []
    samples = adata.obs[sample_key].unique()
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs[sample_mask]
        stats = {"sample": sample, "n_cells": len(sample_data)}
        for metric in qc_metrics:
            stats[f"{metric}_mean"] = sample_data[metric].mean()
            stats[f"{metric}_median"] = sample_data[metric].median()
            stats[f"{metric}_std"] = sample_data[metric].std()
        summary_stats.append(stats)

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(Path(save_dir) / "qc_summary_statistics.csv", index=False)

    n_metrics = len(qc_metrics)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(12, 4 * n_metrics))
    if n_metrics == 1:
        axes = [axes]

    for i, metric in enumerate(qc_metrics):
        ax = axes[i]
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

    if include_before_after and adata_before is not None:
        _plot_before_after_comparison(adata_before, adata, save_dir, sample_key, qc_metrics)

    outlier_cols = [col for col in adata.obs.columns if col.startswith("outlier_")]
    if outlier_cols:
        outlier_summary = []
        for sample in samples:
            sample_mask = adata.obs[sample_key] == sample
            sample_data = adata.obs[sample_mask]
            stats = {"sample": sample, "n_cells": len(sample_data)}

            for col in outlier_cols:
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
    trace_review_summary = qc_trace.get("review_summary", {}).get("data", {})
    trace_doublet_summary = (
        trace_review_summary.get("doublet_evidence_summary", {})
        if isinstance(trace_review_summary, dict)
        else {}
    )
    trace_doublet_benchmark_decision = (
        trace_doublet_summary.get("benchmark_decision", {})
        if isinstance(trace_doublet_summary, dict)
        else {}
    )
    trace_reviewer_table = (
        trace_review_summary.get("qc_reviewer_table", [])
        if isinstance(trace_review_summary, dict)
        else []
    )

    report_summary = {
        "dataset_shape_after": [adata.n_obs, adata.n_vars],
        "dataset_shape_before": (
            [adata_before.n_obs, adata_before.n_vars] if adata_before is not None else None
        ),
        "context": trace_context,
        "recommendation": trace_recommendation,
        "filtering_summary": trace_filtering,
        "tumor_aware_flags": trace_tumor_flags,
        "qc_reviewer_table": trace_reviewer_table,
        "doublet_benchmark_decision": trace_doublet_benchmark_decision,
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

    if trace_doublet_benchmark_decision:
        md_lines.extend(["", "## Doublet Benchmark Decision", ""])
        for key in [
            "recommended_default_mode",
            "recommended_primary_method",
            "recommended_algorithm_weight",
            "recommended_fusion_method",
            "review_required",
            "risk_note",
        ]:
            if key in trace_doublet_benchmark_decision:
                md_lines.append(f"- **{key}**: {trace_doublet_benchmark_decision[key]}")

    if trace_reviewer_table:
        md_lines.extend(
            [
                "",
                "## QC Reviewer Table",
                "",
                "| Item | Category | Recommended | Applied | Source | Confidence | Affected cells | Review | Risk note |",
                "| --- | --- | --- | --- | --- | --- | ---: | --- | --- |",
            ]
        )
        reviewer_rows = (
            trace_reviewer_table.values()
            if isinstance(trace_reviewer_table, dict)
            else trace_reviewer_table
        )
        for row in list(reviewer_rows)[:20]:
            if not isinstance(row, dict):
                continue
            md_lines.append(
                "| "
                + " | ".join(
                    [
                        str(row.get("item", "")),
                        str(row.get("category", "")),
                        str(row.get("recommended_value", "")),
                        str(row.get("applied_value", "")),
                        str(row.get("source", "")),
                        str(row.get("confidence", "")),
                        str(row.get("affected_cells", "")),
                        str(row.get("review_required", "")),
                        str(row.get("biological_risk_note", "")).replace("|", "/"),
                    ]
                )
                + " |"
            )

    (Path(save_dir) / "qc_summary.md").write_text("\n".join(md_lines))

    try:
        generate_qc_html_report(
            adata,
            output_path=str(Path(save_dir) / "qc_report.html"),
            adata_before=adata_before,
            title="scLucid Quality Control Report",
        )
    except Exception as exc:
        log.warning(f"Enhanced QC HTML report generation skipped: {exc}")

    log.info("QC report generation completed")


class EnhancedQCReport:
    """
    Generate comprehensive HTML reports for QC analysis.

    Features:
    - Publication-ready formatting
    - Embedded interactive visualizations
    - Automatic recommendations
    - Downloadable tables
    - Statistical summaries
    """

    def __init__(
        self,
        adata: AnnData,
        adata_before: Optional[AnnData] = None,
    ):
        """
        Initialize the report generator.

        Args:
            adata: AnnData after QC filtering
            adata_before: AnnData before QC (for comparison)
        """
        self.adata = adata
        self.adata_before = adata_before or adata

    def generate_html_report(
        self,
        output_path: str,
        title: str = "Quality Control Report",
        author: Optional[str] = None,
        include_plots: bool = True,
        include_recommendations: bool = True,
    ):
        """
        Generate a comprehensive HTML report.

        Args:
            output_path: Path to save the HTML report
            title: Report title
            author: Report author
            include_plots: Whether to include plots
            include_recommendations: Whether to include recommendations
        """
        # Gather report data
        report_data = self._gather_report_data(include_plots=include_plots)

        # Generate HTML
        html_content = self._generate_html(
            title=title,
            author=author,
            data=report_data,
            include_plots=include_plots,
            include_recommendations=include_recommendations,
        )

        # Write to file
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, "w", encoding="utf-8") as f:
            f.write(html_content)

        log.info(f"QC report saved to {output_path}")

    def _gather_report_data(self, include_plots: bool = True) -> Dict[str, Any]:
        """Gather all data for the report."""
        qc_trace = self._get_qc_trace()
        data = {
            "trace": qc_trace,
            "metadata": self._get_metadata(),
            "summary": self._get_summary_statistics(),
            "metrics": self._get_metrics_summary(),
            "filtering": self._get_filtering_summary(),
            "plots": self._get_plot_data() if include_plots else {},
            "recommendations": self._get_recommendations(),
            "tables": self._get_table_data(),
        }
        return data

    def _get_qc_trace(self) -> Dict[str, Any]:
        """Safely read the unified QC trace."""
        return self.adata.uns.get("sclucid", {}).get("qc", {})

    @staticmethod
    def _extract_data(value: Any, default: Any) -> Any:
        if isinstance(value, dict) and "data" in value:
            return value.get("data", default)
        return value if value is not None else default

    def _get_retention_summary(self) -> Dict[str, Any]:
        """Compute canonical retention summary, then enrich with QC trace details."""
        qc_trace = self._get_qc_trace()
        context = self._extract_data(qc_trace.get("context"), {})
        filtering_results = self._extract_data(qc_trace.get("filtering_summary"), {})

        retention = compute_retention_metrics(
            self.adata_before,
            self.adata,
            sample_key=context.get("sample_key"),
        )
        retention["criteria_used"] = filtering_results.get("criteria_used", [])
        retention["combination_logic"] = filtering_results.get("combination_logic")
        retention["criteria_counts"] = filtering_results.get("criteria_counts", {})
        retention["removal_rate"] = 1.0 - float(retention.get("retention_rate", 0.0))
        return retention

    def _get_metadata(self) -> Dict[str, Any]:
        """Get report metadata."""
        qc_trace = self._get_qc_trace()
        context = self._extract_data(qc_trace.get("context"), {})
        recommendation = self._extract_data(qc_trace.get("recommendation"), {})
        retention = self._get_retention_summary()
        return {
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "n_cells_after": self.adata.n_obs,
            "n_cells_before": self.adata_before.n_obs,
            "n_genes": self.adata.n_vars,
            "retention_rate": retention.get("retention_rate", 0.0),
            "sample_key": context.get("sample_key"),
            "n_samples": context.get("n_samples"),
            "threshold_mode": context.get("threshold_mode"),
            "tissue_type": context.get("tissue_type"),
            "strategy": recommendation.get("overall_strategy"),
            "overall_confidence": recommendation.get("overall_confidence"),
        }

    def _get_summary_statistics(self) -> Dict[str, Any]:
        """Get summary statistics."""
        stats = {
            "total_cells": self.adata.n_obs,
            "total_genes": self.adata.n_vars,
        }

        # QC metric statistics
        qc_metrics = [
            "n_genes_by_counts",
            "total_counts",
            "log1p_total_counts",
            "pct_counts_mt",
        ]

        for metric in qc_metrics:
            if metric in self.adata.obs:
                values = self.adata.obs[metric].values
                if len(values) == 0:
                    continue
                stats[metric] = {
                    "mean": float(np.mean(values)),
                    "median": float(np.median(values)),
                    "std": float(np.std(values)),
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                }

        return stats

    def _get_metrics_summary(self) -> List[Dict[str, Any]]:
        """Get per-metric summary."""
        summary = []

        for metric in self.adata.obs.columns:
            if not metric.startswith(("log1p_", "pct_", "n_", "total_")):
                continue

            values = self.adata.obs[metric].values
            if len(values) == 0:
                continue
            summary.append(
                {
                    "name": metric,
                    "mean": float(np.mean(values)),
                    "median": float(np.median(values)),
                    "std": float(np.std(values)),
                }
            )

        return summary

    def _get_filtering_summary(self) -> Dict[str, Any]:
        """Get filtering summary."""
        retention = self._get_retention_summary()
        return {
            "n_before": int(retention.get("initial_cells", self.adata_before.n_obs)),
            "n_after": int(retention.get("final_cells", self.adata.n_obs)),
            "n_removed": int(retention.get("removed_cells", 0)),
            "retention_rate": float(retention.get("retention_rate", 0.0)),
            "removal_rate": float(retention.get("removal_rate", 0.0)),
            "criteria_used": retention.get("criteria_used", []),
            "combination_logic": retention.get("combination_logic"),
            "criteria_counts": retention.get("criteria_counts", {}),
        }

    def _get_recommendations(self) -> List[Dict[str, str]]:
        """Generate QC recommendations."""
        recommendations = []
        qc_trace = self._get_qc_trace()
        recommendation_trace = self._extract_data(qc_trace.get("recommendation"), {})
        warnings = self._extract_data(qc_trace.get("warnings"), [])
        tumor_flags = self._extract_data(qc_trace.get("tumor_aware_flags"), {})

        # Check filtering rate
        retention_rate = float(self._get_retention_summary().get("retention_rate", 0.0))

        if retention_rate > 0.95:
            recommendations.append(
                {
                    "category": "Filtering",
                    "severity": "warning",
                    "message": "Very few cells were filtered. Consider if QC thresholds are appropriate.",
                }
            )
        elif retention_rate < 0.5:
            recommendations.append(
                {
                    "category": "Filtering",
                    "severity": "warning",
                    "message": "More than 50% of cells were removed. Review QC parameters.",
                }
            )

        # Check mitochondrial percentage
        if "pct_counts_mt" in self.adata.obs:
            mt_values = self.adata.obs["pct_counts_mt"].values
            mean_mt = np.mean(mt_values)

            if mean_mt > 20:
                recommendations.append(
                    {
                        "category": "Mitochondrial",
                        "severity": "info",
                        "message": f"High mean MT% ({mean_mt:.1f}%). Possible cell stress or dying cells.",
                    }
                )

        # Check gene counts
        if "n_genes_by_counts" in self.adata.obs:
            gene_counts = self.adata.obs["n_genes_by_counts"].values

            if np.median(gene_counts) < 500:
                recommendations.append(
                    {
                        "category": "Gene Detection",
                        "severity": "warning",
                        "message": f"Low median gene count ({np.median(gene_counts):.0}). Possible low-quality data.",
                    }
                )

        for concern in recommendation_trace.get("concerns", []):
            recommendations.append(
                {
                    "category": "Recommendation Engine",
                    "severity": "warning",
                    "message": concern,
                }
            )

        for consideration in recommendation_trace.get("tumor_specific_considerations", []):
            recommendations.append(
                {
                    "category": "Tumor Context",
                    "severity": "info",
                    "message": consideration,
                }
            )

        if tumor_flags.get("tumor_aware_enabled"):
            recommendations.append(
                {
                    "category": "Tumor-aware QC",
                    "severity": "info",
                    "message": tumor_flags.get(
                        "note",
                        "Tumor-aware QC is enabled; flagged populations should be reviewed before hard filtering.",
                    ),
                }
            )

        for warning in warnings:
            recommendations.append(
                {
                    "category": "Workflow Warning",
                    "severity": "warning",
                    "message": warning,
                }
            )

        return recommendations

    def _get_plot_data(self) -> Dict[str, Any]:
        """Get data for plots."""
        plots = {}

        # Prepare data for violin plots
        for metric in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
            if metric in self.adata.obs:
                plots[metric] = {
                    "values": self.adata.obs[metric].tolist(),
                    "name": metric,
                }

        return plots

    def _get_table_data(self) -> Dict[str, Any]:
        """Get data for tables."""
        tables = {}
        qc_trace = self._get_qc_trace()
        recommendation = self._extract_data(qc_trace.get("recommendation"), {})
        filtering = self._extract_data(qc_trace.get("filtering_summary"), {})
        sample_thresholds = self._extract_data(qc_trace.get("sample_thresholds"), {})
        overrides = self._extract_data(qc_trace.get("user_overrides"), {})

        # Summary statistics table
        summary_data = []
        for metric, stats in self._get_summary_statistics().items():
            if isinstance(stats, dict):
                summary_data.append(
                    {
                        "Metric": metric,
                        "Mean": f"{stats.get('mean', 0):.2f}",
                        "Median": f"{stats.get('median', 0):.2f}",
                        "Std": f"{stats.get('std', 0):.2f}",
                    }
                )

        tables["summary"] = summary_data

        # QC metrics table
        qc_metrics_data = []
        for item in self._get_metrics_summary():
            qc_metrics_data.append(
                {
                    "Metric": item["name"],
                    "Mean": f"{item['mean']:.2f}",
                    "Median": f"{item['median']:.2f}",
                    "Std": f"{item['std']:.2f}",
                }
            )

        tables["metrics"] = qc_metrics_data

        recommendation_rows = []
        for name in ["min_genes", "n_counts", "max_mt_percent", "doublet_threshold"]:
            rec = recommendation.get(name)
            if not isinstance(rec, dict):
                continue
            recommendation_rows.append(
                {
                    "Parameter": name,
                    "Threshold": f"{rec.get('threshold')}",
                    "CI": f"{rec.get('ci_lower')} - {rec.get('ci_upper')}",
                    "Confidence": f"{float(rec.get('confidence', 0.0)):.2f}",
                    "Method": rec.get("method", ""),
                }
            )
        tables["recommendations"] = recommendation_rows

        filtering_rows = []
        for criterion, count in filtering.get("criteria_counts", {}).items():
            filtering_rows.append(
                {
                    "Criterion": criterion,
                    "Flagged Cells": str(count),
                }
            )
        tables["filtering"] = filtering_rows

        threshold_rows = []
        for sample_id, metrics in sample_thresholds.items():
            for metric_name, bounds in metrics.items():
                threshold_rows.append(
                    {
                        "Sample": sample_id,
                        "Metric": metric_name,
                        "Lower": f"{bounds.get('lower')}",
                        "Upper": f"{bounds.get('upper')}",
                    }
                )
        tables["sample_thresholds"] = threshold_rows

        override_rows = []
        for name, values in overrides.items():
            override_rows.append(
                {
                    "Parameter": name,
                    "Recommended": f"{values.get('recommended')}",
                    "User Config": f"{values.get('actual')}",
                }
            )
        tables["overrides"] = override_rows

        return tables

    def _generate_html(
        self,
        title: str,
        author: Optional[str],
        data: Dict[str, Any],
        include_plots: bool,
        include_recommendations: bool,
    ) -> str:
        """Generate complete HTML report."""
        metadata = data["metadata"]
        filtering = data["filtering"]
        recommendations = data.get("recommendations", [])
        tables = data.get("tables", {})
        trace = data.get("trace", {})
        warnings = self._extract_data(trace.get("warnings"), [])

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Georgia, serif;
            line-height: 1.6;
            color: #1f2a2e;
            background:
                radial-gradient(circle at top left, rgba(186, 218, 209, 0.35), transparent 30%),
                linear-gradient(180deg, #f4f1e8 0%, #eef3f2 100%);
        }}

        .container {{
            max-width: 1240px;
            margin: 0 auto;
            padding: 20px;
        }}

        .header {{
            background: linear-gradient(135deg, #274c4d 0%, #486b5c 100%);
            color: #f6f3ec;
            padding: 40px;
            border-radius: 18px;
            margin-bottom: 30px;
            box-shadow: 0 18px 40px rgba(39, 76, 77, 0.18);
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header .meta {{
            font-size: 1.1em;
            opacity: 0.9;
        }}

        .section {{
            background: rgba(255, 255, 255, 0.92);
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 18px;
            box-shadow: 0 10px 24px rgba(32, 43, 45, 0.08);
            border: 1px solid rgba(39, 76, 77, 0.08);
        }}

        .section h2 {{
            color: #274c4d;
            border-bottom: 3px solid #88a99b;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}

        .section h3 {{
            color: #6a4e33;
            margin-top: 25px;
            margin-bottom: 15px;
        }}

        .stat-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}

        .stat-card {{
            background: linear-gradient(135deg, #faf7f0 0%, #e8f0ed 100%);
            color: #274c4d;
            padding: 20px;
            border-radius: 14px;
            text-align: center;
            border: 1px solid rgba(39, 76, 77, 0.08);
        }}

        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            margin-bottom: 5px;
        }}

        .stat-card .label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}

        table th {{
            background: #274c4d;
            color: #f9f7f0;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}

        table td {{
            padding: 12px;
            border-bottom: 1px solid #ddd;
        }}

        table tr:hover {{
            background: #f5f5f5;
        }}

        .recommendation {{
            padding: 15px;
            margin: 10px 0;
            border-radius: 5px;
            border-left: 4px solid #ccc;
        }}

        .recommendation.warning {{
            background: #fff2e0;
            border-left-color: #c17b2f;
        }}

        .recommendation.info {{
            background: #e6f1ef;
            border-left-color: #3f7d7a;
        }}

        .recommendation.success {{
            background: #d4edda;
            border-left-color: #28a745;
        }}

        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: 600;
        }}

        .badge.warning {{
            background: #c17b2f;
            color: #fff;
        }}

        .badge.info {{
            background: #3f7d7a;
            color: #fff;
        }}

        .badge.success {{
            background: #28a745;
            color: #fff;
        }}

        .footer {{
            text-align: center;
            padding: 20px;
            color: #5f6c70;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="meta">
                Generated: {metadata['date']}<br>
                {f"Author: {author}" if author else ""}
            </div>
        </div>

        <div class="section">
            <h2>Executive Summary</h2>
            <div class="stat-grid">
                <div class="stat-card">
                    <div class="value">{metadata['n_cells_after']:,}</div>
                    <div class="label">Cells After QC</div>
                </div>
                <div class="stat-card">
                    <div class="value">{metadata['n_genes']:,}</div>
                    <div class="label">Total Genes</div>
                </div>
                <div class="stat-card">
                    <div class="value">{filtering['retention_rate']:.1%}</div>
                    <div class="label">Retention Rate</div>
                </div>
                <div class="stat-card">
                    <div class="value">{filtering['n_removed']:,}</div>
                    <div class="label">Cells Removed</div>
                </div>
            </div>
            <p><strong>QC strategy</strong>: {escape(str(metadata.get('strategy') or 'not available'))} |
               <strong>Threshold mode</strong>: {escape(str(metadata.get('threshold_mode') or 'not available'))} |
               <strong>Tissue type</strong>: {escape(str(metadata.get('tissue_type') or 'not available'))}</p>
        </div>

        <div class="section">
            <h2>Detailed Metrics</h2>
"""

        # Add summary table
        if "summary" in tables:
            html += self._generate_table("Summary Statistics", tables["summary"])

        # Add metrics table
        if "metrics" in tables:
            html += self._generate_table("QC Metrics", tables["metrics"])
        if tables.get("recommendations"):
            html += self._generate_table("Recommended Thresholds", tables["recommendations"])
        if tables.get("filtering"):
            html += self._generate_table("Filtering Criteria Summary", tables["filtering"])
        if tables.get("sample_thresholds"):
            html += self._generate_table("Per-sample Thresholds", tables["sample_thresholds"])
        if tables.get("overrides"):
            html += self._generate_table("User Overrides", tables["overrides"])

        # Add recommendations
        html += """
        </div>

        <div class="section">
            <h2>Interpretation And Warnings</h2>
"""

        if recommendations:
            for rec in recommendations:
                severity_class = rec["severity"]
                html += f"""
            <div class="recommendation {severity_class}">
                <strong>{escape(rec['category'])}:</strong> <span class="badge {severity_class}">{escape(rec['severity'].upper())}</span><br>
                {escape(rec['message'])}
            </div>
"""
        else:
            html += "<p>No specific recommendations. Your data looks good!</p>"

        if warnings:
            html += "<h3>Workflow warnings</h3><ul>"
            for warning in warnings:
                html += f"<li>{escape(str(warning))}</li>"
            html += "</ul>"

        html += f"""
        </div>

        <div class="footer">
            <p>Generated by scLucid QC</p>
            <p>Report generated on {metadata['date']}</p>
        </div>
    </div>

    <script>
        // Add any interactive JavaScript here
        document.addEventListener('DOMContentLoaded', function() {{
            console.log('QC Report loaded');
        }});
    </script>
</body>
</html>
"""

        return html

    def _generate_table(self, title: str, data: List[Dict[str, str]]) -> str:
        """Generate HTML table from data."""
        if not data:
            return ""

        html = f"<h3>{title}</h3>\n<table>\n"

        # Header
        html += "<thead><tr>"
        for key in data[0].keys():
            html += f"<th>{key}</th>"
        html += "</tr></thead>\n"

        # Body
        html += "<tbody>"
        for row in data:
            html += "<tr>"
            for value in row.values():
                html += f"<td>{value}</td>"
            html += "</tr>\n"
        html += "</tbody>\n"

        html += "</table>\n"

        return html


def generate_qc_html_report(
    adata: AnnData,
    output_path: str,
    adata_before: Optional[AnnData] = None,
    title: str = "Quality Control Report",
    author: Optional[str] = None,
) -> str:
    """
    Convenience function to generate QC HTML report.

    Args:
        adata: AnnData after QC
        output_path: Path to save report
        adata_before: AnnData before QC
        title: Report title
        author: Report author

    Returns:
        Path to generated report
    """
    reporter = EnhancedQCReport(adata, adata_before)
    reporter.generate_html_report(
        output_path=output_path,
        title=title,
        author=author,
    )

    return output_path


class InteractiveReportGenerator:
    """
    Generate interactive HTML reports with embedded JavaScript.

    Creates reports with interactive plots using Plotly.js.
    """

    def __init__(self, adata: AnnData):
        """
        Initialize interactive report generator.

        Args:
            adata: AnnData object with QC results
        """
        self.adata = adata

    def generate_interactive_html(
        self,
        output_path: str,
        title: str = "Interactive QC Report",
    ):
        """
        Generate interactive HTML report with embedded Plotly charts.

        Args:
            output_path: Path to save report
            title: Report title
        """
        if go is None:
            raise ImportError("Plotly is required for interactive reports")

        # Create plots
        plots = self._create_interactive_plots()

        # Generate HTML with embedded plots
        html = self._generate_interactive_html(title, plots)

        # Write to file
        with open(output_path, "w") as f:
            f.write(html)

        log.info(f"Interactive report saved to {output_path}")

    def _create_interactive_plots(self) -> Dict[str, "go.Figure"]:
        """Create interactive Plotly figures."""
        plots = {}

        # Violin plot for key metrics
        metrics_to_plot = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        available_metrics = [m for m in metrics_to_plot if m in self.adata.obs]

        if available_metrics:
            from plotly.subplots import make_subplots

            fig = make_subplots(
                rows=1,
                cols=len(available_metrics),
                subplot_titles=available_metrics,
            )

            for i, metric in enumerate(available_metrics):
                fig.add_trace(
                    go.Violin(
                        y=self.adata.obs[metric],
                        name=metric,
                        box_visible=True,
                        meanline_visible=True,
                    ),
                    row=1,
                    col=i + 1,
                )

            fig.update_layout(
                title_text="QC Metrics Distribution",
                showlegend=False,
                height=400,
            )

            plots["violin"] = fig

        return plots

    def _generate_interactive_html(
        self,
        title: str,
        plots: Dict[str, "go.Figure"],
    ) -> str:
        """Generate HTML with embedded Plotly charts."""
        import plotly.io as pio

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f8f9fa;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        .plot-container {{
            background: white;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <p>Interactive Quality Control Report</p>
        </div>
"""

        # Add plots
        for plot_name, fig in plots.items():
            plot_html = pio.to_html(fig, full_html=False, include_plotlyjs=False)
            html += f"""
        <div class="plot-container">
            <h3>{plot_name.replace('_', ' ').title()}</h3>
            {plot_html}
        </div>
"""

        html += """
    </div>

    <script>
        // Initialize plots
        console.log('Interactive QC Report loaded');
    </script>
</body>
</html>
"""

        return html
