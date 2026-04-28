"""Benchmark metrics and report templates for QC recommendation evaluation."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
from anndata import AnnData
from scipy import sparse as sp

QC_BENCHMARK_SCHEMA_VERSION = "1.1"

PBMC_MARKER_SETS: dict[str, list[str]] = {
    "T_cells": ["CD3D", "CD3E", "TRAC", "IL7R"],
    "B_cells": ["MS4A1", "CD79A", "CD74"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1"],
    "Monocytes": ["LYZ", "S100A8", "S100A9", "FCGR3A"],
    "Dendritic": ["FCER1A", "CST3"],
}

TISSUE_MARKER_SETS: dict[str, list[str]] = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Immune": ["PTPRC", "LST1", "TYROBP"],
}

TUMOR_MARKER_SETS: dict[str, list[str]] = {
    "Epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA"],
    "EMT": ["VIM", "ZEB1", "SNAI2"],
    "Immune_microenvironment": ["PTPRC", "CD3D", "LYZ"],
    "Stromal_microenvironment": ["COL1A1", "COL1A2", "DCN"],
}

BENCHMARK_PROFILES: dict[str, dict[str, Any]] = {
    "pbmc": {
        "label": "PBMC / immune suspension",
        "min_retention": 0.55,
        "max_retention": 0.98,
        "min_marker_fidelity": 0.60,
        "marker_sets": PBMC_MARKER_SETS,
    },
    "tissue": {
        "label": "Primary tissue",
        "min_retention": 0.45,
        "max_retention": 0.98,
        "min_marker_fidelity": 0.55,
        "marker_sets": TISSUE_MARKER_SETS,
    },
    "tumor": {
        "label": "Tumor tissue",
        "min_retention": 0.35,
        "max_retention": 0.99,
        "min_marker_fidelity": 0.50,
        "marker_sets": TUMOR_MARKER_SETS,
    },
    "cell_line": {
        "label": "Cell line",
        "min_retention": 0.65,
        "max_retention": 0.995,
        "min_marker_fidelity": 0.50,
        "marker_sets": {},
    },
}


def infer_qc_benchmark_profile(
    *,
    dataset_type: str | None = None,
    tissue_type: str | None = None,
    tissue: str | None = None,
) -> str:
    """Infer the benchmark profile used to score QC output."""
    text = " ".join(str(x).lower() for x in [dataset_type, tissue_type, tissue] if x)
    if "tumor" in text or "cancer" in text or "malignan" in text:
        return "tumor"
    if "pbmc" in text or "blood" in text or "immune" in text:
        return "pbmc"
    if "cell_line" in text or "cell line" in text or "cellline" in text:
        return "cell_line"
    return "tissue"


def default_marker_sets_for_profile(profile: str) -> dict[str, list[str]]:
    """Return default marker sets for a benchmark profile."""
    return dict(BENCHMARK_PROFILES.get(profile, BENCHMARK_PROFILES["tissue"])["marker_sets"])


def detect_marker_sets_from_var(adata: AnnData) -> dict[str, list[str]]:
    """Detect synthetic/curated marker sets stored as boolean columns in ``adata.var``."""
    marker_sets: dict[str, list[str]] = {}
    for column in adata.var.columns:
        if not str(column).startswith("marker_"):
            continue
        values = np.asarray(adata.var[column]).astype(bool)
        if values.sum() == 0:
            continue
        marker_sets[str(column).replace("marker_", "")] = adata.var_names[values].tolist()
    return marker_sets


def compute_retention_metrics(
    adata_before: AnnData,
    adata_after: AnnData,
    *,
    sample_key: str | None = None,
    cell_type_key: str | None = None,
) -> dict[str, Any]:
    """Compute overall and stratified cell-retention metrics."""
    before_obs = adata_before.obs
    after_index = set(adata_after.obs_names)
    n_before = int(adata_before.n_obs)
    n_after = int(adata_after.n_obs)

    result: dict[str, Any] = {
        "initial_cells": n_before,
        "final_cells": n_after,
        "removed_cells": max(0, n_before - n_after),
        "retention_rate": float(n_after / n_before) if n_before else 0.0,
        "per_sample": {},
        "per_cell_type": {},
    }

    def _stratified_retention(key: str) -> dict[str, dict[str, Any]]:
        output: dict[str, dict[str, Any]] = {}
        if key not in before_obs:
            return output
        for group, index in before_obs.groupby(key, observed=False).groups.items():
            names = set(index)
            retained = len(names & after_index)
            total = len(names)
            output[str(group)] = {
                "initial_cells": int(total),
                "final_cells": int(retained),
                "retention_rate": float(retained / total) if total else 0.0,
            }
        return output

    if sample_key:
        result["per_sample"] = _stratified_retention(sample_key)
    if cell_type_key:
        result["per_cell_type"] = _stratified_retention(cell_type_key)
    return result


def compute_marker_fidelity(
    adata_before: AnnData,
    adata_after: AnnData,
    *,
    marker_sets: Mapping[str, list[str]] | None = None,
) -> dict[str, Any]:
    """Measure how well known marker signal is preserved after QC filtering."""
    markers = dict(marker_sets or {})
    if not markers:
        markers = detect_marker_sets_from_var(adata_before)

    per_set: dict[str, dict[str, Any]] = {}
    scores: list[float] = []
    before_names = {gene.upper(): gene for gene in adata_before.var_names}
    after_names = {gene.upper(): gene for gene in adata_after.var_names}

    for set_name, genes in markers.items():
        available_before = [before_names[g.upper()] for g in genes if g.upper() in before_names]
        available_after = [after_names[g.upper()] for g in genes if g.upper() in after_names]
        common = sorted(set(available_before) & set(available_after))
        if not common:
            per_set[set_name] = {
                "available_markers": [],
                "n_markers": 0,
                "status": "not_available",
                "fidelity_score": None,
            }
            continue

        before_expr = _mean_marker_expression(adata_before, common)
        after_expr = _mean_marker_expression(adata_after, common)
        before_detected = _marker_detection_fraction(adata_before, common)
        after_detected = _marker_detection_fraction(adata_after, common)
        expression_ratio = _safe_ratio(after_expr, before_expr)
        detection_ratio = _safe_ratio(after_detected, before_detected)
        score = _bounded_preservation_score(expression_ratio, detection_ratio)
        scores.append(score)

        per_set[set_name] = {
            "available_markers": common,
            "n_markers": len(common),
            "mean_expression_before": before_expr,
            "mean_expression_after": after_expr,
            "mean_expression_ratio": expression_ratio,
            "detection_fraction_before": before_detected,
            "detection_fraction_after": after_detected,
            "detection_fraction_delta": after_detected - before_detected,
            "fidelity_score": score,
            "status": "ok",
        }

    return {
        "available": bool(scores),
        "overall_marker_fidelity": float(np.mean(scores)) if scores else None,
        "n_marker_sets_available": len(scores),
        "per_marker_set": per_set,
    }


def evaluate_qc_benchmark(
    adata_before: AnnData,
    adata_after: AnnData,
    *,
    dataset_type: str | None = None,
    tissue_type: str | None = None,
    tissue: str | None = None,
    sample_key: str | None = None,
    cell_type_key: str | None = None,
    marker_sets: Mapping[str, list[str]] | None = None,
) -> dict[str, Any]:
    """Evaluate QC output against profile-specific retention and marker criteria."""
    profile = infer_qc_benchmark_profile(
        dataset_type=dataset_type,
        tissue_type=tissue_type,
        tissue=tissue,
    )
    profile_spec = BENCHMARK_PROFILES.get(profile, BENCHMARK_PROFILES["tissue"])
    markers = dict(marker_sets or profile_spec.get("marker_sets", {}))
    if not markers:
        markers = detect_marker_sets_from_var(adata_before)

    retention = compute_retention_metrics(
        adata_before,
        adata_after,
        sample_key=sample_key if sample_key in adata_before.obs else None,
        cell_type_key=cell_type_key if cell_type_key in adata_before.obs else None,
    )
    marker_fidelity = compute_marker_fidelity(
        adata_before,
        adata_after,
        marker_sets=markers,
    )

    checks = _evaluate_benchmark_checks(retention, marker_fidelity, profile_spec)
    assessment = build_qc_benchmark_assessment(
        checks=checks,
        retention=retention,
        marker_fidelity=marker_fidelity,
        profile=profile,
        profile_spec=profile_spec,
    )
    return {
        "schema_version": QC_BENCHMARK_SCHEMA_VERSION,
        "profile": profile,
        "profile_label": profile_spec["label"],
        "criteria": {
            "min_retention": profile_spec["min_retention"],
            "max_retention": profile_spec["max_retention"],
            "min_marker_fidelity": profile_spec["min_marker_fidelity"],
        },
        "retention": retention,
        "marker_fidelity": marker_fidelity,
        "checks": checks,
        "assessment": assessment,
        "status": assessment["status"],
    }


def build_qc_benchmark_assessment(
    *,
    checks: list[Mapping[str, Any]],
    retention: Mapping[str, Any],
    marker_fidelity: Mapping[str, Any],
    profile: str,
    profile_spec: Mapping[str, Any],
) -> dict[str, Any]:
    """Convert benchmark checks into a publication-facing QC assessment."""
    failed = [check for check in checks if not check.get("passed")]
    critical = [check for check in failed if check.get("severity") == "critical"]
    high = [check for check in failed if check.get("severity") == "high"]

    if critical:
        status = "fail"
        risk_level = "critical"
    elif high:
        status = "review_required"
        risk_level = "high"
    elif failed:
        status = "review_required"
        risk_level = "moderate"
    else:
        status = "pass"
        risk_level = "low"

    recommendations = [
        {
            "priority": _priority_for_severity(str(check.get("severity", "moderate"))),
            "action": str(check.get("recommendation")),
            "rationale": str(check.get("interpretation") or ""),
            "evidence_key": f"benchmark_summary.checks.{check.get('name')}",
        }
        for check in failed
        if check.get("recommendation")
    ]
    if not recommendations and status == "pass":
        recommendations.append(
            {
                "priority": "optional",
                "action": "Archive QC benchmark outputs with the analysis record.",
                "rationale": "All configured benchmark checks passed for the inferred profile.",
                "evidence_key": "benchmark_summary.assessment",
            }
        )

    return _json_ready(
        {
            "status": status,
            "risk_level": risk_level,
            "review_required": status != "pass",
            "summary": _benchmark_summary_sentence(
                status=status,
                profile_label=str(profile_spec.get("label", profile)),
                retention_rate=retention.get("retention_rate"),
                marker_fidelity=marker_fidelity.get("overall_marker_fidelity"),
                failed_checks=len(failed),
            ),
            "reasons": [
                str(check.get("interpretation") or check.get("name")) for check in failed
            ],
            "recommendations": recommendations,
            "profile_assumptions": {
                "profile": profile,
                "label": profile_spec.get("label", profile),
                "min_retention": profile_spec.get("min_retention"),
                "max_retention": profile_spec.get("max_retention"),
                "min_marker_fidelity": profile_spec.get("min_marker_fidelity"),
            },
        }
    )


def export_qc_benchmark_report(
    benchmark: Mapping[str, Any],
    output_dir: str | Path,
    *,
    prefix: str = "qc_benchmark",
) -> dict[str, str]:
    """Export benchmark results as JSON and Markdown."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / f"{prefix}.json"
    md_path = out / f"{prefix}.md"
    json_path.write_text(json.dumps(benchmark, indent=2, default=str), encoding="utf-8")
    md_path.write_text(render_qc_benchmark_markdown(benchmark), encoding="utf-8")
    return {"json": str(json_path), "markdown": str(md_path)}


def render_qc_benchmark_markdown(benchmark: Mapping[str, Any]) -> str:
    """Render a human-reviewable QC benchmark summary."""
    retention = benchmark.get("retention", {})
    marker = benchmark.get("marker_fidelity", {})
    assessment = benchmark.get("assessment", {})
    lines = [
        "# QC Benchmark Summary",
        "",
        f"- **Profile**: {benchmark.get('profile_label')} ({benchmark.get('profile')})",
        f"- **Status**: {benchmark.get('status')}",
        f"- **Risk level**: {assessment.get('risk_level', 'unknown')}",
        f"- **Retention rate**: {_format_percent(retention.get('retention_rate'))}",
        f"- **Marker fidelity**: {_format_float(marker.get('overall_marker_fidelity'))}",
        "",
        "## Assessment",
        "",
        str(assessment.get("summary", "No benchmark assessment was generated.")),
        "",
        "### Recommended Actions",
        "",
        "| Priority | Action | Rationale |",
        "|----------|--------|-----------|",
    ]
    for item in assessment.get("recommendations", []):
        lines.append(
            f"| {item.get('priority')} | {item.get('action')} | {item.get('rationale')} |"
        )
    if not assessment.get("recommendations"):
        lines.append("| optional | No benchmark action required. | All checks passed or unavailable. |")

    lines.extend(
        [
            "",
            "## Checks",
            "",
            "| Check | Passed | Severity | Value | Threshold | Interpretation |",
            "|-------|--------|----------|-------|-----------|----------------|",
        ]
    )
    for check in benchmark.get("checks", []):
        lines.append(
            f"| {check.get('name')} | {check.get('passed')} | "
            f"{check.get('severity', 'review')} | {check.get('value')} | "
            f"{check.get('threshold')} | {check.get('interpretation', '')} |"
        )

    lines.extend(["", "## Marker Sets", "", "| Marker set | Markers | Fidelity | Status |"])
    lines.append("|------------|---------|----------|--------|")
    for name, payload in marker.get("per_marker_set", {}).items():
        lines.append(
            f"| {name} | {payload.get('n_markers', 0)} | "
            f"{_format_float(payload.get('fidelity_score'))} | {payload.get('status')} |"
        )
    return "\n".join(lines)


def render_qc_benchmark_compact_markdown(benchmark: Mapping[str, Any]) -> str:
    """Render a compact benchmark block for embedding in larger reports."""
    assessment = benchmark.get("assessment", {})
    retention = benchmark.get("retention", {})
    marker = benchmark.get("marker_fidelity", {})
    lines = [
        "### QC Benchmark",
        "",
        f"- Profile: {benchmark.get('profile_label')} ({benchmark.get('profile')})",
        f"- Status: {benchmark.get('status')}",
        f"- Risk: {assessment.get('risk_level', 'unknown')}",
        f"- Retention: {_format_percent(retention.get('retention_rate'))}",
        f"- Marker fidelity: {_format_float(marker.get('overall_marker_fidelity'))}",
        f"- Summary: {assessment.get('summary', 'No assessment generated.')}",
        "",
        "Action items:",
        "",
    ]
    for item in assessment.get("recommendations", []):
        lines.append(f"- [{item.get('priority')}] {item.get('action')}")
    if not assessment.get("recommendations"):
        lines.append("- [optional] No benchmark action required.")
    return "\n".join(lines)


def _evaluate_benchmark_checks(
    retention: Mapping[str, Any],
    marker_fidelity: Mapping[str, Any],
    profile_spec: Mapping[str, Any],
) -> list[dict[str, Any]]:
    retention_rate = float(retention.get("retention_rate", 0.0))
    fidelity = marker_fidelity.get("overall_marker_fidelity")
    checks = [
        {
            "name": "minimum_retention",
            "passed": retention_rate >= profile_spec["min_retention"],
            "value": retention_rate,
            "threshold": profile_spec["min_retention"],
            "severity": "critical" if retention_rate < 0.05 else "high",
            "interpretation": (
                "Cell retention is below the profile-specific lower bound; QC may be "
                "over-filtering or input quality may be poor."
            ),
            "recommendation": (
                "Review min_genes/min_counts/mitochondrial and doublet filters; inspect "
                "per-sample retention before preprocessing."
            ),
        },
        {
            "name": "maximum_retention",
            "passed": retention_rate <= profile_spec["max_retention"],
            "value": retention_rate,
            "threshold": profile_spec["max_retention"],
            "severity": "moderate",
            "interpretation": (
                "Cell retention is above the profile-specific upper bound; QC may be too "
                "permissive or filters may not have been applied."
            ),
            "recommendation": (
                "Confirm low-quality, high-mitochondrial, and doublet flags are active "
                "and biologically appropriate for this dataset."
            ),
        },
    ]
    checks.extend(_stratified_retention_checks(retention, profile_spec))
    if fidelity is not None:
        fidelity = float(fidelity)
        checks.append(
            {
                "name": "marker_fidelity",
                "passed": fidelity >= profile_spec["min_marker_fidelity"],
                "value": fidelity,
                "threshold": profile_spec["min_marker_fidelity"],
                "severity": "high",
                "interpretation": (
                    "Known marker signal is poorly preserved after QC, suggesting "
                    "selective loss of biological populations or excessive filtering."
                ),
                "recommendation": (
                    "Inspect marker-set retention and compare filtered versus unfiltered "
                    "embeddings or cell-type summaries."
                ),
            }
        )
    else:
        checks.append(
            {
                "name": "marker_fidelity",
                "passed": True,
                "value": None,
                "threshold": profile_spec["min_marker_fidelity"],
                "severity": "informational",
                "interpretation": "No configured marker set was available in this dataset.",
                "recommendation": (
                    "Provide tissue- or study-specific marker sets for stronger QC benchmarking."
                ),
                "note": "No configured marker set was available in this dataset.",
            }
        )
    return checks


def _stratified_retention_checks(
    retention: Mapping[str, Any],
    profile_spec: Mapping[str, Any],
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    min_retention = float(profile_spec["min_retention"])
    for key, label in [
        ("per_sample", "sample"),
        ("per_cell_type", "cell_type"),
    ]:
        strata = retention.get(key, {})
        if not isinstance(strata, Mapping) or len(strata) < 2:
            continue
        rates = [
            float(payload.get("retention_rate", 0.0))
            for payload in strata.values()
            if isinstance(payload, Mapping)
        ]
        if len(rates) < 2:
            continue
        min_rate = min(rates)
        max_rate = max(rates)
        spread = max_rate - min_rate
        label_text = label.replace("_", " ")
        checks.append(
            {
                "name": f"minimum_{label}_retention",
                "passed": min_rate >= min_retention * 0.75,
                "value": min_rate,
                "threshold": min_retention * 0.75,
                "severity": "high" if label == "sample" else "moderate",
                "interpretation": (
                    f"At least one {label_text} has much lower retention than expected, "
                    "which can bias downstream comparisons."
                ),
                "recommendation": (
                    f"Review {label_text}-level QC distributions and consider sample-aware "
                    "thresholds or documenting biological justification."
                ),
            }
        )
        checks.append(
            {
                "name": f"{label}_retention_spread",
                "passed": spread <= 0.35,
                "value": spread,
                "threshold": 0.35,
                "severity": "moderate",
                "interpretation": (
                    f"Retention varies substantially across {label_text} strata, which may "
                    "confound downstream differential analyses."
                ),
                "recommendation": (
                    f"Inspect retained/removed fractions by {label_text} before interpreting "
                    "downstream abundance or expression shifts."
                ),
            }
        )
    return checks


def _mean_marker_expression(adata: AnnData, genes: list[str]) -> float:
    if adata.n_obs == 0 or not genes:
        return 0.0
    X = adata[:, genes].X
    return float(X.mean() if not sp.issparse(X) else X.mean())


def _marker_detection_fraction(adata: AnnData, genes: list[str]) -> float:
    if adata.n_obs == 0 or not genes:
        return 0.0
    X = adata[:, genes].X
    detected = np.asarray((X > 0).sum(axis=1)).ravel() if sp.issparse(X) else (X > 0).sum(axis=1)
    return float(np.mean(detected > 0))


def _safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator <= 0:
        return None if numerator <= 0 else 1.0
    return float(numerator / denominator)


def _bounded_preservation_score(*ratios: float | None) -> float:
    scores = []
    for ratio in ratios:
        if ratio is None:
            continue
        if ratio <= 0:
            scores.append(0.0)
        else:
            scores.append(float(min(ratio, 1.0 / ratio, 1.0)))
    return float(np.mean(scores)) if scores else 1.0


def _format_percent(value: Any) -> str:
    return "N/A" if value is None else f"{float(value):.1%}"


def _format_float(value: Any) -> str:
    return "N/A" if value is None else f"{float(value):.3f}"


def _priority_for_severity(severity: str) -> str:
    if severity == "critical":
        return "blocking"
    if severity == "high":
        return "required"
    if severity == "moderate":
        return "review"
    return "optional"


def _benchmark_summary_sentence(
    *,
    status: str,
    profile_label: str,
    retention_rate: Any,
    marker_fidelity: Any,
    failed_checks: int,
) -> str:
    retention_text = _format_percent(retention_rate)
    marker_text = _format_float(marker_fidelity)
    if status == "pass":
        return (
            f"QC output passed the {profile_label} benchmark with {retention_text} "
            f"cell retention and marker fidelity {marker_text}."
        )
    if status == "fail":
        return (
            f"QC output failed the {profile_label} benchmark; {failed_checks} critical "
            "or high-risk check(s) require correction before downstream analysis."
        )
    return (
        f"QC output requires review under the {profile_label} benchmark; "
        f"{failed_checks} check(s) flagged potential retention or marker-fidelity issues."
    )


def _json_ready(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [_json_ready(v) for v in value]
    if hasattr(value, "item"):
        return value.item()
    return value
