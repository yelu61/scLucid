"""QC trace and review-summary schema helpers.

This module keeps the workflow-facing QC audit contract independent from the
execution code. It turns recommendations, applied config, and filtering output
into a compact machine-readable decision table and evidence chain.
"""

from __future__ import annotations

from collections.abc import Mapping
from enum import Enum
from typing import Any

from anndata import AnnData

from ..utils.evidence import (
    DecisionRecord,
    EvidenceBundle,
    EvidenceItem,
    ReviewAction,
    model_to_dict,
)

QC_TRACE_SCHEMA_VERSION = "1.0"

QC_REQUIRED_REVIEW_SECTIONS = {
    "recommendation_summary",
    "recommended_threshold_summary",
    "applied_threshold_summary",
    "user_override_summary",
    "sample_threshold_summary",
    "tumor_aware_summary",
    "filtering_summary",
    "warnings",
    "decision_table",
    "evidence_chain",
    "execution_trace",
    "output_health",
    "downstream_preprocess_recommendations",
    "qc_readiness",
    "review_action_items",
    "reproducibility_manifest",
    "evidence_bundle",
}

QC_REQUIRED_OBS_METRICS = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
]

_PARAMETER_SPECS = [
    {
        "parameter": "min_genes",
        "recommended_key": "min_genes",
        "applied_path": ("marking_config", "thresholds", "min_genes"),
        "obs_metric": "n_genes_by_counts",
        "filtering_flag": "outlier_min_genes",
        "direction": "lower_bound",
    },
    {
        "parameter": "max_genes",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "max_genes"),
        "obs_metric": "n_genes_by_counts",
        "filtering_flag": "outlier_max_genes",
        "direction": "upper_bound",
    },
    {
        "parameter": "n_counts",
        "recommended_key": "n_counts",
        "applied_path": ("marking_config", "thresholds", "min_counts"),
        "obs_metric": "total_counts",
        "filtering_flag": "outlier_min_counts",
        "direction": "lower_bound",
    },
    {
        "parameter": "max_counts",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "max_counts"),
        "obs_metric": "total_counts",
        "filtering_flag": "outlier_max_counts",
        "direction": "upper_bound",
    },
    {
        "parameter": "max_mt_percent",
        "recommended_key": "max_mt_percent",
        "applied_path": ("marking_config", "thresholds", "pc_mt"),
        "obs_metric": "pct_counts_mt",
        "filtering_flag": "outlier_mt",
        "direction": "upper_bound",
    },
    {
        "parameter": "max_hb_percent",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "pc_hb"),
        "obs_metric": "pct_counts_hb",
        "filtering_flag": "outlier_hb",
        "direction": "upper_bound",
    },
    {
        "parameter": "doublet_threshold",
        "recommended_key": "doublet_threshold",
        "applied_path": ("doublet_config", "score_threshold"),
        "obs_metric": "doublet_score",
        "filtering_flag": "predicted_doublet",
        "direction": "upper_bound",
    },
    {
        "parameter": "nmads",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "nmads"),
        "obs_metric": "qc_metric_distribution",
        "filtering_flag": "outlier_qc_metrics",
        "direction": "mad_outlier",
    },
]


def _json_safe(value: Any) -> Any:
    """Convert common scientific/Python objects into JSON-safe values."""
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def _to_dict(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, dict):
        return _json_safe(value)
    if hasattr(value, "to_dict"):
        return _json_safe(value.to_dict())
    if hasattr(value, "model_dump"):
        return _json_safe(value.model_dump())
    return {}


def _get_nested(data: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = data
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _recommendation_value(rec_dict: Mapping[str, Any], key: str | None, field: str) -> Any:
    if key is None:
        return None
    rec = rec_dict.get(key)
    if not isinstance(rec, Mapping):
        return None
    return rec.get(field)


def _decision_source(
    *,
    parameter: str,
    applied_value: Any,
    recommended_value: Any,
    user_overrides: Mapping[str, Any],
) -> str:
    if parameter in user_overrides:
        return "user_override"
    if recommended_value is not None and applied_value == recommended_value:
        return "recommendation"
    if applied_value is None:
        return "disabled_or_not_available"
    return "default_or_config"


def build_qc_decision_table(
    config: Any,
    original_config: Any,
    recommendation: Any,
    user_overrides: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build a machine-readable table of QC threshold decisions."""
    applied_config = _to_dict(config)
    original_config_dict = _to_dict(original_config)
    rec_dict = _to_dict(recommendation)
    overrides = dict(user_overrides or {})
    criteria = set(
        _get_nested(applied_config, ("filter_config", "criteria_to_filter")) or []
    )
    rows: list[dict[str, Any]] = []

    for spec in _PARAMETER_SPECS:
        recommended_key = spec["recommended_key"]
        applied_value = _get_nested(applied_config, spec["applied_path"])
        user_value = _get_nested(original_config_dict, spec["applied_path"])
        recommended_value = _recommendation_value(rec_dict, recommended_key, "threshold")
        row = {
            "parameter": spec["parameter"],
            "obs_metric": spec["obs_metric"],
            "filtering_flag": spec["filtering_flag"],
            "direction": spec["direction"],
            "recommended": recommended_value,
            "applied": applied_value,
            "user_provided": user_value,
            "source": _decision_source(
                parameter=recommended_key or spec["parameter"],
                applied_value=applied_value,
                recommended_value=recommended_value,
                user_overrides=overrides,
            ),
            "recommendation_method": _recommendation_value(rec_dict, recommended_key, "method"),
            "confidence": _recommendation_value(rec_dict, recommended_key, "confidence"),
            "ci_lower": _recommendation_value(rec_dict, recommended_key, "ci_lower"),
            "ci_upper": _recommendation_value(rec_dict, recommended_key, "ci_upper"),
            "evidence": _recommendation_value(rec_dict, recommended_key, "evidence") or {},
            "is_filtering_enabled": spec["filtering_flag"] in criteria,
        }
        rows.append(_json_safe(row))

    return rows


def build_qc_execution_trace(
    *,
    context: Mapping[str, Any],
    recommendation: Any,
    sample_thresholds: Mapping[str, Any],
    warnings: list[str],
    steps_executed: list[str] | None = None,
) -> dict[str, Any]:
    """Build the high-level execution trace for QC review."""
    rec_dict = _to_dict(recommendation)
    return {
        "qc_schema_version": QC_TRACE_SCHEMA_VERSION,
        "steps_executed": list(steps_executed or []),
        "sample_key": context.get("sample_key"),
        "threshold_mode": context.get("threshold_mode"),
        "n_samples": context.get("n_samples"),
        "tissue_type": context.get("tissue_type"),
        "use_recommendations": context.get("use_recommendations"),
        "recommendation_available": bool(rec_dict),
        "tumor_aware_enabled": _is_tumor_context(context.get("tissue_type")),
        "sample_thresholds_computed": bool(sample_thresholds),
        "warnings_count": len(warnings),
    }


def build_qc_recommended_threshold_summary(
    *,
    recommendation: Any,
    decision_table: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Summarize data-driven threshold recommendations in a stable schema."""
    rec_dict = _to_dict(recommendation)
    parameters: dict[str, Any] = {}
    unavailable: list[str] = []

    for row in decision_table:
        parameter = row.get("parameter")
        if not parameter:
            continue
        recommended = row.get("recommended")
        if recommended is None:
            unavailable.append(str(parameter))
            continue
        parameters[str(parameter)] = {
            "recommended": recommended,
            "applied": row.get("applied"),
            "source": row.get("source"),
            "method": row.get("recommendation_method"),
            "confidence": row.get("confidence"),
            "ci_lower": row.get("ci_lower"),
            "ci_upper": row.get("ci_upper"),
            "evidence": row.get("evidence", {}),
            "filtering_flag": row.get("filtering_flag"),
            "is_filtering_enabled": row.get("is_filtering_enabled"),
        }

    return _json_safe(
        {
            "available": bool(rec_dict),
            "overall_strategy": rec_dict.get("overall_strategy"),
            "overall_confidence": rec_dict.get("overall_confidence"),
            "data_quality_score": rec_dict.get("data_quality_score"),
            "parameters": parameters,
            "unavailable_parameters": unavailable,
        }
    )


def build_qc_output_health(
    adata: AnnData,
    filtering_summary: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Summarize whether QC output is usable for downstream workflow steps."""
    fs = dict(filtering_summary or {})
    initial_cells = fs.get("initial_cells")
    final_cells = fs.get("final_cells", adata.n_obs)
    missing_metrics = [metric for metric in QC_REQUIRED_OBS_METRICS if metric not in adata.obs]
    retention_fraction = None
    if initial_cells not in (None, 0) and final_cells is not None:
        retention_fraction = float(final_cells) / float(initial_cells)

    issues: list[str] = []
    if adata.n_obs == 0:
        issues.append("QC output contains zero cells.")
    if missing_metrics:
        issues.append(f"Missing required QC obs metrics: {', '.join(missing_metrics)}.")
    if retention_fraction is not None and retention_fraction < 0.05:
        issues.append("QC retained fewer than 5% of input cells; thresholds should be reviewed.")

    return {
        "status": "review_required" if issues else "ok",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "initial_cells": initial_cells,
        "final_cells": final_cells,
        "retention_fraction": retention_fraction,
        "missing_required_obs_metrics": missing_metrics,
        "issues": issues,
    }


def build_downstream_preprocess_recommendations(
    *,
    adata: AnnData,
    context: Mapping[str, Any],
    sample_thresholds: Mapping[str, Any],
    filtering_summary: Mapping[str, Any],
    output_health: Mapping[str, Any],
    decision_table: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Recommend safe preprocessing choices based on QC decisions and output health."""
    is_tumor = _is_tumor_context(context.get("tissue_type"))
    n_samples = int(context.get("n_samples") or 1)
    sample_key = context.get("sample_key")
    retention = output_health.get("retention_fraction")
    has_counts_layer = "counts" in adata.layers
    recommendations: list[dict[str, Any]] = []
    blockers = list(output_health.get("issues", []))

    def add(
        *,
        target: str,
        recommendation: str,
        priority: str,
        rationale: str,
        suggested_config: Mapping[str, Any] | None = None,
    ) -> None:
        recommendations.append(
            {
                "target": target,
                "recommendation": recommendation,
                "priority": priority,
                "rationale": rationale,
                "suggested_config": dict(suggested_config or {}),
            }
        )

    add(
        target="counts_layer",
        recommendation=(
            "Use adata.layers['counts'] as the preprocessing input."
            if has_counts_layer
            else "Create or preserve adata.layers['counts'] before normalization when raw counts are available."
        ),
        priority="required" if not has_counts_layer else "recommended",
        rationale="Preprocess needs an auditable raw-count source after QC filtering.",
        suggested_config={"layer": "counts" if has_counts_layer else None},
    )

    add(
        target="normalization",
        recommendation="Run library-size normalization followed by log1p transformation.",
        priority="required",
        rationale="QC has filtered cells but does not normalize expression for HVG, PCA, or clustering.",
        suggested_config={"normalize_total": True, "log1p": True},
    )

    if n_samples > 1:
        add(
            target="batch_aware_hvg",
            recommendation=f"Use sample-aware HVG selection with batch_key='{sample_key}'.",
            priority="recommended",
            rationale=(
                "Multiple samples were detected; sample-aware HVG selection reduces sample-specific "
                "technical dominance while preserving shared biology."
            ),
            suggested_config={"batch_key": sample_key, "sample_thresholds_available": bool(sample_thresholds)},
        )

    if is_tumor:
        add(
            target="tumor_preservation",
            recommendation="Avoid automatic mitochondrial regression or hard MT-based removal before tumor-state review.",
            priority="review",
            rationale=(
                "Tumor-aware QC was active; high mitochondrial signal can reflect malignant state, stress, "
                "hypoxia, or tissue dissociation rather than pure low quality."
            ),
            suggested_config={"regress_out": [], "review_mt_programs": True},
        )
    else:
        mt_row = next(
            (row for row in decision_table if row.get("parameter") == "max_mt_percent"),
            {},
        )
        if mt_row.get("applied") is not None:
            add(
                target="mitochondrial_covariate",
                recommendation="Consider pct_counts_mt as a covariate only after checking biological relevance.",
                priority="optional",
                rationale="QC recorded an MT threshold; downstream regression should remain an explicit choice.",
                suggested_config={"candidate_covariates": ["pct_counts_mt"]},
            )

    if isinstance(retention, (int, float)) and retention < 0.5:
        add(
            target="retention_review",
            recommendation="Inspect QC plots and threshold decisions before continuing to preprocessing.",
            priority="review",
            rationale=f"QC retained {retention:.1%} of cells, which may indicate over-filtering or poor input quality.",
            suggested_config={"review_before_preprocess": True},
        )

    return _json_safe(
        {
            "ready_for_preprocess": output_health.get("status") == "ok",
            "status": "review_required" if blockers else "ready",
            "blockers": blockers,
            "recommendations": recommendations,
            "input_assumptions": {
                "has_counts_layer": has_counts_layer,
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_samples": n_samples,
                "sample_key": sample_key,
                "tumor_aware": is_tumor,
            },
            "filtering_context": {
                "initial_cells": filtering_summary.get("initial_cells"),
                "final_cells": filtering_summary.get("final_cells"),
                "retention_fraction": retention,
                "sample_specific_thresholds_available": bool(sample_thresholds),
            },
        }
    )


def build_qc_readiness_assessment(
    *,
    output_health: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    benchmark_summary: Mapping[str, Any] | None,
    tumor_aware_summary: Mapping[str, Any] | None,
    warnings: list[str],
    decision_table: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Assess whether QC output is ready for downstream analysis."""
    blockers = list(output_health.get("issues", []))
    blockers.extend(downstream_recommendations.get("blockers", []))
    blockers = list(dict.fromkeys(str(item) for item in blockers))

    review_reasons: list[str] = []
    if warnings:
        review_reasons.append(f"{len(warnings)} workflow warning(s) were recorded.")

    tumor_warnings = list((tumor_aware_summary or {}).get("warnings", []))
    review_reasons.extend(tumor_warnings)

    user_overrides = [
        row["parameter"] for row in decision_table if row.get("source") == "user_override"
    ]
    if user_overrides:
        review_reasons.append(
            "User overrides should be documented before publication: "
            + ", ".join(user_overrides)
        )

    benchmark_status = None
    if benchmark_summary:
        benchmark_status = benchmark_summary.get("status")
        if benchmark_status and benchmark_status != "pass":
            review_reasons.append(f"QC benchmark status is {benchmark_status}.")
        assessment = benchmark_summary.get("assessment", {})
        if isinstance(assessment, Mapping):
            for reason in assessment.get("reasons", []):
                review_reasons.append(str(reason))
            if assessment.get("status") == "fail":
                blockers.append(str(assessment.get("summary", "QC benchmark failed.")))

    downstream_status = downstream_recommendations.get("status")
    if downstream_status == "review_required":
        review_reasons.append("Downstream preprocessing recommendations require review.")

    if blockers:
        status = "blocked"
    elif review_reasons:
        status = "review_required"
    else:
        status = "ready"

    score = 100
    score -= min(60, 30 * len(blockers))
    score -= min(30, 8 * len(review_reasons))
    score = max(0, score)

    if status == "ready":
        verdict = "QC output is ready for preprocessing and downstream analysis."
    elif status == "review_required":
        verdict = "QC output can proceed after the listed review items are checked."
    else:
        verdict = "QC output should not proceed until blocking issues are resolved."

    return _json_safe(
        {
            "status": status,
            "score": score,
            "verdict": verdict,
            "blockers": blockers,
            "review_reasons": review_reasons,
            "benchmark_status": benchmark_status,
            "output_health_status": output_health.get("status"),
            "downstream_status": downstream_status,
        }
    )


def build_qc_review_action_items(
    *,
    readiness: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    tumor_aware_summary: Mapping[str, Any] | None,
    benchmark_summary: Mapping[str, Any] | None,
    decision_table: list[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Create human-readable QC review actions from trace evidence."""
    actions: list[dict[str, Any]] = []

    def add(
        *,
        priority: str,
        action: str,
        rationale: str,
        evidence_key: str,
    ) -> None:
        actions.append(
            {
                "priority": priority,
                "action": action,
                "rationale": rationale,
                "evidence_key": evidence_key,
            }
        )

    for blocker in readiness.get("blockers", []):
        add(
            priority="blocking",
            action="Resolve QC output health issue before preprocessing.",
            rationale=str(blocker),
            evidence_key="output_health.issues",
        )

    for item in downstream_recommendations.get("recommendations", []):
        if item.get("priority") in {"required", "review"}:
            add(
                priority=item.get("priority", "review"),
                action=item.get("recommendation", "Review downstream preprocessing choice."),
                rationale=item.get("rationale", ""),
                evidence_key=f"downstream_preprocess_recommendations.{item.get('target')}",
            )

    for warning in (tumor_aware_summary or {}).get("warnings", []):
        add(
            priority="review",
            action="Document tumor-aware QC handling in methods or supplementary QC.",
            rationale=str(warning),
            evidence_key="tumor_aware_summary.warnings",
        )

    overridden = [
        row for row in decision_table if row.get("source") == "user_override"
    ]
    for row in overridden:
        add(
            priority="review",
            action=f"Justify user override for {row.get('parameter')}.",
            rationale=(
                f"Recommended={row.get('recommended')}, applied={row.get('applied')}, "
                f"method={row.get('recommendation_method')}."
            ),
            evidence_key="decision_table",
        )

    if benchmark_summary and benchmark_summary.get("status") != "pass":
        assessment = benchmark_summary.get("assessment", {})
        assessment_actions = assessment.get("recommendations", []) if isinstance(assessment, Mapping) else []
        for item in assessment_actions:
            if not isinstance(item, Mapping):
                continue
            add(
                priority=item.get("priority", "review"),
                action=item.get(
                    "action",
                    "Inspect QC benchmark checks before finalizing downstream analysis.",
                ),
                rationale=item.get("rationale", f"Benchmark status is {benchmark_summary.get('status')}."),
                evidence_key=item.get("evidence_key", "benchmark_summary.checks"),
            )
        if not assessment_actions:
            add(
                priority="review",
                action="Inspect QC benchmark checks before finalizing downstream analysis.",
                rationale=f"Benchmark status is {benchmark_summary.get('status')}.",
                evidence_key="benchmark_summary.checks",
            )

    if not actions:
        add(
            priority="optional",
            action="Archive the QC review summary with downstream analysis outputs.",
            rationale="No blocking or mandatory review items were detected.",
            evidence_key="review_summary",
        )

    priority_order = {"blocking": 0, "required": 1, "review": 2, "optional": 3}
    actions.sort(key=lambda item: priority_order.get(str(item.get("priority")), 9))
    return _json_safe(actions)


def build_qc_reproducibility_manifest(
    *,
    adata: AnnData,
    config: Any,
    original_config: Any,
    context: Mapping[str, Any],
    steps_executed: list[str] | None,
    decision_table: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Record reproducibility-critical state for the QC run."""
    required_obs_present = {
        metric: metric in adata.obs for metric in QC_REQUIRED_OBS_METRICS
    }
    applied_thresholds = {
        row.get("parameter"): row.get("applied") for row in decision_table
    }
    threshold_sources = {
        row.get("parameter"): row.get("source") for row in decision_table
    }

    return _json_safe(
        {
            "schema_version": QC_TRACE_SCHEMA_VERSION,
            "workflow": "run_standard_qc",
            "storage_path": 'adata.uns["sclucid"]["qc"]["review_summary"]["data"]',
            "steps_executed": list(steps_executed or []),
            "data_shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
            "layers_present": sorted(str(key) for key in adata.layers.keys()),
            "required_obs_metrics_present": required_obs_present,
            "context": dict(context),
            "applied_thresholds": applied_thresholds,
            "threshold_sources": threshold_sources,
            "applied_config": _to_dict(config),
            "original_config": _to_dict(original_config),
        }
    )


def build_qc_evidence_chain(
    *,
    recommendation: Any,
    sample_thresholds: Mapping[str, Any],
    filtering_summary: Mapping[str, Any],
    output_health: Mapping[str, Any],
    decision_table: list[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Build a compact ordered evidence chain for QC decisions."""
    rec_dict = _to_dict(recommendation)
    overrides = [row["parameter"] for row in decision_table if row.get("source") == "user_override"]
    return [
        {
            "stage": "recommendation",
            "available": bool(rec_dict),
            "strategy": rec_dict.get("overall_strategy"),
            "confidence": rec_dict.get("overall_confidence"),
            "data_quality_score": rec_dict.get("data_quality_score"),
            "concerns": rec_dict.get("concerns", []),
        },
        {
            "stage": "threshold_application",
            "n_decisions": len(decision_table),
            "user_overrides": overrides,
            "recommendation_driven": [
                row["parameter"] for row in decision_table if row.get("source") == "recommendation"
            ],
        },
        {
            "stage": "sample_thresholds",
            "computed": bool(sample_thresholds),
            "n_samples_with_thresholds": len(sample_thresholds),
        },
        {
            "stage": "filtering",
            "initial_cells": filtering_summary.get("initial_cells"),
            "final_cells": filtering_summary.get("final_cells"),
            "removed_cells": filtering_summary.get("removed_cells"),
            "removed_fraction": filtering_summary.get("removed_fraction"),
            "criteria_used": filtering_summary.get("criteria_used", []),
        },
        {
            "stage": "output_health",
            "status": output_health.get("status"),
            "issues": output_health.get("issues", []),
        },
    ]


def _evidence_source_for_stage(stage: Any) -> str:
    mapping = {
        "recommendation": "recommendation",
        "threshold_application": "metric",
        "sample_thresholds": "metric",
        "filtering": "metric",
        "output_health": "output_health",
    }
    return mapping.get(str(stage), "metric")


def build_qc_evidence_bundle(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Convert QC-specific review fields into the shared EvidenceBundle schema."""
    decisions: list[DecisionRecord] = []
    for row in summary.get("decision_table", []):
        if not isinstance(row, Mapping):
            continue
        evidence = []
        row_evidence = row.get("evidence")
        if row_evidence:
            evidence.append(
                EvidenceItem(
                    source="recommendation",
                    name=f"{row.get('parameter')}_recommendation_evidence",
                    value=row_evidence,
                    confidence=row.get("confidence"),
                    rationale=f"Evidence attached to {row.get('parameter')} recommendation.",
                    related_keys=["decision_table"],
                )
            )
        decisions.append(
            DecisionRecord(
                parameter=str(row.get("parameter")),
                recommended=row.get("recommended"),
                applied=row.get("applied"),
                source=str(row.get("source") or "unknown"),
                confidence=row.get("confidence"),
                evidence=evidence,
                user_override=row.get("source") == "user_override",
                downstream_impact=(
                    f"Controls {row.get('filtering_flag')}"
                    if row.get("filtering_flag")
                    else None
                ),
            )
        )

    evidence_chain: list[EvidenceItem] = []
    for item in summary.get("evidence_chain", []):
        if not isinstance(item, Mapping):
            continue
        stage = item.get("stage", "unknown")
        confidence = item.get("confidence")
        evidence_chain.append(
            EvidenceItem(
                source=_evidence_source_for_stage(stage),
                name=str(stage),
                value=dict(item),
                confidence=confidence if isinstance(confidence, (int, float)) else None,
                rationale=f"QC evidence stage: {stage}.",
                related_keys=["evidence_chain"],
            )
        )

    for issue in summary.get("output_health", {}).get("issues", []):
        evidence_chain.append(
            EvidenceItem(
                source="output_health",
                name="output_health_issue",
                value=issue,
                rationale="Output health issue requiring review.",
                limitations=[str(issue)],
                related_keys=["output_health.issues"],
            )
        )

    benchmark_summary = summary.get("benchmark_summary", {})
    if isinstance(benchmark_summary, Mapping):
        assessment = benchmark_summary.get("assessment", {})
        evidence_chain.append(
            EvidenceItem(
                source="benchmark",
                name="qc_benchmark_assessment",
                value={
                    "status": benchmark_summary.get("status"),
                    "profile": benchmark_summary.get("profile"),
                    "risk_level": assessment.get("risk_level") if isinstance(assessment, Mapping) else None,
                    "summary": assessment.get("summary") if isinstance(assessment, Mapping) else None,
                },
                rationale="Profile-aware benchmark assessment for QC output.",
                limitations=[
                    "Benchmark thresholds are heuristic and should be interpreted with dataset context."
                ],
                related_keys=["benchmark_summary.assessment"],
            )
        )

    action_items = [
        ReviewAction(
            priority=item.get("priority", "review"),
            action=str(item.get("action", "")),
            rationale=str(item.get("rationale", "")),
            evidence_keys=[str(item.get("evidence_key"))] if item.get("evidence_key") else [],
        )
        for item in summary.get("review_action_items", [])
        if isinstance(item, Mapping)
    ]

    readiness = summary.get("qc_readiness", {})
    confidence = None
    if isinstance(readiness, Mapping) and isinstance(readiness.get("score"), (int, float)):
        confidence = max(0.0, min(1.0, float(readiness["score"]) / 100.0))

    bundle = EvidenceBundle(
        module="qc",
        stage="run_standard_qc",
        status=str(readiness.get("status", "unknown")) if isinstance(readiness, Mapping) else "unknown",
        confidence=confidence,
        context=dict(summary.get("execution_trace", {})),
        decisions=decisions,
        evidence_chain=evidence_chain,
        action_items=action_items,
        reproducibility=dict(summary.get("reproducibility_manifest", {})),
        related_review_keys=[
            "decision_table",
            "evidence_chain",
            "qc_readiness",
            "review_action_items",
            "reproducibility_manifest",
            "benchmark_summary",
        ],
    )
    return model_to_dict(bundle)


def enrich_qc_review_summary(
    summary: dict[str, Any],
    *,
    adata: AnnData,
    config: Any,
    original_config: Any,
    recommendation: Any,
    sample_thresholds: Mapping[str, Any],
    filtering_summary: Mapping[str, Any],
    warnings: list[str],
    context: Mapping[str, Any],
    steps_executed: list[str] | None = None,
) -> dict[str, Any]:
    """Add benchmark-grade QC audit fields to the review summary."""
    user_overrides = summary.get("user_override_summary", {}).get("details", {})
    decision_table = build_qc_decision_table(
        config,
        original_config,
        recommendation,
        user_overrides=user_overrides,
    )
    output_health = build_qc_output_health(adata, filtering_summary)
    benchmark_summary = summary.get("benchmark_summary", {})
    summary["qc_schema_version"] = QC_TRACE_SCHEMA_VERSION
    summary["decision_table"] = decision_table
    summary["recommended_threshold_summary"] = build_qc_recommended_threshold_summary(
        recommendation=recommendation,
        decision_table=decision_table,
    )
    summary["execution_trace"] = build_qc_execution_trace(
        context=context,
        recommendation=recommendation,
        sample_thresholds=sample_thresholds,
        warnings=warnings,
        steps_executed=steps_executed,
    )
    summary["output_health"] = output_health
    summary["evidence_chain"] = build_qc_evidence_chain(
        recommendation=recommendation,
        sample_thresholds=sample_thresholds,
        filtering_summary=filtering_summary,
        output_health=output_health,
        decision_table=decision_table,
    )
    summary["required_obs_metrics"] = list(QC_REQUIRED_OBS_METRICS)
    downstream_recommendations = build_downstream_preprocess_recommendations(
        adata=adata,
        context=context,
        sample_thresholds=sample_thresholds,
        filtering_summary=filtering_summary,
        output_health=output_health,
        decision_table=decision_table,
    )
    summary["downstream_preprocess_recommendations"] = downstream_recommendations
    readiness = build_qc_readiness_assessment(
        output_health=output_health,
        downstream_recommendations=downstream_recommendations,
        benchmark_summary=benchmark_summary,
        tumor_aware_summary=summary.get("tumor_aware_summary", {}),
        warnings=warnings,
        decision_table=decision_table,
    )
    summary["qc_readiness"] = readiness
    summary["review_action_items"] = build_qc_review_action_items(
        readiness=readiness,
        downstream_recommendations=downstream_recommendations,
        tumor_aware_summary=summary.get("tumor_aware_summary", {}),
        benchmark_summary=benchmark_summary,
        decision_table=decision_table,
    )
    summary["reproducibility_manifest"] = build_qc_reproducibility_manifest(
        adata=adata,
        config=config,
        original_config=original_config,
        context=context,
        steps_executed=steps_executed,
        decision_table=decision_table,
    )
    summary["evidence_bundle"] = build_qc_evidence_bundle(summary)
    return _json_safe(summary)


def validate_qc_review_summary(
    summary: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """Validate QC-specific review-summary sections."""
    errors: list[str] = []
    missing = sorted(QC_REQUIRED_REVIEW_SECTIONS - set(summary.keys()))
    if missing:
        errors.append(f"QC review summary missing required sections: {missing}")
    if not isinstance(summary.get("decision_table"), list):
        errors.append("QC review summary field 'decision_table' must be a list.")
    if not isinstance(summary.get("evidence_chain"), list):
        errors.append("QC review summary field 'evidence_chain' must be a list.")
    execution_trace = summary.get("execution_trace")
    if not isinstance(execution_trace, Mapping):
        errors.append("QC review summary field 'execution_trace' must be a mapping.")
    elif execution_trace.get("qc_schema_version") != QC_TRACE_SCHEMA_VERSION:
        errors.append("QC execution trace has an unsupported schema version.")
    output_health = summary.get("output_health")
    if not isinstance(output_health, Mapping):
        errors.append("QC review summary field 'output_health' must be a mapping.")
    readiness = summary.get("qc_readiness")
    if not isinstance(readiness, Mapping):
        errors.append("QC review summary field 'qc_readiness' must be a mapping.")
    actions = summary.get("review_action_items")
    if not isinstance(actions, list):
        errors.append("QC review summary field 'review_action_items' must be a list.")
    manifest = summary.get("reproducibility_manifest")
    if not isinstance(manifest, Mapping):
        errors.append("QC review summary field 'reproducibility_manifest' must be a mapping.")
    bundle = summary.get("evidence_bundle")
    if not isinstance(bundle, Mapping):
        errors.append("QC review summary field 'evidence_bundle' must be a mapping.")
    elif bundle.get("module") != "qc":
        errors.append("QC evidence_bundle.module must be 'qc'.")

    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def _is_tumor_context(tissue_type: Any) -> bool:
    if not tissue_type:
        return False
    tissue_text = str(tissue_type).lower()
    return "tumor" in tissue_text or "cancer" in tissue_text
