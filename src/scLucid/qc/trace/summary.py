"""QC review summary, maturity, and validation helpers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from anndata import AnnData

from scLucid.utils.contracts import _review_payload

from ...utils.context import is_tumor_context as _shared_is_tumor_context
from ...utils.evidence import (
    DecisionRecord,
    EvidenceBundle,
    EvidenceItem,
    ReviewAction,
    model_to_dict,
)
from . import (
    QC_MODULE_MATURITY_SCHEMA_VERSION,
    QC_REQUIRED_OBS_METRICS,
    QC_REQUIRED_REVIEW_SECTIONS,
    QC_TRACE_SCHEMA_VERSION,
    _json_safe,
)
from .review import (
    build_ambient_evidence_summary,
    build_doublet_evidence_summary,
    build_downstream_preprocess_recommendations,
    build_post_annotation_qc_review,
    build_qc_benchmark_scorecard,
    build_qc_decision_table,
    build_qc_evidence_chain,
    build_qc_execution_trace,
    build_qc_filtering_policy_summary,
    build_qc_handoff_readiness,
    build_qc_output_health,
    build_qc_policy_flow,
    build_qc_readiness_assessment,
    build_qc_recommended_threshold_summary,
    build_qc_reproducibility_manifest,
    build_qc_retention_audit_summary,
    build_qc_review_action_items,
    build_qc_reviewer_table,
    enrich_qc_decision_table_for_review,
    get_qc_module_contract,
)

def _evidence_source_for_stage(stage: Any) -> str:
    mapping = {
        "recommendation": "recommendation",
        "threshold_application": "metric",
        "sample_thresholds": "metric",
        "filtering": "metric",
        "filtering_policy": "contract",
        "retention_audit": "metric",
        "ambient_evidence": "metric",
        "post_annotation_qc_review": "context",
        "qc_benchmark_scorecard": "benchmark",
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

    handoff = summary.get("qc_handoff_readiness", {})
    if isinstance(handoff, Mapping):
        evidence_chain.append(
            EvidenceItem(
                source="contract",
                name="qc_handoff_readiness",
                value=handoff,
                rationale=(
                    "Declares the QC-to-preprocess handoff, including recommended counts layer, "
                    "cell decision columns, review/sensitivity cells, and downstream safety."
                ),
                limitations=list(str(item) for item in handoff.get("warnings", [])),
                related_keys=["qc_handoff_readiness", "downstream_preprocess_recommendations"],
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
            "qc_reviewer_table",
            "evidence_chain",
            "qc_handoff_readiness",
            "qc_readiness",
            "review_action_items",
            "reproducibility_manifest",
            "qc_filtering_policy_summary",
            "qc_retention_audit_summary",
            "ambient_evidence_summary",
            "post_annotation_qc_review",
            "qc_benchmark_scorecard",
            "benchmark_summary",
        ],
    )
    return model_to_dict(bundle)


def build_qc_module_maturity_assessment(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Assess whether a QC review summary satisfies the benchmark module contract."""
    payload = _review_payload(summary)
    required_sections = set(QC_REQUIRED_REVIEW_SECTIONS)
    # ``module_maturity`` is produced by this function, so it cannot be required
    # while the assessment is being built.
    required_sections.discard("module_maturity")
    missing_sections = sorted(required_sections - set(payload.keys()))
    decision_table = payload.get("decision_table")
    evidence_bundle = payload.get("evidence_bundle")
    readiness = payload.get("qc_readiness", {})
    handoff = payload.get("qc_handoff_readiness", {})
    output_health = payload.get("output_health", {})
    manifest = payload.get("reproducibility_manifest", {})

    issues: list[str] = []
    if missing_sections:
        issues.append("Missing required QC review sections: " + ", ".join(missing_sections))
    if not isinstance(decision_table, list) or not decision_table:
        issues.append("QC decision_table must be a non-empty list.")
    if not isinstance(evidence_bundle, Mapping) or evidence_bundle.get("module") != "qc":
        issues.append("QC evidence_bundle must be present and identify module='qc'.")
    if not isinstance(readiness, Mapping) or "status" not in readiness:
        issues.append("QC readiness assessment must be present.")
    if not isinstance(handoff, Mapping) or "status" not in handoff:
        issues.append("QC handoff readiness assessment must be present.")
    if not isinstance(output_health, Mapping) or "status" not in output_health:
        issues.append("QC output_health summary must be present.")
    if not isinstance(manifest, Mapping) or manifest.get("workflow") != "run_standard_qc":
        issues.append("QC reproducibility_manifest must identify workflow='run_standard_qc'.")

    review_required = []
    if isinstance(readiness, Mapping) and readiness.get("status") != "ready":
        review_required.append(f"qc_readiness.status={readiness.get('status')}")
    if isinstance(handoff, Mapping) and handoff.get("status") != "ready":
        review_required.append(f"qc_handoff_readiness.status={handoff.get('status')}")
    if isinstance(output_health, Mapping) and output_health.get("status") != "ok":
        review_required.append(f"output_health.status={output_health.get('status')}")

    if issues:
        status = "incomplete"
    elif review_required:
        status = "review_required"
    else:
        status = "complete"

    return _json_safe(
        {
            "schema_version": QC_MODULE_MATURITY_SCHEMA_VERSION,
            "module": "qc",
            "status": status,
            "status_scope": "review_contract_completeness_only",
            "scientific_validation_status": "REVIEW",
            "core_position": "withheld_until_locked_acceptance_passes",
            "superiority_claim": "unsupported_pending_blinded_validation",
            "issues": issues,
            "review_required": review_required,
            "contract": get_qc_module_contract(),
            "summary": (
                "QC review summary satisfies the benchmark module contract."
                if status == "complete"
                else "QC review summary is present but requires review."
                if status == "review_required"
                else "QC review summary does not satisfy the benchmark module contract."
            ),
        }
    )


def summarize_qc_review_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Return a compact user-facing summary of the QC review bundle."""
    payload = _review_payload(summary)
    readiness = payload.get("qc_readiness", {}) if isinstance(payload, Mapping) else {}
    handoff = payload.get("qc_handoff_readiness", {}) if isinstance(payload, Mapping) else {}
    filtering = payload.get("filtering_summary", {}) if isinstance(payload, Mapping) else {}
    recommendation = (
        payload.get("recommended_threshold_summary", {}) if isinstance(payload, Mapping) else {}
    )
    maturity = payload.get("module_maturity", {}) if isinstance(payload, Mapping) else {}
    decision_table = payload.get("decision_table", []) if isinstance(payload, Mapping) else []
    action_items = payload.get("review_action_items", []) if isinstance(payload, Mapping) else []
    doublet_summary = (
        payload.get("doublet_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    )
    doublet_predictions = (
        doublet_summary.get("predictions", {}) if isinstance(doublet_summary, Mapping) else {}
    )
    final_doublets = (
        doublet_predictions.get("predicted_doublet", {})
        if isinstance(doublet_predictions, Mapping)
        else {}
    )
    doublet_benchmark_decision = (
        doublet_summary.get("benchmark_decision", {})
        if isinstance(doublet_summary, Mapping)
        else {}
    )
    filtering_policy = (
        payload.get("qc_filtering_policy_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    retention_audit = (
        payload.get("qc_retention_audit_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    ambient_evidence = (
        payload.get("ambient_evidence_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    post_annotation = (
        payload.get("post_annotation_qc_review", {})
        if isinstance(payload, Mapping)
        else {}
    )
    scorecard = (
        payload.get("qc_benchmark_scorecard", {})
        if isinstance(payload, Mapping)
        else {}
    )
    benchmark = payload.get("benchmark_summary", {}) if isinstance(payload, Mapping) else {}
    benchmark_assessment = (
        benchmark.get("assessment", {}) if isinstance(benchmark, Mapping) else {}
    )
    benchmark_guide = (
        benchmark_assessment.get("interpretation_guide", {})
        if isinstance(benchmark_assessment, Mapping)
        else {}
    )

    applied_thresholds = {
        row.get("parameter"): row.get("applied")
        for row in decision_table
        if isinstance(row, Mapping) and row.get("parameter")
    }
    threshold_sources = {
        row.get("parameter"): row.get("source")
        for row in decision_table
        if isinstance(row, Mapping) and row.get("parameter")
    }

    return _json_safe(
        {
            "module": "qc",
            "maturity_status": maturity.get("status"),
            "readiness_status": readiness.get("status"),
            "readiness_score": readiness.get("score"),
            "qc_handoff_status": handoff.get("status") if isinstance(handoff, Mapping) else None,
            "ready_for_preprocess": (
                handoff.get("ready_for_preprocess") if isinstance(handoff, Mapping) else None
            ),
            "recommended_preprocess_counts_layer": (
                handoff.get("recommended_preprocess_counts_layer")
                if isinstance(handoff, Mapping)
                else None
            ),
            "review_required_cell_fraction": (
                handoff.get("handoff_cell_fractions", {}).get("review_required")
                if isinstance(handoff, Mapping)
                and isinstance(handoff.get("handoff_cell_fractions"), Mapping)
                else None
            ),
            "sensitivity_cell_fraction": (
                handoff.get("handoff_cell_fractions", {}).get("sensitivity_only")
                if isinstance(handoff, Mapping)
                and isinstance(handoff.get("handoff_cell_fractions"), Mapping)
                else None
            ),
            "verdict": readiness.get("verdict"),
            "recommendation_available": recommendation.get("available"),
            "overall_confidence": recommendation.get("overall_confidence"),
            "initial_cells": filtering.get("initial_cells"),
            "final_cells": filtering.get("final_cells"),
            "removed_fraction": filtering.get("removed_fraction"),
            "qc_decision_filter_mode": (
                filtering_policy.get("qc_decision_filter_mode")
                if isinstance(filtering_policy, Mapping)
                else None
            ),
            "final_filter_basis": (
                filtering_policy.get("final_filter_basis")
                if isinstance(filtering_policy, Mapping)
                else None
            ),
            "filtering_policy_review_required": (
                filtering_policy.get("review_required")
                if isinstance(filtering_policy, Mapping)
                else None
            ),
            "recommended_reviewer_first_path": (
                filtering_policy.get("recommended_reviewer_first_path")
                if isinstance(filtering_policy, Mapping)
                else None
            ),
            "retention_audit_available": (
                retention_audit.get("available")
                if isinstance(retention_audit, Mapping)
                else None
            ),
            "retention_audit_review_required": (
                retention_audit.get("retention_review_required")
                if isinstance(retention_audit, Mapping)
                else None
            ),
            "retention_audit_group_keys": (
                retention_audit.get("group_keys_reviewed", {})
                if isinstance(retention_audit, Mapping)
                else {}
            ),
            "ambient_status": (
                ambient_evidence.get("status")
                if isinstance(ambient_evidence, Mapping)
                else None
            ),
            "ambient_risk_level": (
                ambient_evidence.get("risk_level")
                if isinstance(ambient_evidence, Mapping)
                else None
            ),
            "ambient_review_required": (
                ambient_evidence.get("review_required")
                if isinstance(ambient_evidence, Mapping)
                else None
            ),
            "cell_probability_available": (
                ambient_evidence.get("cell_probability", {}).get("available")
                if isinstance(ambient_evidence, Mapping)
                and isinstance(ambient_evidence.get("cell_probability"), Mapping)
                else None
            ),
            "post_annotation_qc_available": (
                post_annotation.get("available")
                if isinstance(post_annotation, Mapping)
                else None
            ),
            "post_annotation_qc_review_required": (
                post_annotation.get("review_required")
                if isinstance(post_annotation, Mapping)
                else None
            ),
            "qc_benchmark_scorecard_status": (
                scorecard.get("status") if isinstance(scorecard, Mapping) else None
            ),
            "qc_benchmark_scorecard_score": (
                scorecard.get("score") if isinstance(scorecard, Mapping) else None
            ),
            "benchmark_status": benchmark.get("status") if isinstance(benchmark, Mapping) else None,
            "benchmark_risk_level": (
                benchmark_assessment.get("risk_level")
                if isinstance(benchmark_assessment, Mapping)
                else None
            ),
            "benchmark_summary": (
                benchmark_assessment.get("summary")
                if isinstance(benchmark_assessment, Mapping)
                else None
            ),
            "benchmark_main_risk": (
                benchmark_guide.get("main_risk")
                if isinstance(benchmark_guide, Mapping)
                else None
            ),
            "benchmark_next_step": (
                benchmark_guide.get("next_step")
                if isinstance(benchmark_guide, Mapping)
                else None
            ),
            "doublet_status": doublet_summary.get("status") if isinstance(doublet_summary, Mapping) else None,
            "predicted_doublets": final_doublets.get("count") if isinstance(final_doublets, Mapping) else None,
            "predicted_doublet_fraction": (
                final_doublets.get("fraction") if isinstance(final_doublets, Mapping) else None
            ),
            "doublet_recommended_default_mode": (
                doublet_benchmark_decision.get("recommended_default_mode")
                if isinstance(doublet_benchmark_decision, Mapping)
                else None
            ),
            "doublet_recommended_primary_method": (
                doublet_benchmark_decision.get("recommended_primary_method")
                if isinstance(doublet_benchmark_decision, Mapping)
                else None
            ),
            "doublet_recommended_algorithm_weight": (
                doublet_benchmark_decision.get("recommended_algorithm_weight")
                if isinstance(doublet_benchmark_decision, Mapping)
                else None
            ),
            "applied_thresholds": applied_thresholds,
            "threshold_sources": threshold_sources,
            "top_review_action": (
                action_items[0].get("action")
                if isinstance(action_items, list)
                and action_items
                and isinstance(action_items[0], Mapping)
                else None
            ),
            "n_review_action_items": len(action_items) if isinstance(action_items, list) else None,
        }
    )


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
    adata_before_filtering: AnnData | None = None,
) -> dict[str, Any]:
    """Add benchmark-grade QC audit fields to the review summary."""
    user_overrides = summary.get("user_override_summary", {}).get("details", {})
    decision_table = build_qc_decision_table(
        config,
        original_config,
        recommendation,
        user_overrides=user_overrides,
    )
    decision_table = enrich_qc_decision_table_for_review(
        decision_table,
        filtering_summary=filtering_summary,
        context=context,
    )
    output_health = build_qc_output_health(adata, filtering_summary)
    benchmark_summary = summary.get("benchmark_summary", {})
    doublet_evidence_summary = build_doublet_evidence_summary(adata)
    ambient_evidence_summary = build_ambient_evidence_summary(adata)
    post_annotation_review = build_post_annotation_qc_review(
        adata,
        sample_key=str(context.get("sample_key")) if context.get("sample_key") else None,
    )
    summary["qc_schema_version"] = QC_TRACE_SCHEMA_VERSION
    summary["decision_table"] = decision_table
    summary["threshold_reviewer_table"] = decision_table
    summary["doublet_evidence_summary"] = doublet_evidence_summary
    summary["ambient_evidence_summary"] = ambient_evidence_summary
    summary["post_annotation_qc_review"] = post_annotation_review
    qc_decision_summary = (
        adata.uns.get("sclucid", {}).get("qc", {}).get("qc_decision_summary", {})
    )
    if isinstance(qc_decision_summary, Mapping):
        summary["qc_decision_summary"] = dict(qc_decision_summary)
    filtering_policy_summary = build_qc_filtering_policy_summary(
        config=config,
        filtering_summary=filtering_summary,
        qc_decision_summary=summary.get("qc_decision_summary", {}),
    )
    summary["qc_filtering_policy_summary"] = filtering_policy_summary
    retention_audit_summary = build_qc_retention_audit_summary(
        adata_before_filtering=adata_before_filtering,
        adata_after_filtering=adata,
        context=context,
        filtering_summary=filtering_summary,
    )
    summary["qc_retention_audit_summary"] = retention_audit_summary
    qc_benchmark_scorecard = build_qc_benchmark_scorecard(
        benchmark_summary=benchmark_summary,
        doublet_evidence_summary=doublet_evidence_summary,
        ambient_evidence_summary=ambient_evidence_summary,
        retention_audit_summary=retention_audit_summary,
        post_annotation_qc_review=post_annotation_review,
    )
    summary["qc_benchmark_scorecard"] = qc_benchmark_scorecard
    summary["qc_reviewer_table"] = build_qc_reviewer_table(
        adata,
        decision_table=decision_table,
        qc_decision_summary=summary.get("qc_decision_summary", {}),
    )
    summary["policy_flow"] = build_qc_policy_flow(
        decision_table=decision_table,
        filtering_summary=filtering_summary,
        recommendation=recommendation,
        sample_thresholds=sample_thresholds,
    )
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
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
        ambient_evidence_summary=ambient_evidence_summary,
        post_annotation_qc_review=post_annotation_review,
        qc_benchmark_scorecard=qc_benchmark_scorecard,
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
    summary["qc_handoff_readiness"] = build_qc_handoff_readiness(
        adata=adata,
        context=context,
        filtering_summary=filtering_summary,
        output_health=output_health,
        downstream_recommendations=downstream_recommendations,
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
        doublet_evidence_summary=doublet_evidence_summary,
        ambient_evidence_summary=ambient_evidence_summary,
        tumor_aware_summary=summary.get("tumor_aware_summary", {}),
        post_annotation_qc_review=post_annotation_review,
    )
    readiness = build_qc_readiness_assessment(
        output_health=output_health,
        downstream_recommendations=downstream_recommendations,
        benchmark_summary=benchmark_summary,
        tumor_aware_summary=summary.get("tumor_aware_summary", {}),
        doublet_evidence_summary=doublet_evidence_summary,
        ambient_evidence_summary=ambient_evidence_summary,
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
        post_annotation_qc_review=post_annotation_review,
        qc_benchmark_scorecard=qc_benchmark_scorecard,
        warnings=warnings,
        decision_table=decision_table,
    )
    summary["qc_readiness"] = readiness
    summary["review_action_items"] = build_qc_review_action_items(
        readiness=readiness,
        downstream_recommendations=downstream_recommendations,
        tumor_aware_summary=summary.get("tumor_aware_summary", {}),
        benchmark_summary=benchmark_summary,
        doublet_evidence_summary=doublet_evidence_summary,
        ambient_evidence_summary=ambient_evidence_summary,
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
        post_annotation_qc_review=post_annotation_review,
        qc_benchmark_scorecard=qc_benchmark_scorecard,
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
    summary["module_maturity"] = build_qc_module_maturity_assessment(summary)
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
    if not isinstance(summary.get("decision_table"), (list, dict)):
        errors.append("QC review summary field 'decision_table' must be a list or dict.")
    else:
        decision_rows = summary.get("decision_table")
        if isinstance(decision_rows, Mapping):
            decision_rows = decision_rows.values()
        required_decision_columns = {
            "parameter",
            "recommended",
            "applied",
            "source",
            "confidence",
            "evidence",
            "review_required",
            "affected_cells",
            "biological_guardrail",
            "risk_note",
        }
        for idx, row in enumerate(decision_rows):
            if not isinstance(row, Mapping):
                errors.append(f"QC decision_table row {idx} must be a mapping.")
                continue
            missing_columns = sorted(required_decision_columns - set(row.keys()))
            if missing_columns:
                errors.append(f"QC decision_table row {idx} missing columns: {missing_columns}")
    if not isinstance(summary.get("evidence_chain"), (list, dict)):
        errors.append("QC review summary field 'evidence_chain' must be a list or dict.")
    execution_trace = summary.get("execution_trace")
    if not isinstance(execution_trace, Mapping):
        errors.append("QC review summary field 'execution_trace' must be a mapping.")
    elif execution_trace.get("qc_schema_version") != QC_TRACE_SCHEMA_VERSION:
        errors.append("QC execution trace has an unsupported schema version.")
    output_health = summary.get("output_health")
    if not isinstance(output_health, Mapping):
        errors.append("QC review summary field 'output_health' must be a mapping.")
    doublet_summary = summary.get("doublet_evidence_summary")
    if not isinstance(doublet_summary, Mapping):
        errors.append("QC review summary field 'doublet_evidence_summary' must be a mapping.")
    readiness = summary.get("qc_readiness")
    if not isinstance(readiness, Mapping):
        errors.append("QC review summary field 'qc_readiness' must be a mapping.")
    handoff = summary.get("qc_handoff_readiness")
    if not isinstance(handoff, Mapping):
        errors.append("QC review summary field 'qc_handoff_readiness' must be a mapping.")
    else:
        required_handoff_fields = {
            "status",
            "ready_for_preprocess",
            "recommended_preprocess_counts_layer",
            "counts_layer_present",
            "filtering_decision_columns",
            "handoff_cell_counts",
            "handoff_cell_fractions",
            "safe_to_continue",
            "required_downstream_handling",
        }
        missing_handoff = sorted(required_handoff_fields - set(handoff.keys()))
        if missing_handoff:
            errors.append(f"QC qc_handoff_readiness missing fields: {missing_handoff}")
    actions = summary.get("review_action_items")
    if not isinstance(actions, (list, dict)):
        errors.append("QC review summary field 'review_action_items' must be a list or dict.")
    manifest = summary.get("reproducibility_manifest")
    if not isinstance(manifest, Mapping):
        errors.append("QC review summary field 'reproducibility_manifest' must be a mapping.")
    bundle = summary.get("evidence_bundle")
    if not isinstance(bundle, Mapping):
        errors.append("QC review summary field 'evidence_bundle' must be a mapping.")
    elif bundle.get("module") != "qc":
        errors.append("QC evidence_bundle.module must be 'qc'.")
    maturity = summary.get("module_maturity")
    if not isinstance(maturity, Mapping):
        errors.append("QC review summary field 'module_maturity' must be a mapping.")
    elif maturity.get("module") != "qc":
        errors.append("QC module_maturity.module must be 'qc'.")

    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def validate_qc_module_completeness(
    adata: AnnData,
    *,
    require_ready: bool = False,
    raise_on_error: bool = False,
) -> dict[str, Any]:
    """Validate that an AnnData object contains a benchmark-grade QC result."""
    issues: list[str] = []
    warnings: list[str] = []

    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    if not isinstance(qc_ns, Mapping):
        issues.append('Missing or invalid adata.uns["sclucid"]["qc"] namespace.')
        qc_ns = {}

    review_summary = qc_ns.get("review_summary")
    payload = _review_payload(review_summary) if isinstance(review_summary, Mapping) else {}
    if not payload:
        issues.append('Missing adata.uns["sclucid"]["qc"]["review_summary"].')
        maturity = build_qc_module_maturity_assessment({})
    else:
        validation_errors = validate_qc_review_summary(payload)
        issues.extend(validation_errors)
        maturity = build_qc_module_maturity_assessment(payload)
        if maturity.get("status") == "incomplete":
            issues.extend(maturity.get("issues", []))
        elif maturity.get("status") == "review_required":
            warnings.extend(maturity.get("review_required", []))

    missing_metrics = [metric for metric in QC_REQUIRED_OBS_METRICS if metric not in adata.obs]
    if missing_metrics:
        issues.append("Missing required QC obs metrics: " + ", ".join(missing_metrics))

    readiness = payload.get("qc_readiness", {}) if isinstance(payload, Mapping) else {}
    if require_ready and readiness.get("status") != "ready":
        issues.append(f"QC readiness is {readiness.get('status')!r}, expected 'ready'.")

    result = {
        "schema_version": QC_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "qc",
        "valid": len(issues) == 0,
        "status": "valid" if not issues else "invalid",
        "issues": list(dict.fromkeys(str(item) for item in issues)),
        "warnings": list(dict.fromkeys(str(item) for item in warnings)),
        "maturity": maturity,
        "summary": summarize_qc_review_summary(payload) if payload else {},
    }

    if result["issues"] and raise_on_error:
        raise ValueError("; ".join(result["issues"]))
    return _json_safe(result)


def _is_tumor_context(tissue_type: Any) -> bool:
    return _shared_is_tumor_context(tissue_type)
