"""QC trace review builders and validation helpers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd
from anndata import AnnData

from ...utils.context import is_tumor_context as _shared_is_tumor_context
from . import (
    QC_EXPECTED_ARTIFACTS,
    QC_MODULE_MATURITY_SCHEMA_VERSION,
    QC_REQUIRED_OBS_METRICS,
    QC_REQUIRED_REVIEW_SECTIONS,
    QC_STABLE_ENTRYPOINTS,
    QC_TRACE_SCHEMA_VERSION,
    _PARAMETER_SPECS,
    _decision_source,
    _get_nested,
    _json_safe,
    _recommendation_value,
    _to_dict,
)


def _is_tumor_context(tissue_type: Any) -> bool:
    return _shared_is_tumor_context(tissue_type)


def get_qc_module_contract() -> dict[str, Any]:
    """Return the frozen QC module maturity contract."""
    return {
        "schema_version": QC_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "qc",
        "stable_entrypoints": list(QC_STABLE_ENTRYPOINTS),
        "required_obs_metrics": list(QC_REQUIRED_OBS_METRICS),
        "required_review_sections": sorted(QC_REQUIRED_REVIEW_SECTIONS),
        "expected_sidecar_artifacts": list(QC_EXPECTED_ARTIFACTS),
        "canonical_namespace": 'adata.uns["sclucid"]["qc"]',
        "readiness_key": "qc_readiness",
        "handoff_key": "qc_handoff_readiness",
        "decision_table_key": "decision_table",
        "evidence_bundle_key": "evidence_bundle",
    }


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
            "metric": spec["obs_metric"],
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


def enrich_qc_decision_table_for_review(
    decision_table: list[Mapping[str, Any]],
    *,
    filtering_summary: Mapping[str, Any],
    context: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Add reviewer-facing impact and risk columns to QC threshold rows."""
    criteria_counts = filtering_summary.get("criteria_counts", {})
    if not isinstance(criteria_counts, Mapping):
        criteria_counts = {}
    review_counts = filtering_summary.get("review_criteria_counts", {})
    if not isinstance(review_counts, Mapping):
        review_counts = {}
    criteria_used = set(filtering_summary.get("criteria_used", []) or [])
    is_tumor = _is_tumor_context(context.get("tissue_type"))
    enriched: list[dict[str, Any]] = []
    for row in decision_table:
        item = dict(row)
        flag = item.get("filtering_flag")
        affected_cells = int(
            criteria_counts.get(flag, review_counts.get(flag, 0)) or 0
        ) if flag else 0
        filtering_enabled = bool(item.get("is_filtering_enabled") or flag in criteria_used)
        source = item.get("source")
        parameter = item.get("parameter")
        risk_note = ""
        biological_guardrail = ""
        review_required = False
        if source == "user_override":
            review_required = True
            risk_note = "User override applied; document why the configured value differs from the recommendation."
        elif is_tumor and parameter == "max_mt_percent":
            review_required = True
            biological_guardrail = (
                "Preserve high-mt malignant/stress/program signal until marker/program evidence is reviewed."
            )
            risk_note = (
                "Tumor context: high mitochondrial fraction can mark stress or malignant biology; "
                "review marker/program retention before hard deletion."
            )
        elif filtering_enabled and affected_cells > 0 and item.get("recommended") is None:
            review_required = True
            risk_note = "Filtering is enabled without an available recommendation; verify the configured threshold."
        elif not filtering_enabled and item.get("applied") is not None:
            risk_note = "Threshold is recorded but its filtering flag is not active."

        item["affected_cells"] = affected_cells
        item["affected_fraction"] = (
            float(item["affected_cells"] / max(float(filtering_summary.get("initial_cells") or 0), 1.0))
            if filtering_summary.get("initial_cells") is not None
            else None
        )
        item["review_required"] = review_required
        item["biological_guardrail"] = biological_guardrail
        item["risk_note"] = risk_note
        item["evidence"] = item.get("evidence") or {}
        enriched.append(_json_safe(item))
    return enriched


def build_qc_reviewer_table(
    adata: AnnData,
    *,
    decision_table: list[Mapping[str, Any]],
    qc_decision_summary: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build one reviewer-facing QC table across thresholds and evidence flags."""
    rows: list[dict[str, Any]] = []
    for row in decision_table:
        rows.append(
            _json_safe(
                {
                    "item": row.get("parameter"),
                    "category": "threshold",
                    "metric": row.get("metric") or row.get("obs_metric"),
                    "recommended_value": row.get("recommended"),
                    "applied_value": row.get("applied"),
                    "source": row.get("source"),
                    "confidence": row.get("confidence"),
                    "affected_cells": int(row.get("affected_cells") or 0),
                    "affected_fraction": row.get("affected_fraction"),
                    "biological_risk_note": row.get("risk_note")
                    or row.get("biological_guardrail")
                    or "",
                    "review_required": bool(row.get("review_required")),
                    "decision_column": row.get("filtering_flag"),
                }
            )
        )

    evidence_summary = {}
    if isinstance(qc_decision_summary, Mapping):
        evidence_summary = qc_decision_summary.get("evidence_summary", {}) or {}
    if not isinstance(evidence_summary, Mapping):
        evidence_summary = {}

    evidence_specs = [
        (
            "ambient_risk",
            "contamination",
            "ambient_fraction/ambient_score",
            "review_high_ambient_signal",
            "Ambient RNA evidence should be reviewed before downstream interpretation; prefer corrected counts when available.",
        ),
        (
            "hemoglobin_contamination",
            "contamination",
            "pct_counts_hb/hemoglobin_score",
            "review_high_hemoglobin_signal",
            "Hemoglobin signal can indicate RBC contamination or fragile sample composition; avoid deleting rare biology without marker review.",
        ),
        (
            "platelet_contamination",
            "contamination",
            "platelet_score",
            "review_high_platelet_signal",
            "Platelet markers can reflect contamination or platelet-bearing immune interactions; inspect clusters and samples.",
        ),
        (
            "stress_high",
            "stress",
            "stress_score",
            "mark_for_review_or_sensitivity",
            "Dissociation/stress programs can be technical or biological; prefer sensitivity analysis after annotation.",
        ),
        (
            "apoptosis_high",
            "stress",
            "apoptosis_score",
            "joint_evidence_required_for_removal",
            "Apoptosis-like signal supports low-quality removal only when combined with other QC failures.",
        ),
        (
            "predicted_doublet",
            "doublet",
            "combined_doublet_score/predicted_doublet",
            "sample_aware_doublet_review",
            "Doublet calls should be reviewed around cell-type boundaries and with sample-level expected rates.",
        ),
        (
            "qc_remove",
            "decision",
            "qc_decision",
            "remove_when_decision_remove",
            "Final removal decision from the unified QC decision engine.",
        ),
    ]
    confidence = (
        pd.to_numeric(adata.obs["qc_confidence"], errors="coerce")
        if "qc_confidence" in adata.obs
        else pd.Series(dtype=float)
    )
    for column, category, metric, recommended, risk_note in evidence_specs:
        if column in adata.obs:
            mask = adata.obs[column].fillna(False).astype(bool)
            affected = int(mask.sum())
        else:
            affected = int(evidence_summary.get(column, 0) or 0)
            mask = pd.Series(False, index=adata.obs_names)
        mean_confidence = None
        if affected and not confidence.empty:
            mean_confidence = float(confidence.loc[mask].mean())
        rows.append(
            _json_safe(
                {
                    "item": column,
                    "category": category,
                    "metric": metric,
                    "recommended_value": recommended,
                    "applied_value": "flagged" if affected else "not_flagged",
                    "source": "qc_decision_engine",
                    "confidence": mean_confidence,
                    "affected_cells": affected,
                    "affected_fraction": affected / max(float(adata.n_obs), 1.0),
                    "biological_risk_note": risk_note,
                    "review_required": bool(
                        affected and category in {"contamination", "stress", "doublet"}
                    ),
                    "decision_column": column,
                }
            )
        )
    return rows


def build_qc_policy_flow(
    *,
    decision_table: list[Mapping[str, Any]],
    filtering_summary: Mapping[str, Any],
    recommendation: Any,
    sample_thresholds: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Describe the canonical QC decision flow in reviewable stages."""
    rec_dict = _to_dict(recommendation)
    filtering_applied = filtering_summary.get("initial_cells") is not None
    return _json_safe(
        [
            {
                "stage": "profile_dataset",
                "status": "complete",
                "evidence": "QC metrics and dataset context were collected.",
            },
            {
                "stage": "propose_candidate_thresholds",
                "status": "complete" if rec_dict or sample_thresholds else "not_available",
                "evidence": "Intelligent and/or sample-aware threshold candidates are available."
                if rec_dict or sample_thresholds
                else "No intelligent or sample-aware threshold candidates were generated.",
            },
            {
                "stage": "score_biological_risk",
                "status": "review_required"
                if any(row.get("biological_guardrail") for row in decision_table)
                else "complete",
                "evidence": "Tumor-aware biological guardrails are present."
                if any(row.get("biological_guardrail") for row in decision_table)
                else "No tumor-specific biological guardrail was triggered.",
            },
            {
                "stage": "choose_recommend_policy",
                "status": "complete",
                "evidence": "Applied thresholds and their sources are recorded in decision_table.",
            },
            {
                "stage": "emit_reviewer_table",
                "status": "complete" if decision_table else "missing",
                "evidence": f"{len(decision_table)} threshold decision row(s) emitted.",
            },
            {
                "stage": "optionally_apply",
                "status": "complete" if filtering_applied else "not_applied",
                "evidence": "Filtering summary is available."
                if filtering_applied
                else "Policy was recommended without applying filters.",
            },
        ]
    )


def build_qc_filtering_policy_summary(
    *,
    config: Any,
    filtering_summary: Mapping[str, Any],
    qc_decision_summary: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Summarize whether final filtering used legacy thresholds or qc_decision."""
    mode = str(getattr(config, "qc_decision_filter_mode", "off"))
    decision_engine_enabled = bool(getattr(config, "run_decision_engine", False))
    filtering_applied = filtering_summary.get("initial_cells") is not None
    criteria = list(filtering_summary.get("criteria_used", []) or [])
    decision_counts = {}
    if isinstance(qc_decision_summary, Mapping):
        decision_counts = qc_decision_summary.get("decision_counts", {}) or {}

    if mode == "replace":
        basis = "qc_decision_remove"
        recommendation = "reviewer_first_decision_filtering"
        review_required = False
        risk_note = "Final filtering used cells marked qc_decision == 'remove'."
    elif mode == "append":
        basis = "legacy_thresholds_plus_qc_remove"
        recommendation = "hybrid_filtering_review"
        review_required = True
        risk_note = (
            "Filtering appended qc_remove to legacy criteria; audit retention by sample, "
            "condition, and annotation before downstream interpretation."
        )
    else:
        basis = "legacy_threshold_filtering"
        recommendation = "prefer_run_iterative_qc_or_qc_decision_replace_for_reviewer_first_qc"
        review_required = bool(decision_engine_enabled and filtering_applied)
        risk_note = (
            "QC decision labels were built but final filtering used legacy filter_config "
            "criteria. This is supported for compatibility, but reviewer-first workflows "
            "should prefer run_iterative_qc(..., final_filter_policy='decision_remove') "
            "or qc_decision_filter_mode='replace'."
        )

    return _json_safe(
        {
            "schema_version": "qc_filtering_policy_summary_v1",
            "qc_decision_filter_mode": mode,
            "decision_engine_enabled": decision_engine_enabled,
            "filtering_applied": filtering_applied,
            "final_filter_basis": basis,
            "criteria_used": criteria,
            "decision_counts": decision_counts,
            "recommended_reviewer_first_path": recommendation,
            "review_required": review_required,
            "risk_note": risk_note,
        }
    )


def _first_existing_obs_key(adata: AnnData, candidates: tuple[str, ...]) -> str | None:
    for key in candidates:
        if key in adata.obs:
            return key
    return None


def _build_retention_table(
    before: AnnData,
    retained_mask: pd.Series,
    *,
    key: str,
    label: str,
    min_group_size: int = 5,
) -> list[dict[str, Any]]:
    if key not in before.obs:
        return []
    rows: list[dict[str, Any]] = []
    group_values = before.obs[key].astype(str).fillna("NA")
    for value, idx in group_values.groupby(group_values, observed=False).groups.items():
        idx_list = list(idx)
        before_n = len(idx_list)
        if before_n < min_group_size:
            continue
        after_n = int(retained_mask.loc[idx_list].sum())
        removed = before_n - after_n
        retention_rate = float(after_n / before_n) if before_n else None
        rows.append(
            {
                "scope": label,
                "key": key,
                "group": str(value),
                "before": int(before_n),
                "after": int(after_n),
                "removed": int(removed),
                "retention_rate": retention_rate,
                "removed_fraction": float(removed / before_n) if before_n else None,
                "review_flag": (
                    "low_retention"
                    if retention_rate is not None and retention_rate < 0.5
                    else "high_loss"
                    if retention_rate is not None and retention_rate < 0.75
                    else "ok"
                ),
            }
        )
    rows.sort(key=lambda row: (row["retention_rate"] is None, row["retention_rate"] or 0.0))
    return _json_safe(rows)


def build_qc_retention_audit_summary(
    *,
    adata_before_filtering: AnnData | None,
    adata_after_filtering: AnnData,
    context: Mapping[str, Any],
    filtering_summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit filtering retention by sample, condition, annotation, and cluster."""
    if adata_before_filtering is None:
        return {
            "schema_version": "qc_retention_audit_summary_v1",
            "available": False,
            "filtering_applied": filtering_summary.get("initial_cells") is not None,
            "retention_review_required": False,
            "reason": "Pre-filter AnnData was not provided; stratified retention audit is unavailable.",
            "tables": {},
            "flags": [],
        }

    before = adata_before_filtering
    retained_index = before.obs_names.isin(adata_after_filtering.obs_names)
    retained_mask = pd.Series(retained_index, index=before.obs_names)
    sample_key = str(context.get("sample_key")) if context.get("sample_key") else None
    condition_key = _first_existing_obs_key(
        before,
        (
            "condition",
            "group",
            "treatment",
            "disease",
            "status",
            "phenotype",
            "timepoint",
        ),
    )
    annotation_key = _first_existing_obs_key(
        before,
        (
            "cell_type",
            "cell_type_auto",
            "celltype",
            "celltype_lineage",
            "celltype_lineage_auto",
            "annotation",
            "predicted_cell_type",
        ),
    )
    cluster_key = _first_existing_obs_key(
        before,
        (
            "leiden_clusters",
            "leiden",
            "louvain",
            "cluster",
            "seurat_clusters",
        ),
    )
    key_specs = [
        ("sample", sample_key if sample_key in before.obs else None),
        ("condition", condition_key),
        ("annotation", annotation_key),
        ("cluster", cluster_key),
    ]
    tables: dict[str, list[dict[str, Any]]] = {}
    flags: list[dict[str, Any]] = []
    for label, key in key_specs:
        if not key:
            continue
        table = _build_retention_table(before, retained_mask, key=key, label=label)
        tables[label] = table
        if not table:
            continue
        rates = [row["retention_rate"] for row in table if row.get("retention_rate") is not None]
        if rates and min(rates) < 0.5:
            flags.append(
                {
                    "scope": label,
                    "flag": "low_group_retention",
                    "key": key,
                    "min_retention_rate": min(rates),
                    "review_required": True,
                }
            )
        if rates and max(rates) - min(rates) > 0.35:
            flags.append(
                {
                    "scope": label,
                    "flag": "retention_imbalance",
                    "key": key,
                    "retention_spread": max(rates) - min(rates),
                    "review_required": True,
                }
            )

    filtering_applied = filtering_summary.get("initial_cells") is not None
    mandatory_review = bool(filtering_applied)
    review_required = mandatory_review or any(flag.get("review_required") for flag in flags)
    return _json_safe(
        {
            "schema_version": "qc_retention_audit_summary_v1",
            "available": True,
            "filtering_applied": filtering_applied,
            "retention_review_required": review_required,
            "mandatory_review": mandatory_review,
            "mandatory_review_note": (
                "Review retention by sample, condition, annotation, and cluster before "
                "treating filtered data as final."
            )
            if mandatory_review
            else "Filtering was not applied in this run.",
            "overall": {
                "initial_cells": int(before.n_obs),
                "final_cells": int(adata_after_filtering.n_obs),
                "removed_cells": int(before.n_obs - adata_after_filtering.n_obs),
                "retention_rate": float(adata_after_filtering.n_obs / max(before.n_obs, 1)),
            },
            "group_keys_reviewed": {
                label: key for label, key in key_specs if key is not None
            },
            "tables": tables,
            "flags": flags,
        }
    )


def build_doublet_evidence_summary(adata: AnnData) -> dict[str, Any]:
    """Summarize doublet evidence for the QC review contract.

    This is intentionally report-only: it reads predictions, scores, and stored
    method metadata but never changes doublet calls or filtering decisions.
    """
    prediction_cols = [
        "predicted_doublet",
        "algorithm_predicted_doublet",
        "scrublet_predicted",
        "scanpy_scrublet_predicted",
        "scdblfinder_predicted",
        "doubletdetection_predicted",
        "heuristic_predicted",
        "external_doublet_evidence",
    ]
    score_cols = [
        "combined_doublet_score",
        "algorithm_doublet_score",
        "doublet_score",
        "scrublet_score",
        "scanpy_scrublet_score",
        "scdblfinder_score",
        "heuristic_confidence_score",
        "heterotypic_doublet_risk",
        "homotypic_doublet_risk",
    ]

    def _benchmark_decision_summary(evidence: Mapping[str, Any]) -> dict[str, Any]:
        recommendations = evidence.get("algorithm_weight_recommendations", [])
        if isinstance(recommendations, Mapping):
            recommendations = list(recommendations.values())
        if not isinstance(recommendations, list):
            recommendations = []

        normalized_recommendations: list[dict[str, Any]] = []
        preferred = None
        for item in recommendations:
            if not isinstance(item, Mapping):
                continue
            record = {
                "base_method": item.get("base_method"),
                "recommended_default_mode": item.get("recommended_default_mode"),
                "recommended_algorithm_weight": item.get("recommended_algorithm_weight"),
                "recommended_method": item.get("recommended_method"),
                "f1_delta_vs_algorithm_only": item.get("f1_delta_vs_algorithm_only"),
                "recall_delta_vs_algorithm_only": item.get("recall_delta_vs_algorithm_only"),
                "precision_delta_vs_algorithm_only": item.get(
                    "precision_delta_vs_algorithm_only"
                ),
                "review_required": bool(item.get("review_required", False)),
                "risk_note": item.get("risk_note", ""),
            }
            normalized_recommendations.append(record)
            mode = str(record.get("recommended_default_mode") or "")
            if preferred is None and mode.startswith("algorithm_only"):
                preferred = record
            elif preferred is None and record.get("base_method"):
                preferred = record

        parity = evidence.get("python_r_parity", {})
        parity_review_required = (
            bool(parity.get("review_required", False)) if isinstance(parity, Mapping) else False
        )
        threshold_reviews = evidence.get("threshold_calibration_review", [])
        if isinstance(threshold_reviews, Mapping):
            threshold_reviews = list(threshold_reviews.values())
        threshold_review_count = len(threshold_reviews) if isinstance(threshold_reviews, list) else 0

        review_required = bool(
            parity_review_required
            or threshold_review_count > 0
            or any(row.get("review_required") for row in normalized_recommendations)
        )
        decision = {
            "schema_version": evidence.get("schema_version"),
            "best_method": evidence.get("best_method"),
            "best_method_f1": evidence.get("best_method_f1"),
            "best_method_auc": evidence.get("best_method_auc"),
            "recommended_default_mode": (
                preferred.get("recommended_default_mode") if preferred else None
            ),
            "recommended_primary_method": preferred.get("base_method") if preferred else None,
            "recommended_algorithm_weight": (
                preferred.get("recommended_algorithm_weight") if preferred else None
            ),
            "recommended_fusion_method": (
                preferred.get("recommended_method") if preferred else None
            ),
            "algorithm_weight_recommendations": normalized_recommendations,
            "threshold_calibration_review_count": threshold_review_count,
            "python_r_parity_review_required": parity_review_required,
            "review_required": review_required,
        }
        if preferred and preferred.get("risk_note"):
            decision["risk_note"] = preferred.get("risk_note")
        return _json_safe(decision)

    predictions: dict[str, dict[str, Any]] = {}
    for col in prediction_cols:
        if col not in adata.obs:
            continue
        values = adata.obs[col]
        if values.dtype == bool or set(pd.Series(values).dropna().unique()).issubset({0, 1}):
            bool_values = values.astype(bool)
            count = int(bool_values.sum())
            predictions[col] = {
                "count": count,
                "fraction": float(count / max(adata.n_obs, 1)),
                "percent": float(count / max(adata.n_obs, 1) * 100.0),
            }

    scores: dict[str, dict[str, Any]] = {}
    for col in score_cols:
        if col not in adata.obs:
            continue
        vals = pd.to_numeric(adata.obs[col], errors="coerce").dropna()
        if vals.empty:
            continue
        scores[col] = {
            "median": float(vals.median()),
            "p90": float(vals.quantile(0.90)),
            "p95": float(vals.quantile(0.95)),
            "max": float(vals.max()),
        }

    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    doublet_params = qc_ns.get("doublet_params", {}) if isinstance(qc_ns, Mapping) else {}
    doublet_params = _json_safe(doublet_params) if isinstance(doublet_params, Mapping) else {}
    benchmark_evidence = (
        qc_ns.get("doublet_benchmark_evidence", {}) if isinstance(qc_ns, Mapping) else {}
    )
    benchmark_evidence = (
        _json_safe(benchmark_evidence) if isinstance(benchmark_evidence, Mapping) else {}
    )
    benchmark_decision = (
        _benchmark_decision_summary(benchmark_evidence) if benchmark_evidence else {}
    )
    risk_decomposition = (
        doublet_params.get("risk_decomposition", {}) if isinstance(doublet_params, Mapping) else {}
    )
    external_evidence = (
        doublet_params.get("external_doublet_evidence", {})
        if isinstance(doublet_params, Mapping)
        else {}
    )

    final = predictions.get("predicted_doublet", {})
    predicted_fraction = final.get("fraction")
    review_required = False
    notes: list[str] = []
    if not predictions and not scores:
        status = "not_run"
        notes.append("No doublet prediction or score columns were found.")
    else:
        status = "available"
        if predicted_fraction is not None and predicted_fraction > 0.20:
            review_required = True
            notes.append("Predicted doublet fraction exceeds 20%; inspect sample-level rates.")
        if external_evidence:
            notes.append("External doublet evidence is present; compare overlap with algorithmic calls.")
        if risk_decomposition:
            notes.append("Heterotypic/homotypic risk decomposition is available.")
        if benchmark_evidence:
            notes.append("External doublet benchmark evidence is attached to the QC report.")
        if benchmark_decision.get("recommended_default_mode"):
            notes.append(
                "Benchmark-supported doublet mode: "
                f"{benchmark_decision['recommended_default_mode']}."
            )
        if benchmark_decision.get("review_required"):
            review_required = True
            notes.append("Doublet benchmark evidence indicates review is required.")

    return _json_safe(
        {
            "status": status,
            "n_cells": int(adata.n_obs),
            "predictions": predictions,
            "scores": scores,
            "risk_decomposition": risk_decomposition,
            "external_evidence": external_evidence,
            "benchmark_evidence": benchmark_evidence,
            "benchmark_decision": benchmark_decision,
            "method_metadata_keys": sorted(doublet_params.keys()) if isinstance(doublet_params, Mapping) else [],
            "review_required": review_required,
            "notes": notes,
        }
    )


def _numeric_obs_summary(adata: AnnData, column: str) -> dict[str, Any]:
    if column not in adata.obs:
        return {"available": False}
    values = pd.to_numeric(adata.obs[column], errors="coerce").dropna()
    if values.empty:
        return {"available": False}
    return {
        "available": True,
        "median": float(values.median()),
        "p90": float(values.quantile(0.90)),
        "p95": float(values.quantile(0.95)),
        "max": float(values.max()),
        "n_measured": int(values.shape[0]),
    }


def build_ambient_evidence_summary(adata: AnnData) -> dict[str, Any]:
    """Summarize ambient RNA evidence, correction status, and canonical fields."""
    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    ambient_summary = qc_ns.get("ambient_rna_summary", {}) if isinstance(qc_ns, Mapping) else {}
    empty_summary = qc_ns.get("empty_droplet_summary", {}) if isinstance(qc_ns, Mapping) else {}
    correction_summary = (
        qc_ns.get("ambient_correction_summary", {}) if isinstance(qc_ns, Mapping) else {}
    )
    correction_status = (
        qc_ns.get("ambient_correction_status", {}) if isinstance(qc_ns, Mapping) else {}
    )
    layer_contract = (
        qc_ns.get("ambient_layer_contract", {}) if isinstance(qc_ns, Mapping) else {}
    )
    probability_schema = (
        qc_ns.get("qc_probability_schema", {}) if isinstance(qc_ns, Mapping) else {}
    )

    ambient_risk_count = (
        int(adata.obs["ambient_risk"].fillna(False).astype(bool).sum())
        if "ambient_risk" in adata.obs
        else 0
    )
    corrected = bool(
        correction_summary.get("corrected", False)
        or correction_status.get("corrected", False)
        or layer_contract.get("corrected_layer_present", False)
    )
    risk_level = ambient_summary.get("risk_level", "unknown") if isinstance(ambient_summary, Mapping) else "unknown"
    review_required = bool(
        risk_level in {"moderate", "high", "unknown"}
        or ambient_risk_count > 0
        or layer_contract.get("review_required", False)
    )
    if corrected and risk_level == "low" and ambient_risk_count == 0:
        review_required = False

    notes: list[str] = []
    if risk_level in {"moderate", "high"}:
        notes.append(f"Ambient RNA diagnostic risk is {risk_level}.")
    if not corrected:
        notes.append("Ambient correction was not applied or not registered.")
    if layer_contract:
        recommended = layer_contract.get("recommended_preprocess_counts_layer")
        if recommended:
            notes.append(f"Recommended preprocessing counts layer: {recommended}.")
    if empty_summary and not empty_summary.get("available", False):
        notes.append(
            "Empty-droplet diagnostic is unavailable or diagnostic-only; cell_probability may be missing."
        )

    return _json_safe(
        {
            "schema_version": "ambient_evidence_summary_v1",
            "status": "available" if ambient_summary or layer_contract else "not_run",
            "risk_level": risk_level,
            "risk_score": ambient_summary.get("risk_score") if isinstance(ambient_summary, Mapping) else None,
            "ambient_risk_cells": ambient_risk_count,
            "ambient_fraction": _numeric_obs_summary(adata, "ambient_fraction"),
            "ambient_score": _numeric_obs_summary(adata, "ambient_score"),
            "cell_probability": _numeric_obs_summary(adata, "cell_probability"),
            "empty_droplet_probability": _numeric_obs_summary(
                adata, "empty_droplet_probability"
            ),
            "correction": correction_summary or correction_status,
            "empty_droplet_summary": empty_summary,
            "layer_contract": layer_contract,
            "probability_schema": probability_schema,
            "review_required": review_required,
            "notes": notes,
        }
    )


def build_post_annotation_qc_review(
    adata: AnnData,
    *,
    sample_key: str | None = None,
    cell_type_key: str | None = None,
    min_cells: int = 10,
) -> dict[str, Any]:
    """Review QC risk after annotation when cell-type labels are available.

    This does not filter cells.  It identifies whether retained review signals
    are sample-driven, cell-type-specific, or broadly concordant, so downstream
    analysis can decide between keep, remove, sensitivity analysis, or explicit
    modeling.
    """
    if cell_type_key is None:
        cell_type_key = _first_existing_obs_key(
            adata,
            (
                "cell_type",
                "cell_type_auto",
                "celltype",
                "celltype_lineage",
                "celltype_lineage_auto",
                "annotation",
                "predicted_cell_type",
            ),
        )
    if sample_key is None:
        sample_key = _first_existing_obs_key(
            adata,
            ("sampleID", "sample", "Sample", "orig.ident", "batch", "Batch"),
        )

    review_columns = [
        "stress_high",
        "ambient_risk",
        "predicted_doublet",
        "qc_high_mt",
        "apoptosis_high",
        "hemoglobin_contamination",
        "platelet_contamination",
    ]
    available = [col for col in review_columns if col in adata.obs]
    if not available:
        return {
            "schema_version": "post_annotation_qc_review_v1",
            "available": False,
            "reason": "No QC review columns are available.",
            "review_required": False,
            "recommended_actions": [],
            "tables": {},
        }

    tables: dict[str, list[dict[str, Any]]] = {}
    recommended_actions: list[dict[str, Any]] = []
    for column in available:
        mask = adata.obs[column].fillna(False).astype(bool)
        total = int(mask.sum())
        if total == 0:
            continue
        action = "sensitivity_only" if column in {"stress_high", "qc_high_mt", "ambient_risk"} else "review"
        if column == "predicted_doublet":
            action = "review_boundary_or_remove_high_confidence"
        recommended_actions.append(
            {
                "signal": column,
                "affected_cells": total,
                "affected_fraction": float(total / max(adata.n_obs, 1)),
                "recommended_action": action,
                "rationale": (
                    "Review after annotation before irreversible exclusion; this signal can be technical or biological."
                ),
            }
        )

        if cell_type_key and cell_type_key in adata.obs:
            rows = []
            for value, group in adata.obs.groupby(cell_type_key, observed=True):
                n = int(group.shape[0])
                if n < min_cells:
                    continue
                affected = int(group[column].fillna(False).astype(bool).sum())
                rows.append(
                    {
                        "signal": column,
                        "scope": "cell_type",
                        "key": cell_type_key,
                        "group": str(value),
                        "n_cells": n,
                        "affected_cells": affected,
                        "affected_fraction": float(affected / max(n, 1)),
                    }
                )
            if rows:
                tables.setdefault("cell_type", []).extend(rows)

        if sample_key and sample_key in adata.obs:
            rows = []
            for value, group in adata.obs.groupby(sample_key, observed=True):
                n = int(group.shape[0])
                if n < min_cells:
                    continue
                affected = int(group[column].fillna(False).astype(bool).sum())
                rows.append(
                    {
                        "signal": column,
                        "scope": "sample",
                        "key": sample_key,
                        "group": str(value),
                        "n_cells": n,
                        "affected_cells": affected,
                        "affected_fraction": float(affected / max(n, 1)),
                    }
                )
            if rows:
                tables.setdefault("sample", []).extend(rows)

    flags: list[dict[str, Any]] = []
    for scope, rows in tables.items():
        by_signal: dict[str, list[Mapping[str, Any]]] = {}
        for row in rows:
            by_signal.setdefault(str(row.get("signal")), []).append(row)
        for signal, signal_rows in by_signal.items():
            fractions = [float(row.get("affected_fraction") or 0.0) for row in signal_rows]
            if fractions and max(fractions) >= 0.5:
                flags.append(
                    {
                        "scope": scope,
                        "signal": signal,
                        "flag": "signal_enriched_group",
                        "max_fraction": max(fractions),
                        "review_required": True,
                    }
                )
            if len(fractions) >= 2 and max(fractions) - min(fractions) >= 0.35:
                flags.append(
                    {
                        "scope": scope,
                        "signal": signal,
                        "flag": "signal_imbalance",
                        "spread": max(fractions) - min(fractions),
                        "review_required": True,
                    }
                )

    return _json_safe(
        {
            "schema_version": "post_annotation_qc_review_v1",
            "available": True,
            "cell_type_key": cell_type_key,
            "sample_key": sample_key,
            "review_columns": available,
            "review_required": bool(recommended_actions or flags),
            "recommended_actions": recommended_actions,
            "tables": tables,
            "flags": flags,
            "note": (
                "Post-annotation QC review is advisory. It should guide sensitivity analyses "
                "or explicit covariate modeling rather than automatic deletion."
            ),
        }
    )


def build_qc_benchmark_scorecard(
    *,
    benchmark_summary: Mapping[str, Any] | None,
    doublet_evidence_summary: Mapping[str, Any] | None,
    ambient_evidence_summary: Mapping[str, Any] | None,
    retention_audit_summary: Mapping[str, Any] | None,
    post_annotation_qc_review: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Create a compact scorecard spanning the seven benchmark-QC pillars."""
    rows: list[dict[str, Any]] = []

    def add(claim: str, status: str, evidence: str, next_action: str, score: float) -> None:
        rows.append(
            {
                "claim": claim,
                "status": status,
                "evidence": evidence,
                "next_action": next_action,
                "score": float(max(0.0, min(1.0, score))),
            }
        )

    filtering_basis = None
    if isinstance(benchmark_summary, Mapping):
        filtering_basis = benchmark_summary.get("status")
    add(
        "threshold_and_retention_audit",
        "pass" if filtering_basis == "pass" else "review_required",
        f"benchmark_status={filtering_basis}",
        "Inspect benchmark_summary and retention audit before publication.",
        1.0 if filtering_basis == "pass" else 0.65,
    )

    doublet_status = doublet_evidence_summary.get("status") if isinstance(doublet_evidence_summary, Mapping) else "not_run"
    add(
        "sample_aware_doublet_evidence",
        "review_required" if doublet_evidence_summary and doublet_evidence_summary.get("review_required") else ("pass" if doublet_status == "available" else "missing"),
        f"doublet_status={doublet_status}",
        "Use sample-aware expected-rate calibration and inspect boundary cells.",
        0.85 if doublet_status == "available" else 0.35,
    )

    ambient_status = ambient_evidence_summary.get("status") if isinstance(ambient_evidence_summary, Mapping) else "not_run"
    ambient_review = bool(ambient_evidence_summary.get("review_required")) if isinstance(ambient_evidence_summary, Mapping) else True
    add(
        "ambient_rna_contract",
        "review_required" if ambient_review else "pass",
        f"ambient_status={ambient_status}; risk={ambient_evidence_summary.get('risk_level') if isinstance(ambient_evidence_summary, Mapping) else None}",
        "Register external CellBender/SoupX/DecontX evidence when ambient risk is moderate/high.",
        0.75 if ambient_status == "available" else 0.40,
    )

    retention_review = bool(retention_audit_summary.get("retention_review_required")) if isinstance(retention_audit_summary, Mapping) else True
    add(
        "stratified_retention_fairness",
        "review_required" if retention_review else "pass",
        f"retention_available={retention_audit_summary.get('available') if isinstance(retention_audit_summary, Mapping) else False}",
        "Review retention by sample/condition/annotation/cluster.",
        0.75 if retention_audit_summary and retention_audit_summary.get("available") else 0.30,
    )

    post_available = bool(post_annotation_qc_review.get("available")) if isinstance(post_annotation_qc_review, Mapping) else False
    add(
        "post_annotation_sensitivity_review",
        "review_required" if post_annotation_qc_review and post_annotation_qc_review.get("review_required") else ("pass" if post_available else "not_applicable"),
        f"post_annotation_available={post_available}",
        "After annotation, decide keep/remove/sensitivity/model for stress, ambient, MT, and doublet signals.",
        0.80 if post_available else 0.55,
    )

    mean_score = sum(row["score"] for row in rows) / max(len(rows), 1)
    status = (
        "pass"
        if all(row["status"] == "pass" for row in rows)
        else "review_required"
        if any(row["status"] in {"review_required", "missing"} for row in rows)
        else "partial"
    )
    return _json_safe(
        {
            "schema_version": "qc_benchmark_scorecard_v1",
            "status": status,
            "score": mean_score,
            "rows": rows,
            "note": "Scorecard summarizes evidence maturity; it is not a substitute for dataset-specific review.",
        }
    )


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
    ambient_contract = (
        adata.uns.get("sclucid", {}).get("qc", {}).get("ambient_layer_contract", {})
    )
    ambient_counts_layer = ambient_contract.get("recommended_preprocess_counts_layer")
    ambient_corrected_present = bool(ambient_contract.get("corrected_layer_present"))
    recommended_counts_layer = (
        ambient_counts_layer
        if ambient_corrected_present and ambient_counts_layer
        else ("counts" if has_counts_layer else None)
    )

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
            f"Use adata.layers['{recommended_counts_layer}'] as the preprocessing input."
            if recommended_counts_layer
            else "Create or preserve adata.layers['counts'] before normalization when raw counts are available."
        ),
        priority="required" if recommended_counts_layer is None else "recommended",
        rationale=(
            "Ambient RNA correction produced a reviewer-visible corrected counts layer."
            if ambient_corrected_present
            else "Preprocess needs an auditable raw-count source after QC filtering."
        ),
        suggested_config={"layer": recommended_counts_layer},
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
                "ambient_layer_contract": ambient_contract or None,
            },
            "filtering_context": {
                "initial_cells": filtering_summary.get("initial_cells"),
                "final_cells": filtering_summary.get("final_cells"),
                "retention_fraction": retention,
                "sample_specific_thresholds_available": bool(sample_thresholds),
            },
        }
    )


def build_qc_handoff_readiness(
    *,
    adata: AnnData,
    context: Mapping[str, Any],
    filtering_summary: Mapping[str, Any],
    output_health: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    filtering_policy_summary: Mapping[str, Any] | None = None,
    retention_audit_summary: Mapping[str, Any] | None = None,
    doublet_evidence_summary: Mapping[str, Any] | None = None,
    ambient_evidence_summary: Mapping[str, Any] | None = None,
    tumor_aware_summary: Mapping[str, Any] | None = None,
    post_annotation_qc_review: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Declare the QC-to-preprocess handoff contract for downstream use."""
    blockers = list(output_health.get("issues", []))
    blockers.extend(downstream_recommendations.get("blockers", []))
    blockers = list(dict.fromkeys(str(item) for item in blockers))
    review_items: list[str] = []
    warnings: list[str] = []

    downstream_inputs = downstream_recommendations.get("input_assumptions", {})
    filtering_context = downstream_recommendations.get("filtering_context", {})
    recommended_counts_layer = None
    for item in downstream_recommendations.get("recommendations", []):
        if isinstance(item, Mapping) and item.get("target") == "counts_layer":
            suggested = item.get("suggested_config", {})
            if isinstance(suggested, Mapping):
                recommended_counts_layer = suggested.get("layer")
            break
    if recommended_counts_layer is None and "counts" in adata.layers:
        recommended_counts_layer = "counts"

    obs = adata.obs

    def _bool_count(column: str) -> int:
        if column not in obs:
            return 0
        return int(obs[column].fillna(False).astype(bool).sum())

    def _fraction(count: int) -> float | None:
        return float(count) / float(adata.n_obs) if adata.n_obs else None

    qc_decision_counts = (
        obs["qc_decision"].astype(str).value_counts().to_dict()
        if "qc_decision" in obs
        else {}
    )
    remove_count = _bool_count("qc_remove")
    review_count = _bool_count("qc_review_required")
    sensitivity_count = int(qc_decision_counts.get("sensitivity_only", 0))
    doublet_like_count = _bool_count("predicted_doublet")
    ambient_risk_count = _bool_count("ambient_risk")
    stress_high_count = _bool_count("stress_high")
    high_mt_count = _bool_count("qc_high_mt") or _bool_count("outlier_mt")

    if not recommended_counts_layer:
        blockers.append("No auditable raw counts layer was available for preprocessing handoff.")
    elif recommended_counts_layer not in adata.layers:
        blockers.append(
            f"Recommended preprocess counts layer {recommended_counts_layer!r} is not present in adata.layers."
        )

    policy = filtering_policy_summary or {}
    if policy.get("review_required"):
        review_items.append(str(policy.get("risk_note")))
    retention = output_health.get("retention_fraction")
    if isinstance(retention, (int, float)) and retention < 0.5:
        review_items.append(
            f"QC retained {retention:.1%} of cells; inspect sample/cell-type retention before preprocessing."
        )
    retention_audit = retention_audit_summary or {}
    if retention_audit.get("retention_review_required"):
        review_items.append(
            str(
                retention_audit.get(
                    "mandatory_review_note",
                    "Retention audit requires review before downstream interpretation.",
                )
            )
        )
    if (doublet_evidence_summary or {}).get("review_required"):
        review_items.append("Doublet evidence requires review before irreversible exclusion.")
    if (ambient_evidence_summary or {}).get("review_required"):
        review_items.append("Ambient RNA evidence requires review before choosing preprocess counts.")
    if (post_annotation_qc_review or {}).get("review_required"):
        review_items.append("Post-annotation QC review recommends sensitivity/review handling.")

    tumor_aware = bool((tumor_aware_summary or {}).get("enabled")) or _is_tumor_context(
        context.get("tissue_type")
    )
    if tumor_aware:
        tumor_warnings = list((tumor_aware_summary or {}).get("warnings", []))
        review_items.extend(str(item) for item in tumor_warnings)
        warnings.append(
            "Tumor-aware QC should keep high-MT, stress-high, or cycling malignant-like states in review/sensitivity records unless clearly technical."
        )

    if review_count or sensitivity_count:
        warnings.append(
            "Review and sensitivity_only cells should be tracked through preprocessing outputs for sensitivity analysis, not silently dropped from interpretation."
        )

    status = "blocked" if blockers else ("review_required" if review_items or warnings else "ready")
    return _json_safe(
        {
            "schema_version": "qc_handoff_readiness_v1",
            "status": status,
            "ready_for_preprocess": not blockers,
            "recommended_preprocess_counts_layer": recommended_counts_layer,
            "counts_layer_present": bool(
                recommended_counts_layer and recommended_counts_layer in adata.layers
            ),
            "filtering_decision_columns": {
                "qc_decision": "qc_decision" in obs,
                "qc_remove": "qc_remove" in obs,
                "qc_review_required": "qc_review_required" in obs,
            },
            "cell_decision_counts": qc_decision_counts,
            "handoff_cell_counts": {
                "n_cells": int(adata.n_obs),
                "remove": remove_count,
                "review_required": review_count,
                "sensitivity_only": sensitivity_count,
                "doublet_like": doublet_like_count,
                "ambient_risk": ambient_risk_count,
                "stress_high": stress_high_count,
                "high_mitochondrial": high_mt_count,
            },
            "handoff_cell_fractions": {
                "remove": _fraction(remove_count),
                "review_required": _fraction(review_count),
                "sensitivity_only": _fraction(sensitivity_count),
                "doublet_like": _fraction(doublet_like_count),
                "ambient_risk": _fraction(ambient_risk_count),
                "stress_high": _fraction(stress_high_count),
                "high_mitochondrial": _fraction(high_mt_count),
            },
            "safe_to_continue": {
                "preprocess": not blockers,
                "graph_analysis": not blockers,
                "marker_annotation": not blockers,
                "condition_de": not blockers and not review_items,
                "tumor_state_interpretation": not blockers and not tumor_aware,
            },
            "required_downstream_handling": [
                "Use the recommended counts layer for normalization and layer contracts.",
                "Carry qc_decision/qc_remove/qc_review_required columns into downstream AnnData objects.",
                "Use review/sensitivity_only cells for sensitivity analysis when biological fragility is plausible.",
            ],
            "input_assumptions": {
                "sample_key": context.get("sample_key"),
                "n_samples": context.get("n_samples"),
                "tumor_aware": tumor_aware,
                "downstream_ready_for_preprocess": downstream_recommendations.get(
                    "ready_for_preprocess"
                ),
                "filtering_context": filtering_context,
                "downstream_input_assumptions": downstream_inputs,
            },
            "blockers": blockers,
            "review_items": list(dict.fromkeys(str(item) for item in review_items if item)),
            "warnings": warnings,
        }
    )


def build_qc_readiness_assessment(
    *,
    output_health: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    benchmark_summary: Mapping[str, Any] | None,
    tumor_aware_summary: Mapping[str, Any] | None,
    doublet_evidence_summary: Mapping[str, Any] | None = None,
    ambient_evidence_summary: Mapping[str, Any] | None = None,
    filtering_policy_summary: Mapping[str, Any] | None = None,
    retention_audit_summary: Mapping[str, Any] | None = None,
    post_annotation_qc_review: Mapping[str, Any] | None = None,
    qc_benchmark_scorecard: Mapping[str, Any] | None = None,
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

    if doublet_evidence_summary and doublet_evidence_summary.get("review_required"):
        doublet_notes = doublet_evidence_summary.get("notes", [])
        review_reasons.append(
            "Doublet evidence requires review"
            + (": " + "; ".join(str(note) for note in doublet_notes) if doublet_notes else ".")
        )
    if ambient_evidence_summary and ambient_evidence_summary.get("review_required"):
        ambient_notes = ambient_evidence_summary.get("notes", [])
        review_reasons.append(
            "Ambient RNA evidence requires review"
            + (": " + "; ".join(str(note) for note in ambient_notes) if ambient_notes else ".")
        )
    if filtering_policy_summary and filtering_policy_summary.get("review_required"):
        review_reasons.append(str(filtering_policy_summary.get("risk_note")))
    if retention_audit_summary and retention_audit_summary.get("retention_review_required"):
        review_reasons.append(
            str(
                retention_audit_summary.get(
                    "mandatory_review_note",
                    "Retention audit requires review before downstream interpretation.",
                )
            )
        )
        for flag in retention_audit_summary.get("flags", []):
            if isinstance(flag, Mapping):
                review_reasons.append(
                    f"Retention audit {flag.get('scope')}:{flag.get('flag')}"
                )
    if post_annotation_qc_review and post_annotation_qc_review.get("review_required"):
        review_reasons.append(
            "Post-annotation QC review recommends sensitivity/review handling for retained QC signals."
        )
    if qc_benchmark_scorecard and qc_benchmark_scorecard.get("status") != "pass":
        review_reasons.append(
            f"QC benchmark scorecard status is {qc_benchmark_scorecard.get('status')}."
        )

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
    doublet_evidence_summary: Mapping[str, Any] | None,
    ambient_evidence_summary: Mapping[str, Any] | None = None,
    filtering_policy_summary: Mapping[str, Any] | None = None,
    retention_audit_summary: Mapping[str, Any] | None = None,
    post_annotation_qc_review: Mapping[str, Any] | None = None,
    qc_benchmark_scorecard: Mapping[str, Any] | None = None,
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

    if doublet_evidence_summary:
        doublet_status = doublet_evidence_summary.get("status")
        doublet_notes = doublet_evidence_summary.get("notes", [])
        if doublet_evidence_summary.get("review_required"):
            add(
                priority="review",
                action="Inspect doublet evidence before finalizing QC filtering.",
                rationale="; ".join(str(note) for note in doublet_notes)
                or "Doublet evidence summary requires review.",
                evidence_key="doublet_evidence_summary",
            )
        elif doublet_status == "not_run":
            add(
                priority="optional",
                action="Run doublet detection when the experiment has meaningful multiplet risk.",
                rationale="No doublet prediction or score columns were found in this QC result.",
                evidence_key="doublet_evidence_summary",
            )

    if ambient_evidence_summary:
        if ambient_evidence_summary.get("review_required"):
            add(
                priority="review",
                action="Review ambient RNA evidence and counts-layer contract.",
                rationale="; ".join(str(note) for note in ambient_evidence_summary.get("notes", []))
                or "Ambient evidence summary requires review.",
                evidence_key="ambient_evidence_summary",
            )
        if not ambient_evidence_summary.get("cell_probability", {}).get("available", False):
            add(
                priority="optional",
                action="Register CellBender/EmptyDrops-style cell probabilities when available.",
                rationale="cell_probability is part of the canonical QC schema but no measured values were detected.",
                evidence_key="ambient_evidence_summary.cell_probability",
            )

    if filtering_policy_summary and filtering_policy_summary.get("review_required"):
        add(
            priority="review",
            action="Review the final QC filtering basis before preprocessing.",
            rationale=str(filtering_policy_summary.get("risk_note", "")),
            evidence_key="qc_filtering_policy_summary",
        )
    if retention_audit_summary and retention_audit_summary.get("retention_review_required"):
        add(
            priority="review",
            action="Review QC retention by sample, condition, annotation, and cluster.",
            rationale=str(
                retention_audit_summary.get(
                    "mandatory_review_note",
                    "Stratified retention should be checked before downstream interpretation.",
                )
            ),
            evidence_key="qc_retention_audit_summary",
        )
    if post_annotation_qc_review and post_annotation_qc_review.get("review_required"):
        add(
            priority="review",
            action="Run or inspect post-annotation QC sensitivity decisions.",
            rationale=(
                "Stress, ambient, MT, apoptosis, contamination, or doublet-like retained cells "
                "should be interpreted after cell-type annotation."
            ),
            evidence_key="post_annotation_qc_review",
        )
    if qc_benchmark_scorecard and qc_benchmark_scorecard.get("status") != "pass":
        add(
            priority="review",
            action="Inspect QC benchmark scorecard before claiming benchmark-module readiness.",
            rationale=f"Scorecard status is {qc_benchmark_scorecard.get('status')}.",
            evidence_key="qc_benchmark_scorecard",
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
    filtering_policy_summary: Mapping[str, Any],
    retention_audit_summary: Mapping[str, Any],
    ambient_evidence_summary: Mapping[str, Any],
    post_annotation_qc_review: Mapping[str, Any],
    qc_benchmark_scorecard: Mapping[str, Any],
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
            "stage": "filtering_policy",
            "final_filter_basis": filtering_policy_summary.get("final_filter_basis"),
            "qc_decision_filter_mode": filtering_policy_summary.get(
                "qc_decision_filter_mode"
            ),
            "review_required": filtering_policy_summary.get("review_required"),
        },
        {
            "stage": "retention_audit",
            "available": retention_audit_summary.get("available"),
            "review_required": retention_audit_summary.get("retention_review_required"),
            "group_keys_reviewed": retention_audit_summary.get("group_keys_reviewed", {}),
            "flags": retention_audit_summary.get("flags", []),
        },
        {
            "stage": "ambient_evidence",
            "status": ambient_evidence_summary.get("status"),
            "risk_level": ambient_evidence_summary.get("risk_level"),
            "review_required": ambient_evidence_summary.get("review_required"),
        },
        {
            "stage": "post_annotation_qc_review",
            "available": post_annotation_qc_review.get("available"),
            "review_required": post_annotation_qc_review.get("review_required"),
            "cell_type_key": post_annotation_qc_review.get("cell_type_key"),
        },
        {
            "stage": "qc_benchmark_scorecard",
            "status": qc_benchmark_scorecard.get("status"),
            "score": qc_benchmark_scorecard.get("score"),
        },
        {
            "stage": "output_health",
            "status": output_health.get("status"),
            "issues": output_health.get("issues", []),
        },
    ]
