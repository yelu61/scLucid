"""QC trace and review-summary schema helpers.

This module keeps the workflow-facing QC audit contract independent from the
execution code. It turns recommendations, applied config, and filtering output
into a compact machine-readable decision table and evidence chain.
"""

from __future__ import annotations

from collections.abc import Mapping
from enum import Enum
from typing import Any

import pandas as pd
from anndata import AnnData

from scLucid.utils.contracts import _review_payload

from ..utils.context import is_tumor_context as _shared_is_tumor_context
from ..utils.evidence import (
    DecisionRecord,
    EvidenceBundle,
    EvidenceItem,
    ReviewAction,
    model_to_dict,
)

QC_TRACE_SCHEMA_VERSION = "1.0"
QC_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

QC_REQUIRED_REVIEW_SECTIONS = {
    "recommendation_summary",
    "recommended_threshold_summary",
    "applied_threshold_summary",
    "user_override_summary",
    "sample_threshold_summary",
    "tumor_aware_summary",
    "filtering_summary",
    "qc_filtering_policy_summary",
    "qc_retention_audit_summary",
    "doublet_evidence_summary",
    "warnings",
    "decision_table",
    "qc_reviewer_table",
    "evidence_chain",
    "execution_trace",
    "output_health",
    "downstream_preprocess_recommendations",
    "qc_readiness",
    "review_action_items",
    "reproducibility_manifest",
    "evidence_bundle",
    "module_maturity",
}

QC_REQUIRED_OBS_METRICS = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
]

QC_STABLE_ENTRYPOINTS = (
    "scLucid.qc.run_qc",
    "scLucid.qc.recommend_qc_policy",
    "scLucid.qc.apply_qc_policy",
    "scLucid.qc.run_standard_qc",
    "scLucid.qc.calculate_qc_metric",
    "scLucid.qc.recommend_intelligent_qc",
    "scLucid.qc.run_qc_threshold_decision",
    "scLucid.qc.build_qc_decisions",
    "scLucid.qc.filter_cells",
)

QC_EXPECTED_ARTIFACTS = (
    "qc_review_summary.json",
    "qc_review_summary.md",
    "qc_benchmark.json",
    "qc_benchmark.md",
)

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


def build_qc_readiness_assessment(
    *,
    output_health: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    benchmark_summary: Mapping[str, Any] | None,
    tumor_aware_summary: Mapping[str, Any] | None,
    doublet_evidence_summary: Mapping[str, Any] | None = None,
    filtering_policy_summary: Mapping[str, Any] | None = None,
    retention_audit_summary: Mapping[str, Any] | None = None,
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
    filtering_policy_summary: Mapping[str, Any] | None = None,
    retention_audit_summary: Mapping[str, Any] | None = None,
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
        "filtering_policy": "contract",
        "retention_audit": "metric",
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
            "qc_reviewer_table",
            "evidence_chain",
            "qc_readiness",
            "review_action_items",
            "reproducibility_manifest",
            "qc_filtering_policy_summary",
            "qc_retention_audit_summary",
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
    if not isinstance(output_health, Mapping) or "status" not in output_health:
        issues.append("QC output_health summary must be present.")
    if not isinstance(manifest, Mapping) or manifest.get("workflow") != "run_standard_qc":
        issues.append("QC reproducibility_manifest must identify workflow='run_standard_qc'.")

    review_required = []
    if isinstance(readiness, Mapping) and readiness.get("status") != "ready":
        review_required.append(f"qc_readiness.status={readiness.get('status')}")
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
    summary["qc_schema_version"] = QC_TRACE_SCHEMA_VERSION
    summary["decision_table"] = decision_table
    summary["threshold_reviewer_table"] = decision_table
    summary["doublet_evidence_summary"] = doublet_evidence_summary
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
        doublet_evidence_summary=doublet_evidence_summary,
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
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
        filtering_policy_summary=filtering_policy_summary,
        retention_audit_summary=retention_audit_summary,
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
