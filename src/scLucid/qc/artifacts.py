"""QC artifact contract and cleanup helpers.

The QC workflow intentionally creates several classes of outputs:
columns used by filtering, review-only evidence columns, optional preprocessing
layers, and temporary intermediate scores.  Keeping that contract explicit
prevents diagnostic evidence from being mistaken for a hard filtering rule.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Literal, Optional

from anndata import AnnData
import pandas as pd

FILTERING_REQUIRED_OBS_COLUMNS = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_hb",
    "predicted_doublet",
    "combined_doublet_score",
    "doublet_score",
    "ambient_score",
    "ambient_fraction",
    "cell_probability",
    "empty_droplet_probability",
    "outlier_min_genes",
    "outlier_max_genes",
    "outlier_min_counts",
    "outlier_max_counts",
    "outlier_mt",
    "outlier_hb",
    "outlier_qc_metrics",
    "outlier_count",
    "qc_decision",
    "qc_remove",
    "qc_reason",
    "qc_confidence",
]

REVIEW_OBS_COLUMNS = [
    "pct_counts_ribo",
    "pct_counts_in_top_20_genes",
    "pct_counts_in_top_50_genes",
    "pct_counts_in_top_100_genes",
    "algorithm_doublet_score",
    "algorithm_predicted_doublet",
    "heuristic_confidence_score",
    "heuristic_predicted",
    "heterotypic_doublet_risk",
    "homotypic_doublet_risk",
    "external_doublet_evidence",
    "qc_low_counts",
    "qc_low_genes",
    "qc_high_mt",
    "qc_low_complexity",
    "qc_high_hb",
    "platelet_contamination",
    "hemoglobin_contamination",
    "ambient_risk",
    "stress_high",
    "apoptosis_high",
    "qc_review_required",
    "qc_biological_risk_note",
    "qc_phase",
    "S_score",
    "G2M_score",
    "phase",
]

OPTIONAL_PREPROCESS_LAYERS = [
    "counts",
    "ambient_corrected_counts",
]

INTERMEDIATE_OBS_PATTERNS = [
    "scrublet_",
    "scanpy_scrublet_",
    "solo_",
    "doubletdetection_",
    "scdblfinder_",
    "expected_total_doublet_rate",
    "expected_heterotypic_doublet_rate",
    "expected_homotypic_doublet_rate",
]

QC_ARTIFACT_CONTRACT_SCHEMA_VERSION = "qc_artifact_contract_v1"

QC_ARTIFACT_CONTRACT: Dict[str, Any] = {
    "schema_version": QC_ARTIFACT_CONTRACT_SCHEMA_VERSION,
    "decision_flow": [
        "threshold_recommendation",
        "threshold_decision",
        "mark_evidence",
        "qc_decision",
        "filter_cells",
        "benchmark_review",
    ],
    "filtering_required_obs_columns": FILTERING_REQUIRED_OBS_COLUMNS,
    "review_obs_columns": REVIEW_OBS_COLUMNS,
    "optional_preprocess_layers": OPTIONAL_PREPROCESS_LAYERS,
    "intermediate_obs_patterns": INTERMEDIATE_OBS_PATTERNS,
    "semantics": {
        "filtering_required": "Columns consumed by filtering or final QC decisions.",
        "review": "Columns intended for diagnostics, plots, or manual review.",
        "optional_preprocess_layer": "Layers that may be selected by preprocessing.",
        "intermediate": "Temporary or method-specific values safe to remove after final QC review.",
        "threshold_recommendation": "Candidate thresholds from statistical or intelligent recommenders.",
        "threshold_decision": "Resolved thresholds selected for marking, including source policy.",
        "mark_evidence": "Boolean evidence columns created from thresholds and other QC signals.",
        "qc_decision": "Cell-level quality label derived from combined evidence.",
        "filter_cells": "Irreversible/subsetting step that consumes labels or evidence columns.",
        "benchmark_review": "Post-filter retention and biological-fidelity checks.",
    },
}


def get_qc_artifact_contract() -> Dict[str, Any]:
    """Return a copy of the QC artifact contract."""
    return dict(QC_ARTIFACT_CONTRACT)


def record_qc_artifact_contract(adata: AnnData) -> Dict[str, Any]:
    """Store the QC artifact contract in ``adata.uns`` and return it."""
    contract = get_qc_artifact_contract()
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["artifact_contract"] = contract
    return contract


def _qc_namespace(adata: AnnData) -> Dict[str, Any]:
    """Return the scLucid QC namespace, creating it if needed."""
    return adata.uns.setdefault("sclucid", {}).setdefault("qc", {})


def _safe_count_bool_obs(adata: AnnData, columns: Iterable[str]) -> Dict[str, int]:
    """Count truthy values for existing boolean-like obs columns."""
    counts: Dict[str, int] = {}
    for col in columns:
        if col in adata.obs:
            counts[col] = int(pd.Series(adata.obs[col]).fillna(False).astype(bool).sum())
    return counts


def record_threshold_recommendation(
    adata: AnnData,
    *,
    source: str,
    payload: Dict[str, Any],
    append: bool = True,
) -> Dict[str, Any]:
    """Record candidate threshold recommendations before final resolution."""
    record = {
        "schema_version": "qc_threshold_recommendation_v1",
        "source": source,
        "payload": payload,
    }
    qc_ns = _qc_namespace(adata)
    if append:
        qc_ns.setdefault("threshold_recommendations", []).append(record)
    else:
        qc_ns["threshold_recommendations"] = [record]
    record_qc_artifact_contract(adata)
    return record


def record_threshold_decision(
    adata: AnnData,
    *,
    resolved_thresholds: Dict[str, Any],
    policy: str,
    sources: Dict[str, Any],
    sample_key: Optional[str] = None,
    sample_thresholds: Optional[Dict[str, Any]] = None,
    notes: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Record the final threshold decision used to mark QC evidence."""
    record = {
        "schema_version": "qc_threshold_decision_v1",
        "resolved_thresholds": resolved_thresholds,
        "threshold_policy": policy,
        "sources": sources,
        "sample_key": sample_key,
        "sample_thresholds": sample_thresholds or {},
        "notes": list(notes or []),
        "decision_role": "Input to QC evidence marking; not itself a cell-level label.",
    }
    _qc_namespace(adata)["threshold_decision"] = record
    record_qc_artifact_contract(adata)
    return record


def record_mark_evidence(
    adata: AnnData,
    *,
    evidence_columns: Iterable[str],
    thresholds: Dict[str, Any],
    sample_key: Optional[str] = None,
    sample_thresholds: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Record boolean evidence columns created by QC marking."""
    evidence_columns = [col for col in evidence_columns if col in adata.obs]
    record = {
        "schema_version": "qc_mark_evidence_v1",
        "sample_key": sample_key,
        "evidence_columns": evidence_columns,
        "evidence_counts": _safe_count_bool_obs(adata, evidence_columns),
        "thresholds": thresholds,
        "sample_thresholds": sample_thresholds or {},
        "note": "Marking creates evidence columns. Removal is handled by qc_decision/filter_cells.",
    }
    _qc_namespace(adata)["mark_evidence"] = record
    record_qc_artifact_contract(adata)
    return record


def record_qc_decision_artifact(
    adata: AnnData,
    *,
    summary: Dict[str, Any],
    evidence_columns: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Record cell-level QC decision metadata in the artifact flow."""
    record = {
        "schema_version": "qc_decision_artifact_v1",
        "summary": summary,
        "decision_columns": [
            col
            for col in [
                "qc_decision",
                "qc_remove",
                "qc_reason",
                "qc_confidence",
                "qc_review_required",
            ]
            if col in adata.obs
        ],
        "evidence_columns": list(evidence_columns or []),
        "decision_counts": (
            adata.obs["qc_decision"].astype(str).value_counts().to_dict()
            if "qc_decision" in adata.obs
            else {}
        ),
    }
    _qc_namespace(adata)["qc_decision_artifact"] = record
    record_qc_artifact_contract(adata)
    return record


def record_filter_result(
    adata: AnnData,
    *,
    stats: Dict[str, Any],
    result_key: str = "filtering_results",
) -> Dict[str, Any]:
    """Record filtering/subsetting results under the unified artifact flow."""
    record = {
        "schema_version": "qc_filter_result_v1",
        **stats,
        "filter_role": "Subsetting step that consumes evidence/decision columns.",
    }
    qc_ns = _qc_namespace(adata)
    qc_ns[result_key] = record
    qc_ns["filter_result"] = record
    record_qc_artifact_contract(adata)
    return record


def record_benchmark_review(
    adata: AnnData,
    *,
    benchmark_summary: Dict[str, Any],
) -> Dict[str, Any]:
    """Record post-filter benchmark/review output in the artifact flow."""
    record = {
        "schema_version": "qc_benchmark_review_v1",
        "benchmark_summary": benchmark_summary,
        "review_role": "Post-filter retention and biological-fidelity review.",
    }
    _qc_namespace(adata)["benchmark_review"] = record
    record_qc_artifact_contract(adata)
    return record


def _matches_any_prefix(name: str, prefixes: Iterable[str]) -> bool:
    """Return True if ``name`` starts with any configured intermediate prefix."""
    return any(name.startswith(prefix) for prefix in prefixes)


def cleanup_qc_intermediates(
    adata: AnnData,
    *,
    mode: Literal["review", "minimal"] = "review",
    drop_obsm: bool = False,
    extra_keep_obs: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Drop selected intermediate QC artifacts after final QC decisions.

    ``mode='review'`` keeps filtering and review evidence columns, removing only
    method-specific temporary columns. ``mode='minimal'`` keeps only filtering
    and final decision columns. The function never removes layers or ``.uns``
    records; it reports what was dropped so the cleanup itself is auditable.
    """
    contract = get_qc_artifact_contract()
    keep = set(contract["filtering_required_obs_columns"])
    if mode == "review":
        keep.update(contract["review_obs_columns"])
    keep.update(extra_keep_obs or [])

    existing = set(adata.obs.columns)
    prefixes = tuple(contract["intermediate_obs_patterns"])
    if mode == "minimal":
        candidates = existing - keep
    else:
        candidates = {col for col in existing if _matches_any_prefix(col, prefixes) and col not in keep}

    dropped_obs = sorted(candidates)
    if dropped_obs:
        adata.obs.drop(columns=dropped_obs, inplace=True)

    dropped_obsm = []
    if drop_obsm and "lineage_module_scores" in adata.obsm:
        del adata.obsm["lineage_module_scores"]
        dropped_obsm.append("lineage_module_scores")

    summary = {
        "schema_version": "qc_artifact_cleanup_v1",
        "mode": mode,
        "dropped_obs_columns": dropped_obs,
        "dropped_obsm_keys": dropped_obsm,
        "contract_schema_version": contract["schema_version"],
    }
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["artifact_cleanup"] = summary
    record_qc_artifact_contract(adata)
    return summary


__all__ = [
    "QC_ARTIFACT_CONTRACT_SCHEMA_VERSION",
    "QC_ARTIFACT_CONTRACT",
    "FILTERING_REQUIRED_OBS_COLUMNS",
    "REVIEW_OBS_COLUMNS",
    "OPTIONAL_PREPROCESS_LAYERS",
    "INTERMEDIATE_OBS_PATTERNS",
    "get_qc_artifact_contract",
    "record_qc_artifact_contract",
    "record_threshold_recommendation",
    "record_threshold_decision",
    "record_mark_evidence",
    "record_qc_decision_artifact",
    "record_filter_result",
    "record_benchmark_review",
    "cleanup_qc_intermediates",
]
