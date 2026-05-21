"""Review-summary enrichment for the analysis workflow."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd
from anndata import AnnData

from ..utils import sanitize_for_hdf5
from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict

ANALYSIS_TRACE_SCHEMA_VERSION = "1.0"
ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

ANALYSIS_REQUIRED_REVIEW_SECTIONS = {
    "analysis_schema_version",
    "preprocess_input_context",
    "clustering_evidence_summary",
    "annotation_evidence_summary",
    "annotation_consensus_summary",
    "malignancy_interpretation_summary",
    "analysis_readiness",
    "review_action_items",
    "evidence_bundle",
    "module_maturity",
}

ANALYSIS_STABLE_ENTRYPOINTS = (
    "scLucid.analysis.run_standard_analysis",
    "scLucid.analysis.run_clustering_review",
    "scLucid.analysis.run_annotation_evidence",
    "scLucid.analysis.build_annotation_consensus",
    "scLucid.analysis.run_malignancy_interpretation",
    "scLucid.analysis.cluster_cells",
    "scLucid.analysis.find_markers",
    "scLucid.analysis.run_annotation",
)

ANALYSIS_EXPECTED_OUTPUTS = (
    "adata.obs['leiden_clusters']",
    "adata.uns['rank_genes_groups']",
    "adata.uns['sclucid']['analysis']['annotation']['annotation_review_table']",
    "adata.uns['sclucid']['analysis']['malignancy']['malignancy_interpretation_summary']",
    "adata.uns['sclucid']['analysis']['review_summary']",
)


def enrich_analysis_review_summary(
    summary: dict[str, Any],
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    cluster_key: str,
) -> dict[str, Any]:
    """Add benchmark-grade review fields to an analysis summary."""
    summary = dict(summary)
    preprocess_context = build_preprocess_input_context(adata)
    clustering = build_clustering_evidence_summary(adata, cluster_key)
    annotation = build_annotation_evidence_summary(adata, config=config)
    consensus = build_annotation_consensus_summary(adata, config=config)
    malignancy = build_malignancy_interpretation_summary(adata, config=config)
    readiness = build_analysis_readiness_assessment(
        adata=adata,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
        preprocess_context=preprocess_context,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        malignancy_summary=malignancy,
    )
    actions = build_analysis_review_action_items(
        readiness=readiness,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        malignancy_summary=malignancy,
    )

    summary["analysis_schema_version"] = ANALYSIS_TRACE_SCHEMA_VERSION
    summary["preprocess_input_context"] = preprocess_context
    summary["clustering_evidence_summary"] = clustering
    summary["annotation_evidence_summary"] = annotation
    summary["annotation_consensus_summary"] = consensus
    summary["malignancy_interpretation_summary"] = malignancy
    summary["analysis_readiness"] = readiness
    summary["review_action_items"] = actions
    summary["evidence_bundle"] = build_analysis_evidence_bundle(summary)
    summary["module_maturity"] = build_analysis_module_maturity_assessment(summary)
    return _json_safe(summary)


def get_analysis_module_contract() -> dict[str, Any]:
    """Return the frozen analysis module maturity contract."""
    return {
        "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "analysis",
        "stable_entrypoints": list(ANALYSIS_STABLE_ENTRYPOINTS),
        "required_review_sections": sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS),
        "expected_outputs": list(ANALYSIS_EXPECTED_OUTPUTS),
        "canonical_namespace": 'adata.uns["sclucid"]["analysis"]',
        "readiness_key": "analysis_readiness",
        "clustering_evidence_key": "clustering_evidence_summary",
        "annotation_evidence_key": "annotation_evidence_summary",
        "annotation_consensus_key": "annotation_consensus_summary",
        "malignancy_interpretation_key": "malignancy_interpretation_summary",
        "preprocess_input_key": "preprocess_input_context",
    }


def build_preprocess_input_context(adata: AnnData) -> dict[str, Any]:
    """Summarize the preprocessing state consumed by analysis."""
    pp_ns = adata.uns.get("sclucid", {}).get("preprocess", {})
    review = _review_payload(pp_ns.get("review_summary", {})) if isinstance(pp_ns, Mapping) else {}
    readiness = review.get("preprocess_readiness", {}) if isinstance(review, Mapping) else {}
    maturity = review.get("module_maturity", {}) if isinstance(review, Mapping) else {}

    return _json_safe(
        {
            "available": bool(review),
            "preprocess_readiness_status": (
                readiness.get("status") if isinstance(readiness, Mapping) else None
            ),
            "preprocess_readiness_score": (
                readiness.get("score") if isinstance(readiness, Mapping) else None
            ),
            "preprocess_maturity_status": (
                maturity.get("status") if isinstance(maturity, Mapping) else None
            ),
            "normalized_layer_present": "normalized" in adata.layers,
            "pca_present": "X_pca" in adata.obsm,
            "neighbors_present": "neighbors" in adata.uns,
            "umap_present": "X_umap" in adata.obsm,
            "hvg_present": "highly_variable" in adata.var,
        }
    )


def build_clustering_evidence_summary(adata: AnnData, cluster_key: str) -> dict[str, Any]:
    """Summarize clustering output and optional resolution-review evidence."""
    clustering_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("clustering", {})
    review_summary = (
        clustering_ns.get("clustering_review_summary", {})
        if isinstance(clustering_ns, Mapping)
        else {}
    )
    cluster_series = adata.obs[cluster_key].astype(str) if cluster_key in adata.obs else None
    counts = cluster_series.value_counts() if cluster_series is not None else pd.Series(dtype=int)
    review_required = []
    if cluster_key not in adata.obs:
        review_required.append(f"cluster_key_missing:{cluster_key}")
    if isinstance(review_summary, Mapping) and review_summary.get("review_required_clusters"):
        review_required.append("resolution_review_flagged_clusters")

    return _json_safe(
        {
            "cluster_key": cluster_key,
            "available": cluster_key in adata.obs,
            "n_clusters": int(counts.shape[0]) if not counts.empty else 0,
            "min_cluster_size": int(counts.min()) if not counts.empty else None,
            "median_cluster_size": float(counts.median()) if not counts.empty else None,
            "resolution_review_available": bool(review_summary),
            "recommended_resolution": review_summary.get("recommended_resolution")
            if isinstance(review_summary, Mapping)
            else None,
            "recommended_cluster_key": review_summary.get("recommended_cluster_key")
            if isinstance(review_summary, Mapping)
            else None,
            "recommendation_rationale": review_summary.get("rationale")
            if isinstance(review_summary, Mapping)
            else "",
            "review_required": review_required,
        }
    )


def build_annotation_evidence_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize available annotation evidence tables."""
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    if not isinstance(annotation_ns, Mapping):
        annotation_ns = {}
    review_table = annotation_ns.get("annotation_review_table")
    marker_evidence = annotation_ns.get("marker_annotation_evidence")
    llm_bundle = annotation_ns.get("llm_annotation_bundle")
    annotation_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )

    n_review_rows = int(review_table.shape[0]) if hasattr(review_table, "shape") else 0
    n_needs_review = 0
    if isinstance(review_table, pd.DataFrame) and "needs_review" in review_table.columns:
        n_needs_review = int(review_table["needs_review"].fillna(True).astype(bool).sum())

    return _json_safe(
        {
            "annotation_key": annotation_key,
            "annotation_obs_present": annotation_key in adata.obs,
            "review_table_available": isinstance(review_table, pd.DataFrame),
            "review_table_rows": n_review_rows,
            "needs_review_clusters": n_needs_review,
            "marker_evidence_available": isinstance(marker_evidence, pd.DataFrame),
            "llm_bundle_available": isinstance(llm_bundle, Mapping),
            "evidence_methods": list(getattr(config, "annotation_methods", ()) or ()),
            "annotation_level": getattr(config, "annotation_level", None),
        }
    )


def build_annotation_consensus_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize final consensus labels applied to cells."""
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    review_table = (
        annotation_ns.get("annotation_review_table") if isinstance(annotation_ns, Mapping) else None
    )
    final_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )
    confidence_key = f"{final_key}_confidence"
    status_key = f"{final_key}_status"

    confidence = None
    if confidence_key in adata.obs:
        values = pd.to_numeric(adata.obs[confidence_key], errors="coerce")
        confidence = float(values.mean()) if values.notna().any() else None

    return _json_safe(
        {
            "final_key": final_key,
            "final_obs_present": final_key in adata.obs,
            "n_final_labels": int(adata.obs[final_key].nunique()) if final_key in adata.obs else 0,
            "mean_confidence": confidence,
            "status_key": status_key if status_key in adata.obs else None,
            "needs_review_cells": (
                int((adata.obs[status_key].astype(str) == "needs_review").sum())
                if status_key in adata.obs
                else 0
            ),
            "review_table_available": isinstance(review_table, pd.DataFrame),
        }
    )


def build_malignancy_interpretation_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize optional malignancy interpretation evidence."""
    malignancy_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("malignancy", {})
    stored = (
        malignancy_ns.get("malignancy_interpretation_summary", {})
        if isinstance(malignancy_ns, Mapping)
        else {}
    )
    call_key = getattr(config, "malignancy_key_added", "malignancy_call")
    score_key = getattr(config, "malignancy_score_key", "malignancy_interpretation_score")
    if isinstance(stored, Mapping) and stored:
        summary = dict(stored)
    else:
        summary = {
            "available": call_key in adata.obs,
            "call_key": call_key,
            "score_key": score_key,
            "review_required": False,
            "evidence_sources": [],
        }
    if call_key in adata.obs:
        calls = adata.obs[call_key].astype(str)
        summary.update(
            {
                "available": True,
                "n_malignant": int((calls == "malignant").sum()),
                "n_suspect_malignant": int((calls == "suspect_malignant").sum()),
                "n_non_malignant": int((calls == "non_malignant").sum()),
                "n_unresolved": int((calls == "unresolved").sum()),
            }
        )
    if score_key in adata.obs:
        scores = pd.to_numeric(adata.obs[score_key], errors="coerce")
        summary["mean_score"] = float(scores.mean()) if scores.notna().any() else None
    summary["enabled"] = bool(getattr(config, "run_malignancy_interpretation", False))
    return _json_safe(summary)


def build_analysis_readiness_assessment(
    *,
    adata: AnnData,
    successful_steps: list[str],
    cluster_key: str,
    preprocess_context: Mapping[str, Any],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    malignancy_summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Assess whether analysis outputs are ready for downstream interpretation."""
    score = 1.0
    reasons: list[str] = []
    blockers: list[str] = []

    if not preprocess_context.get("pca_present"):
        blockers.append("preprocess_input_context.pca_present=False")
        score -= 0.35
    if cluster_key not in adata.obs:
        blockers.append(f"cluster_key_missing:{cluster_key}")
        score -= 0.35
    if "markers" in successful_steps and "rank_genes_groups" not in adata.uns:
        reasons.append("markers_step_ran_but_rank_genes_groups_missing")
        score -= 0.15
    if annotation_summary.get("needs_review_clusters", 0):
        reasons.append("annotation_review_required_clusters_present")
        score -= 0.15
    if "annotation" in successful_steps and not consensus_summary.get("final_obs_present"):
        reasons.append("annotation_consensus_not_applied")
        score -= 0.15
    if "malignancy_interpretation" in successful_steps:
        if not malignancy_summary.get("available"):
            reasons.append("malignancy_interpretation_missing")
            score -= 0.10
        elif malignancy_summary.get("review_required"):
            reasons.append("malignancy_interpretation_review_required")
            score -= 0.05

    score = float(max(0.0, min(1.0, score)))
    if blockers:
        status = "blocked"
    elif reasons:
        status = "review_required"
    else:
        status = "ready"

    return _json_safe(
        {
            "status": status,
            "score": score,
            "blockers": blockers,
            "review_reasons": reasons,
            "ready_for_downstream": status in {"ready", "review_required"},
        }
    )


def build_analysis_review_action_items(
    *,
    readiness: Mapping[str, Any],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    malignancy_summary: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Create human-review action items from analysis evidence."""
    actions: list[ReviewAction] = []
    if readiness.get("blockers"):
        actions.append(
            ReviewAction(
                priority="blocking",
                action=(
                    "Resolve missing preprocessing or clustering outputs before downstream "
                    "analysis."
                ),
                rationale="Analysis readiness has blocking failures.",
                evidence_keys=["analysis_readiness", "preprocess_input_context"],
            )
        )
    if clustering_summary.get("review_required"):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review clustering resolution and flagged clusters before final annotation.",
                rationale="Resolution evidence marked one or more clustering concerns.",
                evidence_keys=["clustering_evidence_summary"],
            )
        )
    if annotation_summary.get("needs_review_clusters", 0):
        actions.append(
            ReviewAction(
                priority="review",
                action="Manually inspect annotation clusters marked as needing review.",
                rationale="Consensus annotation found weak or conflicting evidence.",
                evidence_keys=["annotation_evidence_summary", "annotation_consensus_summary"],
            )
        )
    if not consensus_summary.get("final_obs_present"):
        actions.append(
            ReviewAction(
                priority="optional",
                action="Apply consensus labels when a final cell-type column is required.",
                rationale=(
                    "Annotation evidence exists but no final consensus obs column was detected."
                ),
                evidence_keys=["annotation_consensus_summary"],
            )
        )
    if malignancy_summary.get("enabled") and malignancy_summary.get("review_required"):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review malignant/suspect/unresolved calls before tumor downstream analysis.",
                rationale=(
                    "Malignancy interpretation combines annotation, marker, CNV, and signature "
                    "evidence and found calls that require human confirmation."
                ),
                evidence_keys=["malignancy_interpretation_summary"],
            )
        )
    return [model_to_dict(action) for action in actions]


def build_analysis_evidence_bundle(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Build the common evidence bundle for analysis."""
    readiness = summary.get("analysis_readiness", {}) if isinstance(summary, Mapping) else {}
    clustering = (
        summary.get("clustering_evidence_summary", {}) if isinstance(summary, Mapping) else {}
    )
    annotation = (
        summary.get("annotation_evidence_summary", {}) if isinstance(summary, Mapping) else {}
    )
    malignancy = (
        summary.get("malignancy_interpretation_summary", {})
        if isinstance(summary, Mapping)
        else {}
    )

    evidence_chain = [
        EvidenceItem(
            source="context",
            name="preprocess_input_context",
            value=summary.get("preprocess_input_context", {}),
            confidence=None,
            rationale="Analysis depends on preprocessing outputs such as PCA and neighbors.",
            related_keys=["preprocess_input_context"],
        ),
        EvidenceItem(
            source="metric",
            name="clustering_evidence",
            value=clustering,
            confidence=None,
            rationale="Cluster interpretability gates marker discovery and annotation.",
            related_keys=["clustering_evidence_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="annotation_evidence",
            value=annotation,
            confidence=None,
            rationale="Annotation evidence summarizes marker/reference/data-driven agreement.",
            related_keys=["annotation_evidence_summary", "annotation_consensus_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="malignancy_interpretation",
            value=malignancy,
            confidence=None,
            rationale=(
                "Malignancy interpretation summarizes optional tumor-context evidence for "
                "downstream tumor analysis."
            ),
            related_keys=["malignancy_interpretation_summary"],
        ),
    ]
    bundle = EvidenceBundle(
        module="analysis",
        stage="run_standard_analysis",
        status=str(readiness.get("status", "unknown")),
        confidence=readiness.get("score") if isinstance(readiness, Mapping) else None,
        context={"workflow_name": summary.get("workflow_name", "standard")},
        evidence_chain=evidence_chain,
        action_items=[
            ReviewAction(**item)
            for item in summary.get("review_action_items", [])
            if isinstance(item, Mapping)
        ],
        related_review_keys=sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS),
    )
    return _json_safe(model_to_dict(bundle))


def build_analysis_module_maturity_assessment(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Assess whether an analysis review summary satisfies the module contract."""
    payload = _review_payload(summary)
    required_sections = set(ANALYSIS_REQUIRED_REVIEW_SECTIONS)
    required_sections.discard("module_maturity")
    missing = sorted(required_sections - set(payload.keys()))
    issues = [f"missing_required_section:{key}" for key in missing]
    review_required = []

    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    if readiness.get("status") == "blocked":
        issues.extend(readiness.get("blockers", []))
    elif readiness.get("status") == "review_required":
        review_required.extend(readiness.get("review_reasons", []))

    if issues:
        status = "incomplete"
    elif review_required:
        status = "review_required"
    else:
        status = "complete"

    return _json_safe(
        {
            "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
            "module": "analysis",
            "status": status,
            "issues": issues,
            "review_required": review_required,
            "contract": get_analysis_module_contract(),
            "summary": (
                "Analysis review summary satisfies the benchmark module contract."
                if status == "complete"
                else "Analysis review summary is present but requires review."
                if status == "review_required"
                else "Analysis review summary does not satisfy the benchmark module contract."
            ),
        }
    )


def summarize_analysis_review_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Return a compact user-facing summary of analysis output."""
    payload = _review_payload(summary)
    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    maturity = payload.get("module_maturity", {}) if isinstance(payload, Mapping) else {}
    clustering = (
        payload.get("clustering_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    )
    annotation = (
        payload.get("annotation_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    )
    consensus = (
        payload.get("annotation_consensus_summary", {}) if isinstance(payload, Mapping) else {}
    )
    malignancy = (
        payload.get("malignancy_interpretation_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    return _json_safe(
        {
            "module": "analysis",
            "maturity_status": maturity.get("status"),
            "readiness_status": readiness.get("status"),
            "readiness_score": readiness.get("score"),
            "cluster_key": clustering.get("cluster_key"),
            "n_clusters": clustering.get("n_clusters"),
            "recommended_resolution": clustering.get("recommended_resolution"),
            "annotation_key": annotation.get("annotation_key"),
            "review_table_rows": annotation.get("review_table_rows"),
            "needs_review_clusters": annotation.get("needs_review_clusters"),
            "final_key": consensus.get("final_key"),
            "n_final_labels": consensus.get("n_final_labels"),
            "mean_confidence": consensus.get("mean_confidence"),
            "malignancy_enabled": malignancy.get("enabled"),
            "n_malignant": malignancy.get("n_malignant"),
            "n_suspect_malignant": malignancy.get("n_suspect_malignant"),
        }
    )


def validate_analysis_review_summary(
    summary: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """Validate analysis-specific review-summary sections."""
    payload = _review_payload(summary)
    errors: list[str] = []
    missing = sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS - set(payload.keys()))
    if missing:
        errors.append(f"Analysis review summary missing required sections: {missing}")
    bundle = payload.get("evidence_bundle")
    if not isinstance(bundle, Mapping):
        errors.append("Analysis review summary field 'evidence_bundle' must be a mapping.")
    elif bundle.get("module") != "analysis":
        errors.append("Analysis evidence_bundle.module must be 'analysis'.")
    maturity = payload.get("module_maturity")
    if not isinstance(maturity, Mapping):
        errors.append("Analysis review summary field 'module_maturity' must be a mapping.")
    elif maturity.get("module") != "analysis":
        errors.append("Analysis module_maturity.module must be 'analysis'.")
    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def validate_analysis_module_completeness(
    adata: AnnData,
    *,
    require_ready: bool = False,
    raise_on_error: bool = False,
) -> dict[str, Any]:
    """Validate that an AnnData object contains a benchmark-grade analysis result."""
    issues: list[str] = []
    warnings: list[str] = []
    analysis_ns = adata.uns.get("sclucid", {}).get("analysis", {})
    if not isinstance(analysis_ns, Mapping):
        issues.append('Missing or invalid adata.uns["sclucid"]["analysis"] namespace.')
        analysis_ns = {}

    review_summary = analysis_ns.get("review_summary")
    payload = _review_payload(review_summary) if isinstance(review_summary, Mapping) else {}
    if not payload:
        issues.append('Missing adata.uns["sclucid"]["analysis"]["review_summary"].')
        maturity = build_analysis_module_maturity_assessment({})
    else:
        issues.extend(validate_analysis_review_summary(payload))
        maturity = build_analysis_module_maturity_assessment(payload)
        if maturity.get("status") == "incomplete":
            issues.extend(maturity.get("issues", []))
        elif maturity.get("status") == "review_required":
            warnings.extend(maturity.get("review_required", []))

    if not any(key in adata.obs for key in ("leiden_clusters", "leiden")):
        issues.append("Missing canonical analysis cluster column: 'leiden_clusters' or 'leiden'.")

    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    if require_ready and readiness.get("status") != "ready":
        issues.append(f"Analysis readiness is {readiness.get('status')!r}, expected 'ready'.")

    result = {
        "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "analysis",
        "valid": len(issues) == 0,
        "status": "valid" if not issues else "invalid",
        "issues": list(dict.fromkeys(str(item) for item in issues)),
        "warnings": list(dict.fromkeys(str(item) for item in warnings)),
        "maturity": maturity,
        "summary": summarize_analysis_review_summary(payload) if payload else {},
    }
    if result["issues"] and raise_on_error:
        raise ValueError("; ".join(result["issues"]))
    return _json_safe(result)


def _review_payload(summary: Mapping[str, Any]) -> Mapping[str, Any]:
    """Return the canonical review payload from flat or mirrored summaries."""
    if not isinstance(summary, Mapping):
        return {}
    data = summary.get("data")
    if isinstance(data, Mapping):
        return data
    return summary


def _json_safe(value: Any) -> Any:
    return sanitize_for_hdf5(value)


__all__ = [
    "ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION",
    "ANALYSIS_REQUIRED_REVIEW_SECTIONS",
    "ANALYSIS_TRACE_SCHEMA_VERSION",
    "build_analysis_module_maturity_assessment",
    "build_malignancy_interpretation_summary",
    "enrich_analysis_review_summary",
    "get_analysis_module_contract",
    "summarize_analysis_review_summary",
    "validate_analysis_module_completeness",
    "validate_analysis_review_summary",
]
