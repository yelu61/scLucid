"""Review-summary enrichment for the preprocessing workflow."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from anndata import AnnData

from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict

PREPROCESS_TRACE_SCHEMA_VERSION = "1.0"

PREPROCESS_REQUIRED_REVIEW_SECTIONS = {
    "preprocess_schema_version",
    "applied_parameter_summary",
    "layer_transition_summary",
    "tumor_aware_batch_correction_warnings",
    "hvg_selection_evidence_summary",
    "downstream_analysis_recommendations",
    "preprocess_readiness",
    "review_action_items",
    "evidence_bundle",
}


def enrich_preprocessing_review_summary(
    summary: dict[str, Any],
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    tissue_type: str,
    keep_intermediate_layers: bool,
) -> dict[str, Any]:
    """Add benchmark-grade review fields to a preprocessing summary."""
    summary = dict(summary)
    hvg_summary = build_hvg_selection_evidence_summary(adata, config, successful_steps)
    tumor_warnings = build_tumor_aware_batch_correction_warnings(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        tissue_type=tissue_type,
    )
    downstream = build_downstream_analysis_recommendations(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
    )
    readiness = build_preprocess_readiness_assessment(
        adata=adata,
        downstream_recommendations=downstream,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
    )
    actions = build_preprocess_review_action_items(
        readiness=readiness,
        downstream_recommendations=downstream,
        tumor_warnings=tumor_warnings,
        hvg_summary=hvg_summary,
    )

    summary["preprocess_schema_version"] = PREPROCESS_TRACE_SCHEMA_VERSION
    summary["applied_parameter_summary"] = build_applied_parameter_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
    )
    summary["layer_transition_summary"] = build_layer_transition_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        keep_intermediate_layers=keep_intermediate_layers,
    )
    summary["tumor_aware_batch_correction_warnings"] = tumor_warnings
    summary["hvg_selection_evidence_summary"] = hvg_summary
    summary["downstream_analysis_recommendations"] = downstream
    summary["preprocess_readiness"] = readiness
    summary["review_action_items"] = actions
    summary["evidence_bundle"] = build_preprocess_evidence_bundle(summary)
    return _json_safe(summary)


def build_applied_parameter_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Summarize the effective preprocessing parameters that were applied."""
    hvg_meta = _preprocess_namespace(adata).get("hvg", {})
    integration = _preprocess_namespace(adata).get("integration", {}).get("workflow", {})
    return _json_safe(
        {
            "normalization": {
                "executed": "normalization" in successful_steps,
                "method": config.normalization.method,
                "target_sum": config.normalization.target_sum,
                "input_layer": config.counts_layer,
                "output_layer": config.normalized_layer,
                "update_X": config.normalization.update_X,
            },
            "hvg_selection": {
                "executed": "hvg_selection" in successful_steps,
                "method": config.hvg.method,
                "flavor": config.hvg.flavor,
                "requested_n_top_genes": config.hvg.n_top_genes,
                "batch_key": config.hvg.batch_key,
                "sample_key": config.hvg.sample_key,
                "input_layer": hvg_meta.get("input_layer"),
                "output_key": hvg_meta.get("output_key"),
            },
            "regression": {
                "executed": "regression" in successful_steps and bool(config.scaling.vars_to_regress),
                "vars_to_regress": list(config.scaling.vars_to_regress or []),
                "input_layer": config.normalized_layer,
                "output_layer": config.regressed_layer,
            },
            "scaling": {
                "executed": "scaling" in successful_steps,
                "method": config.scaling.scale_method,
                "max_value": config.scaling.max_value,
                "output_layer": config.scaled_layer,
            },
            "pca": {
                "executed": "pca" in successful_steps,
                "requested_n_pcs": config.graph.n_pcs,
                "actual_n_pcs": int(adata.obsm["X_pca"].shape[1]) if "X_pca" in adata.obsm else None,
            },
            "batch_correction": {
                "executed": "batch_correction" in successful_steps,
                "applied": bool(integration),
                "method": integration.get("method", config.integration.method),
                "batch_key": integration.get("batch_key", config.integration.batch_key),
                "use_rep": integration.get("use_rep", config.integration.use_rep),
                "output_key": integration.get("output_key", config.integration.output_key),
            },
            "neighbors_umap": {
                "executed": "neighbors_umap" in successful_steps,
                "requested_n_neighbors": config.graph.n_neighbors,
                "requested_n_pcs": config.graph.n_pcs,
                "umap_computed": "X_umap" in adata.obsm,
            },
        }
    )


def build_layer_transition_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    keep_intermediate_layers: bool,
) -> dict[str, Any]:
    """Describe how expression data moved across layers and embeddings."""
    transitions = [
        {
            "step": "normalization",
            "executed": "normalization" in successful_steps,
            "input": f"layers['{config.counts_layer}'] or X",
            "output": f"layers['{config.normalized_layer}']",
            "output_present": config.normalized_layer in adata.layers,
        },
        {
            "step": "set_raw",
            "executed": "set_raw" in successful_steps,
            "input": f"layers['{config.normalized_layer}']",
            "output": "raw",
            "output_present": adata.raw is not None,
        },
        {
            "step": "regression",
            "executed": "regression" in successful_steps and bool(config.scaling.vars_to_regress),
            "input": f"layers['{config.normalized_layer}']",
            "output": f"layers['{config.regressed_layer}']",
            "output_present": config.regressed_layer in adata.layers,
        },
        {
            "step": "scaling",
            "executed": "scaling" in successful_steps,
            "input": f"layers['{config.regressed_layer}'] or layers['{config.normalized_layer}']",
            "output": f"layers['{config.scaled_layer}']",
            "output_present": config.scaled_layer in adata.layers,
        },
        {
            "step": "pca",
            "executed": "pca" in successful_steps,
            "input": f"layers['{config.scaled_layer}'] or layers['{config.regressed_layer}'] or layers['{config.normalized_layer}']",
            "output": "obsm['X_pca']",
            "output_present": "X_pca" in adata.obsm,
        },
        {
            "step": "batch_correction",
            "executed": "batch_correction" in successful_steps,
            "input": f"obsm['{config.integration.use_rep}']",
            "output": f"obsm['{config.integration.output_key or f'X_{config.integration.method}'}']",
            "output_present": (config.integration.output_key or f"X_{config.integration.method}") in adata.obsm
            if config.integration.method
            else False,
        },
    ]
    return _json_safe(
        {
            "keep_intermediate_layers": keep_intermediate_layers,
            "layers_present": sorted(str(layer) for layer in adata.layers.keys()),
            "obsm_present": sorted(str(key) for key in adata.obsm.keys()),
            "raw_present": adata.raw is not None,
            "transitions": transitions,
        }
    )


def build_tumor_aware_batch_correction_warnings(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    tissue_type: str,
) -> dict[str, Any]:
    """Warn when tumor heterogeneity could be affected by batch correction."""
    is_tumor = _is_tumor_context(tissue_type)
    batch_key = config.integration.batch_key
    batch_applied = "batch_correction" in successful_steps and bool(
        config.integration.method and batch_key
    )
    n_batches = _n_batches(adata, batch_key)
    warnings: list[str] = []
    if is_tumor and batch_applied:
        warnings.append(
            "Batch correction is enabled in a tumor context; review whether malignant-state, clone, patient, or microenvironment heterogeneity is being over-corrected."
        )
    if is_tumor and batch_applied and n_batches and n_batches > 1:
        warnings.append(
            f"Tumor data include {n_batches} batch/sample groups for batch key {batch_key!r}; inspect embeddings before and after correction."
        )
    return {
        "enabled": is_tumor,
        "tissue_type": tissue_type,
        "batch_correction_applied": batch_applied,
        "method": config.integration.method if batch_applied else None,
        "batch_key": batch_key if batch_applied else None,
        "n_batches": n_batches,
        "warnings": warnings,
    }


def build_hvg_selection_evidence_summary(
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Summarize evidence supporting HVG selection."""
    hvg_meta = _preprocess_namespace(adata).get("hvg", {})
    output_key = hvg_meta.get("output_key")
    n_hvg = hvg_meta.get("n_hvg")
    if n_hvg is None and output_key in adata.var:
        n_hvg = int(adata.var[output_key].sum())
    selected_fraction = float(n_hvg / adata.n_vars) if n_hvg is not None and adata.n_vars else None
    status = "not_run"
    warnings: list[str] = []
    if "hvg_selection" in successful_steps:
        status = "ok"
        if not output_key or output_key not in adata.var:
            status = "review_required"
            warnings.append("HVG selection ran but no output key was found in adata.var.")
        elif not n_hvg:
            status = "review_required"
            warnings.append("HVG selection produced zero HVGs.")
        elif selected_fraction is not None and selected_fraction < 0.02:
            status = "review_required"
            warnings.append("Very small HVG fraction selected; downstream PCA/clustering may be unstable.")
    return _json_safe(
        {
            "status": status,
            "executed": "hvg_selection" in successful_steps,
            "output_key": output_key,
            "input_layer": hvg_meta.get("input_layer"),
            "method": hvg_meta.get("method", config.hvg.method),
            "flavor": config.hvg.flavor,
            "requested_n_top_genes": config.hvg.n_top_genes,
            "n_hvg_selected": n_hvg,
            "n_input_genes": int(adata.n_vars),
            "selected_fraction": selected_fraction,
            "input_stats": hvg_meta.get("input_stats", {}),
            "excluded_gene_types": hvg_meta.get("excluded_gene_types", {}),
            "warnings": warnings,
        }
    )


def build_downstream_analysis_recommendations(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    hvg_summary: Mapping[str, Any],
    tumor_warnings: Mapping[str, Any],
) -> dict[str, Any]:
    """Generate next-step recommendations after preprocessing."""
    blockers: list[str] = []
    recommendations: list[dict[str, Any]] = []
    if "X_pca" not in adata.obsm:
        blockers.append("PCA embedding obsm['X_pca'] is missing.")
    if "normalization" in successful_steps and config.normalized_layer not in adata.layers:
        blockers.append(f"Normalized layer {config.normalized_layer!r} is missing.")
    if hvg_summary.get("status") == "review_required":
        recommendations.append(
            {
                "target": "hvg_selection",
                "priority": "review",
                "recommendation": "Review HVG selection before clustering or trajectory analysis.",
                "rationale": "; ".join(hvg_summary.get("warnings", [])),
            }
        )
    integration_key = config.integration.output_key or f"X_{config.integration.method}"
    if config.integration.method and integration_key in adata.obsm:
        recommendations.append(
            {
                "target": "batch_corrected_embedding",
                "priority": "required",
                "recommendation": f"Use obsm['{integration_key}'] or document why raw PCA is preferred.",
                "rationale": "Batch correction output is available and should be considered for graph construction and clustering.",
            }
        )
    elif "X_pca" in adata.obsm:
        recommendations.append(
            {
                "target": "embedding",
                "priority": "required",
                "recommendation": "Use obsm['X_pca'] as the primary representation for downstream graph construction.",
                "rationale": "PCA is available and no batch-corrected embedding was detected.",
            }
        )
    if tumor_warnings.get("warnings"):
        recommendations.append(
            {
                "target": "tumor_batch_review",
                "priority": "review",
                "recommendation": "Compare tumor embeddings before and after batch correction.",
                "rationale": "Tumor-aware preprocessing warnings were generated.",
            }
        )
    status = "blocked" if blockers else ("review_required" if any(r["priority"] == "review" for r in recommendations) else "ready")
    return {
        "ready_for_analysis": not blockers,
        "status": status,
        "blockers": blockers,
        "recommendations": recommendations,
    }


def build_preprocess_readiness_assessment(
    *,
    adata: AnnData,
    downstream_recommendations: Mapping[str, Any],
    hvg_summary: Mapping[str, Any],
    tumor_warnings: Mapping[str, Any],
) -> dict[str, Any]:
    """Assess whether preprocessing output is ready for analysis."""
    blockers = list(downstream_recommendations.get("blockers", []))
    review_reasons = list(hvg_summary.get("warnings", [])) + list(tumor_warnings.get("warnings", []))
    if adata.n_obs == 0 or adata.n_vars == 0:
        blockers.append("Preprocessed AnnData is empty.")
    if blockers:
        status = "blocked"
    elif review_reasons:
        status = "review_required"
    else:
        status = "ready"
    score = max(0, 100 - min(60, 30 * len(blockers)) - min(30, 8 * len(review_reasons)))
    return {
        "status": status,
        "score": score,
        "blockers": blockers,
        "review_reasons": review_reasons,
        "output_health": {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "has_normalized_layer": any(str(layer).endswith("normalized") for layer in adata.layers.keys()),
            "has_pca": "X_pca" in adata.obsm,
            "has_umap": "X_umap" in adata.obsm,
        },
    }


def build_preprocess_review_action_items(
    *,
    readiness: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    tumor_warnings: Mapping[str, Any],
    hvg_summary: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Create prioritized preprocessing review actions."""
    actions: list[dict[str, Any]] = []
    for blocker in readiness.get("blockers", []):
        actions.append(
            {
                "priority": "blocking",
                "action": "Resolve preprocessing output issue before analysis.",
                "rationale": str(blocker),
                "evidence_key": "preprocess_readiness.blockers",
            }
        )
    for item in downstream_recommendations.get("recommendations", []):
        if item.get("priority") in {"required", "review"}:
            actions.append(
                {
                    "priority": item.get("priority", "review"),
                    "action": item.get("recommendation", "Review preprocessing output."),
                    "rationale": item.get("rationale", ""),
                    "evidence_key": f"downstream_analysis_recommendations.{item.get('target')}",
                }
            )
    for warning in tumor_warnings.get("warnings", []):
        actions.append(
            {
                "priority": "review",
                "action": "Document tumor-aware batch-correction decision.",
                "rationale": str(warning),
                "evidence_key": "tumor_aware_batch_correction_warnings.warnings",
            }
        )
    if hvg_summary.get("status") == "review_required":
        actions.append(
            {
                "priority": "review",
                "action": "Inspect HVG evidence summary before downstream analysis.",
                "rationale": "; ".join(hvg_summary.get("warnings", [])),
                "evidence_key": "hvg_selection_evidence_summary",
            }
        )
    if not actions:
        actions.append(
            {
                "priority": "optional",
                "action": "Archive preprocessing review summary with analysis outputs.",
                "rationale": "No blocking or required preprocessing review items were detected.",
                "evidence_key": "review_summary",
            }
        )
    priority_order = {"blocking": 0, "required": 1, "review": 2, "optional": 3}
    actions.sort(key=lambda item: priority_order.get(str(item.get("priority")), 9))
    return _json_safe(actions)


def build_preprocess_evidence_bundle(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Convert preprocessing review fields into the shared EvidenceBundle schema."""
    readiness = summary.get("preprocess_readiness", {})
    confidence = None
    if isinstance(readiness, Mapping) and isinstance(readiness.get("score"), (int, float)):
        confidence = max(0.0, min(1.0, float(readiness["score"]) / 100.0))
    evidence_items = [
        EvidenceItem(
            source="metric",
            name="layer_transition_summary",
            value=summary.get("layer_transition_summary", {}),
            rationale="Tracks expression-layer and embedding transitions across preprocessing.",
            related_keys=["layer_transition_summary"],
        ),
        EvidenceItem(
            source="metric",
            name="hvg_selection_evidence_summary",
            value=summary.get("hvg_selection_evidence_summary", {}),
            rationale="Summarizes HVG selection method, input layer, selected genes, and warnings.",
            related_keys=["hvg_selection_evidence_summary"],
        ),
        EvidenceItem(
            source="warning",
            name="tumor_aware_batch_correction_warnings",
            value=summary.get("tumor_aware_batch_correction_warnings", {}),
            rationale="Flags tumor-specific risks around batch correction.",
            related_keys=["tumor_aware_batch_correction_warnings"],
        ),
        EvidenceItem(
            source="downstream",
            name="downstream_analysis_recommendations",
            value=summary.get("downstream_analysis_recommendations", {}),
            rationale="Preprocessing-to-analysis handoff guidance.",
            related_keys=["downstream_analysis_recommendations"],
        ),
    ]
    actions = [
        ReviewAction(
            priority=item.get("priority", "review"),
            action=str(item.get("action", "")),
            rationale=str(item.get("rationale", "")),
            evidence_keys=[str(item.get("evidence_key"))] if item.get("evidence_key") else [],
        )
        for item in summary.get("review_action_items", [])
        if isinstance(item, Mapping)
    ]
    bundle = EvidenceBundle(
        module="preprocess",
        stage="run_preprocessing",
        status=str(readiness.get("status", "unknown")) if isinstance(readiness, Mapping) else "unknown",
        confidence=confidence,
        context={
            "steps_executed": list(summary.get("steps_executed", [])),
            "data_shape": dict(summary.get("data_shape", {})),
        },
        evidence_chain=evidence_items,
        action_items=actions,
        reproducibility={
            "workflow": "run_preprocessing",
            "applied_parameters": summary.get("applied_parameter_summary", {}),
        },
        related_review_keys=[
            "applied_parameter_summary",
            "layer_transition_summary",
            "hvg_selection_evidence_summary",
            "tumor_aware_batch_correction_warnings",
            "downstream_analysis_recommendations",
            "preprocess_readiness",
            "review_action_items",
        ],
    )
    return model_to_dict(bundle)


def validate_preprocessing_review_summary(
    summary: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """Validate preprocessing-specific review-summary sections."""
    errors: list[str] = []
    missing = sorted(PREPROCESS_REQUIRED_REVIEW_SECTIONS - set(summary.keys()))
    if missing:
        errors.append(f"Preprocessing review summary missing required sections: {missing}")
    bundle = summary.get("evidence_bundle")
    if not isinstance(bundle, Mapping):
        errors.append("Preprocessing review summary field 'evidence_bundle' must be a mapping.")
    elif bundle.get("module") != "preprocess":
        errors.append("Preprocessing evidence_bundle.module must be 'preprocess'.")
    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def _preprocess_namespace(adata: AnnData) -> dict[str, Any]:
    return adata.uns.get("sclucid", {}).get("preprocess", {})


def _is_tumor_context(tissue_type: Any) -> bool:
    text = str(tissue_type or "").lower()
    return any(token in text for token in ["tumor", "cancer", "malignan"])


def _n_batches(adata: AnnData, batch_key: Any) -> int | None:
    if isinstance(batch_key, str) and batch_key in adata.obs:
        return int(adata.obs[batch_key].nunique())
    if isinstance(batch_key, list):
        return max((int(adata.obs[key].nunique()) for key in batch_key if key in adata.obs), default=0)
    return None


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


__all__ = [
    "PREPROCESS_REQUIRED_REVIEW_SECTIONS",
    "PREPROCESS_TRACE_SCHEMA_VERSION",
    "build_applied_parameter_summary",
    "build_downstream_analysis_recommendations",
    "build_hvg_selection_evidence_summary",
    "build_layer_transition_summary",
    "build_preprocess_evidence_bundle",
    "build_preprocess_readiness_assessment",
    "build_preprocess_review_action_items",
    "build_tumor_aware_batch_correction_warnings",
    "enrich_preprocessing_review_summary",
    "validate_preprocessing_review_summary",
]
