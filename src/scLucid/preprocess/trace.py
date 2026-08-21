"""Review-summary enrichment for the preprocessing workflow."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from anndata import AnnData

from scLucid.utils.contracts import _review_payload

from ..utils.context import is_tumor_context as _shared_is_tumor_context
from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict

PREPROCESS_TRACE_SCHEMA_VERSION = "1.0"
PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

PREPROCESS_REQUIRED_REVIEW_SECTIONS = {
    "preprocess_schema_version",
    "applied_parameter_summary",
    "normalization_decision_policy",
    "preprocess_layer_contract",
    "layer_transition_summary",
    "layer_transition_table",
    "preprocess_decision_summary",
    "preprocess_reviewer_table",
    "preprocess_method_semantics",
    "step_evidence_summary",
    "tumor_aware_batch_correction_warnings",
    "hvg_selection_evidence_summary",
    "downstream_analysis_recommendations",
    "analysis_handoff_readiness",
    "preprocess_readiness",
    "review_action_items",
    "evidence_bundle",
    "qc_input_context",
    "module_maturity",
}

PREPROCESS_STABLE_ENTRYPOINTS = (
    "scLucid.preprocess.run_preprocessing",
    "scLucid.preprocess.normalize_data",
    "scLucid.preprocess.find_hvgs",
    "scLucid.preprocess.scale_data",
    "scLucid.preprocess.batch_correction",
)

PREPROCESS_EXPECTED_OUTPUTS = (
    'adata.uns["sclucid"]["preprocess"]["gene_filtering"]',
    "adata.layers['normalized']",
    "adata.var['highly_variable']",
    "adata.obsm['X_pca']",
    "adata.obsm['X_umap']",
    'adata.uns["sclucid"]["preprocess"]["review_summary"]',
)


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
    qc_input_context = build_qc_input_context(adata)
    hvg_summary = build_hvg_selection_evidence_summary(adata, config, successful_steps)
    tumor_warnings = build_tumor_aware_batch_correction_warnings(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        tissue_type=tissue_type,
    )
    summary["preprocess_schema_version"] = PREPROCESS_TRACE_SCHEMA_VERSION
    summary["qc_input_context"] = qc_input_context
    summary["applied_parameter_summary"] = build_applied_parameter_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
    )
    normalization_policy = build_normalization_decision_policy(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
    )
    summary["normalization_decision_policy"] = normalization_policy
    layer_transition_summary = build_layer_transition_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        keep_intermediate_layers=keep_intermediate_layers,
    )
    layer_transition_table = build_layer_transition_table(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        keep_intermediate_layers=keep_intermediate_layers,
    )
    summary["layer_transition_summary"] = layer_transition_summary
    summary["layer_transition_table"] = layer_transition_table
    summary["preprocess_layer_contract"] = build_preprocess_layer_contract(
        adata=adata,
        config=config,
        layer_transition_table=layer_transition_table,
        keep_intermediate_layers=keep_intermediate_layers,
    )
    summary["step_evidence_summary"] = build_step_evidence_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
    )
    method_semantics = build_preprocess_method_semantics(adata)
    summary["preprocess_method_semantics"] = method_semantics
    summary["tumor_aware_batch_correction_warnings"] = tumor_warnings
    summary["hvg_selection_evidence_summary"] = hvg_summary
    downstream = build_downstream_analysis_recommendations(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
    )
    summary["downstream_analysis_recommendations"] = downstream
    handoff = build_analysis_handoff_readiness(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        downstream_recommendations=downstream,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
    )
    summary["analysis_handoff_readiness"] = handoff
    decision_summary = build_preprocess_decision_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
        downstream_recommendations=downstream,
        layer_contract=summary["preprocess_layer_contract"],
        step_evidence=summary["step_evidence_summary"],
        normalization_policy=normalization_policy,
        method_semantics=method_semantics,
    )
    summary["preprocess_decision_summary"] = decision_summary
    summary["preprocess_reviewer_table"] = build_preprocess_reviewer_table(
        decision_summary,
        method_semantics=method_semantics,
    )
    readiness = build_preprocess_readiness_assessment(
        adata=adata,
        downstream_recommendations=downstream,
        hvg_summary=hvg_summary,
        tumor_warnings=tumor_warnings,
    )
    summary["preprocess_readiness"] = readiness
    actions = build_preprocess_review_action_items(
        readiness=readiness,
        downstream_recommendations=downstream,
        tumor_warnings=tumor_warnings,
        hvg_summary=hvg_summary,
    )
    summary["review_action_items"] = actions
    summary["evidence_bundle"] = build_preprocess_evidence_bundle(summary)
    summary["module_maturity"] = build_preprocess_module_maturity_assessment(summary)
    return _json_safe(summary)


def get_preprocess_module_contract() -> dict[str, Any]:
    """Return the frozen preprocessing module maturity contract."""
    return {
        "schema_version": PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "preprocess",
        "stable_entrypoints": list(PREPROCESS_STABLE_ENTRYPOINTS),
        "required_review_sections": sorted(PREPROCESS_REQUIRED_REVIEW_SECTIONS),
        "expected_outputs": list(PREPROCESS_EXPECTED_OUTPUTS),
        "canonical_namespace": 'adata.uns["sclucid"]["preprocess"]',
        "readiness_key": "preprocess_readiness",
        "layer_contract_key": "preprocess_layer_contract",
        "normalization_policy_key": "normalization_decision_policy",
        "layer_transition_key": "layer_transition_summary",
        "layer_transition_table_key": "layer_transition_table",
        "decision_summary_key": "preprocess_decision_summary",
        "reviewer_table_key": "preprocess_reviewer_table",
        "method_semantics_key": "preprocess_method_semantics",
        "step_evidence_key": "step_evidence_summary",
        "qc_input_key": "qc_input_context",
        "analysis_handoff_key": "analysis_handoff_readiness",
    }


def build_qc_input_context(adata: AnnData) -> dict[str, Any]:
    """Summarize the QC state consumed by preprocessing."""
    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    review = _review_payload(qc_ns.get("review_summary", {})) if isinstance(qc_ns, Mapping) else {}
    readiness = review.get("qc_readiness", {}) if isinstance(review, Mapping) else {}
    filtering = review.get("filtering_summary", {}) if isinstance(review, Mapping) else {}
    maturity = review.get("module_maturity", {}) if isinstance(review, Mapping) else {}

    return _json_safe(
        {
            "available": bool(review),
            "qc_readiness_status": readiness.get("status") if isinstance(readiness, Mapping) else None,
            "qc_readiness_score": readiness.get("score") if isinstance(readiness, Mapping) else None,
            "qc_maturity_status": maturity.get("status") if isinstance(maturity, Mapping) else None,
            "initial_cells": filtering.get("initial_cells") if isinstance(filtering, Mapping) else None,
            "post_qc_cells": int(adata.n_obs),
            "counts_layer_present": "counts" in adata.layers,
            "required_obs_metrics_present": {
                key: key in adata.obs
                for key in ("n_genes_by_counts", "total_counts", "pct_counts_mt")
            },
            "review_reasons": (
                readiness.get("review_reasons", [])
                if isinstance(readiness, Mapping)
                else []
            ),
        }
    )


def build_applied_parameter_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Summarize the effective preprocessing parameters that were applied."""
    hvg_meta = _preprocess_namespace(adata).get("hvg", {})
    norm_meta = _preprocess_namespace(adata).get("normalization", {})
    integration = _preprocess_namespace(adata).get("integration", {}).get("workflow", {})
    return _json_safe(
        {
            "normalization": {
                "executed": "normalization" in successful_steps,
                "method": config.normalization.method,
                "log_transformed": norm_meta.get("log_transformed")
                if isinstance(norm_meta, Mapping)
                else None,
                "target_sum": config.normalization.target_sum,
                "input_layer": config.counts_layer,
                "output_layer": config.normalized_layer,
                "update_X": config.normalization.update_X,
            },
            "gene_filtering": {
                "executed": "gene_filtering" in successful_steps,
                "min_cells_per_gene": getattr(config, "min_cells_per_gene", None),
                "metadata": _preprocess_namespace(adata).get("gene_filtering", {}),
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
                "zero_center": config.scaling.scale_method == "zscore",
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
            "step": "gene_filtering",
            "executed": "gene_filtering" in successful_steps,
            "input": f"layers['{config.counts_layer}'] or X",
            "output": "var subset",
            "output_present": "gene_filtering" in _preprocess_namespace(adata),
        },
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


def build_layer_transition_table(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    keep_intermediate_layers: bool,
) -> list[dict[str, Any]]:
    """Build a reviewer-facing table of layer and embedding semantics."""
    has_counts = config.counts_layer in adata.layers
    normalized_present = config.normalized_layer in adata.layers
    regressed_present = config.regressed_layer in adata.layers
    scaled_present = config.scaled_layer in adata.layers
    integration_key = config.integration.output_key or (
        f"X_{config.integration.method}" if config.integration.method else None
    )

    rows = [
        {
            "step": "input",
            "executed": True,
            "input_layer": f"layers['{config.counts_layer}']" if has_counts else "adata.X",
            "output_slot": "adata.X",
            "adata_X_semantics_before": "raw_counts_expected",
            "adata_X_semantics_after": "raw_counts_or_workflow_input",
            "raw_semantics": "unset_or_source_defined",
            "output_present": True,
            "review_required": not has_counts,
            "risk_note": ""
            if has_counts
            else f"layers['{config.counts_layer}'] missing; preprocessing fell back to adata.X.",
        },
        {
            "step": "normalization",
            "executed": "normalization" in successful_steps,
            "input_layer": f"layers['{config.counts_layer}'] or adata.X",
            "output_slot": f"layers['{config.normalized_layer}']; adata.X if update_X",
            "adata_X_semantics_before": "raw_counts_or_workflow_input",
            "adata_X_semantics_after": "normalized_counts"
            if config.normalization.update_X and "normalization" in successful_steps
            else "unchanged",
            "raw_semantics": "unchanged",
            "output_present": normalized_present,
            "review_required": "normalization" in successful_steps and not normalized_present,
            "risk_note": ""
            if normalized_present or "normalization" not in successful_steps
            else f"Normalization ran but layers['{config.normalized_layer}'] is missing.",
        },
        {
            "step": "set_raw",
            "executed": "set_raw" in successful_steps,
            "input_layer": f"layers['{config.normalized_layer}']",
            "output_slot": "adata.raw",
            "adata_X_semantics_before": "normalized_counts_or_current_X",
            "adata_X_semantics_after": "unchanged",
            "raw_semantics": "normalized_full_gene_snapshot" if adata.raw is not None else "unset",
            "output_present": adata.raw is not None,
            "review_required": "set_raw" in successful_steps and adata.raw is None,
            "risk_note": ""
            if adata.raw is not None or "set_raw" not in successful_steps
            else "set_raw ran but adata.raw is missing.",
        },
        {
            "step": "hvg_selection",
            "executed": "hvg_selection" in successful_steps,
            "input_layer": f"layers['{config.normalized_layer}']",
            "output_slot": "var['highly_variable'] or configured hvg key",
            "adata_X_semantics_before": "normalized_counts",
            "adata_X_semantics_after": "unchanged",
            "raw_semantics": "unchanged",
            "output_present": "highly_variable" in adata.var
            or bool(_preprocess_namespace(adata).get("hvg", {}).get("output_key") in adata.var),
            "review_required": False,
            "risk_note": "",
        },
        {
            "step": "scaling",
            "executed": "scaling" in successful_steps,
            "input_layer": f"layers['{config.regressed_layer}'] or layers['{config.normalized_layer}']",
            "output_slot": f"layers['{config.scaled_layer}']",
            "adata_X_semantics_before": "normalized_or_hvg_subset",
            "adata_X_semantics_after": "scaled_or_unchanged_depending_config",
            "raw_semantics": "unchanged",
            "output_present": scaled_present,
            "review_required": "scaling" in successful_steps and not scaled_present,
            "risk_note": ""
            if scaled_present or "scaling" not in successful_steps
            else f"Scaling ran but layers['{config.scaled_layer}'] is missing.",
        },
        {
            "step": "pca",
            "executed": "pca" in successful_steps,
            "input_layer": f"layers['{config.scaled_layer}'] or layers['{config.regressed_layer}'] or layers['{config.normalized_layer}']",
            "output_slot": "obsm['X_pca']",
            "adata_X_semantics_before": "unchanged",
            "adata_X_semantics_after": "unchanged",
            "raw_semantics": "unchanged",
            "output_present": "X_pca" in adata.obsm,
            "review_required": "pca" in successful_steps and "X_pca" not in adata.obsm,
            "risk_note": "" if "X_pca" in adata.obsm or "pca" not in successful_steps else "PCA ran but obsm['X_pca'] is missing.",
        },
        {
            "step": "batch_correction",
            "executed": "batch_correction" in successful_steps,
            "input_layer": f"obsm['{config.integration.use_rep}']",
            "output_slot": f"obsm['{integration_key}']" if integration_key else "",
            "adata_X_semantics_before": "unchanged",
            "adata_X_semantics_after": "unchanged",
            "raw_semantics": "unchanged",
            "output_present": bool(integration_key and integration_key in adata.obsm),
            "review_required": bool("batch_correction" in successful_steps and integration_key and integration_key not in adata.obsm),
            "risk_note": ""
            if not ("batch_correction" in successful_steps and integration_key and integration_key not in adata.obsm)
            else f"Batch correction ran but obsm['{integration_key}'] is missing.",
        },
        {
            "step": "neighbors_umap",
            "executed": "neighbors_umap" in successful_steps,
            "input_layer": "obsm['X_pca'] or corrected representation",
            "output_slot": "uns['neighbors']; obsp connectivities/distances; obsm['X_umap']",
            "adata_X_semantics_before": "unchanged",
            "adata_X_semantics_after": "unchanged",
            "raw_semantics": "unchanged",
            "output_present": "neighbors" in adata.uns and "X_umap" in adata.obsm,
            "review_required": "neighbors_umap" in successful_steps
            and not ("neighbors" in adata.uns and "X_umap" in adata.obsm),
            "risk_note": ""
            if "neighbors_umap" not in successful_steps or ("neighbors" in adata.uns and "X_umap" in adata.obsm)
            else "Graph/UMAP step ran but neighbors or UMAP output is missing.",
        },
    ]
    canonical_order = {
        "input": (0, "counts"),
        "normalization": (1, "normalized"),
        "set_raw": (2, "raw"),
        "hvg_selection": (3, "HVG"),
        "scaling": (4, "scaled"),
        "pca": (5, "PCA"),
        "batch_correction": (6, "batch_correction"),
        "neighbors_umap": (7, "graph"),
    }
    for row in rows:
        order, stage = canonical_order.get(str(row.get("step")), (99, str(row.get("step"))))
        row["canonical_order"] = order
        row["canonical_stage"] = stage
        row["canonical_flow"] = "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph"
        row["keep_intermediate_layers"] = keep_intermediate_layers
        row["layers_present"] = sorted(str(layer) for layer in adata.layers.keys())
        row["obsm_present"] = sorted(str(key) for key in adata.obsm.keys())
        row["regressed_layer_present"] = regressed_present
    return _json_safe(rows)


def build_preprocess_layer_contract(
    *,
    adata: AnnData,
    config: Any,
    layer_transition_table: list[Mapping[str, Any]],
    keep_intermediate_layers: bool,
) -> dict[str, Any]:
    """Build the canonical preprocessing layer contract for reviewers."""
    qc_ambient_contract = (
        adata.uns.get("sclucid", {}).get("qc", {}).get("ambient_layer_contract", {})
    )
    ambient_recommended = (
        qc_ambient_contract.get("recommended_preprocess_counts_layer")
        if isinstance(qc_ambient_contract, Mapping)
        else None
    )
    counts_layer = (
        ambient_recommended
        if ambient_recommended and ambient_recommended in adata.layers
        else config.counts_layer
        if config.counts_layer in adata.layers
        else None
    )
    row_by_stage = {
        str(row.get("canonical_stage")): row
        for row in layer_transition_table
        if isinstance(row, Mapping)
    }
    canonical_stages = [
        ("counts", "input", True),
        ("normalized", "normalization", True),
        ("raw", "set_raw", True),
        ("HVG", "hvg_selection", True),
        ("scaled", "scaling", True),
        ("PCA", "pca", True),
        ("graph", "neighbors_umap", bool(getattr(config, "run_neighbors", True))),
    ]
    stages: list[dict[str, Any]] = []
    review_flags: list[str] = []
    for stage, step, expected in canonical_stages:
        row = row_by_stage.get(stage, {})
        executed = bool(row.get("executed", False))
        output_present = bool(row.get("output_present", False))
        review_required = bool(row.get("review_required", False))
        if review_required:
            review_flags.append(f"{step}: {row.get('risk_note') or 'review required'}")
        status = (
            "complete"
            if output_present and (executed or step == "input")
            else "skipped"
            if not expected and not executed
            else "review_required"
            if review_required
            else "not_run"
            if not executed
            else "missing_output"
        )
        stages.append(
            {
                "stage": stage,
                "step": step,
                "expected": expected,
                "status": status,
                "input_layer": row.get("input_layer"),
                "output_slot": row.get("output_slot"),
                "output_present": output_present,
                "adata_X_semantics_after": row.get("adata_X_semantics_after"),
                "raw_semantics": row.get("raw_semantics"),
                "risk_note": row.get("risk_note", ""),
            }
        )

    missing_required = [
        item["stage"]
        for item in stages
        if item["expected"] and item["status"] in {"missing_output", "review_required"}
    ]
    return _json_safe(
        {
            "schema_version": "preprocess_layer_contract_v1",
            "canonical_flow": "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph",
            "recommended_counts_layer": counts_layer,
            "counts_layer_source": "ambient_layer_contract"
            if counts_layer and counts_layer == ambient_recommended
            else "workflow_config"
            if counts_layer
            else "adata.X_fallback",
            "normalized_layer": config.normalized_layer,
            "raw_source_layer": config.normalized_layer,
            "hvg_input_layer": config.normalized_layer,
            "scaled_layer": config.scaled_layer,
            "pca_key": "X_pca",
            "graph_keys": ["neighbors", "connectivities", "distances", "X_umap"],
            "keep_intermediate_layers": keep_intermediate_layers,
            "stage_contracts": stages,
            "layers_present": sorted(str(layer) for layer in adata.layers.keys()),
            "obsm_present": sorted(str(key) for key in adata.obsm.keys()),
            "raw_present": adata.raw is not None,
            "review_required": bool(review_flags or missing_required),
            "review_flags": review_flags,
            "missing_required_stages": missing_required,
            "ambient_layer_contract": qc_ambient_contract or None,
        }
    )


def build_normalization_decision_policy(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Describe why the normalization method/input layer is acceptable or reviewable."""
    norm_meta = _preprocess_namespace(adata).get("normalization", {})
    params = norm_meta.get("params", {}) if isinstance(norm_meta, Mapping) else {}
    input_stats = norm_meta.get("input_stats", {}) if isinstance(norm_meta, Mapping) else {}
    output_stats = norm_meta.get("output_stats", {}) if isinstance(norm_meta, Mapping) else {}
    ambient_contract = (
        adata.uns.get("sclucid", {}).get("qc", {}).get("ambient_layer_contract", {})
    )
    ambient_recommended = (
        ambient_contract.get("recommended_preprocess_counts_layer")
        if isinstance(ambient_contract, Mapping)
        else None
    )

    method = str(params.get("method", config.normalization.method))
    input_layer = str(norm_meta.get("input_layer", config.normalization.input_layer))
    zero_frac = input_stats.get("zero_frac") if isinstance(input_stats, Mapping) else None
    min_val = input_stats.get("min") if isinstance(input_stats, Mapping) else None
    review_reasons: list[str] = []
    recommended_method = "standard"
    recommended_input_layer = ambient_recommended or config.counts_layer
    source = "workflow_default"

    if ambient_recommended:
        source = "qc_ambient_layer_contract"
        if input_layer != ambient_recommended:
            review_reasons.append(
                f"QC recommends ambient-corrected counts layer {ambient_recommended!r}, but normalization used {input_layer!r}."
            )

    if isinstance(min_val, (int, float)) and min_val < 0:
        review_reasons.append("Normalization input contains negative values; count-based normalization may be inappropriate.")
    if isinstance(zero_frac, (int, float)) and zero_frac < 0.05:
        review_reasons.append("Normalization input has very few zeros; it may already be transformed.")

    if method == "pearson_residuals":
        recommended_method = "pearson_residuals_for_feature_selection_or_PCA"
        review_reasons.append(
            "Pearson residuals are useful for feature selection/PCA but should not replace count-level inputs for DE."
        )
    elif method == "scran":
        recommended_method = "scran_when_size_factor_estimation_is_needed"
    elif method == "quality_aware":
        recommended_method = "quality_aware_when_QC_metrics_drive_depth_or_composition"
        source = "quality_aware_policy"
    elif method == "clr":
        recommended_method = "clr_for_compositional_or_ADT_like_data"
        review_reasons.append("CLR is task-specific; verify it is appropriate for RNA expression workflows.")

    if "normalization" in successful_steps and config.normalized_layer not in adata.layers:
        review_reasons.append(f"Normalization ran but layer {config.normalized_layer!r} is missing.")

    status = "review_required" if review_reasons else "ok"
    confidence = 0.65 if review_reasons else 0.9
    return _json_safe(
        {
            "schema_version": "normalization_decision_policy_v1",
            "status": status,
            "executed": "normalization" in successful_steps,
            "recommended_method": recommended_method,
            "applied_method": method,
            "transformation_type": norm_meta.get("transformation_type"),
            "model_type": norm_meta.get("model_type"),
            "claim_level": norm_meta.get("claim_level"),
            "review_note": norm_meta.get("review_note"),
            "recommended_input_layer": recommended_input_layer,
            "applied_input_layer": input_layer,
            "source": source,
            "confidence": confidence,
            "review_required": bool(review_reasons),
            "review_reasons": review_reasons,
            "input_stats": input_stats,
            "output_stats": output_stats,
            "ambient_layer_contract": ambient_contract or None,
            "downstream_note": (
                "Use normalized/log expression for PCA/visualization/marker display; "
                "use counts or pseudobulk counts for formal condition DE."
            ),
        }
    )


def build_preprocess_method_semantics(adata: AnnData) -> dict[str, Any]:
    """Collect scientific claim levels for heuristic and model-based preprocess steps."""
    pp_ns = _preprocess_namespace(adata)
    rows: list[dict[str, Any]] = []

    def _add_row(
        *,
        step: str,
        source_key: str,
        meta: Mapping[str, Any] | None,
        method: Any = None,
        model_type_key: str = "model_type",
        claim_level_key: str = "claim_level",
        review_note_key: str = "review_note",
    ) -> None:
        if not isinstance(meta, Mapping):
            return
        model_type = meta.get(model_type_key)
        claim_level = meta.get(claim_level_key)
        review_note = meta.get(review_note_key)
        if model_type is None and claim_level is None and review_note is None:
            return
        claim_text = str(claim_level or "")
        rows.append(
            {
                "step": step,
                "source_key": source_key,
                "method": method if method is not None else meta.get("method"),
                "model_type": model_type,
                "claim_level": claim_level,
                "review_note": review_note,
                "review_required": any(
                    token in claim_text
                    for token in (
                        "heuristic",
                        "experimental",
                        "approximation",
                        "compositional",
                    )
                ),
            }
        )

    norm_meta = pp_ns.get("normalization", {})
    _add_row(step="normalization", source_key="normalization", meta=norm_meta)
    adaptive_meta = pp_ns.get("adaptive_normalization", {})
    _add_row(step="normalization", source_key="adaptive_normalization", meta=adaptive_meta)
    quality_policy = pp_ns.get("quality_aware_normalization_policy", {})
    _add_row(
        step="normalization",
        source_key="quality_aware_normalization_policy",
        meta=quality_policy,
        method="quality_aware",
    )
    _add_row(step="scaling", source_key="scaling", meta=pp_ns.get("scaling", {}))
    _add_row(
        step="pca",
        source_key="pca_n_pcs_selection",
        meta=pp_ns.get("pca_n_pcs_selection", {}),
        method=pp_ns.get("pca_n_pcs_selection", {}).get("method")
        if isinstance(pp_ns.get("pca_n_pcs_selection", {}), Mapping)
        else None,
    )
    _add_row(
        step="neighbors_umap",
        source_key="neighbors_optimization",
        meta=pp_ns.get("neighbors_optimization", {}),
        method="silhouette_grid_search",
    )

    integration_eval = pp_ns.get("integration", {}).get("evaluation", {})
    _add_row(
        step="batch_correction",
        source_key="integration.evaluation.kbet",
        meta=integration_eval,
        method="kbet",
        model_type_key="kbet_model_type",
        claim_level_key="kbet_claim_level",
        review_note_key="kbet_review_note",
    )
    _add_row(
        step="batch_correction",
        source_key="integration.evaluation.batch_asw",
        meta=integration_eval,
        method="batch_asw",
        model_type_key="batch_asw_model_type",
        claim_level_key="batch_asw_claim_level",
        review_note_key="batch_asw_review_note",
    )

    claim_counts: dict[str, int] = {}
    for row in rows:
        claim = str(row.get("claim_level") or "unspecified")
        claim_counts[claim] = claim_counts.get(claim, 0) + 1

    return _json_safe(
        {
            "schema_version": "preprocess_method_semantics_v1",
            "status": "review_required"
            if any(row.get("review_required") for row in rows)
            else "ok"
            if rows
            else "not_available",
            "rows": rows,
            "claim_level_counts": claim_counts,
            "review_required_steps": sorted(
                {str(row.get("step")) for row in rows if row.get("review_required")}
            ),
            "interpretation": (
                "These rows define the boundary between standard preprocessing, "
                "heuristic recommendations, approximations, and experimental transforms. "
                "Review-required rows should be interpreted as audit prompts, not automatic optima."
            ),
        }
    )


def build_preprocess_decision_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    normalization_policy: Mapping[str, Any],
    hvg_summary: Mapping[str, Any],
    tumor_warnings: Mapping[str, Any],
    downstream_recommendations: Mapping[str, Any],
    layer_contract: Mapping[str, Any],
    step_evidence: Mapping[str, Any],
    method_semantics: Mapping[str, Any],
) -> dict[str, Any]:
    """Build step-level preprocessing decisions and scientific guardrails."""

    def _step_status_from_evidence(step: str) -> str:
        steps = step_evidence.get("steps", []) if isinstance(step_evidence, Mapping) else []
        if isinstance(steps, Mapping):
            steps = steps.values()
        for item in steps:
            if isinstance(item, Mapping) and item.get("step") == step:
                return str(item.get("status", "unknown"))
        return "unknown"

    semantics_rows = (
        method_semantics.get("rows", []) if isinstance(method_semantics, Mapping) else []
    )
    semantics_by_step: dict[str, list[Mapping[str, Any]]] = {}
    if isinstance(semantics_rows, Mapping):
        semantics_rows = list(semantics_rows.values())
    for row in semantics_rows:
        if isinstance(row, Mapping):
            semantics_by_step.setdefault(str(row.get("step")), []).append(row)

    def _decision_from_status(
        *,
        step: str,
        applied: Any,
        recommended: Any,
        source: str,
        confidence: float,
        risk_note: str = "",
        review_required: bool = False,
        blocker: bool = False,
        downstream_target: str | None = None,
    ) -> dict[str, Any]:
        status = _step_status_from_evidence(step)
        if blocker or status == "missing_output":
            decision = "blocked"
            review_required = True
        elif review_required:
            decision = "review"
        elif status == "skipped":
            decision = "skipped"
        elif status == "not_run":
            decision = "not_run"
        else:
            decision = "use"
        step_semantics = semantics_by_step.get(step, [])
        return {
            "step": step,
            "decision": decision,
            "recommended": recommended,
            "applied": applied,
            "source": source,
            "confidence": confidence,
            "review_required": bool(review_required),
            "risk_note": risk_note,
            "status": status,
            "downstream_target": downstream_target,
            "method_semantics": list(step_semantics),
            "claim_levels": [
                row.get("claim_level") for row in step_semantics if row.get("claim_level")
            ],
            "semantic_review_required": any(row.get("review_required") for row in step_semantics),
        }

    normalized_present = config.normalized_layer in adata.layers
    pca_present = "X_pca" in adata.obsm
    graph_present = "neighbors" in adata.uns and "X_umap" in adata.obsm
    integration_key = config.integration.output_key or (
        f"X_{config.integration.method}" if config.integration.method else None
    )
    integration_present = bool(integration_key and integration_key in adata.obsm)
    sample_depth = downstream_recommendations.get("sample_depth_diagnostic", {})
    cell_cycle = downstream_recommendations.get("cell_cycle_regression_diagnostic", {})
    blockers = list(downstream_recommendations.get("blockers", []))

    hvg_review = hvg_summary.get("status") == "review_required"
    integration_review = bool(tumor_warnings.get("warnings")) or (
        isinstance(sample_depth, Mapping) and sample_depth.get("status") == "review_required"
    )
    regression_vars = list(getattr(config.scaling, "vars_to_regress", []) or [])
    regression_review = bool(regression_vars) or (
        isinstance(cell_cycle, Mapping) and cell_cycle.get("status") == "review_required"
    )

    primary_rep = integration_key if integration_present else "X_pca" if pca_present else None
    decisions = [
        _decision_from_status(
            step="normalization",
            applied=config.normalization.method,
            recommended=normalization_policy.get(
                "recommended_method",
                "standard_log_normalize_by_default; consider scran/pearson_residuals by dataset/task",
            ),
            source=normalization_policy.get("source", "workflow_config"),
            confidence=float(normalization_policy.get("confidence", 0.85 if normalized_present else 0.2)),
            review_required=bool(normalization_policy.get("review_required")) or not normalized_present,
            blocker=not normalized_present and "normalization" in successful_steps,
            risk_note=(
                "Normalization output is missing."
                if not normalized_present and "normalization" in successful_steps
                else "; ".join(str(item) for item in normalization_policy.get("review_reasons", []))
                or normalization_policy.get(
                    "downstream_note",
                    "Normalization choice should match downstream task; DE should use counts/pseudobulk rather than scaled data.",
                )
            ),
            downstream_target=config.normalized_layer,
        ),
        _decision_from_status(
            step="hvg_selection",
            applied={
                "method": hvg_summary.get("method"),
                "flavor": hvg_summary.get("flavor"),
                "n_hvg_selected": hvg_summary.get("n_hvg_selected"),
            },
            recommended="batch-aware HVG with protected marker/program genes when biology could be variance-sparse",
            source="hvg_selection_evidence_summary",
            confidence=0.65 if hvg_review else 0.9,
            review_required=hvg_review,
            risk_note="; ".join(str(item) for item in hvg_summary.get("warnings", []))
            or "Review marker/program preservation when studying rare states, tumor heterogeneity, or immune activation.",
            downstream_target=hvg_summary.get("output_key") or "highly_variable",
        ),
        _decision_from_status(
            step="regression",
            applied=regression_vars,
            recommended="off_by_default; enable only after confounding review",
            source="workflow_config",
            confidence=0.55 if regression_review else 0.9,
            review_required=regression_review,
            risk_note=(
                "Regression can remove biological gradients; document why these covariates are technical."
                if regression_vars
                else str(cell_cycle.get("message", ""))
                if isinstance(cell_cycle, Mapping) and cell_cycle.get("status") == "review_required"
                else ""
            ),
            downstream_target=config.regressed_layer if regression_vars else config.normalized_layer,
        ),
        _decision_from_status(
            step="scaling",
            applied=config.scaling.scale_method,
            recommended="zscore_for_PCA; retain normalized/raw layers for expression interpretation",
            source="workflow_config",
            confidence=0.85 if config.scaled_layer in adata.layers else 0.4,
            review_required="scaling" in successful_steps and config.scaled_layer not in adata.layers,
            risk_note=(
                "Scaled values are for PCA/graph, not expression-level biological interpretation."
            ),
            downstream_target=config.scaled_layer,
        ),
        _decision_from_status(
            step="pca",
            applied={"requested_n_pcs": config.graph.n_pcs, "present": pca_present},
            recommended="use PCA before graph; review n_pcs with variance/marker preservation diagnostics",
            source="workflow_config",
            confidence=0.9 if pca_present else 0.2,
            review_required=not pca_present and "pca" in successful_steps,
            blocker=not pca_present and ("pca" in successful_steps or "PCA embedding" in " ".join(blockers)),
            risk_note="PCA is the primary handoff to graph/integration; missing PCA blocks standard analysis."
            if not pca_present
            else "",
            downstream_target="X_pca",
        ),
        _decision_from_status(
            step="batch_correction",
            applied={
                "method": config.integration.method,
                "batch_key": config.integration.batch_key,
                "output_key": integration_key if integration_present else None,
            },
            recommended=(
                "use_corrected_representation_after_review"
                if integration_present
                else "keep_unintegrated_X_pca_unless_batch_effect_is_demonstrated"
            ),
            source="integration_guardrail",
            confidence=0.55 if integration_review else 0.8,
            review_required=integration_review,
            risk_note="; ".join(str(item) for item in tumor_warnings.get("warnings", []))
            or (
                sample_depth.get("message", "")
                if isinstance(sample_depth, Mapping)
                and sample_depth.get("status") == "review_required"
                else "Integration is optional; compare integrated and unintegrated embeddings when batch correction is used."
            ),
            downstream_target=integration_key if integration_present else "X_pca",
        ),
        _decision_from_status(
            step="neighbors_umap",
            applied={"n_neighbors": config.graph.n_neighbors, "n_pcs": config.graph.n_pcs},
            recommended=f"use obsm['{primary_rep}'] for graph" if primary_rep else "run PCA before graph",
            source="downstream_analysis_recommendations",
            confidence=0.85 if graph_present else 0.6 if not getattr(config, "run_neighbors", True) else 0.25,
            review_required=("neighbors_umap" in successful_steps and not graph_present),
            blocker=("neighbors_umap" in successful_steps and not graph_present),
            risk_note="Graph/UMAP missing after requested graph step."
            if "neighbors_umap" in successful_steps and not graph_present
            else "If integration was applied, verify graph uses the intended representation.",
            downstream_target=primary_rep,
        ),
    ]

    decision_counts: dict[str, int] = {}
    for item in decisions:
        decision_counts[item["decision"]] = decision_counts.get(item["decision"], 0) + 1

    return _json_safe(
        {
            "schema_version": "preprocess_decision_summary_v1",
            "canonical_flow": layer_contract.get(
                "canonical_flow",
                "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph",
            ),
            "decisions": decisions,
            "decision_counts": decision_counts,
            "review_required_steps": [
                item["step"] for item in decisions if item.get("review_required")
            ],
            "primary_downstream_representation": primary_rep,
            "risk_note": (
                "Preprocess decisions are step/representation-level guardrails. "
                "They do not change data automatically; they define what should be reviewed before analysis."
            ),
        }
    )


def build_preprocess_reviewer_table(
    decision_summary: Mapping[str, Any],
    *,
    method_semantics: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build a single reviewer table for preprocessing choices and risks."""
    decisions = decision_summary.get("decisions", []) if isinstance(decision_summary, Mapping) else []
    if isinstance(decisions, Mapping):
        decisions = decisions.values()
    semantics_rows = (
        method_semantics.get("rows", [])
        if isinstance(method_semantics, Mapping)
        else []
    )
    if isinstance(semantics_rows, Mapping):
        semantics_rows = semantics_rows.values()
    semantics_by_step: dict[str, list[Mapping[str, Any]]] = {}
    for row in semantics_rows:
        if isinstance(row, Mapping):
            semantics_by_step.setdefault(str(row.get("step")), []).append(row)
    rows = []
    for item in decisions:
        if not isinstance(item, Mapping):
            continue
        step_semantics = semantics_by_step.get(str(item.get("step")), [])
        rows.append(
            _json_safe(
                {
                    "item": item.get("step"),
                    "category": "preprocess_step",
                    "recommended_value": item.get("recommended"),
                    "applied_value": item.get("applied"),
                    "source": item.get("source"),
                    "confidence": item.get("confidence"),
                    "affected_representation": item.get("downstream_target"),
                    "preprocess_decision": item.get("decision"),
                    "review_required": bool(item.get("review_required")),
                    "biological_risk_note": item.get("risk_note", ""),
                    "status": item.get("status"),
                    "claim_levels": item.get("claim_levels", []),
                    "model_types": [
                        row.get("model_type")
                        for row in step_semantics
                        if row.get("model_type")
                    ],
                    "semantic_review_required": bool(
                        item.get("semantic_review_required")
                    ),
                    "scientific_semantics": step_semantics,
                }
            )
        )
    return rows


def build_step_evidence_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Build auditable evidence records for each preprocessing step."""
    hvg_meta = _preprocess_namespace(adata).get("hvg", {})
    integration = _preprocess_namespace(adata).get("integration", {}).get("workflow", {})
    hvg_key = hvg_meta.get("output_key", "highly_variable")
    gene_filtering_meta = _preprocess_namespace(adata).get("gene_filtering", {})
    n_hvg = hvg_meta.get("n_hvg")
    if n_hvg is None and hvg_key in adata.var:
        n_hvg = int(adata.var[hvg_key].sum())
    pca_variance_top3 = None
    if "pca" in adata.uns and isinstance(adata.uns["pca"], Mapping):
        variance_ratio = adata.uns["pca"].get("variance_ratio")
        if variance_ratio is not None:
            pca_variance_top3 = [round(float(value), 4) for value in list(variance_ratio)[:3]]

    steps = [
        {
            "step": "gene_filtering",
            "status": _step_status(
                "gene_filtering",
                successful_steps,
                output_present="gene_filtering" in _preprocess_namespace(adata),
                skipped=not getattr(config, "run_gene_filtering", True),
            ),
            "input": {
                "layer": gene_filtering_meta.get("source", config.counts_layer),
                "layer_present": config.counts_layer in adata.layers,
            },
            "output": {
                "n_genes_before": gene_filtering_meta.get("initial_genes"),
                "n_genes_after": gene_filtering_meta.get("final_genes", int(adata.n_vars)),
                "n_genes_removed": gene_filtering_meta.get("removed_genes"),
                "skipped": gene_filtering_meta.get("skipped", False),
            },
            "parameters": {
                "min_cells_per_gene": getattr(config, "min_cells_per_gene", None),
            },
            "audit_fields": [
                "applied_parameter_summary.gene_filtering",
                "layer_transition_summary.transitions.gene_filtering",
            ],
            "review_flags": _missing_output_flags(
                "gene_filtering",
                successful_steps,
                {"uns['gene_filtering']": "gene_filtering" in _preprocess_namespace(adata)},
            ),
        },
        {
            "step": "normalization",
            "status": _step_status(
                "normalization",
                successful_steps,
                output_present=config.normalized_layer in adata.layers,
            ),
            "input": {
                "layer": config.counts_layer,
                "fallback": "X",
                "layer_present": config.counts_layer in adata.layers,
            },
            "output": {
                "layer": config.normalized_layer,
                "present": config.normalized_layer in adata.layers,
                "shape": _layer_shape(adata, config.normalized_layer),
            },
            "parameters": {
                "method": config.normalization.method,
                "target_sum": config.normalization.target_sum,
                "update_X": config.normalization.update_X,
            },
            "audit_fields": [
                "applied_parameter_summary.normalization",
                "layer_transition_summary.transitions.normalization",
            ],
            "review_flags": _missing_output_flags(
                "normalization",
                successful_steps,
                {f"layers['{config.normalized_layer}']": config.normalized_layer in adata.layers},
            ),
        },
        {
            "step": "set_raw",
            "status": _step_status("set_raw", successful_steps, output_present=adata.raw is not None),
            "input": {"layer": config.normalized_layer},
            "output": {"slot": "raw", "present": adata.raw is not None},
            "parameters": {"source_layer": config.normalized_layer},
            "audit_fields": ["layer_transition_summary.transitions.set_raw"],
            "review_flags": _missing_output_flags(
                "set_raw",
                successful_steps,
                {"raw": adata.raw is not None},
            ),
        },
        {
            "step": "regression",
            "status": _step_status(
                "regression",
                successful_steps,
                output_present=(not config.scaling.vars_to_regress)
                or config.regressed_layer in adata.layers,
                skipped=not bool(config.scaling.vars_to_regress),
            ),
            "input": {"layer": config.normalized_layer},
            "output": {
                "layer": config.regressed_layer,
                "present": config.regressed_layer in adata.layers,
                "shape": _layer_shape(adata, config.regressed_layer),
            },
            "parameters": {"vars_to_regress": list(config.scaling.vars_to_regress or [])},
            "audit_fields": ["applied_parameter_summary.regression"],
            "review_flags": _missing_output_flags(
                "regression",
                successful_steps,
                {f"layers['{config.regressed_layer}']": config.regressed_layer in adata.layers},
            )
            if config.scaling.vars_to_regress
            else [],
        },
        {
            "step": "hvg_selection",
            "status": _step_status(
                "hvg_selection",
                successful_steps,
                output_present=hvg_key in adata.var,
            ),
            "input": {
                "layer": hvg_meta.get("input_layer"),
                "layer_present": hvg_meta.get("input_layer") in adata.layers
                if hvg_meta.get("input_layer")
                else None,
            },
            "output": {
                "var_key": hvg_key,
                "present": hvg_key in adata.var,
                "n_hvg_selected": n_hvg,
                "n_input_genes": int(adata.n_vars),
            },
            "parameters": {
                "method": hvg_meta.get("method", config.hvg.method),
                "flavor": config.hvg.flavor,
                "requested_n_top_genes": config.hvg.n_top_genes,
                "batch_key": config.hvg.batch_key,
                "sample_key": config.hvg.sample_key,
            },
            "audit_fields": [
                "applied_parameter_summary.hvg_selection",
                "hvg_selection_evidence_summary",
            ],
            "review_flags": _missing_output_flags(
                "hvg_selection",
                successful_steps,
                {f"var['{hvg_key}']": hvg_key in adata.var},
            )
            + ([] if n_hvg else ["HVG selection produced no selected genes."] if "hvg_selection" in successful_steps else []),
        },
        {
            "step": "subset_hvg",
            "status": _step_status(
                "subset_hvg",
                successful_steps,
                output_present=True,
            ),
            "input": {"var_key": hvg_key},
            "output": {"n_vars_after_subset": int(adata.n_vars)},
            "parameters": {"mode": "direct", "keep_raw": False},
            "audit_fields": ["hvg_selection_evidence_summary.n_hvg_selected"],
            "review_flags": [],
        },
        {
            "step": "scaling",
            "status": _step_status(
                "scaling",
                successful_steps,
                output_present=config.scaled_layer in adata.layers,
            ),
            "input": {
                "preferred_layers": [config.regressed_layer, config.normalized_layer],
            },
            "output": {
                "layer": config.scaled_layer,
                "present": config.scaled_layer in adata.layers,
                "shape": _layer_shape(adata, config.scaled_layer),
            },
            "parameters": {
                "method": config.scaling.scale_method,
                "max_value": config.scaling.max_value,
            },
            "audit_fields": [
                "applied_parameter_summary.scaling",
                "layer_transition_summary.transitions.scaling",
            ],
            "review_flags": _missing_output_flags(
                "scaling",
                successful_steps,
                {f"layers['{config.scaled_layer}']": config.scaled_layer in adata.layers},
            ),
        },
        {
            "step": "pca",
            "status": _step_status("pca", successful_steps, output_present="X_pca" in adata.obsm),
            "input": {
                "preferred_layers": [
                    config.scaled_layer,
                    config.regressed_layer,
                    config.normalized_layer,
                ],
            },
            "output": {
                "obsm_key": "X_pca",
                "present": "X_pca" in adata.obsm,
                "shape": _obsm_shape(adata, "X_pca"),
                "variance_explained_top3": pca_variance_top3,
            },
            "parameters": {
                "requested_n_pcs": config.graph.n_pcs,
                "actual_n_pcs": int(adata.obsm["X_pca"].shape[1]) if "X_pca" in adata.obsm else None,
            },
            "audit_fields": [
                "applied_parameter_summary.pca",
                "layer_transition_summary.transitions.pca",
            ],
            "review_flags": _missing_output_flags(
                "pca",
                successful_steps,
                {"obsm['X_pca']": "X_pca" in adata.obsm},
            ),
        },
        {
            "step": "batch_correction",
            "status": _step_status(
                "batch_correction",
                successful_steps,
                output_present=bool(integration)
                or not (config.integration.method and config.integration.batch_key),
                skipped=not bool(config.integration.method and config.integration.batch_key),
            ),
            "input": {"obsm_key": config.integration.use_rep, "batch_key": config.integration.batch_key},
            "output": {
                "obsm_key": integration.get(
                    "output_key",
                    config.integration.output_key or f"X_{config.integration.method}",
                )
                if config.integration.method
                else None,
                "present": (
                    integration.get(
                        "output_key",
                        config.integration.output_key or f"X_{config.integration.method}",
                    )
                    in adata.obsm
                )
                if config.integration.method
                else False,
            },
            "parameters": {
                "method": integration.get("method", config.integration.method),
                "batch_key": integration.get("batch_key", config.integration.batch_key),
                "use_rep": integration.get("use_rep", config.integration.use_rep),
            },
            "audit_fields": [
                "applied_parameter_summary.batch_correction",
                "tumor_aware_batch_correction_warnings",
            ],
            "review_flags": [],
        },
        {
            "step": "neighbors_umap",
            "status": _step_status(
                "neighbors_umap",
                successful_steps,
                output_present="neighbors" in adata.uns and "X_umap" in adata.obsm,
            ),
            "input": {"obsm_key": "X_pca"},
            "output": {
                "neighbors_present": "neighbors" in adata.uns,
                "umap_present": "X_umap" in adata.obsm,
                "umap_shape": _obsm_shape(adata, "X_umap"),
            },
            "parameters": {
                "requested_n_neighbors": config.graph.n_neighbors,
                "requested_n_pcs": config.graph.n_pcs,
                "effective_n_pcs": min(config.graph.n_pcs, adata.obsm["X_pca"].shape[1])
                if "X_pca" in adata.obsm
                else None,
                "effective_n_neighbors": min(config.graph.n_neighbors, max(2, adata.n_obs - 1)),
            },
            "audit_fields": [
                "applied_parameter_summary.neighbors_umap",
                "downstream_analysis_recommendations",
            ],
            "review_flags": _missing_output_flags(
                "neighbors_umap",
                successful_steps,
                {"neighbors": "neighbors" in adata.uns, "obsm['X_umap']": "X_umap" in adata.obsm},
            ),
        },
    ]
    status_counts: dict[str, int] = {}
    for item in steps:
        status_counts[item["status"]] = status_counts.get(item["status"], 0) + 1
    return _json_safe(
        {
            "schema_version": PREPROCESS_TRACE_SCHEMA_VERSION,
            "steps": steps,
            "status_counts": status_counts,
            "review_required_steps": [
                item["step"] for item in steps if item["status"] in {"missing_output", "review_required"}
            ],
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
    auto_decide = bool(getattr(config.integration, "auto_decide", False))
    evaluate = bool(getattr(config.integration, "evaluate", False))
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
    if batch_applied and not auto_decide:
        warnings.append(
            "Batch correction was applied without IntegrationConfig.auto_decide=True; document the batch effect evidence and biological-preservation check."
        )
    if batch_applied and not evaluate:
        warnings.append(
            "Batch correction was applied without IntegrationConfig.evaluate=True; compare integrated and unintegrated representations before interpreting clusters."
        )
    return {
        "enabled": is_tumor,
        "tissue_type": tissue_type,
        "batch_correction_applied": batch_applied,
        "method": config.integration.method if batch_applied else None,
        "batch_key": batch_key if batch_applied else None,
        "auto_decide": auto_decide if batch_applied else None,
        "evaluate": evaluate if batch_applied else None,
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
    excluded_gene_types = hvg_meta.get("excluded_gene_types", {})
    excluded_total = (
        sum(int(v) for v in excluded_gene_types.values() if isinstance(v, (int, float)))
        if isinstance(excluded_gene_types, Mapping)
        else 0
    )
    pre_exclusion_total = (int(n_hvg) + excluded_total) if n_hvg is not None else None
    excluded_fraction = (
        float(excluded_total / pre_exclusion_total)
        if pre_exclusion_total and pre_exclusion_total > 0
        else None
    )
    if "hvg_selection" in successful_steps and excluded_fraction is not None and excluded_fraction > 0.3:
        status = "review_required"
        warnings.append(
            "More than 30% of initially selected HVGs overlapped excluded technical gene types; inspect HVG/PCA diagnostics before analysis."
        )
    return _json_safe(
        {
            "status": status,
            "executed": "hvg_selection" in successful_steps,
            "output_key": output_key,
            "input_layer": hvg_meta.get("input_layer"),
            "method": hvg_meta.get("method", config.hvg.method),
            "flavor": config.hvg.flavor,
            "requested_n_top_genes": config.hvg.n_top_genes,
            "method_report": hvg_meta.get("method_report", {}),
            "n_hvg_selected": n_hvg,
            "n_input_genes": int(adata.n_vars),
            "selected_fraction": selected_fraction,
            "input_stats": hvg_meta.get("input_stats", {}),
            "excluded_gene_types": excluded_gene_types,
            "excluded_gene_type_total": excluded_total,
            "excluded_gene_type_fraction": excluded_fraction,
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
    sample_depth = build_sample_depth_diagnostic(adata, config)
    cell_cycle = build_cell_cycle_regression_diagnostic(adata, config)
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
    if sample_depth.get("status") == "review_required":
        recommendations.append(
            {
                "target": "sample_depth",
                "priority": "review",
                "recommendation": "Review sample-level sequencing depth differences before choosing integration/downsampling.",
                "rationale": sample_depth.get("message", ""),
            }
        )
    if cell_cycle.get("status") == "review_required":
        recommendations.append(
            {
                "target": "cell_cycle_regression",
                "priority": "review",
                "recommendation": "Review whether cell-cycle regression is biologically appropriate before enabling it.",
                "rationale": cell_cycle.get("message", ""),
            }
        )
    status = "blocked" if blockers else ("review_required" if any(r["priority"] == "review" for r in recommendations) else "ready")
    return {
        "ready_for_analysis": not blockers,
        "status": status,
        "blockers": blockers,
        "sample_depth_diagnostic": sample_depth,
        "cell_cycle_regression_diagnostic": cell_cycle,
        "recommendations": recommendations,
    }


def build_analysis_handoff_readiness(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    downstream_recommendations: Mapping[str, Any],
    hvg_summary: Mapping[str, Any],
    tumor_warnings: Mapping[str, Any],
) -> dict[str, Any]:
    """Declare which preprocessing outputs are safe for each analysis use."""
    blockers = list(downstream_recommendations.get("blockers", []))
    review_items: list[str] = []
    warnings: list[str] = []

    normalized_layer = getattr(config, "normalized_layer", "normalized")
    scaled_layer = getattr(config, "scaled_layer", "scaled")
    hvg_key = hvg_summary.get("output_key") or "highly_variable"

    integrated_embeddings = [
        str(key)
        for key in adata.obsm.keys()
        if str(key).startswith("X_")
        and str(key) not in {"X_pca", "X_umap"}
        and any(
            token in str(key).lower()
            for token in ("harmony", "bbknn", "scanorama", "scvi", "integrat")
        )
    ]
    configured_integration_key = getattr(getattr(config, "integration", None), "output_key", None)
    if (
        configured_integration_key
        and configured_integration_key in adata.obsm
        and str(configured_integration_key) not in integrated_embeddings
    ):
        integrated_embeddings.append(str(configured_integration_key))

    graph_representation = integrated_embeddings[0] if integrated_embeddings else "X_pca"
    if graph_representation not in adata.obsm:
        graph_representation = None

    expression_for_markers = (
        "adata.raw.X"
        if adata.raw is not None
        else f"adata.layers['{normalized_layer}']"
        if normalized_layer in adata.layers
        else None
    )
    expression_for_de = (
        f"adata.layers['{normalized_layer}']"
        if normalized_layer in adata.layers
        else expression_for_markers
    )

    if graph_representation is None:
        blockers.append(
            "No graph-ready embedding was found; expected obsm['X_pca'] or a reviewed integrated embedding."
        )
    if expression_for_markers is None:
        blockers.append("No normalized/raw expression source is available for marker review.")
    if expression_for_de is None:
        blockers.append("No unintegrated normalized expression source is available for DE handoff.")
    if hvg_summary.get("status") == "review_required":
        review_items.extend(str(item) for item in hvg_summary.get("warnings", []))
    if tumor_warnings.get("warnings"):
        review_items.extend(str(item) for item in tumor_warnings.get("warnings", []))
    if integrated_embeddings:
        warnings.append(
            "Integrated embeddings are graph/UMAP evidence only; marker, annotation, and DE claims should use unintegrated normalized/raw expression."
        )
    scaling_config = getattr(config, "scaling", None)
    if "regression" in successful_steps and getattr(scaling_config, "vars_to_regress", None):
        review_items.append(
            "Regression was configured; confirm regressed covariates are technical rather than biological signals."
        )

    status = "blocked" if blockers else ("review_required" if review_items or warnings else "ready")
    return _json_safe(
        {
            "schema_version": "analysis_handoff_readiness_v1",
            "status": status,
            "ready_for_analysis": not blockers,
            "graph_representation": graph_representation,
            "diagnostic_embedding": "X_pca" if "X_pca" in adata.obsm else None,
            "integrated_embeddings": integrated_embeddings,
            "expression_for_markers": expression_for_markers,
            "expression_for_annotation": expression_for_markers,
            "expression_for_de": expression_for_de,
            "hvg_key": hvg_key if hvg_key in adata.var else None,
            "raw_present": adata.raw is not None,
            "normalized_layer_present": normalized_layer in adata.layers,
            "scaled_layer_present": scaled_layer in adata.layers,
            "safe_uses": {
                "graph_clustering_umap": graph_representation,
                "marker_annotation_review": expression_for_markers,
                "condition_de": expression_for_de,
                "proportion_analysis": "obs cell labels plus sample metadata",
            },
            "unsafe_uses": [
                "Do not use integrated embeddings or scaled expression as the primary matrix for marker or condition-DE claims.",
                "Do not treat UMAP geometry as quantitative biological distance.",
            ],
            "blockers": blockers,
            "review_items": review_items,
            "warnings": warnings,
        }
    )


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
    review_reasons.extend(
        str(item.get("rationale", ""))
        for item in downstream_recommendations.get("recommendations", [])
        if isinstance(item, Mapping) and item.get("priority") == "review" and item.get("rationale")
    )
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


def build_sample_depth_diagnostic(adata: AnnData, config: Any) -> dict[str, Any]:
    """Flag large sample-level sequencing-depth differences for review."""
    candidate_keys = []
    batch_key = getattr(getattr(config, "integration", None), "batch_key", None)
    hvg_sample_key = getattr(getattr(config, "hvg", None), "sample_key", None)
    for key in (batch_key, hvg_sample_key, "sampleID", "sample", "batch"):
        if isinstance(key, str) and key not in candidate_keys:
            candidate_keys.append(key)
    group_key = next((key for key in candidate_keys if key in adata.obs), None)
    if group_key is None or "total_counts" not in adata.obs:
        return {
            "status": "not_available",
            "group_key": group_key,
            "message": "Sample-level depth diagnostic requires a group key and obs['total_counts'].",
        }
    medians = adata.obs.groupby(group_key, observed=False)["total_counts"].median().dropna()
    if len(medians) < 2:
        return {"status": "ok", "group_key": group_key, "n_groups": int(len(medians))}
    min_median = float(medians.min())
    max_median = float(medians.max())
    ratio = max_median / max(min_median, 1.0)
    status = "review_required" if ratio > 2.0 else "ok"
    return _json_safe(
        {
            "status": status,
            "group_key": group_key,
            "n_groups": int(len(medians)),
            "min_median_total_counts": min_median,
            "max_median_total_counts": max_median,
            "max_to_min_ratio": ratio,
            "message": (
                f"Median total_counts differs by {ratio:.2f}x across {group_key!r}; inspect PCA/UMAP and consider downsampling or integration settings."
                if status == "review_required"
                else "No large sample-level sequencing-depth disparity detected."
            ),
        }
    )


def build_cell_cycle_regression_diagnostic(adata: AnnData, config: Any) -> dict[str, Any]:
    """Provide review guidance for cell-cycle regression without enabling it by default."""
    from .scale import diagnose_cell_cycle_regression

    has_scores = {"S_score", "G2M_score"}.issubset(set(adata.obs.columns)) or "phase" in adata.obs
    vars_to_regress = set(getattr(getattr(config, "scaling", None), "vars_to_regress", []) or [])
    regress_in_scale = bool(getattr(getattr(config, "scaling", None), "regress_in_scale", False))
    configured = bool({"S_score", "G2M_score", "phase"} & vars_to_regress) or regress_in_scale
    if not has_scores:
        return {
            "status": "not_available",
            "scores_present": False,
            "regression_configured": configured,
            "message": "Cell-cycle scores were not detected.",
        }
    diagnostic = diagnose_cell_cycle_regression(
        adata,
        condition_key=getattr(config, "condition_key", None),
        batch_key=getattr(getattr(config, "integration", None), "batch_key", None)
        if isinstance(getattr(getattr(config, "integration", None), "batch_key", None), str)
        else None,
        cell_type_key="cell_type_final" if "cell_type_final" in adata.obs else None,
        record=False,
    )
    diagnostic["regression_configured"] = configured
    diagnostic["message"] = diagnostic.get("recommendation", "")
    if configured and diagnostic.get("status") == "low_risk":
        diagnostic["status"] = "ok"
        diagnostic["message"] = (
            "Cell-cycle covariates are configured for regression; document biological rationale."
        )
    return diagnostic


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
            name="preprocess_layer_contract",
            value=summary.get("preprocess_layer_contract", {}),
            rationale="Freezes the canonical counts-to-graph layer contract used by preprocessing.",
            related_keys=["preprocess_layer_contract", "layer_transition_table"],
        ),
        EvidenceItem(
            source="metric",
            name="layer_transition_summary",
            value=summary.get("layer_transition_summary", {}),
            rationale="Tracks expression-layer and embedding transitions across preprocessing.",
            related_keys=["layer_transition_summary"],
        ),
        EvidenceItem(
            source="recommendation",
            name="normalization_decision_policy",
            value=summary.get("normalization_decision_policy", {}),
            rationale="Explains the normalization method/input-layer choice and downstream expression semantics.",
            related_keys=["normalization_decision_policy"],
        ),
        EvidenceItem(
            source="output_health",
            name="step_evidence_summary",
            value=summary.get("step_evidence_summary", {}),
            rationale="Records status, inputs, outputs, parameters, and review flags for each preprocessing step.",
            related_keys=["step_evidence_summary"],
        ),
        EvidenceItem(
            source="contract",
            name="preprocess_method_semantics",
            value=summary.get("preprocess_method_semantics", {}),
            rationale="Labels preprocessing methods as standard, heuristic, approximate, or experimental for reviewer interpretation.",
            related_keys=["preprocess_method_semantics", "preprocess_reviewer_table"],
        ),
        EvidenceItem(
            source="recommendation",
            name="preprocess_decision_summary",
            value=summary.get("preprocess_decision_summary", {}),
            rationale="Converts preprocessing outputs into step-level use/review/block decisions.",
            related_keys=["preprocess_decision_summary", "preprocess_reviewer_table"],
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
        EvidenceItem(
            source="contract",
            name="analysis_handoff_readiness",
            value=summary.get("analysis_handoff_readiness", {}),
            rationale="Declares which preprocessing outputs are safe for graph, marker, annotation, DE, and tumor-analysis handoff.",
            related_keys=["analysis_handoff_readiness", "preprocess_layer_contract"],
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
            "step_evidence": summary.get("step_evidence_summary", {}),
        },
        related_review_keys=[
            "applied_parameter_summary",
            "normalization_decision_policy",
            "preprocess_layer_contract",
            "layer_transition_summary",
            "preprocess_decision_summary",
            "preprocess_reviewer_table",
            "preprocess_method_semantics",
            "step_evidence_summary",
            "hvg_selection_evidence_summary",
            "tumor_aware_batch_correction_warnings",
            "downstream_analysis_recommendations",
            "analysis_handoff_readiness",
            "preprocess_readiness",
            "review_action_items",
        ],
    )
    return model_to_dict(bundle)


def build_preprocess_module_maturity_assessment(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Assess whether preprocessing satisfies the benchmark module contract."""
    payload = _review_payload(summary)
    required_sections = set(PREPROCESS_REQUIRED_REVIEW_SECTIONS)
    required_sections.discard("module_maturity")
    missing_sections = sorted(required_sections - set(payload.keys()))

    layer_summary = payload.get("layer_transition_summary")
    layer_contract = payload.get("preprocess_layer_contract")
    step_evidence = payload.get("step_evidence_summary")
    method_semantics = payload.get("preprocess_method_semantics")
    decision_summary = payload.get("preprocess_decision_summary")
    reviewer_table = payload.get("preprocess_reviewer_table")
    parameter_summary = payload.get("applied_parameter_summary")
    normalization_policy = payload.get("normalization_decision_policy")
    hvg_summary = payload.get("hvg_selection_evidence_summary")
    readiness = payload.get("preprocess_readiness", {})
    handoff = payload.get("analysis_handoff_readiness", {})
    evidence_bundle = payload.get("evidence_bundle")
    qc_context = payload.get("qc_input_context", {})

    issues: list[str] = []
    if missing_sections:
        issues.append(
            "Missing required preprocessing review sections: " + ", ".join(missing_sections)
        )
    if not isinstance(layer_summary, Mapping):
        issues.append("layer_transition_summary must be present.")
    if not isinstance(layer_contract, Mapping):
        issues.append("preprocess_layer_contract must be present.")
    elif layer_contract.get("canonical_flow") != "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph":
        issues.append("preprocess_layer_contract.canonical_flow is not the canonical counts-to-graph flow.")
    if not isinstance(step_evidence, Mapping):
        issues.append("step_evidence_summary must be present.")
    elif not isinstance(step_evidence.get("steps"), list):
        issues.append("step_evidence_summary.steps must be present.")
    if not isinstance(method_semantics, Mapping):
        issues.append("preprocess_method_semantics must be present.")
    elif not isinstance(method_semantics.get("rows"), (list, dict)):
        issues.append("preprocess_method_semantics.rows must be present.")
    if not isinstance(decision_summary, Mapping):
        issues.append("preprocess_decision_summary must be present.")
    elif not isinstance(decision_summary.get("decisions"), list):
        issues.append("preprocess_decision_summary.decisions must be present.")
    if not isinstance(reviewer_table, list):
        issues.append("preprocess_reviewer_table must be present.")
    if not isinstance(parameter_summary, Mapping):
        issues.append("applied_parameter_summary must be present.")
    if not isinstance(normalization_policy, Mapping):
        issues.append("normalization_decision_policy must be present.")
    if not isinstance(hvg_summary, Mapping):
        issues.append("hvg_selection_evidence_summary must be present.")
    if not isinstance(readiness, Mapping) or "status" not in readiness:
        issues.append("preprocess_readiness assessment must be present.")
    if not isinstance(handoff, Mapping) or "status" not in handoff:
        issues.append("analysis_handoff_readiness assessment must be present.")
    if not isinstance(evidence_bundle, Mapping) or evidence_bundle.get("module") != "preprocess":
        issues.append("evidence_bundle must be present and identify module='preprocess'.")
    if not isinstance(qc_context, Mapping):
        issues.append("qc_input_context must be present.")

    review_required: list[str] = []
    if isinstance(readiness, Mapping) and readiness.get("status") != "ready":
        review_required.append(f"preprocess_readiness.status={readiness.get('status')}")
    if isinstance(handoff, Mapping) and handoff.get("status") != "ready":
        review_required.append(f"analysis_handoff_readiness.status={handoff.get('status')}")
    if isinstance(qc_context, Mapping) and not qc_context.get("available"):
        review_required.append("qc_input_context.available=False")
    if isinstance(qc_context, Mapping) and not qc_context.get("counts_layer_present"):
        review_required.append("qc_input_context.counts_layer_present=False")
    if isinstance(step_evidence, Mapping):
        review_steps = step_evidence.get("review_required_steps", [])
        if review_steps:
            review_required.append("step_evidence_summary.review_required_steps=" + ",".join(map(str, review_steps)))
    if isinstance(decision_summary, Mapping):
        decision_review_steps = decision_summary.get("review_required_steps", [])
        if decision_review_steps:
            review_required.append(
                "preprocess_decision_summary.review_required_steps="
                + ",".join(map(str, decision_review_steps))
            )
    if isinstance(method_semantics, Mapping):
        semantic_review_steps = method_semantics.get("review_required_steps", [])
        if semantic_review_steps:
            review_required.append(
                "preprocess_method_semantics.review_required_steps="
                + ",".join(map(str, semantic_review_steps))
            )
    if isinstance(layer_contract, Mapping) and layer_contract.get("review_required"):
        review_required.append("preprocess_layer_contract.review_required=True")
    if isinstance(normalization_policy, Mapping) and normalization_policy.get("review_required"):
        review_required.append("normalization_decision_policy.review_required=True")

    if issues:
        status = "incomplete"
    elif review_required:
        status = "review_required"
    else:
        status = "complete"

    return _json_safe(
        {
            "schema_version": PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION,
            "module": "preprocess",
            "status": status,
            "status_scope": "review_contract_completeness_only",
            "scientific_validation_status": "REVIEW",
            "core_position": "withheld_until_locked_acceptance_passes",
            "superiority_claim": "unsupported_pending_held_out_validation",
            "issues": issues,
            "review_required": review_required,
            "contract": get_preprocess_module_contract(),
            "summary": (
                "Preprocessing review summary satisfies the benchmark module contract."
                if status == "complete"
                else "Preprocessing review summary is present but requires review."
                if status == "review_required"
                else "Preprocessing review summary does not satisfy the benchmark module contract."
            ),
        }
    )


def summarize_preprocess_review_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Return a compact user-facing summary of preprocessing output."""
    payload = _review_payload(summary)
    readiness = payload.get("preprocess_readiness", {}) if isinstance(payload, Mapping) else {}
    maturity = payload.get("module_maturity", {}) if isinstance(payload, Mapping) else {}
    hvg = payload.get("hvg_selection_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    normalization_policy = (
        payload.get("normalization_decision_policy", {})
        if isinstance(payload, Mapping)
        else {}
    )
    layers = payload.get("layer_transition_summary", {}) if isinstance(payload, Mapping) else {}
    layer_contract = payload.get("preprocess_layer_contract", {}) if isinstance(payload, Mapping) else {}
    step_evidence = payload.get("step_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    decision_summary = (
        payload.get("preprocess_decision_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    method_semantics = (
        payload.get("preprocess_method_semantics", {})
        if isinstance(payload, Mapping)
        else {}
    )
    params = payload.get("applied_parameter_summary", {}) if isinstance(payload, Mapping) else {}
    qc_context = payload.get("qc_input_context", {}) if isinstance(payload, Mapping) else {}
    downstream = (
        payload.get("downstream_analysis_recommendations", {})
        if isinstance(payload, Mapping)
        else {}
    )
    handoff = payload.get("analysis_handoff_readiness", {}) if isinstance(payload, Mapping) else {}

    pca = params.get("pca", {}) if isinstance(params.get("pca"), Mapping) else {}
    graph = params.get("neighbors_umap", {}) if isinstance(params.get("neighbors_umap"), Mapping) else {}

    return _json_safe(
        {
            "module": "preprocess",
            "maturity_status": maturity.get("status"),
            "readiness_status": readiness.get("status"),
            "readiness_score": readiness.get("score"),
            "qc_input_available": qc_context.get("available"),
            "qc_readiness_status": qc_context.get("qc_readiness_status"),
            "layers_present": layers.get("layers_present"),
            "obsm_present": layers.get("obsm_present"),
            "raw_present": layers.get("raw_present"),
            "canonical_layer_flow": layer_contract.get("canonical_flow"),
            "recommended_counts_layer": layer_contract.get("recommended_counts_layer"),
            "normalization_status": normalization_policy.get("status"),
            "normalization_recommended_method": normalization_policy.get("recommended_method"),
            "normalization_applied_method": normalization_policy.get("applied_method"),
            "layer_contract_review_required": layer_contract.get("review_required"),
            "step_status_counts": step_evidence.get("status_counts"),
            "review_required_steps": step_evidence.get("review_required_steps"),
            "preprocess_decision_counts": decision_summary.get("decision_counts"),
            "decision_review_required_steps": decision_summary.get("review_required_steps"),
            "method_semantics_status": method_semantics.get("status"),
            "method_claim_level_counts": method_semantics.get("claim_level_counts"),
            "semantic_review_required_steps": method_semantics.get("review_required_steps"),
            "primary_downstream_representation": decision_summary.get(
                "primary_downstream_representation"
            ),
            "hvg_status": hvg.get("status"),
            "n_hvg_selected": hvg.get("n_hvg_selected"),
            "hvg_input_layer": hvg.get("input_layer"),
            "actual_n_pcs": pca.get("actual_n_pcs"),
            "umap_computed": graph.get("umap_computed"),
            "downstream_status": downstream.get("status"),
            "ready_for_analysis": downstream.get("ready_for_analysis"),
            "analysis_handoff_status": handoff.get("status"),
            "graph_representation": handoff.get("graph_representation"),
            "expression_for_markers": handoff.get("expression_for_markers"),
            "expression_for_de": handoff.get("expression_for_de"),
        }
    )


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
    step_evidence = summary.get("step_evidence_summary")
    normalization_policy = summary.get("normalization_decision_policy")
    if not isinstance(normalization_policy, Mapping):
        errors.append(
            "Preprocessing review summary field 'normalization_decision_policy' must be a mapping."
        )
    if not isinstance(step_evidence, Mapping):
        errors.append("Preprocessing review summary field 'step_evidence_summary' must be a mapping.")
    elif not isinstance(step_evidence.get("steps"), (list, dict)):
        errors.append("Preprocessing step_evidence_summary.steps must be a list or dict.")
    method_semantics = summary.get("preprocess_method_semantics")
    if not isinstance(method_semantics, Mapping):
        errors.append(
            "Preprocessing review summary field 'preprocess_method_semantics' must be a mapping."
        )
    elif not isinstance(method_semantics.get("rows"), (list, dict)):
        errors.append("Preprocessing preprocess_method_semantics.rows must be a list or dict.")
    decision_summary = summary.get("preprocess_decision_summary")
    if not isinstance(decision_summary, Mapping):
        errors.append(
            "Preprocessing review summary field 'preprocess_decision_summary' must be a mapping."
        )
    elif not isinstance(decision_summary.get("decisions"), (list, dict)):
        errors.append("Preprocessing preprocess_decision_summary.decisions must be a list or dict.")
    reviewer_table = summary.get("preprocess_reviewer_table")
    if not isinstance(reviewer_table, list) or not reviewer_table:
        errors.append(
            "Preprocessing review summary field 'preprocess_reviewer_table' must be a non-empty list."
        )
    else:
        required_reviewer_columns = {
            "item",
            "category",
            "recommended_value",
            "applied_value",
            "source",
            "confidence",
            "affected_representation",
            "preprocess_decision",
            "review_required",
            "biological_risk_note",
        }
        for idx, row in enumerate(reviewer_table):
            if not isinstance(row, Mapping):
                errors.append(f"Preprocessing preprocess_reviewer_table row {idx} must be a mapping.")
                continue
            missing_columns = sorted(required_reviewer_columns - set(row.keys()))
            if missing_columns:
                errors.append(
                    f"Preprocessing preprocess_reviewer_table row {idx} missing columns: {missing_columns}"
                )
    layer_table = summary.get("layer_transition_table")
    layer_contract = summary.get("preprocess_layer_contract")
    if not isinstance(layer_contract, Mapping):
        errors.append("Preprocessing review summary field 'preprocess_layer_contract' must be a mapping.")
    else:
        if layer_contract.get("canonical_flow") != "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph":
            errors.append("Preprocessing preprocess_layer_contract.canonical_flow is invalid.")
        if not isinstance(layer_contract.get("stage_contracts"), list) or not layer_contract.get("stage_contracts"):
            errors.append(
                "Preprocessing preprocess_layer_contract.stage_contracts must be a non-empty list."
            )
    if not isinstance(layer_table, list) or not layer_table:
        errors.append("Preprocessing review summary field 'layer_transition_table' must be a non-empty list.")
    else:
        required_columns = {
            "step",
            "executed",
            "input_layer",
            "output_slot",
            "adata_X_semantics_before",
            "adata_X_semantics_after",
            "raw_semantics",
            "output_present",
            "review_required",
            "risk_note",
        }
        for idx, row in enumerate(layer_table):
            if not isinstance(row, Mapping):
                errors.append(f"Preprocessing layer_transition_table row {idx} must be a mapping.")
                continue
            missing_columns = sorted(required_columns - set(row.keys()))
            if missing_columns:
                errors.append(
                    f"Preprocessing layer_transition_table row {idx} missing columns: {missing_columns}"
                )
    maturity = summary.get("module_maturity")
    if not isinstance(maturity, Mapping):
        errors.append("Preprocessing review summary field 'module_maturity' must be a mapping.")
    elif maturity.get("module") != "preprocess":
        errors.append("Preprocessing module_maturity.module must be 'preprocess'.")
    handoff = summary.get("analysis_handoff_readiness")
    if not isinstance(handoff, Mapping):
        errors.append(
            "Preprocessing review summary field 'analysis_handoff_readiness' must be a mapping."
        )
    else:
        required_handoff_fields = {
            "status",
            "ready_for_analysis",
            "graph_representation",
            "expression_for_markers",
            "expression_for_de",
            "safe_uses",
            "unsafe_uses",
        }
        missing_handoff = sorted(required_handoff_fields - set(handoff.keys()))
        if missing_handoff:
            errors.append(
                f"Preprocessing analysis_handoff_readiness missing fields: {missing_handoff}"
            )
    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def validate_preprocess_module_completeness(
    adata: AnnData,
    *,
    require_ready: bool = False,
    raise_on_error: bool = False,
) -> dict[str, Any]:
    """Validate that an AnnData object contains a benchmark-grade preprocessing result."""
    issues: list[str] = []
    warnings: list[str] = []
    pp_ns = adata.uns.get("sclucid", {}).get("preprocess", {})
    if not isinstance(pp_ns, Mapping):
        issues.append('Missing or invalid adata.uns["sclucid"]["preprocess"] namespace.')
        pp_ns = {}

    review_summary = pp_ns.get("review_summary")
    payload = _review_payload(review_summary) if isinstance(review_summary, Mapping) else {}
    if not payload:
        issues.append('Missing adata.uns["sclucid"]["preprocess"]["review_summary"].')
        maturity = build_preprocess_module_maturity_assessment({})
    else:
        issues.extend(validate_preprocessing_review_summary(payload))
        maturity = build_preprocess_module_maturity_assessment(payload)
        if maturity.get("status") == "incomplete":
            issues.extend(maturity.get("issues", []))
        elif maturity.get("status") == "review_required":
            warnings.extend(maturity.get("review_required", []))

    if "normalized" not in adata.layers:
        issues.append("Missing required preprocessing layer: 'normalized'.")
    if "X_pca" not in adata.obsm:
        issues.append("Missing required preprocessing embedding: 'X_pca'.")
    if "highly_variable" not in adata.var:
        warnings.append("Missing canonical HVG column: 'highly_variable'.")

    readiness = payload.get("preprocess_readiness", {}) if isinstance(payload, Mapping) else {}
    if require_ready and readiness.get("status") != "ready":
        issues.append(f"Preprocess readiness is {readiness.get('status')!r}, expected 'ready'.")

    result = {
        "schema_version": PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "preprocess",
        "valid": len(issues) == 0,
        "status": "valid" if not issues else "invalid",
        "issues": list(dict.fromkeys(str(item) for item in issues)),
        "warnings": list(dict.fromkeys(str(item) for item in warnings)),
        "maturity": maturity,
        "summary": summarize_preprocess_review_summary(payload) if payload else {},
    }
    if result["issues"] and raise_on_error:
        raise ValueError("; ".join(result["issues"]))
    return _json_safe(result)


def _preprocess_namespace(adata: AnnData) -> dict[str, Any]:
    return adata.uns.get("sclucid", {}).get("preprocess", {})


def _is_tumor_context(tissue_type: Any) -> bool:
    return _shared_is_tumor_context(tissue_type)


def _n_batches(adata: AnnData, batch_key: Any) -> int | None:
    if isinstance(batch_key, str) and batch_key in adata.obs:
        return int(adata.obs[batch_key].nunique())
    if isinstance(batch_key, list):
        return max((int(adata.obs[key].nunique()) for key in batch_key if key in adata.obs), default=0)
    return None


def _step_status(
    step: str,
    successful_steps: list[str],
    *,
    output_present: bool,
    skipped: bool = False,
) -> str:
    if skipped:
        return "skipped"
    if step not in successful_steps:
        return "not_run"
    return "complete" if output_present else "missing_output"


def _missing_output_flags(
    step: str,
    successful_steps: list[str],
    outputs: Mapping[str, bool],
) -> list[str]:
    if step not in successful_steps:
        return []
    return [f"Expected output missing: {name}." for name, present in outputs.items() if not present]


def _layer_shape(adata: AnnData, layer: str) -> list[int] | None:
    if layer not in adata.layers:
        return None
    return [int(value) for value in adata.layers[layer].shape]


def _obsm_shape(adata: AnnData, key: str) -> list[int] | None:
    if key not in adata.obsm:
        return None
    return [int(value) for value in adata.obsm[key].shape]


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if hasattr(value, "item"):
        try:
            return value.item()
        except ValueError:
            # Multi-element array-like (e.g. numpy array with size > 1); convert to list.
            return _json_safe(np.asarray(value).tolist())
    return value


__all__ = [
    "PREPROCESS_REQUIRED_REVIEW_SECTIONS",
    "PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION",
    "PREPROCESS_STABLE_ENTRYPOINTS",
    "PREPROCESS_TRACE_SCHEMA_VERSION",
    "build_applied_parameter_summary",
    "build_downstream_analysis_recommendations",
    "build_analysis_handoff_readiness",
    "build_hvg_selection_evidence_summary",
    "build_layer_transition_summary",
    "build_layer_transition_table",
    "build_normalization_decision_policy",
    "build_preprocess_decision_summary",
    "build_preprocess_layer_contract",
    "build_preprocess_method_semantics",
    "build_step_evidence_summary",
    "build_preprocess_evidence_bundle",
    "build_preprocess_module_maturity_assessment",
    "build_preprocess_readiness_assessment",
    "build_preprocess_reviewer_table",
    "build_preprocess_review_action_items",
    "build_qc_input_context",
    "build_tumor_aware_batch_correction_warnings",
    "enrich_preprocessing_review_summary",
    "get_preprocess_module_contract",
    "summarize_preprocess_review_summary",
    "validate_preprocess_module_completeness",
    "validate_preprocessing_review_summary",
]
