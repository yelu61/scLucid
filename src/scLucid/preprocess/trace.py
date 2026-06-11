"""Review-summary enrichment for the preprocessing workflow."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from anndata import AnnData

from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict
from scLucid.utils.contracts import _review_payload

PREPROCESS_TRACE_SCHEMA_VERSION = "1.0"
PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

PREPROCESS_REQUIRED_REVIEW_SECTIONS = {
    "preprocess_schema_version",
    "applied_parameter_summary",
    "layer_transition_summary",
    "step_evidence_summary",
    "tumor_aware_batch_correction_warnings",
    "hvg_selection_evidence_summary",
    "downstream_analysis_recommendations",
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
    summary["qc_input_context"] = qc_input_context
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
    summary["step_evidence_summary"] = build_step_evidence_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
    )
    summary["tumor_aware_batch_correction_warnings"] = tumor_warnings
    summary["hvg_selection_evidence_summary"] = hvg_summary
    summary["downstream_analysis_recommendations"] = downstream
    summary["preprocess_readiness"] = readiness
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
        "layer_transition_key": "layer_transition_summary",
        "step_evidence_key": "step_evidence_summary",
        "qc_input_key": "qc_input_context",
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
            name="layer_transition_summary",
            value=summary.get("layer_transition_summary", {}),
            rationale="Tracks expression-layer and embedding transitions across preprocessing.",
            related_keys=["layer_transition_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="step_evidence_summary",
            value=summary.get("step_evidence_summary", {}),
            rationale="Records status, inputs, outputs, parameters, and review flags for each preprocessing step.",
            related_keys=["step_evidence_summary"],
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
            "step_evidence": summary.get("step_evidence_summary", {}),
        },
        related_review_keys=[
            "applied_parameter_summary",
            "layer_transition_summary",
            "step_evidence_summary",
            "hvg_selection_evidence_summary",
            "tumor_aware_batch_correction_warnings",
            "downstream_analysis_recommendations",
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
    step_evidence = payload.get("step_evidence_summary")
    parameter_summary = payload.get("applied_parameter_summary")
    hvg_summary = payload.get("hvg_selection_evidence_summary")
    readiness = payload.get("preprocess_readiness", {})
    evidence_bundle = payload.get("evidence_bundle")
    qc_context = payload.get("qc_input_context", {})

    issues: list[str] = []
    if missing_sections:
        issues.append(
            "Missing required preprocessing review sections: " + ", ".join(missing_sections)
        )
    if not isinstance(layer_summary, Mapping):
        issues.append("layer_transition_summary must be present.")
    if not isinstance(step_evidence, Mapping):
        issues.append("step_evidence_summary must be present.")
    elif not isinstance(step_evidence.get("steps"), list):
        issues.append("step_evidence_summary.steps must be present.")
    if not isinstance(parameter_summary, Mapping):
        issues.append("applied_parameter_summary must be present.")
    if not isinstance(hvg_summary, Mapping):
        issues.append("hvg_selection_evidence_summary must be present.")
    if not isinstance(readiness, Mapping) or "status" not in readiness:
        issues.append("preprocess_readiness assessment must be present.")
    if not isinstance(evidence_bundle, Mapping) or evidence_bundle.get("module") != "preprocess":
        issues.append("evidence_bundle must be present and identify module='preprocess'.")
    if not isinstance(qc_context, Mapping):
        issues.append("qc_input_context must be present.")

    review_required: list[str] = []
    if isinstance(readiness, Mapping) and readiness.get("status") != "ready":
        review_required.append(f"preprocess_readiness.status={readiness.get('status')}")
    if isinstance(qc_context, Mapping) and not qc_context.get("available"):
        review_required.append("qc_input_context.available=False")
    if isinstance(qc_context, Mapping) and not qc_context.get("counts_layer_present"):
        review_required.append("qc_input_context.counts_layer_present=False")
    if isinstance(step_evidence, Mapping):
        review_steps = step_evidence.get("review_required_steps", [])
        if review_steps:
            review_required.append("step_evidence_summary.review_required_steps=" + ",".join(map(str, review_steps)))

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
    layers = payload.get("layer_transition_summary", {}) if isinstance(payload, Mapping) else {}
    step_evidence = payload.get("step_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    params = payload.get("applied_parameter_summary", {}) if isinstance(payload, Mapping) else {}
    qc_context = payload.get("qc_input_context", {}) if isinstance(payload, Mapping) else {}
    downstream = (
        payload.get("downstream_analysis_recommendations", {})
        if isinstance(payload, Mapping)
        else {}
    )

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
            "step_status_counts": step_evidence.get("status_counts"),
            "review_required_steps": step_evidence.get("review_required_steps"),
            "hvg_status": hvg.get("status"),
            "n_hvg_selected": hvg.get("n_hvg_selected"),
            "hvg_input_layer": hvg.get("input_layer"),
            "actual_n_pcs": pca.get("actual_n_pcs"),
            "umap_computed": graph.get("umap_computed"),
            "downstream_status": downstream.get("status"),
            "ready_for_analysis": downstream.get("ready_for_analysis"),
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
    if not isinstance(step_evidence, Mapping):
        errors.append("Preprocessing review summary field 'step_evidence_summary' must be a mapping.")
    elif not isinstance(step_evidence.get("steps"), (list, dict)):
        errors.append("Preprocessing step_evidence_summary.steps must be a list or dict.")
    maturity = summary.get("module_maturity")
    if not isinstance(maturity, Mapping):
        errors.append("Preprocessing review summary field 'module_maturity' must be a mapping.")
    elif maturity.get("module") != "preprocess":
        errors.append("Preprocessing module_maturity.module must be 'preprocess'.")
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
    text = str(tissue_type or "").lower()
    return any(token in text for token in ["tumor", "cancer", "malignan"])


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
    if hasattr(value, "item"):
        return value.item()
    return value


__all__ = [
    "PREPROCESS_REQUIRED_REVIEW_SECTIONS",
    "PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION",
    "PREPROCESS_STABLE_ENTRYPOINTS",
    "PREPROCESS_TRACE_SCHEMA_VERSION",
    "build_applied_parameter_summary",
    "build_downstream_analysis_recommendations",
    "build_hvg_selection_evidence_summary",
    "build_layer_transition_summary",
    "build_step_evidence_summary",
    "build_preprocess_evidence_bundle",
    "build_preprocess_module_maturity_assessment",
    "build_preprocess_readiness_assessment",
    "build_preprocess_review_action_items",
    "build_qc_input_context",
    "build_tumor_aware_batch_correction_warnings",
    "enrich_preprocessing_review_summary",
    "get_preprocess_module_contract",
    "summarize_preprocess_review_summary",
    "validate_preprocess_module_completeness",
    "validate_preprocessing_review_summary",
]
