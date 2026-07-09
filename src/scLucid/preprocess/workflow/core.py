"""
High-level preprocessing workflow functions for single-cell RNA-seq data.

This module provides flexible, memory-efficient preprocessing workflows with
fine-grained step control, backend abstraction, progress tracking, and error recovery.
"""

import logging
import traceback
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Callable, Dict, List, Literal, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from anndata import AnnData

from ...utils import (
    PartialResultManager,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    export_review_summary,
    get_progress_bar,
    normalize_review_summary,
    validate_review_summary_schema,
)
from ...utils.context import is_tumor_context
from ..config import WorkflowConfig
from ..gene_biotype import annotate_gene_biotypes, filter_genes_by_biotype
from ..hvg import evaluate_hvg_stability, find_hvgs, select_and_audit_hvgs, select_hvg_sets
from ..integrate import batch_correction, decide_integration
from ..neighbors import run_embedding_pipeline
from ..normalize import normalize_data
from ..scale import regress_out, scale_data
from ..trace import (
    _json_safe,
    enrich_preprocessing_review_summary,
    validate_preprocessing_review_summary,
)

log = logging.getLogger(__name__)

__all__ = [
    "run_preprocessing",
    "run_iterative_preprocessing",
    "WORKFLOW_STEPS",
    "WorkflowError",
]

# Define workflow steps for flexible execution.
# Order follows sc-best-practices: normalize -> set_raw -> select HVGs on
# normalized data -> optional regression on HVG subset -> scale -> PCA ->
# optional batch correction -> neighbors/UMAP.
WORKFLOW_STEPS = [
    "gene_filtering",
    "normalization",
    "set_raw",
    "hvg_selection",
    "subset_hvg",
    "regression",
    "scaling",
    "pca",
    "batch_correction",
    "neighbors_umap",
]

# Keep for backward compatibility
PartialWorkflowResult = PartialResultManager


def run_preprocessing(
    adata: AnnData,
    config: Optional[WorkflowConfig] = None,
    save_dir: Optional[str] = None,
    force: bool = False,
    *,
    # Backward compatibility
    results_dir: Optional[str] = None,
    # Step control
    steps: Optional[Sequence[str]] = None,
    skip_steps: Optional[Sequence[str]] = None,
    # Memory optimization
    inplace: bool = False,
    keep_intermediate_layers: bool = True,
    # Tumor-aware hint
    tissue_type: str = "unknown",
    # Progress tracking
    show_progress: bool = True,
    progress_desc: str = "Preprocessing",
    # Error recovery
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: Literal["raise", "skip", "save"] = "raise",
    # Resume from partial
    resume_from: Optional[str] = None,
    # Custom processing hooks
    custom_pre_step: Optional[Callable[[AnnData, str], AnnData]] = None,
    custom_post_step: Optional[Callable[[AnnData, str], AnnData]] = None,
) -> AnnData:
    """
    Run the preprocessing workflow with flexible step control, progress tracking, and error recovery.

    This is the main entry point for preprocessing single-cell RNA-seq data.
    It provides a 9-step pipeline that can be customized via configuration and
    runtime parameters.

    Args:
        adata: Input AnnData object. Should have raw counts in layers["counts"] or .X.
        config: Preprocessing configuration. If None, uses default WorkflowConfig().
        save_dir: Directory to save results and plots. If None, uses config.save_dir.
        force: Whether to force recomputation of cached steps.
        results_dir: Deprecated. Use save_dir or config.save_dir instead.
        steps: Specific steps to run (default: all). See WORKFLOW_STEPS for valid names.
        skip_steps: Steps to skip (alternative to specifying 'steps').
        inplace: If True, modify adata in-place to save memory. Use with caution.
        keep_intermediate_layers: If False, delete intermediate layers to save memory.
        show_progress: If True, show progress bar for workflow steps.
        progress_desc: Description for progress bar.
        error_recovery: If True, enable error recovery mode.
        recovery_save_dir: Directory to save partial results on error. Required if error_recovery=True.
        on_error: How to handle errors: "raise" (default), "skip" (skip failed step), or "save" (save partial).
        resume_from: Path to partial results directory to resume from.
        custom_pre_step: Optional callable to run before each step. Signature: (adata, step_name) -> adata.
        custom_post_step: Optional callable to run after each step. Signature: (adata, step_name) -> adata.

    Returns:
        AnnData object with preprocessing completed.

    Raises:
        ValueError: If invalid step names are provided.
        KeyError: If required layers are missing.
        WorkflowError: If a step fails and on_error="raise".

    Example:
        >>> # Standard analysis with progress bar
        >>> adata = run_preprocessing(adata, show_progress=True)

        >>> # Error recovery mode
        >>> adata = run_preprocessing(
        ...     adata,
        ...     error_recovery=True,
        ...     recovery_save_dir="./recovery",
        ...     on_error="save"
        ... )

        >>> # Resume from partial results
        >>> adata = run_preprocessing(
        ...     adata,
        ...     resume_from="./recovery",
        ...     show_progress=True
        ... )

        >>> # Skip regression and use all genes (no HVG subsetting)
        >>> adata = run_preprocessing(
        ...     adata,
        ...     skip_steps=["regression", "subset_hvg"],
        ...     config=WorkflowConfig(quick(n_top_genes=None))
        ... )

        >>> # Custom step: add QC filter between normalization and HVG
        >>> def custom_filter(adata, step_name):
        ...     if step_name == "hvg_selection":
        ...         adata = adata[adata.obs.n_genes > 500].copy()
        ...     return adata
        ...
        >>> adata = run_preprocessing(adata, custom_pre_step=custom_filter)

        >>> # Memory-efficient: inplace modification and cleanup
        >>> adata = run_preprocessing(
        ...     adata,
        ...     inplace=True,
        ...     keep_intermediate_layers=False
        ... )
    """
    active_config = _prepare_workflow_config(config)

    # Validate error recovery settings
    if error_recovery and on_error == "save" and not recovery_save_dir:
        raise ValueError(
            "recovery_save_dir is required when error_recovery=True and on_error='save'"
        )

    # Handle resume from partial results
    completed_steps: List[str] = []
    if resume_from:
        manager = PartialResultManager(resume_from)
        adata, checkpoint, _ = manager.load()
        completed_steps = checkpoint.completed_steps
        log.info(f"Resumed from partial results. Completed steps: {completed_steps}")

    # Validate step names
    steps_to_run = _resolve_steps(steps, skip_steps, active_config)
    invalid_steps = set(steps_to_run) - set(WORKFLOW_STEPS)
    if invalid_steps:
        raise ValueError(
            f"Invalid step names: {invalid_steps}. " f"Valid steps are: {WORKFLOW_STEPS}"
        )

    # Skip already completed steps if resuming
    if resume_from and completed_steps:
        steps_to_run = [s for s in steps_to_run if s not in completed_steps]
        log.info(f"Steps to run: {steps_to_run}")

    log.info("=" * 60)
    log.info("=== Starting Preprocessing Workflow ===")
    log.info("=" * 60)
    log.info(f"Steps to run: {steps_to_run}")
    log.info(f"Inplace mode: {inplace}")
    log.info(f"Keep intermediate layers: {keep_intermediate_layers}")
    log.info(f"Show progress: {show_progress}")
    log.info(f"Error recovery: {error_recovery}")

    # Handle inplace vs copy
    if inplace and not resume_from:
        log.warning(
            "Inplace mode enabled. Original adata WILL be modified. "
            "Ensure you have a backup if needed."
        )
    elif not resume_from:
        adata = adata.copy()

    # Handle save_dir priority: explicit > config > deprecated results_dir
    effective_save_dir = save_dir
    if effective_save_dir is None and active_config.save_dir:
        effective_save_dir = active_config.save_dir
    if results_dir is not None:
        log.warning("results_dir is deprecated. Use save_dir or config.save_dir instead.")
        effective_save_dir = results_dir

    results_path = Path(effective_save_dir) if effective_save_dir else None
    if results_path:
        results_path.mkdir(parents=True, exist_ok=True)

    # Initialize progress bar
    step_iterator = get_progress_bar(
        steps_to_run,
        desc=progress_desc,
        enabled=show_progress,
        total=len(steps_to_run),
        unit="step",
    )

    # Track execution
    current_step = None
    successful_steps: List[str] = []

    try:
        for step_name in step_iterator:
            current_step = step_name

            # --- 1. Gene filtering ---
            if step_name == "gene_filtering":
                adata = _run_step(
                    adata,
                    "gene_filtering",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _run_gene_filtering_step(a, active_config),
                )
                successful_steps.append(step_name)

            # --- 2. Normalization ---
            elif step_name == "normalization":
                adata = _run_step(
                    adata,
                    "normalization",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: normalize_data(
                        a,
                        config=active_config.normalization,
                        force=force,
                        save_dir=str(results_path / "normalization") if results_path else None,
                    ),
                )
                successful_steps.append(step_name)

            # --- 2. Set .raw with normalized data BEFORE regression ---
            elif step_name == "set_raw":
                adata = _run_step(
                    adata,
                    "set_raw",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _run_gene_biotype_filtering_step(
                        _set_raw_layer(a, active_config),
                        active_config,
                        initial_genes_before_step=int(a.n_vars),
                        current_stage="after_raw",
                    ),
                )
                successful_steps.append(step_name)

            # --- 3. Regression (Optional) ---
            elif step_name == "regression":
                if active_config.scaling.vars_to_regress:
                    adata = _run_step(
                        adata,
                        "regression",
                        custom_pre_step,
                        custom_post_step,
                        lambda a: regress_out(
                            a,
                            config=active_config.scaling,
                            input_layer=active_config.normalized_layer,
                            output_layer=active_config.regressed_layer,
                        ),
                    )
                    successful_steps.append(step_name)

                    # Optionally clean up normalized layer to save memory
                    if (
                        not keep_intermediate_layers
                        and active_config.normalized_layer in adata.layers
                    ):
                        del adata.layers[active_config.normalized_layer]
                        log.info(f"Removed intermediate layer: {active_config.normalized_layer}")
                else:
                    log.info("Step: Skipping regression (no vars_to_regress).")
                    successful_steps.append(step_name)

            # --- 4. HVG Selection ---
            elif step_name == "hvg_selection":
                adata = _run_step(
                    adata,
                    "hvg_selection",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: find_hvgs(
                        a,
                        config=active_config.hvg,
                        input_layer=_get_hvg_input_layer(a, active_config),
                        force=force,
                        save_dir=str(results_path / "hvg") if results_path else None,
                    ),
                )
                successful_steps.append(step_name)

            # --- 5. Subset to HVGs ---
            elif step_name == "subset_hvg":
                adata = _run_step(
                    adata,
                    "subset_hvg",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _subset_to_hvgs(a, active_config, keep_intermediate_layers),
                )
                successful_steps.append(step_name)

            # --- 6. Scaling ---
            elif step_name == "scaling":
                adata = _run_step(
                    adata,
                    "scaling",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _run_scaling_step(a, active_config),
                )
                successful_steps.append(step_name)

            # --- 7. PCA ---
            elif step_name == "pca":
                adata = _run_step(
                    adata,
                    "pca",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _run_pca(a, active_config, results_path),
                )
                successful_steps.append(step_name)

            # --- 8. Integration/Batch Correction ---
            elif step_name == "batch_correction":
                if active_config.integration.method and active_config.integration.batch_key:
                    adata = _run_step(
                        adata,
                        "batch_correction",
                        custom_pre_step,
                        custom_post_step,
                        lambda a: batch_correction(
                            a,
                            config=active_config.integration,
                            save_dir=str(results_path / "integration") if results_path else None,
                        ),
                    )
                    successful_steps.append(step_name)
                else:
                    log.info("Step: Skipping batch correction (no method or batch_key).")
                    successful_steps.append(step_name)

            # --- 9. Neighbors & UMAP ---
            elif step_name == "neighbors_umap":
                adata = _run_step(
                    adata,
                    "neighbors_umap",
                    custom_pre_step,
                    custom_post_step,
                    lambda a: _run_neighbors_umap(a, active_config, results_path),
                )
                successful_steps.append(step_name)

    except Exception as e:
        error_msg = f"Workflow failed at step '{current_step}': {str(e)}"
        log.error(error_msg)
        log.error(traceback.format_exc())

        if error_recovery and on_error in ["raise", "save"]:
            # Save partial results
            save_dir = recovery_save_dir or (
                str(results_path / "recovery") if results_path else "./recovery"
            )
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(
                completed_steps=successful_steps,
                failed_step=current_step,
                error_message=str(e),
            )
            manager.save(adata, checkpoint, active_config)

            if on_error == "save":
                log.warning(f"Workflow failed but partial results saved to: {save_dir}")
                log.warning(f"To resume, use: run_preprocessing(adata, resume_from='{save_dir}')")
                return adata

        raise WorkflowError(error_msg, step_name=current_step or "unknown", original_error=e)

    # Store final config
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        UnsKeys.WORKFLOW_CONFIG
    ] = active_config.to_dict()
    adata.uns["sclucid"]["preprocess"][UnsKeys.STEPS_EXECUTED] = successful_steps

    # Tumor-aware preprocessing notes
    if is_tumor_context(tissue_type):
        tumor_notes: Dict[str, Any] = {"tissue_type": tissue_type, "tumor_aware_enabled": True}
        if active_config.integration.method and active_config.integration.batch_key:
            tumor_notes["batch_correction_note"] = (
                "Batch correction is enabled. For tumor data, review whether "
                "inter-tumor heterogeneity is being over-corrected."
            )
        adata.uns["sclucid"]["preprocess"]["tumor_aware_notes"] = tumor_notes
        log.info(f"Tumor-aware preprocessing notes stored: {list(tumor_notes.keys())}")

    # Build and store review summary
    enriched_summary = enrich_preprocessing_review_summary(
        _build_preprocessing_review_summary(adata, active_config, successful_steps, tissue_type),
        adata=adata,
        config=active_config,
        successful_steps=successful_steps,
        tissue_type=tissue_type,
        keep_intermediate_layers=keep_intermediate_layers,
    )
    review_summary = normalize_review_summary(
        enriched_summary,
        module="preprocess",
        workflow_name="standard",
        adata=adata,
        steps_executed=successful_steps,
        config=active_config.to_dict(),
        warnings=(
            enriched_summary.get("preprocess_readiness", {}).get("review_reasons", [])
            if isinstance(enriched_summary.get("preprocess_readiness"), dict)
            else []
        ),
    )
    validate_review_summary_schema(review_summary, module="preprocess", raise_on_error=True)
    validate_preprocessing_review_summary(review_summary, raise_on_error=True)
    adata.uns["sclucid"]["preprocess"][UnsKeys.REVIEW_SUMMARY] = review_summary

    # Export review summary to file if save_dir is configured
    if results_path:
        export_review_summary(
            review_summary,
            save_dir=results_path,
            module="preprocess",
            title="Preprocessing Review Summary",
            adata=adata,
        )

    log.info("=" * 60)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)
    return adata


def run_iterative_preprocessing(
    adata: AnnData,
    config: Optional[WorkflowConfig] = None,
    save_dir: Optional[str] = None,
    *,
    tissue_type: str = "unknown",
    sample_key: str = "sampleID",
    biology_columns: Optional[List[str]] = None,
    condition_key: Optional[str] = None,
    integration_policy: Literal["auto_review", "force", "off"] = "auto_review",
    run_hvg_stability: bool = True,
    hvg_stability_kwargs: Optional[Dict[str, Any]] = None,
    hvg_keys: Optional[List[str]] = None,
    hvg_mode: Literal["auto", "direct", "intersection", "union", "difference"] = "auto",
    run_diagnostic_embedding: bool = True,
    optimize_final_graph: bool = True,
    diagnostic_umap_key: str = "X_umap_pca_diagnostic",
    final_umap_key: Optional[str] = None,
    show_progress: bool = True,
    inplace: bool = False,
    keep_intermediate_layers: bool = True,
) -> AnnData:
    """Run reviewer-first preprocessing for real-project Step1 workflows.

    This entrypoint keeps ``run_preprocessing`` as the step-control/resume API
    while adding review-oriented phases that real projects currently run in
    notebooks: HVG strategy audit, HVG stability, diagnostic pre-integration
    embedding, integration decision, optional correction, and final graph/UMAP.
    """
    active_config = _prepare_workflow_config(config)
    results_path = Path(save_dir or active_config.save_dir) if (save_dir or active_config.save_dir) else None
    if results_path:
        results_path.mkdir(parents=True, exist_ok=True)

    available_steps = _resolve_steps(None, None, active_config)
    bootstrap_steps = [
        step
        for step in [
            "gene_filtering",
            "normalization",
            "set_raw",
            "regression",
            "hvg_selection",
            "subset_hvg",
            "scaling",
            "pca",
        ]
        if step in available_steps
    ]
    result = run_preprocessing(
        adata,
        config=active_config,
        save_dir=str(results_path / "bootstrap") if results_path else None,
        steps=bootstrap_steps,
        inplace=inplace,
        keep_intermediate_layers=keep_intermediate_layers,
        tissue_type=tissue_type,
        show_progress=show_progress,
    )

    pp_ns = result.uns.setdefault("sclucid", {}).setdefault("preprocess", {})
    summary: Dict[str, Any] = {
        "schema_version": "iterative_preprocessing_v1",
        "workflow": "run_iterative_preprocessing",
        "tissue_type": tissue_type,
        "sample_key": sample_key,
        "bootstrap_steps": bootstrap_steps,
        "hvg_strategy": {},
        "hvg_stability": {},
        "diagnostic_embedding": {},
        "integration_decision": {},
        "final_embedding": {},
        "review_action_items": [],
    }

    detected_hvg_keys = hvg_keys or [
        str(col)
        for col in result.var.columns
        if str(col).startswith("highly_variable") and result.var[col].dtype == bool
    ]
    if detected_hvg_keys:
        try:
            _, hvg_audit = select_and_audit_hvgs(
                result,
                hvg_keys=detected_hvg_keys,
                mode=hvg_mode,
                subset=False,
                keep_raw=False,
                evaluate_stability=False,
                plot_venn=False,
                save_dir=str(results_path / "hvg_audit") if results_path else None,
            )
            summary["hvg_strategy"] = hvg_audit
        except Exception as exc:
            summary["hvg_strategy"] = {
                "status": "skipped",
                "reason": f"{type(exc).__name__}: {exc}",
                "hvg_keys": detected_hvg_keys,
            }
            summary["review_action_items"].append(
                {
                    "priority": "review",
                    "action": "Inspect HVG strategy manually.",
                    "rationale": "HVG audit could not be generated.",
                }
            )

    if run_hvg_stability and detected_hvg_keys:
        stability_kwargs = dict(hvg_stability_kwargs or {})
        stability_kwargs.setdefault("plot", False)
        try:
            evaluate_hvg_stability(
                result,
                hvg_key=detected_hvg_keys[0],
                **stability_kwargs,
            )
            summary["hvg_stability"] = pp_ns.get("hvg_stability", {})
        except Exception as exc:
            summary["hvg_stability"] = {
                "status": "skipped",
                "reason": f"{type(exc).__name__}: {exc}",
                "hvg_key": detected_hvg_keys[0],
            }
            summary["review_action_items"].append(
                {
                    "priority": "review",
                    "action": "Review HVG stability outside the workflow.",
                    "rationale": "HVG stability evaluation could not be completed.",
                }
            )

    final_rep = "X_pca"
    n_pcs = int(result.obsm["X_pca"].shape[1]) if "X_pca" in result.obsm else None
    if run_diagnostic_embedding and n_pcs is not None:
        try:
            result, diagnostic_grid = run_embedding_pipeline(
                result,
                use_rep="X_pca",
                optimize=False,
                default_n_neighbors=active_config.graph.n_neighbors,
                default_n_pcs=min(active_config.graph.n_pcs, n_pcs),
                umap_key=diagnostic_umap_key,
            )
            diagnostic_meta = dict(pp_ns.get("embedding_workflow", {}))
            diagnostic_meta["n_grid_rows"] = int(len(diagnostic_grid))
            summary["diagnostic_embedding"] = diagnostic_meta
        except Exception as exc:
            summary["diagnostic_embedding"] = {
                "status": "skipped",
                "reason": f"{type(exc).__name__}: {exc}",
            }

    integration_run = False
    integration_warnings: List[str] = []
    integration_risk = None
    batch_key = active_config.integration.batch_key
    if isinstance(batch_key, list):
        batch_key = batch_key[0] if batch_key else None
    biology_cols = biology_columns if biology_columns is not None else active_config.integration.biology_columns
    integration_mode: bool | Literal["auto"]
    if integration_policy == "force":
        integration_mode = True
    elif integration_policy == "off":
        integration_mode = False
    else:
        integration_mode = "auto"

    if batch_key and "X_pca" in result.obsm:
        integration_run, integration_warnings, integration_risk = decide_integration(
            result,
            batch_key=str(batch_key),
            run_integration=integration_mode,
            biology_columns=biology_cols,
            condition_key=condition_key or active_config.integration.condition_key,
            tumor=is_tumor_context(tissue_type),
            before_rep="X_pca",
        )
    else:
        integration_warnings = ["missing batch_key or X_pca; integration decision skipped"]

    integration_output_key = (
        active_config.integration.output_key
        or f"X_{active_config.integration.method}_{batch_key}"
        if active_config.integration.method and batch_key
        else None
    )
    summary["integration_decision"] = {
        "policy": integration_policy,
        "run_integration": bool(integration_run),
        "batch_key": batch_key,
        "biology_columns": biology_cols or [],
        "warnings": integration_warnings,
        "risk": integration_risk,
        "output_key": integration_output_key if integration_run else None,
    }

    if integration_run and integration_output_key:
        integration_config = active_config.integration.model_copy(
            update={
                "batch_key": batch_key,
                "use_rep": "X_pca",
                "output_key": integration_output_key,
                "auto_decide": False,
                "condition_key": condition_key or active_config.integration.condition_key,
                "biology_columns": biology_cols or [],
                "tumor": is_tumor_context(tissue_type),
                "evaluate": False,
            }
        )
        result = batch_correction(
            result,
            config=integration_config,
            save_dir=str(results_path / "integration") if results_path else None,
            force=True,
        )
        if integration_output_key in result.obsm:
            final_rep = integration_output_key
        else:
            summary["review_action_items"].append(
                {
                    "priority": "review",
                    "action": "Inspect integration output.",
                    "rationale": f"Expected representation {integration_output_key!r} was not created.",
                }
            )

    if final_rep in result.obsm:
        try:
            rep_dims = int(result.obsm[final_rep].shape[1])
            active_final_umap_key = final_umap_key or f"X_umap_{final_rep.replace('X_', '').lower()}"
            result, final_grid = run_embedding_pipeline(
                result,
                use_rep=final_rep,
                optimize=optimize_final_graph,
                default_n_neighbors=active_config.graph.n_neighbors,
                default_n_pcs=min(active_config.graph.n_pcs, rep_dims),
                umap_key=active_final_umap_key,
            )
            final_meta = dict(pp_ns.get("embedding_workflow", {}))
            final_meta["n_grid_rows"] = int(len(final_grid))
            summary["final_embedding"] = final_meta
        except Exception as exc:
            summary["final_embedding"] = {
                "status": "skipped",
                "reason": f"{type(exc).__name__}: {exc}",
                "use_rep": final_rep,
            }

    summary["final_representation"] = final_rep
    summary["layers"] = list(result.layers.keys())
    summary["obsm_keys"] = list(result.obsm.keys())
    summary["raw_present"] = result.raw is not None
    pp_ns["iterative_preprocessing_summary"] = _json_safe(summary)

    review_summary = pp_ns.get(UnsKeys.REVIEW_SUMMARY)
    if isinstance(review_summary, dict):
        payload = review_summary.setdefault("data", {})
        if isinstance(payload, dict):
            payload["iterative_preprocessing_summary"] = _json_safe(summary)

    return result


def _resolve_steps(
    steps: Optional[Sequence[str]],
    skip_steps: Optional[Sequence[str]],
    config: Optional[WorkflowConfig] = None,
) -> List[str]:
    """Resolve which steps to run based on steps and skip_steps parameters."""
    if steps is not None and skip_steps is not None:
        raise ValueError("Cannot specify both 'steps' and 'skip_steps'. Choose one.")

    if steps is not None:
        return list(steps)

    if skip_steps is not None:
        return [s for s in WORKFLOW_STEPS if s not in skip_steps]

    if config is None:
        return WORKFLOW_STEPS

    disabled_steps = set()
    if not config.run_gene_filtering:
        disabled_steps.add("gene_filtering")
    if not config.run_regression:
        disabled_steps.add("regression")
    if not config.run_scaling:
        disabled_steps.add("scaling")
    if not config.run_pca:
        disabled_steps.add("pca")
    if not _integration_step_requested(config):
        disabled_steps.add("batch_correction")
    if not config.run_neighbors:
        disabled_steps.add("neighbors_umap")

    return [step for step in WORKFLOW_STEPS if step not in disabled_steps]


def _integration_step_requested(config: WorkflowConfig) -> bool:
    """Return whether the workflow should execute batch correction."""
    if config.run_integration:
        return bool(config.integration.method and config.integration.batch_key)

    if "run_integration" in getattr(config, "model_fields_set", set()):
        return False

    integration_fields = getattr(config.integration, "model_fields_set", set())
    user_touched_integration = bool({"method", "batch_key"} & set(integration_fields))
    return bool(user_touched_integration and config.integration.method and config.integration.batch_key)


def _prepare_workflow_config(config: Optional[WorkflowConfig]) -> WorkflowConfig:
    """Create an isolated workflow config with authoritative workflow layer names."""
    active_config = WorkflowConfig() if config is None else config.model_copy(deep=True)
    active_config.normalization = active_config.normalization.model_copy(
        update={
            "input_layer": active_config.counts_layer,
            "output_layer": active_config.normalized_layer,
        }
    )
    return active_config


def _run_gene_filtering_step(adata: AnnData, config: WorkflowConfig) -> AnnData:
    """Filter low-detection genes and optionally apply biotype-aware filtering."""
    min_cells = min(int(config.min_cells_per_gene), max(int(adata.n_obs), 1))
    source_name = config.counts_layer if config.counts_layer in adata.layers else "X"
    X = adata.layers[config.counts_layer] if source_name != "X" else adata.X

    initial_genes_before_step = int(adata.n_vars)
    if scipy.sparse.issparse(X):
        detected_cells = np.diff(X.tocsc().indptr)
    else:
        detected_cells = np.asarray(X > 0).sum(axis=0)

    detected_cells = np.asarray(detected_cells).ravel()
    keep_mask = detected_cells >= min_cells
    initial_genes = int(adata.n_vars)
    kept_genes = int(keep_mask.sum())
    removed_genes = initial_genes - kept_genes

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["gene_filtering"] = {
        "source": source_name,
        "min_cells_per_gene": min_cells,
        "initial_genes": initial_genes,
        "removed_genes": removed_genes,
        "final_genes": kept_genes if kept_genes > 0 else initial_genes,
        "skipped": kept_genes == 0,
    }

    if kept_genes == 0:
        log.warning(
            "Gene filtering would remove all genes with min_cells_per_gene=%s; skipping.",
            min_cells,
        )
    elif removed_genes > 0:
        log.info(
            "Gene filtering removed %s/%s genes expressed in fewer than %s cells.",
            removed_genes,
            initial_genes,
            min_cells,
        )
        adata._inplace_subset_var(keep_mask)
    else:
        log.info("Gene filtering kept all %s genes.", initial_genes)

    adata = _run_gene_biotype_filtering_step(
        adata,
        config,
        initial_genes_before_step=initial_genes_before_step,
        current_stage="before_normalization",
    )
    return adata


def _run_gene_biotype_filtering_step(
    adata: AnnData,
    config: WorkflowConfig,
    *,
    initial_genes_before_step: int,
    current_stage: Literal["before_normalization", "after_raw"],
) -> AnnData:
    """Optionally annotate/filter genes by biotype and record audit metadata."""
    biotype_config = config.gene_biotype
    meta: Dict[str, Any] = {
        "enabled": bool(biotype_config.annotate or biotype_config.filter),
        "annotate": bool(biotype_config.annotate),
        "filter": bool(biotype_config.filter),
        "species": biotype_config.species,
        "method": biotype_config.method,
        "filter_stage": biotype_config.filter_stage,
        "current_stage": current_stage,
        "initial_genes_before_gene_filtering": int(initial_genes_before_step),
        "initial_genes": int(adata.n_vars),
        "final_genes": int(adata.n_vars),
        "removed_genes": 0,
        "skipped": not bool(biotype_config.annotate or biotype_config.filter),
    }

    if not biotype_config.annotate and not biotype_config.filter:
        adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
            "gene_biotype_filtering"
        ] = meta
        return adata

    try:
        if biotype_config.annotate and (
            biotype_config.overwrite
            or "biotype" not in adata.var.columns
            or "biotype_category" not in adata.var.columns
        ):
            biotype_df = None
            if biotype_config.method == "custom" and biotype_config.custom_biotype_path:
                custom_path = Path(biotype_config.custom_biotype_path)
                sep = "\t" if custom_path.suffix.lower() in {".tsv", ".txt"} else ","
                biotype_df = pd.read_csv(custom_path, sep=sep)
                rename_map = {
                    "external_gene_name": "gene_name",
                    "gene_biotype": "biotype",
                    "ensembl_gene_id": "gene_id",
                    "chromosome_name": "chromosome",
                    "start_position": "start",
                    "end_position": "end",
                }
                biotype_df = biotype_df.rename(
                    columns={k: v for k, v in rename_map.items() if k in biotype_df.columns}
                )

            annotate_gene_biotypes(
                adata,
                species=biotype_config.species,
                method="reference" if biotype_config.method == "ensembl" else biotype_config.method,
                biotype_df=biotype_df,
                fuzzy_match=biotype_config.fuzzy_match,
                overwrite=biotype_config.overwrite,
                allow_download=biotype_config.allow_download,
                cache_dir=biotype_config.cache_dir,
                prefer_bundled=biotype_config.prefer_bundled,
            )

        if "biotype" in adata.var.columns and "biotype_category" in adata.var.columns:
            meta["annotation_rate"] = float(
                adata.var["biotype"].notna().sum() / max(float(adata.n_vars), 1.0)
            )
            meta["biotype_categories"] = {
                str(k): int(v)
                for k, v in adata.var["biotype_category"].value_counts(dropna=False).items()
            }

        should_filter_now = bool(
            biotype_config.filter and biotype_config.filter_stage == current_stage
        )
        if biotype_config.filter and not should_filter_now:
            meta["status"] = "annotated_pending_filter"
        elif should_filter_now:
            before = int(adata.n_vars)
            filtered = filter_genes_by_biotype(
                adata,
                keep_biotypes=biotype_config.keep_biotypes,
                use_recommended=biotype_config.use_recommended,
                copy=True,
            )
            if filtered is not None:
                adata = filtered
            meta["strategy"] = (
                f"custom: {', '.join(biotype_config.keep_biotypes)}"
                if biotype_config.keep_biotypes
                else "recommended biotypes"
            )
            meta["final_genes"] = int(adata.n_vars)
            meta["removed_genes"] = int(before - adata.n_vars)
            meta["status"] = "completed"

        meta["skipped"] = False
        meta.setdefault("status", "completed")
    except Exception as exc:
        meta.update(
            {
                "status": "failed",
                "skipped": True,
                "error_type": type(exc).__name__,
                "message": str(exc),
            }
        )
        log.warning("Gene biotype annotation/filtering failed: %s", exc)
        if biotype_config.fail_on_error:
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
                "gene_biotype_filtering"
            ] = meta
            raise

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "gene_biotype_filtering"
    ] = meta
    return adata


def _run_step(
    adata: AnnData,
    step_name: str,
    custom_pre: Optional[Callable],
    custom_post: Optional[Callable],
    step_func: Callable[[AnnData], AnnData],
) -> AnnData:
    """Execute a single workflow step with optional custom hooks."""
    log.info(f"Step: {step_name}")

    # Pre-step hook
    if custom_pre:
        adata = custom_pre(adata, step_name)

    # Execute main step
    adata = step_func(adata)

    # Post-step hook
    if custom_post:
        adata = custom_post(adata, step_name)

    return adata


def _get_hvg_input_layer(adata: AnnData, config: WorkflowConfig) -> str:
    """Choose the best available layer for HVG calculation."""
    if getattr(config.hvg, "method", None) == "deviance" and config.counts_layer in adata.layers:
        return config.counts_layer
    if config.normalized_layer in adata.layers:
        return config.normalized_layer
    if config.regressed_layer in adata.layers:
        log.warning(
            f"Normalized layer '{config.normalized_layer}' not found; falling back to '{config.regressed_layer}' for HVG selection."
        )
        return config.regressed_layer
    return config.normalized_layer


def _get_preferred_expression_layer(
    adata: AnnData,
    preferred_layers: Sequence[str],
) -> Optional[str]:
    """Return the first available layer from a preference-ordered list."""
    for layer in preferred_layers:
        if layer in adata.layers:
            return layer
    return None


def _run_scaling_step(adata: AnnData, config: WorkflowConfig) -> AnnData:
    """Scale the appropriate workflow layer without re-running inline regression."""
    scaling_input_layer = _get_preferred_expression_layer(
        adata,
        [config.regressed_layer, config.normalized_layer],
    )
    if scaling_input_layer is not None:
        adata.X = adata.layers[scaling_input_layer].copy()

    scaling_config = config.scaling.model_copy(update={"regress_in_scale": False})
    return scale_data(adata, config=scaling_config, output_layer=config.scaled_layer)


def _set_raw_layer(adata: AnnData, config: WorkflowConfig) -> AnnData:
    """Set .raw with normalized data before any regression."""
    log.info(f"Storing data from layer '{config.normalized_layer}' into .raw")
    if config.normalized_layer not in adata.layers:
        raise KeyError(
            f"Layer '{config.normalized_layer}' not found. "
            "Normalization step may have failed or was skipped."
        )
    adata.raw = AnnData(
        X=adata.layers[config.normalized_layer].copy(),
        var=adata.var.copy(),
        obs=adata.obs.copy(),
    )
    return adata


def _subset_to_hvgs(
    adata: AnnData,
    config: WorkflowConfig,
    keep_intermediate_layers: bool,
) -> AnnData:
    """Subset data to HVGs."""
    log.info("Subsetting data to final HVG set")

    # Get HVG key from previous step
    hvg_info = adata.uns.get("sclucid", {}).get("preprocess", {}).get("hvg", {})
    hvg_key = hvg_info.get("output_key", "highly_variable")

    adata = select_hvg_sets(
        adata,
        hvg_keys=[hvg_key],
        mode="direct",
        subset=True,
        keep_raw=False,  # .raw is already correctly set
    )

    # Optionally clean up intermediate layers
    if not keep_intermediate_layers:
        layers_to_clean = [
            config.normalized_layer,
            config.regressed_layer if config.scaling.vars_to_regress else None,
        ]
        for layer in layers_to_clean:
            if layer and layer in adata.layers and layer != config.scaled_layer:
                del adata.layers[layer]
                log.info(f"Removed intermediate layer after subsetting: {layer}")

    return adata


def _select_n_pcs(
    variance_ratio: np.ndarray,
    method: str = "elbow",
    cumulative_threshold: float = 0.9,
) -> int:
    """Automatically select number of PCs.

    Parameters
    ----------
    variance_ratio : np.ndarray
        PCA variance ratio array.
    method : str
        "elbow" (kneedle) or "cumulative".
    cumulative_threshold : float
        For cumulative method, fraction of variance to explain.

    Returns:
    -------
    int
        Recommended number of PCs.
    """
    variance_ratio = np.asarray(variance_ratio, dtype=float)
    n_pcs = len(variance_ratio)
    if n_pcs == 0:
        raise ValueError("variance_ratio must contain at least one component.")
    if n_pcs < 3:
        return n_pcs

    def _bounded(value: int, floor: int = 10) -> int:
        return min(max(floor, value), n_pcs)

    if method == "cumulative":
        cumsum = np.cumsum(variance_ratio)
        hit = np.flatnonzero(cumsum >= cumulative_threshold)
        n_selected = int(hit[0]) + 1 if len(hit) else n_pcs
        return _bounded(n_selected)

    elif method == "elbow":
        # Guard: if variance ratio is effectively constant, kneedle is unreliable.
        # Fall back to cumulative method in that case.
        if np.std(variance_ratio) < 1e-6 or not np.all(np.isfinite(variance_ratio)):
            log.warning(
                "Variance ratio is constant or contains NaN/Inf. "
                "Falling back to cumulative method for n_pcs selection."
            )
            cumsum = np.cumsum(variance_ratio)
            hit = np.flatnonzero(cumsum >= cumulative_threshold)
            n_selected = int(hit[0]) + 1 if len(hit) else n_pcs
            return _bounded(n_selected)

        # Kneedle algorithm on log-transformed variance ratio
        x = np.arange(1, n_pcs + 1)
        y = np.log(variance_ratio + 1e-10)

        # Normalize to unit square
        x_norm = (x - x.min()) / (x.max() - x.min())
        y_norm = (y - y.min()) / (y.max() - y.min())

        # Compute distances from the line connecting first and last points
        line = y_norm[0] + (y_norm[-1] - y_norm[0]) * x_norm
        distances = line - y_norm

        # Find the knee (max distance)
        # Exclude the first and last 5% to avoid edge effects
        margin = max(1, int(0.05 * n_pcs))
        if margin * 2 >= n_pcs:
            return max(1, min(n_pcs, int(np.argmax(distances)) + 1))
        # If all distances are very small, there's no clear elbow
        if distances[margin:-margin].max() < 0.01:
            log.warning(
                "No clear elbow detected in variance ratio curve. "
                "Falling back to cumulative method for n_pcs selection."
            )
            cumsum = np.cumsum(variance_ratio)
            hit = np.flatnonzero(cumsum >= cumulative_threshold)
            n_selected = int(hit[0]) + 1 if len(hit) else n_pcs
            return _bounded(n_selected)
        knee_idx = margin + int(np.argmax(distances[margin:-margin]))
        return _bounded(knee_idx + 1, floor=15)

    else:
        raise ValueError(f"Unknown n_pcs selection method: {method}")


def _run_pca(
    adata: AnnData,
    config: WorkflowConfig,
    results_path: Optional[Path],
) -> AnnData:
    """Run PCA and optionally save variance plot."""
    pca_input_layer = _get_preferred_expression_layer(
        adata,
        [config.scaled_layer, config.regressed_layer, config.normalized_layer],
    )
    if pca_input_layer is not None:
        log.info(
            f"Setting adata.X to layer '{pca_input_layer}' for PCA. "
            "The original expression matrix has been overwritten."
        )
        adata.X = adata.layers[pca_input_layer].copy()

    max_valid_pcs = max(1, min(adata.n_obs, adata.n_vars) - 1)

    # Auto-select n_pcs if configured
    n_comps = config.graph.n_pcs
    auto_select = getattr(config, "auto_select_n_pcs", False)
    if auto_select and hasattr(config, "n_pcs_selection_method"):
        # Run PCA with max possible components first for auto-selection
        temp_n_comps = min(100, max_valid_pcs)
        log.info(f"Running PCA with {temp_n_comps} components for auto-selection...")
        sc.tl.pca(adata, n_comps=temp_n_comps)
        vr = adata.uns["pca"]["variance_ratio"]
        n_comps = _select_n_pcs(vr, method=config.n_pcs_selection_method)
        log.info(f"Auto-selected {n_comps} PCs ({config.n_pcs_selection_method} method)")
        adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
            "pca_n_pcs_selection"
        ] = {
            "schema_version": "pca_n_pcs_selection_v1",
            "method": config.n_pcs_selection_method,
            "selected_n_pcs": int(n_comps),
            "model_type": "log_variance_ratio_kneedle_heuristic"
            if config.n_pcs_selection_method == "elbow"
            else "cumulative_variance_threshold",
            "claim_level": "heuristic_pca_dimension_recommendation",
            "review_note": (
                "The elbow method uses a log-transformed variance-ratio curve and should "
                "be treated as a heuristic recommendation, not a formal optimum."
                if config.n_pcs_selection_method == "elbow"
                else "Cumulative variance threshold is a rule-of-thumb dimension choice."
            ),
        }
        # Re-run with selected n_comps
        sc.tl.pca(adata, n_comps=n_comps)
    else:
        n_comps = min(n_comps, max_valid_pcs)
        if n_comps != config.graph.n_pcs:
            log.info(
                f"PCA requested {config.graph.n_pcs} components but data supports {n_comps}; clipping to valid range."
            )
        log.info(f"PCA (using {n_comps} components)")
        sc.tl.pca(adata, n_comps=n_comps)

    if results_path:
        try:
            old_figdir = sc.settings.figdir
            sc.settings.figdir = results_path
            try:
                sc.pl.pca_variance_ratio(
                    adata,
                    log=True,
                    save="_variance_ratio.png",
                    show=False,
                )
            finally:
                sc.settings.figdir = old_figdir
        except Exception as e:
            log.warning(f"Could not save PCA variance plot: {e}")

    return adata


def _run_neighbors_umap(
    adata: AnnData,
    config: WorkflowConfig,
    results_path: Optional[Path],
) -> AnnData:
    """Run neighbors and UMAP."""
    batch_key = config.integration.batch_key
    if isinstance(batch_key, list):
        batch_key = batch_key[0] if batch_key else None

    # Use the actual output key recorded by batch_correction when it was run.
    integration_meta = adata.uns.get("sclucid", {}).get("preprocess", {}).get("integration", {})
    workflow_output_key = integration_meta.get("workflow", {}).get("output_key")
    integration_output_key = (
        config.integration.output_key
        or workflow_output_key
        or (f"X_{config.integration.method}_{batch_key}" if config.integration.method and batch_key else None)
    )

    # Prefer the corrected embedding when batch correction was run and produced
    # a valid representation; otherwise fall back to PCA.
    if (
        integration_output_key
        and integration_output_key in adata.obsm
        and workflow_output_key is not None
    ):
        use_rep = integration_output_key
        effective_n_pcs = min(config.graph.n_pcs, adata.obsm[use_rep].shape[1])
    elif "X_pca" in adata.obsm:
        use_rep = "X_pca"
        effective_n_pcs = min(config.graph.n_pcs, adata.obsm["X_pca"].shape[1])
    else:
        raise KeyError(
            "No usable embedding found in adata.obsm. Run PCA or batch correction first."
        )

    effective_n_neighbors = min(config.graph.n_neighbors, max(2, adata.n_obs - 1))

    log.info(f"Neighbors graph and UMAP using {use_rep!r}")
    sc.pp.neighbors(
        adata,
        n_pcs=effective_n_pcs,
        n_neighbors=effective_n_neighbors,
        use_rep=use_rep,
    )
    sc.tl.umap(adata)

    # Record which representation was actually used for the final graph so that
    # downstream analysis and review summaries stay consistent.
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["final_neighbors_rep"] = use_rep

    if results_path:
        try:
            color_vars = [
                v for v in [config.integration.batch_key, "phase"] if v and v in adata.obs.columns
            ]
            if color_vars:
                sc.pl.umap(adata, color=color_vars, save="_final.png", show=False, dpi=300)
                Path("./figures/umap_final.png").rename(results_path / "final_umap.png")
        except Exception:
            log.warning("Could not save final UMAP plot.")

    return adata


def _build_preprocessing_review_summary(
    adata: AnnData,
    config: WorkflowConfig,
    successful_steps: List[str],
    tissue_type: str,
) -> Dict[str, Any]:
    """Build a human-reviewable summary of the preprocessing run."""
    summary: Dict[str, Any] = {
        "steps_executed": successful_steps,
        "data_shape": {
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
        },
    }

    # Layers present
    summary["layers_present"] = [
        layer
        for layer in [
            config.counts_layer,
            config.normalized_layer,
            config.regressed_layer,
            config.scaled_layer,
        ]
        if layer in adata.layers
    ]

    if "gene_filtering" in successful_steps:
        summary["gene_filtering"] = (
            adata.uns.get("sclucid", {}).get("preprocess", {}).get("gene_filtering", {})
        )
    if "gene_biotype_filtering" in adata.uns.get("sclucid", {}).get("preprocess", {}):
        summary["gene_biotype_filtering"] = (
            adata.uns.get("sclucid", {})
            .get("preprocess", {})
            .get("gene_biotype_filtering", {})
        )

    # Normalization
    if "normalization" in successful_steps:
        summary["normalization"] = {
            "method": config.normalization.method,
            "target_sum": config.normalization.target_sum,
        }

    # HVG
    if "hvg_selection" in successful_steps:
        hvg_config = config.hvg
        n_hvgs = None
        if "highly_variable" in adata.var.columns:
            n_hvgs = int(adata.var["highly_variable"].sum())
        summary["hvg"] = {
            "n_top_genes": hvg_config.n_top_genes,
            "n_hvgs_selected": n_hvgs,
            "method": hvg_config.method,
        }

    # PCA
    if "pca" in successful_steps:
        pca_info: Dict[str, Any] = {}
        if "X_pca" in adata.obsm:
            pca_info["n_comps"] = int(adata.obsm["X_pca"].shape[1])
        if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
            vr = adata.uns["pca"]["variance_ratio"]
            pca_info["variance_explained_top3"] = [round(float(v), 4) for v in vr[:3]]
        summary["pca"] = pca_info

    # Batch correction
    if "batch_correction" in successful_steps:
        if config.integration.method and config.integration.batch_key:
            summary["batch_correction"] = {
                "method": config.integration.method,
                "batch_key": config.integration.batch_key,
                "status": "applied",
            }
        else:
            summary["batch_correction"] = {"status": "skipped (no method or batch_key)"}

    # Neighbors & UMAP
    if "neighbors_umap" in successful_steps:
        summary["neighbors"] = {
            "n_neighbors": config.graph.n_neighbors,
            "n_pcs": config.graph.n_pcs,
        }
        if "X_umap" in adata.obsm:
            summary["neighbors"]["umap_computed"] = True

    # Tumor aware
    if is_tumor_context(tissue_type):
        summary["tumor_aware"] = {
            "tissue_type": tissue_type,
            "enabled": True,
        }

    return summary
