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

import scanpy as sc
from anndata import AnnData

from ..utils import (
    PartialResultManager,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    export_review_summary,
    get_progress_bar,
    normalize_review_summary,
    validate_review_summary_schema,
)
from .config import WorkflowConfig
from .hvg import find_hvgs, select_hvg_sets
from .integrate import batch_correction
from .normalize import normalize_data
from .scale import regress_out, scale_data
from .trace import enrich_preprocessing_review_summary, validate_preprocessing_review_summary

log = logging.getLogger(__name__)

__all__ = [
    "run_preprocessing",
    "WORKFLOW_STEPS",
    "WorkflowError",
    "PartialWorkflowResult",
]

# Define workflow steps for flexible execution
WORKFLOW_STEPS = [
    "normalization",
    "set_raw",
    "regression",
    "hvg_selection",
    "subset_hvg",
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

            # --- 1. Normalization ---
            if step_name == "normalization":
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
                    lambda a: _set_raw_layer(a, active_config),
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
    if tissue_type and ("tumor" in tissue_type.lower() or "cancer" in tissue_type.lower()):
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
        )

    log.info("=" * 60)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)
    return adata


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
    if not config.run_regression:
        disabled_steps.add("regression")
    if not config.run_scaling:
        disabled_steps.add("scaling")
    if not config.run_pca:
        disabled_steps.add("pca")
    if not config.run_integration:
        disabled_steps.add("batch_correction")
    if not config.run_neighbors:
        disabled_steps.add("neighbors_umap")

    return [step for step in WORKFLOW_STEPS if step not in disabled_steps]


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
        adata.X = adata.layers[pca_input_layer].copy()

    max_valid_pcs = max(1, min(adata.n_obs, adata.n_vars) - 1)
    n_comps = min(config.graph.n_pcs, max_valid_pcs)
    if n_comps != config.graph.n_pcs:
        log.info(
            f"PCA requested {config.graph.n_pcs} components but data supports {n_comps}; clipping to valid range."
        )
    log.info(f"PCA (using {n_comps} components)")
    sc.tl.pca(adata, n_comps=n_comps)

    if results_path:
        try:
            sc.pl.pca_variance_ratio(adata, log=True, save="_variance_ratio.png", show=False)
            fig_path = Path("./figures/pca_variance_ratio.png")
            if fig_path.exists():
                fig_path.rename(results_path / "pca_variance_ratio.png")
        except Exception as e:
            log.warning(f"Could not save PCA variance plot: {e}")

    return adata


def _run_neighbors_umap(
    adata: AnnData,
    config: WorkflowConfig,
    results_path: Optional[Path],
) -> AnnData:
    """Run neighbors and UMAP."""
    if "X_pca" not in adata.obsm:
        raise KeyError(
            "X_pca not found in adata.obsm. Run PCA first or provide a precomputed PCA embedding."
        )

    effective_n_pcs = min(config.graph.n_pcs, adata.obsm["X_pca"].shape[1])
    effective_n_neighbors = min(config.graph.n_neighbors, max(2, adata.n_obs - 1))

    log.info("Neighbors graph and UMAP")
    sc.pp.neighbors(
        adata,
        n_pcs=effective_n_pcs,
        n_neighbors=effective_n_neighbors,
        use_rep="X_pca",
    )
    sc.tl.umap(adata)

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
    if tissue_type and ("tumor" in tissue_type.lower() or "cancer" in tissue_type.lower()):
        summary["tumor_aware"] = {
            "tissue_type": tissue_type,
            "enabled": True,
        }

    return summary
