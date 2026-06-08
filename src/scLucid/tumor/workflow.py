"""High-level tumor analysis workflow for scLucid.

Provides a dual-layer API:
- run_tumor_analysis: compact high-level entry point for wet-lab users
- run_tumor_analysis_expert: explicit per-stage config override for bioinformatics users
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from anndata import AnnData

from ..utils import (
    ReviewAction,
    StepResult,
    sanitize_for_hdf5,
    step_results_to_storage,
    summarize_step_results,
)
from ..utils.storage import save_result, save_workflow_result
from .config import TumorAnalysisConfig, TumorWorkflowConfig

log = logging.getLogger(__name__)

__all__ = [
    "run_tumor_analysis",
    "run_tumor_analysis_expert",
]


def run_tumor_analysis(
    adata: AnnData,
    config: Optional[TumorWorkflowConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Run the complete tumor analysis workflow.

    Pipeline: QC -> Preprocessing -> Standard Analysis -> Tumor Stage

    Parameters
    ----------
    adata : AnnData
        Input single-cell data (raw counts preferred).
    config : TumorWorkflowConfig, optional
        Workflow configuration. If None, uses ``TumorWorkflowConfig.quick()``.
    **kwargs
        Overrides for top-level config fields (e.g., ``save_dir=..., cancer_type=...``).

    Returns:
    -------
    AnnData
        Annotated data with QC, preprocessing, analysis, and tumor results.
    """
    if config is None:
        config = TumorWorkflowConfig.quick(**kwargs)
    else:
        for key, value in kwargs.items():
            if key in config.model_fields:
                setattr(config, key, value)

    return run_tumor_analysis_expert(
        adata,
        qc_config=config.qc_config,
        preprocess_config=config.preprocess_config,
        analysis_config=config.analysis_config,
        tumor_config=config.tumor_config,
        use_recommendations=config.use_recommendations,
        tissue_type=config.tissue_type,
        batch_key=config.batch_key,
        cancer_type=config.cancer_type,
        save_dir=config.save_dir,
        n_jobs=config.n_jobs,
        random_state=config.random_state,
    )


def run_tumor_analysis_expert(
    adata: AnnData,
    qc_config: Optional["QCWorkflowConfig"] = None,
    preprocess_config: Optional["PreprocessWorkflowConfig"] = None,
    analysis_config: Optional["AnalysisWorkflowConfig"] = None,
    tumor_config: Optional[TumorAnalysisConfig] = None,
    *,
    use_recommendations: bool = True,
    tissue_type: str = "tumor",
    batch_key: Optional[str] = None,
    cancer_type: Optional[str] = None,
    save_dir: Optional[str] = None,
    n_jobs: int = -1,
    random_state: int = 42,
) -> AnnData:
    """
    Expert-layer tumor workflow with explicit per-stage configs.

    Parameters
    ----------
    adata : AnnData
        Input data.
    qc_config : Any, optional
        Explicit QC workflow config. If None and recommendations enabled, uses recommended defaults.
    preprocess_config : Any, optional
        Explicit preprocessing config.
    analysis_config : Any, optional
        Explicit analysis workflow config.
    tumor_config : TumorAnalysisConfig, optional
        Tumor-specific stage config.
    use_recommendations : bool
        Whether to run the recommendation engine before executing stages.
    tissue_type : str
        Tissue type hint passed to recommenders and QC/preprocessing.
    batch_key : str, optional
        Batch key for recommendations and integration.
    cancer_type : str, optional
        Cancer type for marker loading and tumor-aware logic.
    save_dir : str, optional
        Root directory for saved outputs.
    n_jobs : int
        Number of parallel jobs.
    random_state : int
        Random seed.

    Returns:
    -------
    AnnData
        Annotated data with full workflow results and execution trace.
    """
    adata = adata.copy()

    if tumor_config is None:
        tumor_config = TumorAnalysisConfig()

    # Resolve default configs if not provided
    if qc_config is None:
        from ..qc.config import QCWorkflowConfig

        qc_config = QCWorkflowConfig()
    if preprocess_config is None:
        from ..preprocess.config import WorkflowConfig as PreprocessWorkflowConfig

        preprocess_config = PreprocessWorkflowConfig()
    if analysis_config is None:
        from ..analysis.config import AnalysisWorkflowConfig

        analysis_config = AnalysisWorkflowConfig()

    # Ensure save_dir propagates
    if save_dir:
        for cfg in (qc_config, preprocess_config, analysis_config):
            if hasattr(cfg, "save_dir") and cfg.save_dir is None:
                (
                    object.__setattr__(cfg, "save_dir", save_dir)
                    if hasattr(cfg, "model_config")
                    else setattr(cfg, "save_dir", save_dir)
                )

    recommendations = None
    warnings_list: List[str] = []
    steps_executed: List[str] = []
    step_results: List[StepResult] = []

    # --- Recommendations ---
    if use_recommendations:
        log.info("Running recommendation engine...")
        try:
            from ..recommendation.config import RecommendationConfig
            from ..recommendation.engine import RecommendationEngine

            rec_modules = ["qc", "preprocess", "clustering", "annotation", "tumor"]
            rec_config = RecommendationConfig(
                modules=rec_modules,
            )
            engine = RecommendationEngine(config=rec_config)
            recommendations = engine.recommend(
                adata,
                batch_key=batch_key,
                tissue_type=tissue_type,
                cancer_type=cancer_type,
                plot=False,
                save_dir=Path(save_dir) if save_dir else None,
            )
            log.info(
                f"Recommendations generated with overall confidence: {recommendations.overall_confidence:.2f}"
            )

            # Apply recommended configs where not explicitly overridden by user
            if recommendations.get_section("qc") is not None:
                qc_config = _apply_qc_recommendations(qc_config, recommendations.get_section("qc"))
            if recommendations.get_section("preprocess") is not None:
                preprocess_config = _apply_preprocess_recommendations(
                    preprocess_config, recommendations.get_section("preprocess")
                )
            if recommendations.get_section("clustering") is not None:
                analysis_config = _apply_clustering_recommendations(
                    analysis_config, recommendations.get_section("clustering")
                )
            if recommendations.get_section("annotation") is not None:
                analysis_config = _apply_annotation_recommendations(
                    analysis_config, recommendations.get_section("annotation")
                )
            if recommendations.get_section("tumor") is not None:
                tumor_config = _apply_tumor_recommendations(
                    tumor_config, recommendations.get_section("tumor")
                )
        except Exception as exc:
            log.warning(f"Recommendation engine failed: {exc}. Proceeding with default configs.")
            warnings_list.append(f"recommendation_engine_failed: {exc}")

    # --- Stage 1: QC ---
    log.info("=" * 60)
    log.info("=== Starting Tumor Workflow: QC ===")
    log.info("=" * 60)
    try:
        from ..qc.workflow import run_standard_qc

        adata = run_standard_qc(adata, config=qc_config, tissue_type=tissue_type)
        steps_executed.append("qc")
    except Exception as exc:
        log.error(f"QC workflow failed: {exc}")
        warnings_list.append(f"qc_failed: {exc}")
        raise

    # --- Stage 2: Preprocessing ---
    log.info("=" * 60)
    log.info("=== Starting Tumor Workflow: Preprocessing ===")
    log.info("=" * 60)
    try:
        from ..preprocess.workflow import run_preprocessing

        adata = run_preprocessing(adata, config=preprocess_config, tissue_type=tissue_type)
        steps_executed.append("preprocessing")
    except Exception as exc:
        log.error(f"Preprocessing workflow failed: {exc}")
        warnings_list.append(f"preprocessing_failed: {exc}")
        raise

    # --- Stage 3: Standard Analysis ---
    log.info("=" * 60)
    log.info("=== Starting Tumor Workflow: Standard Analysis ===")
    log.info("=" * 60)
    try:
        from ..analysis.workflow import run_standard_analysis

        adata = run_standard_analysis(adata, config=analysis_config)
        steps_executed.append("analysis")
    except Exception as exc:
        log.error(f"Analysis workflow failed: {exc}")
        warnings_list.append(f"analysis_failed: {exc}")
        raise

    # --- Stage 4: Tumor-specific analysis ---
    log.info("=" * 60)
    log.info("=== Starting Tumor Workflow: Tumor Stage ===")
    log.info("=" * 60)
    try:
        adata, tumor_steps, tumor_step_results, tumor_warnings = _run_tumor_stage(
            adata, tumor_config, cancer_type=cancer_type
        )
        steps_executed.extend(tumor_steps)
        step_results.extend(tumor_step_results)
        warnings_list.extend(tumor_warnings)
    except Exception as exc:
        log.error(f"Tumor stage failed: {exc}")
        warnings_list.append(f"tumor_stage_failed: {exc}")
        raise

    # --- Execution Trace ---
    user_overrides = _diff_recommendations(
        recommendations,
        {
            "qc": qc_config,
            "preprocess": preprocess_config,
            "analysis": analysis_config,
            "tumor": tumor_config,
        },
    )

    execution_trace = sanitize_for_hdf5({
        "recommended_params": recommendations.to_dict() if recommendations else None,
        "actual_params": {
            "qc": qc_config.to_dict() if hasattr(qc_config, "to_dict") else {},
            "preprocess": (
                preprocess_config.to_dict() if hasattr(preprocess_config, "to_dict") else {}
            ),
            "analysis": analysis_config.to_dict() if hasattr(analysis_config, "to_dict") else {},
            "tumor": tumor_config.to_dict(),
        },
        "user_overrides": user_overrides,
        "warnings": warnings_list,
        "steps_executed": steps_executed,
        "step_results": step_results_to_storage(step_results),
        "tissue_type": tissue_type,
        "batch_key": batch_key,
        "cancer_type": cancer_type,
    })

    review_summary = _build_tumor_review_summary(
        adata=adata,
        tumor_config=tumor_config,
        step_results=step_results,
        warnings=warnings_list,
        cancer_type=cancer_type,
    )

    save_result(adata, "tumor", "execution_trace", execution_trace)
    save_result(adata, "tumor", "review_summary", sanitize_for_hdf5(review_summary))
    save_result(
        adata,
        "tumor",
        "step_results",
        sanitize_for_hdf5(step_results_to_storage(step_results)),
    )
    save_workflow_result(
        adata,
        module="tumor",
        workflow_name="tumor_analysis",
        steps=steps_executed,
        config={
            "use_recommendations": use_recommendations,
            "tissue_type": tissue_type,
            "batch_key": batch_key,
            "cancer_type": cancer_type,
        },
    )

    log.info("=" * 60)
    log.info("=== Tumor Workflow Complete ===")
    log.info(f"Completed steps: {steps_executed}")
    log.info(f"Tumor readiness: {review_summary.get('readiness', {}).get('status', 'unknown')}")
    log.info("=" * 60)

    return adata


def _run_tumor_stage(
    adata: AnnData,
    config: TumorAnalysisConfig,
    cancer_type: Optional[str] = None,
) -> Tuple[AnnData, List[str], List[StepResult], List[str]]:
    """
    Run tumor-specific analysis steps with structured step results.

    Each step is wrapped in try/except so failures degrade gracefully. Step
    results record ``status``, ``evidence_level``, and structured warnings.

    Returns:
    -------
    tuple of (AnnData, executed_step_names, step_results, warnings)
    """
    executed_steps: List[str] = []
    step_results: List[StepResult] = []
    stage_warnings: List[str] = []

    # Malignancy
    if config.run_malignancy:
        step_name = "malignancy"
        try:
            from .malignancy.classification import classify_malignant_cells
            from .malignancy.scoring import score_malignancy

            log.info("Tumor stage: scoring malignancy")
            adata = score_malignancy(adata, key_added="malignancy", cancer_type=cancer_type)
            executed_steps.append("malignancy_scoring")

            log.info("Tumor stage: classifying malignant cells")
            if config.malignancy_method in ("cnv", "combined"):
                classify_malignant_cells(
                    adata, method=config.malignancy_method, key_added="is_malignant"
                )
            elif config.malignancy_method in ("threshold", "ml"):
                # reference key may be used to subset reference cells
                ref_key = config.malignancy_reference_key
                ref_adata = None
                if ref_key and ref_key in adata.obs.columns:
                    ref_adata = adata[
                        adata.obs[ref_key]
                        .astype(str)
                        .str.lower()
                        .isin({"normal", "healthy", "reference", "immune", "stromal"})
                    ].copy()
                    if ref_adata.n_obs == 0:
                        log.warning(
                            f"No reference cells found via '{ref_key}'. Falling back to unsupervised."
                        )
                        ref_adata = None
                else:
                    log.warning(
                        f"malignancy_method='{config.malignancy_method}' may require reference cells."
                    )
                classify_malignant_cells(
                    adata,
                    method=config.malignancy_method,
                    reference_adata=ref_adata,
                    key_added="is_malignant",
                )
            executed_steps.append("malignancy_classification")
            step_results.append(
                StepResult(
                    name=step_name,
                    status="completed",
                    evidence_level="heuristic",
                    outputs={"method": config.malignancy_method},
                    warnings=[],
                )
            )
        except Exception as exc:
            msg = f"Malignancy analysis failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)
            step_results.append(
                StepResult.from_exception(
                    name=step_name,
                    exc=exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )

    # TME
    if config.run_tme:
        step_name = "tme_deconvolution"
        try:
            from .microenvironment.deconvolution import deconvolve_tme

            log.info("Tumor stage: deconvolving TME")
            cell_type_key = config.tme_cell_type_key
            fallback_used = False
            if cell_type_key not in adata.obs.columns:
                log.warning(f"TME cell type key '{cell_type_key}' not found. Trying 'cell_type'.")
                cell_type_key = "cell_type"
                fallback_used = True
            adata = deconvolve_tme(adata, cell_type_key=cell_type_key, key_added="tme")
            executed_steps.append(step_name)

            proportions = adata.uns.get("tme_proportions", {})
            if hasattr(proportions, "to_dict"):
                proportions = proportions.to_dict()
            step_results.append(
                StepResult(
                    name=step_name,
                    status="completed",
                    evidence_level="heuristic",
                    outputs={
                        "cell_type_key": cell_type_key,
                        "proportions": proportions,
                        "immune_score": float(adata.uns.get("tme_immune_score", 0.0)),
                        "stromal_score": float(adata.uns.get("tme_stromal_score", 0.0)),
                    },
                    warnings=[
                        "TME composition is annotation-derived, not bulk deconvolution."
                    ]
                    + ([f"used fallback cell_type key '{cell_type_key}'"] if fallback_used else []),
                )
            )
        except Exception as exc:
            msg = f"TME deconvolution failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)
            step_results.append(
                StepResult.from_exception(
                    name=step_name,
                    exc=exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )

    # CNV
    if config.run_cnv:
        step_name = "cnv_inference"
        try:
            from .cnv.infercnv import infer_cnv

            log.info("Tumor stage: inferring CNV")
            ref_key = config.cnv_reference_key
            ref_cells = None
            has_ref = False
            if ref_key and ref_key in adata.obs.columns:
                ref_mask = (
                    adata.obs[ref_key]
                    .astype(str)
                    .str.lower()
                    .isin({"normal", "healthy", "reference", "immune", "stromal"})
                )
                has_ref = bool(ref_mask.any())
                if has_ref:
                    ref_cells = adata.obs.loc[ref_mask, ref_key].unique().tolist()

            adata = infer_cnv(
                adata,
                reference_cells=ref_cells,
                reference_key=ref_key if ref_key and ref_key in adata.obs.columns else "cell_type",
                key_added="cnv",
            )
            executed_steps.append(step_name)

            cnv_summary = adata.uns.get("cnv_summary", {})
            has_coords = all(c in adata.var.columns for c in ("chromosome", "start", "end"))
            evidence_level: Any = "validated_core" if (has_ref and has_coords) else "heuristic"

            step_results.append(
                StepResult(
                    name=step_name,
                    status="completed",
                    evidence_level=evidence_level,
                    outputs={
                        "score_key": "cnv_score",
                        "predicted_class_key": "cnv_predicted_class",
                        "n_aneuploid": cnv_summary.get("n_aneuploid"),
                        "n_diploid": cnv_summary.get("n_diploid"),
                        "threshold": cnv_summary.get("threshold"),
                        "mean_cnv_score": cnv_summary.get("mean_cnv_score"),
                        "has_reference_cells": has_ref,
                        "has_genomic_coordinates": has_coords,
                    },
                    warnings=[] if has_ref else ["No reference cells provided; threshold-based calls less reliable."],
                )
            )
        except Exception as exc:
            msg = f"CNV inference failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)
            step_results.append(
                StepResult.from_exception(
                    name=step_name,
                    exc=exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )

    # Malignancy interpretation (run after CNV so CNV scores are available)
    if config.run_malignancy:
        step_name = "malignancy_interpretation"
        try:
            from .malignancy.interpretation import run_malignancy_interpretation

            log.info("Tumor stage: malignancy interpretation")
            annotation_key = config.malignancy_annotation_key
            if annotation_key is None:
                annotation_key = "cell_type_auto" if "cell_type_auto" in adata.obs.columns else "cell_type"
            if annotation_key not in adata.obs.columns:
                annotation_key = (
                    "cell_type" if "cell_type" in adata.obs.columns else annotation_key
                )

            run_malignancy_interpretation(
                adata,
                annotation_key=annotation_key,
                cancer_type=cancer_type,
                run_cnv=config.run_cnv,
                run_malignancy_score=True,
            )
            executed_steps.append(step_name)

            malignancy_ns = (
                adata.uns.get("sclucid", {}).get("analysis", {}).get("malignancy", {})
            )
            summary = malignancy_ns.get("malignancy_interpretation_summary", {})
            evidence_sources = summary.get("evidence_sources", [])
            # Upgrade evidence level if multiple independent evidence sources are available
            evidence_level = "heuristic"
            if len(evidence_sources) >= 3:
                evidence_level = "validated_core"
            step_results.append(
                StepResult(
                    name=step_name,
                    status="completed",
                    evidence_level=evidence_level,
                    outputs={
                        "call_key": "malignancy_call",
                        "score_key": "malignancy_interpretation_score",
                        "n_malignant": summary.get("n_malignant"),
                        "n_suspect_malignant": summary.get("n_suspect_malignant"),
                        "tumor_purity_estimate": summary.get("tumor_purity_estimate"),
                        "evidence_sources": evidence_sources,
                        "threshold": summary.get("threshold"),
                        "suspect_threshold": summary.get("suspect_threshold"),
                    },
                    warnings=[],
                )
            )
        except Exception as exc:
            msg = f"Malignancy interpretation failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)
            step_results.append(
                StepResult.from_exception(
                    name=step_name,
                    exc=exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )

    # Therapy
    if config.run_therapy:
        step_name = "therapy_prediction"
        therapy_results: List[StepResult] = []
        try:
            from .therapy.prediction import predict_therapy_response

            log.info("Tumor stage: predicting therapy response")
            drugs = config.therapy_drugs or ["chemotherapy"]
            for drug in drugs:
                drug_step = f"therapy_prediction_{drug}"
                try:
                    predict_therapy_response(
                        adata,
                        therapy_type=drug,
                        method="signature",
                        key_added=f"therapy_response_{drug}",
                    )
                    executed_steps.append(drug_step)
                    therapy_results.append(
                        StepResult(
                            name=drug_step,
                            status="completed",
                            evidence_level="exploratory",
                            outputs={"drug": drug, "key_added": f"therapy_response_{drug}"},
                            warnings=["Signature-based therapy prediction is exploratory."],
                        )
                    )
                except Exception as drug_exc:
                    msg = f"Therapy prediction failed for {drug}: {drug_exc}"
                    log.warning(f"{msg}. Skipping drug.")
                    stage_warnings.append(msg)
                    therapy_results.append(
                        StepResult.from_exception(
                            name=drug_step,
                            exc=drug_exc,
                            degraded=True,
                            evidence_level="unavailable",
                        )
                    )
            # Aggregate therapy step
            if therapy_results:
                n_completed = sum(1 for r in therapy_results if r.status == "completed")
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed" if n_completed == len(therapy_results) else "degraded",
                        evidence_level="exploratory",
                        outputs={
                            "drugs_requested": drugs,
                            "drugs_completed": n_completed,
                            "drug_results": [r.name for r in therapy_results],
                        },
                        warnings=["Therapy prediction uses built-in signatures only."],
                    )
                )
                step_results.extend(therapy_results)
        except Exception as exc:
            msg = f"Therapy response prediction failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)
            step_results.append(
                StepResult.from_exception(
                    name=step_name,
                    exc=exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )

    return adata, executed_steps, step_results, stage_warnings


def _build_tumor_review_summary(
    *,
    adata: AnnData,
    tumor_config: TumorAnalysisConfig,
    step_results: List[StepResult],
    warnings: List[str],
    cancer_type: Optional[str],
) -> Dict[str, Any]:
    """Build a structured tumor review summary for audit and user review."""
    requested_steps = []
    if tumor_config is not None:
        if tumor_config.run_malignancy:
            requested_steps.append("malignancy")
        if tumor_config.run_tme:
            requested_steps.append("tme_deconvolution")
        if tumor_config.run_cnv:
            requested_steps.append("cnv_inference")
        if tumor_config.run_therapy:
            requested_steps.append("therapy_prediction")

    completed_steps = [r.name for r in step_results if r.status == "completed"]
    skipped_or_failed = [r.name for r in step_results if r.status in ("skipped", "failed")]
    degraded_steps = [r.name for r in step_results if r.status == "degraded"]

    evidence_sources: List[str] = []
    for r in step_results:
        evidence_sources.extend(r.outputs.get("evidence_sources", []))
    evidence_sources = list(dict.fromkeys(evidence_sources))

    evidence_rank = {
        "unavailable": 0,
        "exploratory": 1,
        "heuristic": 2,
        "validated_core": 3,
    }
    completed_levels = [
        r.evidence_level for r in step_results if r.status == "completed"
    ]
    if completed_levels:
        claim_boundary = min(completed_levels, key=lambda level: evidence_rank.get(level, 0))
    else:
        claim_boundary = "unavailable"

    readiness_score = 1.0
    readiness_reasons: List[str] = []
    if degraded_steps:
        readiness_score -= 0.15 * len(degraded_steps)
        readiness_reasons.append(f"degraded_steps: {degraded_steps}")
    if skipped_or_failed:
        readiness_score -= 0.25 * len(skipped_or_failed)
        readiness_reasons.append(f"skipped_or_failed_steps: {skipped_or_failed}")
    if not completed_steps:
        readiness_score = 0.0
        readiness_reasons.append("no_completed_tumor_steps")

    readiness_score = max(0.0, min(1.0, readiness_score))
    if skipped_or_failed:
        status = "blocked" if not completed_steps else "degraded"
    elif degraded_steps:
        status = "degraded"
    elif readiness_score >= 0.85:
        status = "ready"
    else:
        status = "review_required"

    action_items: List[Dict[str, Any]] = []
    if degraded_steps:
        action_items.append(
            ReviewAction(
                priority="review",
                action="Review degraded tumor steps before drawing biological conclusions.",
                rationale=f"Steps completed with warnings or fallbacks: {degraded_steps}.",
                evidence_keys=["step_results"],
            ).model_dump(mode="json")
        )
    if skipped_or_failed:
        action_items.append(
            ReviewAction(
                priority="required" if not completed_steps else "review",
                action="Resolve skipped or failed tumor steps.",
                rationale=f"Steps did not complete: {skipped_or_failed}.",
                evidence_keys=["step_results", "warnings"],
            ).model_dump(mode="json")
        )
    if claim_boundary == "exploratory":
        action_items.append(
            ReviewAction(
                priority="optional",
                action="Treat therapy and evolution outputs as exploratory hypotheses.",
                rationale="No external validation database was used.",
                evidence_keys=["claim_boundary"],
            ).model_dump(mode="json")
        )

    malignancy_step = next(
        (r for r in step_results if r.name == "malignancy_interpretation"), None
    )
    tumor_purity = None
    if malignancy_step:
        tumor_purity = malignancy_step.outputs.get("tumor_purity_estimate")

    return {
        "module": "tumor",
        "workflow_name": "tumor_analysis",
        "schema_version": "1.0",
        "cancer_type": cancer_type,
        "requested_steps": requested_steps,
        "completed_steps": completed_steps,
        "degraded_steps": degraded_steps,
        "skipped_or_failed_steps": skipped_or_failed,
        "step_results": summarize_step_results(step_results),
        "evidence_sources": evidence_sources,
        "claim_boundary": claim_boundary,
        "tumor_purity_estimate": tumor_purity,
        "readiness": {
            "status": status,
            "score": round(readiness_score, 3),
            "reasons": readiness_reasons,
        },
        "warnings": warnings,
        "action_items": action_items,
    }


def _diff_recommendations(
    recommendations: Optional[Any],
    actual_configs: Dict[str, Any],
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Compare recommendation parameters against actual config values.

    Returns a nested dict: {stage: {param_name: {"recommended": ..., "actual": ...}}}
    """
    overrides: Dict[str, Dict[str, Dict[str, Any]]] = {}
    if recommendations is None:
        return overrides

    for stage_name, section in recommendations.sections.items():
        actual_cfg = actual_configs.get(stage_name)
        if actual_cfg is None:
            continue
        actual_dict = actual_cfg.to_dict() if hasattr(actual_cfg, "to_dict") else {}
        stage_overrides: Dict[str, Dict[str, Any]] = {}
        for param in section.parameters:
            rec_value = param.value
            actual_value = actual_dict.get(param.name)
            if actual_value is not None and rec_value != actual_value:
                stage_overrides[param.name] = {
                    "recommended": rec_value,
                    "actual": actual_value,
                }
        if stage_overrides:
            overrides[stage_name] = stage_overrides

    return overrides


def _apply_qc_recommendations(qc_config: Any, section: Any) -> Any:
    """Apply recommended QC thresholds to a QCWorkflowConfig."""
    thresholds = section.to_dict().get("parameters", [])
    threshold_map = {p["name"]: p["value"] for p in thresholds}
    if hasattr(qc_config, "filter_config") and qc_config.filter_config is not None:
        fc = qc_config.filter_config
        for key, value in threshold_map.items():
            if hasattr(fc, key) and value is not None:
                setattr(fc, key, value)
    return qc_config


def _apply_preprocess_recommendations(preprocess_config: Any, section: Any) -> Any:
    """Apply recommended preprocessing parameters."""
    raw = section.raw_result
    if raw is not None and hasattr(raw, "to_config"):
        try:
            return raw.to_config(base_config=preprocess_config)
        except Exception:
            pass
    # Fallback: direct parameter assignment
    for param in section.parameters:
        if hasattr(preprocess_config, param.name) and param.value is not None:
            setattr(preprocess_config, param.name, param.value)
    return preprocess_config


def _apply_clustering_recommendations(analysis_config: Any, section: Any) -> Any:
    """Apply recommended clustering parameters to analysis config."""
    if not hasattr(analysis_config, "clustering") or analysis_config.clustering is None:
        from ..analysis.config import ClusteringConfig

        analysis_config.clustering = ClusteringConfig()
    for param in section.parameters:
        if hasattr(analysis_config.clustering, param.name) and param.value is not None:
            setattr(analysis_config.clustering, param.name, param.value)
    return analysis_config


def _apply_annotation_recommendations(analysis_config: Any, section: Any) -> Any:
    """Apply recommended annotation parameters to analysis config."""
    if not hasattr(analysis_config, "annotation") or analysis_config.annotation is None:
        from ..analysis.config import AnnotationConfig

        analysis_config.annotation = AnnotationConfig()
    raw = section.raw_result
    if isinstance(raw, type(analysis_config.annotation)):
        analysis_config.annotation = raw.model_copy()
    else:
        for param in section.parameters:
            if hasattr(analysis_config.annotation, param.name) and param.value is not None:
                setattr(analysis_config.annotation, param.name, param.value)
        for key in ["cluster_key", "marker_species", "marker_tissue", "key_added"]:
            if key in section.metadata and hasattr(analysis_config.annotation, key):
                setattr(analysis_config.annotation, key, section.metadata[key])
    return analysis_config


def _apply_tumor_recommendations(
    tumor_config: TumorAnalysisConfig, section: Any
) -> TumorAnalysisConfig:
    """Apply recommended tumor parameters to TumorAnalysisConfig."""
    raw = section.raw_result
    if isinstance(raw, type(tumor_config)):
        return raw.model_copy()
    for param in section.parameters:
        if hasattr(tumor_config, param.name) and param.value is not None:
            setattr(tumor_config, param.name, param.value)
    return tumor_config
