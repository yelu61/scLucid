"""
High-level tumor analysis workflow for scLucid.

Provides a dual-layer API:
- run_tumor_analysis: compact high-level entry point for wet-lab users
- run_tumor_analysis_expert: explicit per-stage config override for bioinformatics users
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from anndata import AnnData

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
    qc_config: Optional[Any] = None,
    preprocess_config: Optional[Any] = None,
    analysis_config: Optional[Any] = None,
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
        adata, tumor_steps, tumor_warnings = _run_tumor_stage(adata, tumor_config)
        steps_executed.extend(tumor_steps)
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

    execution_trace = {
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
        "tissue_type": tissue_type,
        "batch_key": batch_key,
        "cancer_type": cancer_type,
    }

    save_result(adata, "tumor", "execution_trace", execution_trace)
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
    log.info("=" * 60)

    return adata


def _run_tumor_stage(
    adata: AnnData,
    config: TumorAnalysisConfig,
) -> tuple[AnnData, list[str], list[str]]:
    """
    Run tumor-specific analysis steps.

    Each step is wrapped in try/except so failures degrade gracefully.

    Returns:
    -------
    tuple of (AnnData, executed_steps, warnings)
    """
    executed_steps: list[str] = []
    stage_warnings: list[str] = []

    # Malignancy
    if config.run_malignancy:
        try:
            from .malignancy.classification import classify_malignant_cells
            from .malignancy.scoring import score_malignancy

            log.info("Tumor stage: scoring malignancy")
            adata = score_malignancy(adata, key_added="malignancy")
            executed_steps.append("malignancy_scoring")

            log.info("Tumor stage: classifying malignant cells")
            if config.malignancy_method == "cnv":
                classify_malignant_cells(adata, method="cnv", key_added="is_malignant")
            elif config.malignancy_method in ("threshold", "ml"):
                # reference key may be used to subset reference cells
                ref_key = config.malignancy_reference_key
                if ref_key and ref_key in adata.obs.columns:
                    # current classify_malignant_cells takes reference_adata, not key
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
                    ref_adata = None
                classify_malignant_cells(
                    adata,
                    method=config.malignancy_method,
                    reference_adata=ref_adata,
                    key_added="is_malignant",
                )
            executed_steps.append("malignancy_classification")
        except Exception as exc:
            msg = f"Malignancy analysis failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)

    # TME
    if config.run_tme:
        try:
            from .microenvironment.deconvolution import deconvolve_tme

            log.info("Tumor stage: deconvolving TME")
            cell_type_key = config.tme_cell_type_key
            if cell_type_key not in adata.obs.columns:
                log.warning(f"TME cell type key '{cell_type_key}' not found. Trying 'cell_type'.")
                cell_type_key = "cell_type"
            adata = deconvolve_tme(adata, cell_type_key=cell_type_key, key_added="tme")
            executed_steps.append("tme_deconvolution")
        except Exception as exc:
            msg = f"TME deconvolution failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)

    # CNV
    if config.run_cnv:
        try:
            from .cnv.infercnv import infer_cnv

            log.info("Tumor stage: inferring CNV")
            ref_key = config.cnv_reference_key
            if ref_key and ref_key in adata.obs.columns:
                ref_cells = (
                    adata.obs[ref_key]
                    .astype(str)
                    .str.lower()
                    .isin({"normal", "healthy", "reference", "immune", "stromal"})
                )
                ref_values = adata.obs[ref_key].unique().tolist() if not ref_cells.any() else None
                if ref_values is None:
                    ref_mask = (
                        adata.obs[ref_key]
                        .astype(str)
                        .str.lower()
                        .isin({"normal", "healthy", "reference", "immune", "stromal"})
                    )
                    ref_values = adata.obs.loc[ref_mask, ref_key].unique().tolist()
                infer_cnv(adata, reference_cells=ref_values, reference_key=ref_key, key_added="cnv")
            else:
                infer_cnv(adata, key_added="cnv")
            executed_steps.append("cnv_inference")
        except Exception as exc:
            msg = f"CNV inference failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)

    # Therapy
    if config.run_therapy:
        try:
            from .therapy.prediction import predict_therapy_response

            log.info("Tumor stage: predicting therapy response")
            drugs = config.therapy_drugs or ["chemotherapy"]
            for drug in drugs:
                try:
                    predict_therapy_response(
                        adata,
                        therapy_type=drug,
                        method="signature",
                        key_added=f"therapy_response_{drug}",
                    )
                    executed_steps.append(f"therapy_prediction_{drug}")
                except Exception as drug_exc:
                    msg = f"Therapy prediction failed for {drug}: {drug_exc}"
                    log.warning(f"{msg}. Skipping drug.")
                    stage_warnings.append(msg)
        except Exception as exc:
            msg = f"Therapy response prediction failed: {exc}"
            log.warning(f"{msg}. Skipping.")
            stage_warnings.append(msg)

    return adata, executed_steps, stage_warnings


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
