"""
High-level analysis workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for comprehensive analysis including:
- Clustering and dimensionality reduction
- Cell type annotation
- Marker gene finding and characterization
- Differential expression analysis

Note: This is a convenience wrapper around individual analysis modules.
For more control, use the individual functions directly.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
from anndata import AnnData

from ..base_config import apply_config_overrides
from ..utils import (
    PartialResultManager,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    export_review_summary,
    get_marker_manager,
    get_progress_bar,
    normalize_review_summary,
    validate_review_summary_schema,
)
from .annotation import build_annotation_consensus, run_annotation, run_annotation_evidence
from .clustering import cluster_cells, run_clustering_review
from .config import AnalysisWorkflowConfig, AnnotationConfig
from .differential_expression import characterize_clusters, find_markers
from .malignancy import run_malignancy_interpretation
from .scoring import score_by_gene_sets
from .trace import enrich_analysis_review_summary, validate_analysis_review_summary

log = logging.getLogger(__name__)

# Define workflow steps for flexible execution
ANALYSIS_WORKFLOW_STEPS = [
    "clustering_review",
    "clustering",
    "markers",
    "annotation",
    "annotation_evidence",
    "annotation_consensus",
    "malignancy_interpretation",
    "characterization",
]

# Keep for backward compatibility
AnalysisWorkflowError = WorkflowError
PartialAnalysisResult = PartialResultManager

__all__ = [
    "run_standard_analysis",
    "run_custom_analysis",
    "compare_clustering_resolutions",
    "AnalysisWorkflowError",
    "PartialAnalysisResult",
    "ANALYSIS_WORKFLOW_STEPS",
]


def _resolve_analysis_steps(
    steps: Optional[List[str]],
    skip_steps: Optional[List[str]],
    config: Optional[AnalysisWorkflowConfig] = None,
    completed_steps: Optional[List[str]] = None,
) -> List[str]:
    """Resolve which analysis steps to run."""
    if steps is not None and skip_steps is not None:
        raise ValueError("Cannot specify both 'steps' and 'skip_steps'. Choose one.")

    if steps is not None:
        resolved = list(steps)
    elif skip_steps is not None:
        resolved = [s for s in ANALYSIS_WORKFLOW_STEPS if s not in skip_steps]
    else:
        # Use config flags to determine default steps
        resolved = []
        if config is not None and getattr(config, "run_clustering_review", False):
            resolved.append("clustering_review")
        if config is None or config.clustering is not None:
            resolved.append("clustering")
        if config is None or getattr(config, "find_markers", True):
            resolved.append("markers")
        if config is None or config.annotation is not None:
            resolved.append("annotation")
            if getattr(config, "run_annotation_evidence", True):
                resolved.append("annotation_evidence")
                if getattr(config, "final_annotation_strategy", "consensus") == "consensus":
                    resolved.append("annotation_consensus")
        if config is not None and getattr(config, "run_malignancy_interpretation", False):
            resolved.append("malignancy_interpretation")
        if config is None or getattr(config, "characterize", True):
            resolved.append("characterization")

    invalid = set(resolved) - set(ANALYSIS_WORKFLOW_STEPS)
    if invalid:
        raise ValueError(
            f"Invalid step names: {invalid}. Valid steps are: {ANALYSIS_WORKFLOW_STEPS}"
        )

    if completed_steps:
        resolved = [s for s in resolved if s not in completed_steps]

    return resolved


def _default_groupby_key(adata: AnnData) -> str:
    """Choose the most likely cluster key for downstream analysis steps."""
    if "leiden_clusters" in adata.obs.columns:
        return "leiden_clusters"
    if "leiden" in adata.obs.columns:
        return "leiden"
    return "leiden_clusters"


def run_standard_analysis(
    adata: AnnData,
    config: Optional[AnalysisWorkflowConfig] = None,
    *,
    show_progress: bool = True,
    # Step control
    steps: Optional[List[str]] = None,
    skip_steps: Optional[List[str]] = None,
    # Error recovery
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: str = "raise",
    resume_from: Optional[str] = None,
    **kwargs,
) -> AnnData:
    """
    Run a standard analysis pipeline from clustering to annotation.

    This workflow executes:
    1. Clustering (with automatic resolution selection if needed)
    2. Marker gene identification
    3. Cell type annotation
    4. Cluster characterization

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Should have preprocessing completed (normalized, HVGs, PCA).
    config : AnalysisWorkflowConfig, optional
        Analysis workflow configuration. If None, uses defaults.
    show_progress : bool, default=True
        Show progress bar for workflow steps.
    steps : list of str, optional
        Specific steps to run. See ``ANALYSIS_WORKFLOW_STEPS`` for valid names.
    skip_steps : list of str, optional
        Steps to skip (alternative to specifying ``steps``).
    error_recovery : bool, default=False
        Enable error recovery mode.
    recovery_save_dir : str, optional
        Directory to save partial results on error.
    on_error : {"raise", "skip", "save"}, default="raise"
        How to handle errors.
    resume_from : str, optional
        Path to partial results directory to resume from.
    **kwargs
        Additional parameters to override config values

    Returns:
    -------
    AnnData
        Annotated data with analysis results stored in:
        - adata.obs['leiden']: Cluster labels
        - adata.obs['cell_type']: Cell type annotations
        - adata.uns['markers']: Marker genes
        - adata.uns['characterization']: Cluster characterization

    Examples:
    --------
    >>> # Use defaults with progress bar
    >>> adata = run_standard_analysis(adata, show_progress=True)

    >>> # Skip marker finding and characterization
    >>> adata = run_standard_analysis(adata, skip_steps=["markers", "characterization"])

    >>> # Run only clustering
    >>> adata = run_standard_analysis(adata, steps=["clustering"])

    >>> # Error recovery mode
    >>> adata = run_standard_analysis(
    ...     adata,
    ...     error_recovery=True,
    ...     recovery_save_dir="./recovery",
    ...     on_error="save"
    ... )

    >>> # Resume from partial results
    >>> adata = run_standard_analysis(
    ...     adata,
    ...     resume_from="./recovery",
    ...     show_progress=True
    ... )

    Notes:
    -----
    This function requires that preprocessing has been completed:
    - Normalization
    - HVG selection
    - PCA
    - Neighborhood graph

    Use `scLucid.preprocess.run_preprocessing` to prepare data.
    """
    if config is None:
        from .config import AnalysisWorkflowConfig as DefaultConfig

        config = apply_config_overrides(DefaultConfig(), **kwargs)
    else:
        config = apply_config_overrides(config, **kwargs)

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

    log.info("=" * 60)
    log.info("=== Starting Standard Analysis Pipeline ===")
    log.info("=" * 60)
    log.info(f"Show progress: {show_progress}")
    log.info(f"Error recovery: {error_recovery}")

    # Resolve steps
    steps_to_run = _resolve_analysis_steps(steps, skip_steps, config, completed_steps)
    log.info(f"Steps to run: {steps_to_run}")

    # Initialize progress bar
    step_iterator = get_progress_bar(
        steps_to_run, desc="Analysis", enabled=show_progress, total=len(steps_to_run), unit="step"
    )

    # Track execution
    current_step = None
    successful_steps: List[str] = []
    markers_df: Optional[pd.DataFrame] = None
    annotation_review_table: Optional[pd.DataFrame] = None

    # Determine cluster key from config or existing analysis state. Downstream
    # marker/characterization steps require this key if clustering is skipped.
    from .config import ClusteringConfig

    if isinstance(config.clustering, ClusteringConfig):
        cluster_config = config.clustering
    elif isinstance(config.clustering, dict):
        cluster_config = ClusteringConfig(**config.clustering)
    else:
        cluster_config = ClusteringConfig()
    cluster_key = cluster_config.key_added or f"{cluster_config.method}_clusters"

    def _require_cluster_key(step: str) -> None:
        if cluster_key not in adata.obs.columns:
            existing_key = _default_groupby_key(adata)
            if existing_key in adata.obs.columns:
                log.info(
                    f"Using existing cluster key '{existing_key}' for step '{step}' "
                    f"instead of configured key '{cluster_key}'."
                )
                return
            raise ValueError(
                f"Step '{step}' requires clustering results, but '{cluster_key}' "
                "is not present in adata.obs. Include the 'clustering' step or set "
                "config.clustering.key_added to an existing cluster column."
            )

    try:
        for step_name in step_iterator:
            current_step = step_name

            # Step 0: Clustering resolution evidence
            if step_name == "clustering_review":
                log.info("Step: Clustering resolution evidence")
                review_df = run_clustering_review(
                    adata,
                    resolutions=getattr(config, "candidate_resolutions", None),
                    method=cluster_config.method
                    if cluster_config.method in {"leiden", "louvain"}
                    else "leiden",
                    use_rep=cluster_config.use_rep,
                    random_state=cluster_config.random_state,
                    de_method=config.de.method if config.de else "wilcoxon",
                )
                clustering_ns = (
                    adata.uns.get("sclucid", {}).get("analysis", {}).get("clustering", {})
                )
                review_summary = clustering_ns.get("clustering_review_summary", {})
                recommended = review_summary.get("recommended_resolution")
                if (
                    getattr(config, "use_recommended_resolution", True)
                    and recommended is not None
                    and cluster_config.method in {"leiden", "louvain"}
                ):
                    cluster_config = cluster_config.model_copy(
                        update={"resolution": float(recommended)}
                    )
                    config.clustering = cluster_config
                    log.info(
                        "  Using recommended first-pass resolution: "
                        f"{cluster_config.resolution:g}"
                    )
                log.info(f"  Reviewed {len(review_df)} clustering resolution candidate(s)")
                successful_steps.append(step_name)

            # Step 1: Clustering
            elif step_name == "clustering":
                log.info("Step: Clustering")
                adata = cluster_cells(adata, cluster_config)
                log.info(f"  Clustering complete: {adata.obs[cluster_key].nunique()} clusters")
                successful_steps.append(step_name)

            # Step 2: Marker genes
            elif step_name == "markers":
                log.info("Step: Finding marker genes")
                from .config import DifferentialConfig

                _require_cluster_key(step_name)
                active_cluster_key = (
                    cluster_key if cluster_key in adata.obs.columns else _default_groupby_key(adata)
                )

                marker_config = DifferentialConfig(
                    groupby=active_cluster_key,
                    method=config.marker_method if hasattr(config, "marker_method") else "wilcoxon",
                )
                markers_df = find_markers(adata, marker_config)
                log.info(f"  Found {len(markers_df)} marker rows")
                successful_steps.append(step_name)

            # Step 3: Annotation
            elif step_name == "annotation":
                log.info("Step: Cell type annotation")

                if isinstance(config.annotation, AnnotationConfig):
                    adata = run_annotation(adata, config=config.annotation)
                elif isinstance(config.annotation, dict):
                    adata = run_annotation(adata, config=AnnotationConfig(**config.annotation))
                else:
                    adata = run_annotation(adata, config=AnnotationConfig())
                n_annotated = (
                    adata.obs["cell_type"].notna().sum() if "cell_type" in adata.obs else 0
                )
                log.info(f"  Annotated {n_annotated}/{len(adata)} cells")
                successful_steps.append(step_name)

            # Step 4: Annotation evidence table
            elif step_name == "annotation_evidence":
                log.info("Step: Annotation evidence")
                _require_cluster_key(step_name)
                active_cluster_key = (
                    cluster_key if cluster_key in adata.obs.columns else _default_groupby_key(adata)
                )
                annotation_config = (
                    config.annotation
                    if isinstance(config.annotation, AnnotationConfig)
                    else AnnotationConfig(**config.annotation)
                    if isinstance(config.annotation, dict)
                    else AnnotationConfig()
                )
                marker_config_path = (
                    annotation_config.lineage_marker_config
                    or annotation_config.subtype_marker_config
                    or annotation_config.state_marker_config
                )
                marker_manager = None
                if marker_config_path:
                    marker_manager = marker_config_path
                elif "marker_manager" in tuple(getattr(config, "annotation_methods", ()) or ()):
                    marker_manager = get_marker_manager(
                        species=annotation_config.marker_species,
                        tissue=annotation_config.marker_tissue,
                        view="lineage_annotation",
                    )
                reference_key = (
                    annotation_config.key_added
                    if annotation_config.key_added in adata.obs.columns
                    else None
                )
                active_annotation_methods = tuple(
                    getattr(config, "annotation_methods", ()) or ()
                )
                if reference_key is None and "celltypist" in active_annotation_methods:
                    for candidate in (
                        "celltypist_majority_voting",
                        "celltypist_predicted_labels",
                    ):
                        if candidate in adata.obs.columns:
                            reference_key = candidate
                            break
                confidence_key = (
                    f"{reference_key}_confidence"
                    if reference_key and f"{reference_key}_confidence" in adata.obs.columns
                    else None
                )
                if reference_key == "celltypist_predicted_labels":
                    confidence_key = (
                        "celltypist_conf_score"
                        if "celltypist_conf_score" in adata.obs.columns
                        else confidence_key
                    )
                llm_annotations = getattr(config, "llm_annotations", None)
                if isinstance(llm_annotations, list):
                    llm_annotations = pd.DataFrame(llm_annotations)
                annotation_review_table = run_annotation_evidence(
                    adata,
                    active_cluster_key,
                    markers_df=markers_df,
                    methods=active_annotation_methods,
                    marker_config=marker_manager,
                    reference_key=reference_key,
                    reference_confidence_key=confidence_key,
                    llm_annotations=llm_annotations,
                )
                log.info(
                    "  Built annotation evidence table with "
                    f"{annotation_review_table.shape[0]} cluster rows"
                )
                successful_steps.append(step_name)

            # Step 5: Annotation consensus application
            elif step_name == "annotation_consensus":
                log.info("Step: Annotation consensus")
                _require_cluster_key(step_name)
                active_cluster_key = (
                    cluster_key if cluster_key in adata.obs.columns else _default_groupby_key(adata)
                )
                annotation_config = (
                    config.annotation
                    if isinstance(config.annotation, AnnotationConfig)
                    else AnnotationConfig(**config.annotation)
                    if isinstance(config.annotation, dict)
                    else AnnotationConfig()
                )
                if annotation_review_table is None:
                    annotation_review_table = (
                        adata.uns.get("sclucid", {})
                        .get("analysis", {})
                        .get("annotation", {})
                        .get("annotation_review_table")
                    )
                build_annotation_consensus(
                    adata,
                    active_cluster_key,
                    annotation_review_table,
                    key_added=annotation_config.key_added,
                    lineage_key=annotation_config.lineage_key,
                )
                log.info(f"  Applied consensus labels to obs['{annotation_config.key_added}']")
                successful_steps.append(step_name)

            # Step 6: Malignancy interpretation
            elif step_name == "malignancy_interpretation":
                log.info("Step: Malignancy interpretation")
                active_cluster_key = (
                    cluster_key if cluster_key in adata.obs.columns else _default_groupby_key(adata)
                )
                annotation_config = (
                    config.annotation
                    if isinstance(config.annotation, AnnotationConfig)
                    else AnnotationConfig(**config.annotation)
                    if isinstance(config.annotation, dict)
                    else AnnotationConfig()
                )
                if annotation_config.key_added not in adata.obs.columns:
                    raise ValueError(
                        "Step 'malignancy_interpretation' requires final annotation "
                        f"obs['{annotation_config.key_added}']. Include annotation_consensus "
                        "or set annotation.key_added to an existing column."
                    )
                malignancy_table = run_malignancy_interpretation(
                    adata,
                    annotation_key=annotation_config.key_added,
                    cluster_key=active_cluster_key
                    if active_cluster_key in adata.obs.columns
                    else None,
                    species=annotation_config.marker_species,
                    cancer_type=getattr(config, "malignancy_cancer_type", None),
                    run_cnv=getattr(config, "run_cnv_for_malignancy", False),
                    cnv_score_key=getattr(config, "malignancy_cnv_score_key", None),
                    reference_labels=getattr(config, "malignancy_reference_labels", None),
                    run_malignancy_score=getattr(config, "run_malignancy_score", True),
                    key_added=getattr(config, "malignancy_key_added", "malignancy_call"),
                    score_key=getattr(
                        config,
                        "malignancy_score_key",
                        "malignancy_interpretation_score",
                    ),
                    threshold=getattr(config, "malignancy_threshold", 0.55),
                    suspect_threshold=getattr(config, "malignancy_suspect_threshold", 0.35),
                )
                log.info(
                    "  Built malignancy interpretation table with "
                    f"{malignancy_table.shape[0]} group rows"
                )
                successful_steps.append(step_name)

            # Step 7: Characterization
            elif step_name == "characterization":
                log.info("Step: Cluster characterization")
                try:
                    _require_cluster_key(step_name)
                    active_cluster_key = (
                        cluster_key
                        if cluster_key in adata.obs.columns
                        else _default_groupby_key(adata)
                    )
                    adata = characterize_clusters(
                        adata,
                        groupby=active_cluster_key,
                    )
                    log.info("  Characterization complete")
                    successful_steps.append(step_name)
                except Exception as e:
                    if on_error == "skip":
                        log.warning(f"  Characterization failed: {e}. Skipping...")
                    else:
                        raise

    except Exception as e:
        error_msg = f"Workflow failed at step '{current_step}': {str(e)}"
        log.error(error_msg)
        import traceback

        log.error(traceback.format_exc())

        if error_recovery and on_error in ["raise", "save"]:
            # Save partial results
            save_dir = recovery_save_dir or "./recovery"
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(
                completed_steps=successful_steps,
                failed_step=current_step,
                error_message=str(e),
            )
            manager.save(adata, checkpoint, config)

            if on_error == "save":
                log.warning(f"Workflow failed but partial results saved to: {save_dir}")
                log.warning(
                    f"To resume, use: run_standard_analysis(adata, resume_from='{save_dir}')"
                )
                return adata

        raise WorkflowError(
            f"[analysis] Workflow failed at step '{current_step}': {e}",
            step_name=current_step or "unknown",
            original_error=e,
        )

    # Store final config
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})[
        UnsKeys.WORKFLOW_CONFIG
    ] = config.to_dict()
    adata.uns["sclucid"]["analysis"][UnsKeys.STEPS_EXECUTED] = successful_steps

    # Build and store review summary
    enriched_summary = enrich_analysis_review_summary(
        _build_analysis_review_summary(adata, config, successful_steps, cluster_key),
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
    )
    review_summary = normalize_review_summary(
        enriched_summary,
        module="analysis",
        workflow_name="standard",
        adata=adata,
        steps_executed=successful_steps,
        config=config.to_dict(),
        warnings=(
            enriched_summary.get("analysis_readiness", {}).get("review_reasons", [])
            if isinstance(enriched_summary.get("analysis_readiness"), dict)
            else []
        ),
    )
    validate_review_summary_schema(review_summary, module="analysis", raise_on_error=True)
    validate_analysis_review_summary(review_summary, raise_on_error=True)
    adata.uns["sclucid"]["analysis"][UnsKeys.REVIEW_SUMMARY] = review_summary

    # Export review summary to file if save_dir is configured
    if config.save_dir:
        export_review_summary(
            review_summary,
            save_dir=config.save_dir,
            module="analysis",
            title="Analysis Review Summary",
            adata=adata,
        )

    log.info("=" * 60)
    log.info("=== Standard Analysis Pipeline Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)
    return adata


def _build_analysis_review_summary(
    adata: AnnData,
    config: AnalysisWorkflowConfig,
    successful_steps: List[str],
    cluster_key: str,
) -> Dict[str, Any]:
    """Build a human-reviewable summary of the analysis run."""
    summary: Dict[str, Any] = {
        "module": "analysis",
        "workflow_name": "standard",
        "steps_executed": successful_steps,
        "cluster_key": cluster_key,
        "warnings": [],
        "artifacts": {},
    }

    # Clustering summary
    if "clustering" in successful_steps and cluster_key in adata.obs.columns:
        n_clusters = adata.obs[cluster_key].nunique()
        summary["clustering"] = {
            "n_clusters": int(n_clusters),
            "method": config.clustering.method if config.clustering else "unknown",
            "resolution": round(config.clustering.resolution, 2) if config.clustering else None,
            "use_rep": config.clustering.use_rep if config.clustering else "unknown",
        }
        clustering_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("clustering", {})
        if "clustering_review_summary" in clustering_ns:
            summary["clustering"]["review"] = clustering_ns["clustering_review_summary"]
            summary["artifacts"][
                "clustering_review"
            ] = 'adata.uns["sclucid"]["analysis"]["clustering"]["clustering_review"]'

    # Marker summary
    if "markers" in successful_steps:
        de_key = "rank_genes_groups"
        if de_key in adata.uns:
            n_groups = len(adata.uns[de_key].get("names", []))
            summary["markers"] = {
                "de_key": de_key,
                "n_groups": n_groups,
                "method": config.de.method if config.de else "wilcoxon",
            }

    # Annotation summary
    if "annotation" in successful_steps:
        annotation_key = config.annotation.key_added if config.annotation else "cell_type_auto"
        if annotation_key in adata.obs.columns:
            n_annotated = adata.obs[annotation_key].notna().sum()
            n_types = adata.obs[annotation_key].nunique()
            summary["annotation"] = {
                "key": annotation_key,
                "n_annotated": int(n_annotated),
                "n_cell_types": int(n_types),
                "method": config.annotation.final_method if config.annotation else "unknown",
            }
            annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
            if "annotation_review_table" in annotation_ns:
                review_df = annotation_ns["annotation_review_table"]
                summary["annotation"]["review_table_rows"] = (
                    int(review_df.shape[0]) if hasattr(review_df, "shape") else None
                )
                summary["artifacts"][
                    "annotation_review_table"
                ] = 'adata.uns["sclucid"]["analysis"]["annotation"]["annotation_review_table"]'
            confidence_key = f"{annotation_key}_confidence"
            if confidence_key in adata.obs.columns:
                low_conf = pd.to_numeric(adata.obs[confidence_key], errors="coerce") < 0.5
                if bool(low_conf.any()):
                    summary["warnings"].append("low_confidence_annotation_cells_present")

    # Characterization summary
    if "characterization" in successful_steps:
        summary["characterization"] = {
            "status": "completed",
            "groupby": cluster_key,
        }

    return summary


def run_custom_analysis(
    adata: AnnData,
    steps: List[str],
    step_configs: Optional[Dict[str, dict]] = None,
    save_dir: Optional[Union[str, Path]] = None,
    *,
    show_progress: bool = True,
) -> AnnData:
    """
    Run a custom analysis pipeline with specified steps.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    steps : list of str
        Analysis steps to run. Available steps:
        - 'resolution': Find optimal clustering resolution
        - 'clustering': Perform clustering
        - 'markers': Find marker genes
        - 'annotation': Annotate cell types
        - 'scoring': Score by gene sets
        - 'characterization': Full cluster characterization
    step_configs : dict, optional
        Configuration for each step (step_name -> config_dict)
    save_dir : str or Path, optional
        Directory to save results
    show_progress : bool, default=True
        Show progress bar for workflow steps.

    Returns:
    -------
    AnnData
        Annotated data with analysis results

    Examples:
    --------
    >>> from scLucid.analysis import run_custom_analysis
    >>>
    >>> # Run only clustering and markers with progress bar
    >>> adata = run_custom_analysis(
    ...     adata,
    ...     steps=['clustering', 'markers'],
    ...     step_configs={
    ...         'clustering': {'resolution': 0.8}
    ...     },
    ...     show_progress=True
    ... )
    """
    if step_configs is None:
        step_configs = {}

    log.info(f"Running custom analysis with {len(steps)} steps...")

    # Initialize progress bar
    step_iterator = get_progress_bar(
        steps, desc="Custom Analysis", enabled=show_progress, total=len(steps), unit="step"
    )

    for i, step in enumerate(step_iterator, 1):
        log.info(f"Step {i}/{len(steps)}: {step}")

        if step in {"resolution", "clustering_review"}:
            config = step_configs.get(step, {})
            run_clustering_review(adata, **config)
            summary = (
                adata.uns.get("sclucid", {})
                .get("analysis", {})
                .get("clustering", {})
                .get("clustering_review_summary", {})
            )
            log.info(f"  Recommended resolution: {summary.get('recommended_resolution')}")

        elif step == "clustering":
            from .config import ClusteringConfig

            config = step_configs.get(step, {})
            if isinstance(config, dict):
                config = ClusteringConfig(**config)
            adata = cluster_cells(adata, config)
            cluster_key = config.key_added or f"{config.method}_clusters"
            log.info(f"  Created {adata.obs[cluster_key].nunique()} clusters")

        elif step == "markers":
            from .config import DifferentialConfig

            config = step_configs.get(step, {})
            if "groupby" not in config:
                config["groupby"] = _default_groupby_key(adata)
            if isinstance(config, dict):
                config = DifferentialConfig(**config)
            markers = find_markers(adata, config)
            log.info(f"  Found markers for {len(markers)} clusters")

        elif step == "annotation":
            config = step_configs.get(step, {})
            adata = run_annotation(adata, **config)
            log.info("  Annotation complete")

        elif step == "scoring":
            config = step_configs.get(step, {})
            if "gene_sets" not in config:
                log.warning("  No gene_sets provided, skipping")
                continue
            adata = score_by_gene_sets(adata, **config)
            log.info("  Scoring complete")

        elif step == "characterization":
            config = step_configs.get(step, {})
            if "groupby" not in config:
                config["groupby"] = _default_groupby_key(adata)
            adata = characterize_clusters(adata, save_path=save_dir, **config)
            log.info("  Characterization complete")

        else:
            log.warning(f"  Unknown step: {step}")

    log.info("Custom analysis complete!")
    return adata


def compare_clustering_resolutions(
    adata: AnnData,
    resolutions: List[float],
    metrics: Optional[List[str]] = None,
    save_path: Optional[Union[str, Path]] = None,
    *,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Compare multiple clustering resolutions.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    resolutions : list of float
        Resolutions to test
    metrics : list of str, optional
        Metrics to compute (default: ['n_clusters', 'silhouette'])
    save_path : str or Path, optional
        Path to save comparison results
    show_progress : bool, default=True
        Show progress bar for resolution testing.

    Returns:
    -------
    DataFrame
        Comparison of metrics across resolutions

    Examples:
    --------
    >>> from scLucid.analysis.workflow import compare_clustering_resolutions
    >>>
    >>> results = compare_clustering_resolutions(
    ...     adata,
    ...     resolutions=[0.4, 0.6, 0.8, 1.0, 1.2],
    ...     show_progress=True
    ... )
    >>> print(results)
    """
    from sklearn.metrics import silhouette_score

    from .config import ClusteringConfig

    if metrics is None:
        metrics = ["n_clusters", "silhouette"]

    results = []

    # Initialize progress bar
    res_iterator = get_progress_bar(
        resolutions,
        desc="Resolution Search",
        enabled=show_progress,
        total=len(resolutions),
        unit="res",
    )

    for res in res_iterator:
        log.info(f"Testing resolution: {res}")

        # Cluster at this resolution
        config = ClusteringConfig(resolution=res)
        adata_temp = cluster_cells(adata.copy(), config)

        # Compute metrics
        result = {"resolution": res}
        cluster_key = config.key_added or f"{config.method}_clusters"
        result["n_clusters"] = adata_temp.obs[cluster_key].nunique()

        if "silhouette" in metrics and adata_temp.obsm.get("X_pca") is not None:
            try:
                score = silhouette_score(
                    adata_temp.obsm["X_pca"], adata_temp.obs[cluster_key].astype(int)
                )
                result["silhouette"] = score
            except Exception as e:
                log.warning(f"  Could not compute silhouette: {e}")

        results.append(result)

    df_results = pd.DataFrame(results)

    if save_path:
        df_results.to_csv(save_path, index=False)
        log.info(f"Results saved to: {save_path}")

    return df_results
