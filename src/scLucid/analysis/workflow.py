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
    get_progress_bar,
    normalize_review_summary,
    validate_review_summary_schema,
)
from .annotation import run_annotation
from .clustering import cluster_cells, find_resolution
from .config import AnalysisWorkflowConfig
from .differential_expression import characterize_clusters, find_markers
from .scoring import score_by_gene_sets

log = logging.getLogger(__name__)

# Define workflow steps for flexible execution
ANALYSIS_WORKFLOW_STEPS = [
    "clustering",
    "markers",
    "annotation",
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
        if config is None or config.clustering is not None:
            resolved.append("clustering")
        if config is None or getattr(config, "find_markers", True):
            resolved.append("markers")
        if config is None or config.annotation is not None:
            resolved.append("annotation")
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

            # Step 1: Clustering
            if step_name == "clustering":
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
                markers = find_markers(adata, marker_config)
                log.info(f"  Found markers for {len(markers)} clusters")
                successful_steps.append(step_name)

            # Step 3: Annotation
            elif step_name == "annotation":
                log.info("Step: Cell type annotation")
                from .config import AnnotationConfig

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

            # Step 4: Characterization
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
    review_summary = normalize_review_summary(
        _build_analysis_review_summary(adata, config, successful_steps, cluster_key),
        module="analysis",
        workflow_name="standard",
        adata=adata,
        steps_executed=successful_steps,
        config=config.to_dict(),
    )
    validate_review_summary_schema(review_summary, module="analysis", raise_on_error=True)
    adata.uns["sclucid"]["analysis"][UnsKeys.REVIEW_SUMMARY] = review_summary

    # Export review summary to file if save_dir is configured
    if config.save_dir:
        export_review_summary(
            review_summary,
            save_dir=config.save_dir,
            module="analysis",
            title="Analysis Review Summary",
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
        "steps_executed": successful_steps,
        "cluster_key": cluster_key,
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

        if step == "resolution":
            config = step_configs.get(step, {})
            eval_df, recommended_res = find_resolution(adata, **config)
            log.info(f"  Optimal resolution: {recommended_res}")

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
