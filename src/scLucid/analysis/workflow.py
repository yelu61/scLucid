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
from typing import Dict, List, Optional, Union

from anndata import AnnData
import pandas as pd

from ..base_config import apply_config_overrides
from ..utils import get_progress_bar, PartialResultManager, WorkflowCheckpoint, WorkflowError
from .config import AnalysisWorkflowConfig
from .clustering import find_resolution, cluster_cells
from .annotation import run_annotation
from .differential_expression import find_markers, characterize_clusters
from .scoring import score_by_gene_sets

log = logging.getLogger(__name__)

# Keep for backward compatibility
AnalysisWorkflowError = WorkflowError
PartialAnalysisResult = PartialResultManager

__all__ = [
    "run_standard_analysis",
    "run_custom_analysis",
    "compare_clustering_resolutions",
    "AnalysisWorkflowError",
    "PartialAnalysisResult",
]


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
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: str = "raise",
    resume_from: Optional[str] = None,
    **kwargs
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

    Returns
    -------
    AnnData
        Annotated data with analysis results stored in:
        - adata.obs['leiden']: Cluster labels
        - adata.obs['cell_type']: Cell type annotations
        - adata.uns['markers']: Marker genes
        - adata.uns['characterization']: Cluster characterization

    Examples
    --------
    >>> # Use defaults with progress bar
    >>> adata = run_standard_analysis(adata, show_progress=True)

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

    Notes
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
        raise ValueError("recovery_save_dir is required when error_recovery=True and on_error='save'")

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

    # Define workflow steps
    workflow_steps = [
        ("clustering", config.clustering is not None),
        ("markers", getattr(config, "find_markers", True)),
        ("annotation", config.annotation is not None),
        ("characterization", getattr(config, "characterize", True)),
    ]

    # Filter to only enabled steps and skip completed ones
    steps_to_run = [
        name for name, enabled in workflow_steps
        if enabled and name not in completed_steps
    ]

    log.info(f"Steps to run: {steps_to_run}")

    # Initialize progress bar
    step_iterator = get_progress_bar(
        steps_to_run, desc="Analysis", enabled=show_progress, total=len(steps_to_run), unit="step"
    )

    # Track execution
    current_step = None
    successful_steps: List[str] = []

    # Determine cluster key from config
    from .config import ClusteringConfig
    cluster_config = config.clustering if isinstance(config.clustering, ClusteringConfig) else ClusteringConfig(**config.clustering)
    cluster_key = cluster_config.key_added or f"{cluster_config.method}_clusters"

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

                marker_config = DifferentialConfig(
                    groupby=cluster_key,
                    method=config.marker_method if hasattr(config, 'marker_method') else "wilcoxon"
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
                n_annotated = adata.obs['cell_type'].notna().sum() if 'cell_type' in adata.obs else 0
                log.info(f"  Annotated {n_annotated}/{len(adata)} cells")
                successful_steps.append(step_name)

            # Step 4: Characterization
            elif step_name == "characterization":
                log.info("Step: Cluster characterization")
                try:
                    adata = characterize_clusters(
                        adata,
                        groupby=cluster_key,
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
                log.warning(f"To resume, use: run_standard_analysis(adata, resume_from='{save_dir}')")
                return adata

        raise WorkflowError(
            f"[analysis] Workflow failed at step '{current_step}': {e}",
            step_name=current_step or "unknown",
            original_error=e,
        )

    # Store final config
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})[
        "workflow_config"
    ] = config.to_dict()
    adata.uns["sclucid"]["analysis"]["steps_executed"] = successful_steps

    log.info("=" * 60)
    log.info("=== Standard Analysis Pipeline Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)
    return adata


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

    Returns
    -------
    AnnData
        Annotated data with analysis results

    Examples
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

        if step == 'resolution':
            config = step_configs.get(step, {})
            eval_df, recommended_res = find_resolution(adata, **config)
            log.info(f"  Optimal resolution: {recommended_res}")

        elif step == 'clustering':
            from .config import ClusteringConfig
            config = step_configs.get(step, {})
            if isinstance(config, dict):
                config = ClusteringConfig(**config)
            adata = cluster_cells(adata, config)
            cluster_key = config.key_added or f"{config.method}_clusters"
            log.info(f"  Created {adata.obs[cluster_key].nunique()} clusters")

        elif step == 'markers':
            from .config import DifferentialConfig
            config = step_configs.get(step, {})
            if 'groupby' not in config:
                config['groupby'] = _default_groupby_key(adata)
            if isinstance(config, dict):
                config = DifferentialConfig(**config)
            markers = find_markers(adata, config)
            log.info(f"  Found markers for {len(markers)} clusters")

        elif step == 'annotation':
            config = step_configs.get(step, {})
            adata = run_annotation(adata, **config)
            log.info("  Annotation complete")

        elif step == 'scoring':
            config = step_configs.get(step, {})
            if 'gene_sets' not in config:
                log.warning("  No gene_sets provided, skipping")
                continue
            adata = score_by_gene_sets(adata, **config)
            log.info("  Scoring complete")

        elif step == 'characterization':
            config = step_configs.get(step, {})
            if 'groupby' not in config:
                config['groupby'] = _default_groupby_key(adata)
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

    Returns
    -------
    DataFrame
        Comparison of metrics across resolutions

    Examples
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
    from .config import ClusteringConfig
    from sklearn.metrics import silhouette_score

    if metrics is None:
        metrics = ['n_clusters', 'silhouette']

    results = []

    # Initialize progress bar
    res_iterator = get_progress_bar(
        resolutions, desc="Resolution Search", enabled=show_progress, total=len(resolutions), unit="res"
    )

    for res in res_iterator:
        log.info(f"Testing resolution: {res}")

        # Cluster at this resolution
        config = ClusteringConfig(resolution=res)
        adata_temp = cluster_cells(adata.copy(), config)

        # Compute metrics
        result = {'resolution': res}
        cluster_key = config.key_added or f"{config.method}_clusters"
        result['n_clusters'] = adata_temp.obs[cluster_key].nunique()

        if 'silhouette' in metrics and adata_temp.obsm.get('X_pca') is not None:
            try:
                score = silhouette_score(
                    adata_temp.obsm['X_pca'],
                    adata_temp.obs[cluster_key].astype(int)
                )
                result['silhouette'] = score
            except Exception as e:
                log.warning(f"  Could not compute silhouette: {e}")

        results.append(result)

    df_results = pd.DataFrame(results)

    if save_path:
        df_results.to_csv(save_path, index=False)
        log.info(f"Results saved to: {save_path}")

    return df_results
