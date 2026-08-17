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

import itertools
import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from anndata import AnnData

from ..base_config import apply_config_overrides
from ..utils import (
    PartialResultManager,
    StepResult,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    export_review_summary,
    get_marker_manager,
    get_progress_bar,
    normalize_review_summary,
    sanitize_for_hdf5,
    step_results_to_storage,
    validate_review_summary_schema,
)
from ..utils.context import AnalysisContext, infer_analysis_context
from .annotation import build_annotation_consensus, run_annotation, run_annotation_evidence
from .clustering import cluster_cells, run_clustering_review
from .config import AnalysisWorkflowConfig, AnnotationConfig, ProportionConfig, PseudobulkDEConfig
from .differential_expression import characterize_clusters, find_markers, run_pseudobulk_de
from .proportion import analyze_celltype_proportion
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
    "proportion",
    "pseudobulk_first",
]

_CUSTOM_STEP_ALIASES = {
    "resolution": "clustering_review",
}

# Keep for backward compatibility
AnalysisWorkflowError = WorkflowError
PartialAnalysisResult = PartialResultManager

__all__ = [
    "run_standard_analysis",
    "run_pseudobulk_first_analysis",
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
    else:
        # Use config flags to determine default steps. ``skip_steps`` filters
        # this configured baseline; it must not activate optional steps whose
        # feature flags are disabled.
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
        if config is not None and getattr(config, "run_proportion", False):
            resolved.append("proportion")
        if config is not None and getattr(config, "pseudobulk_first", False):
            resolved.append("pseudobulk_first")
        if skip_steps is not None:
            resolved = [step for step in resolved if step not in skip_steps]

    resolved = [_CUSTOM_STEP_ALIASES.get(step, step) for step in resolved]
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


def _sync_default_annotation_aliases(
    adata: AnnData,
    *,
    annotation_key: str = "cell_type_auto",
    lineage_key: str = "celltype_lineage_auto",
) -> None:
    """Expose canonical downstream annotation aliases without discarding detail."""
    if annotation_key in adata.obs.columns and "cell_type" not in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs[annotation_key]
    if lineage_key in adata.obs.columns and "celltype_lineage" not in adata.obs.columns:
        adata.obs["celltype_lineage"] = adata.obs[lineage_key]


def _resolve_cell_type_key(adata: AnnData, config: AnalysisWorkflowConfig) -> str:
    """Pick the cell-type column to use for pseudobulk grouping."""
    annotation_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )
    if annotation_key in adata.obs.columns:
        return annotation_key
    if "cell_type" in adata.obs.columns:
        return "cell_type"
    raise ValueError(
        "No cell-type column found for pseudobulk grouping. "
        "Run annotation first or provide a 'cell_type' column."
    )


def _build_context_aware_proportion_config(
    adata: AnnData,
    config: AnalysisWorkflowConfig,
    *,
    context: Optional[Union[AnalysisContext, Dict[str, Any]]] = None,
    **kwargs: Any,
) -> ProportionConfig:
    """Build default proportion settings from the resolved study-design contract.

    Explicit ``ProportionConfig`` objects are handled by the caller and never
    reach this helper.  This path exists only for ``run_proportion=True`` with
    no explicit proportion configuration, where falling back to conventional
    column names would otherwise discard the project context.
    """
    resolved_context = infer_analysis_context(adata, context=context)
    sample_col = (
        kwargs.get("sample_col")
        or resolved_context.sample_key
        or resolved_context.experimental_unit_key
    )
    condition_col = kwargs.get("condition_col") or resolved_context.condition_key
    celltype_col = kwargs.get("celltype_col") or resolved_context.cell_type_key
    if not celltype_col:
        try:
            celltype_col = _resolve_cell_type_key(adata, config)
        except ValueError:
            celltype_col = None

    missing = [
        name
        for name, value in (
            ("sample_key/sample_col", sample_col),
            ("condition_key/condition_col", condition_col),
            ("cell_type_key/celltype_col", celltype_col),
        )
        if not value
    ]
    if missing:
        raise ValueError(
            "Step 'proportion' could not resolve the required metadata fields: "
            f"{', '.join(missing)}. Provide them in ProjectContext, pass the "
            "corresponding *_col arguments, or set config.proportion explicitly."
        )

    pairing_col = kwargs.get("pairing_col") or resolved_context.paired_key
    batch_col = kwargs.get("batch_col") or resolved_context.batch_key
    return ProportionConfig(
        celltype_col=str(celltype_col),
        sample_col=str(sample_col),
        condition_col=str(condition_col),
        pairing_col=pairing_col,
        experimental_unit_col=(
            resolved_context.experimental_unit_key or pairing_col or str(sample_col)
        ),
        batch_col=batch_col,
    )


def _build_pseudobulk_first_review_summary(
    adata: AnnData,
    config: AnalysisWorkflowConfig,
    successful_steps: List[str],
    cluster_key: str,
    key_added: str,
    under_replicated: bool,
    warning_message: str,
) -> Dict[str, Any]:
    """Build a trace.py-compatible review summary for the pseudobulk-first workflow."""
    summary: Dict[str, Any] = {
        "module": "analysis",
        "workflow_name": key_added,
        "steps_executed": list(successful_steps),
        "cluster_key": cluster_key,
        "warnings": [warning_message] if under_replicated else [],
        "artifacts": {
            key_added: f'adata.uns["sclucid"]["analysis"]["{key_added}"]',
        },
    }
    enriched_summary = enrich_analysis_review_summary(
        summary,
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
    )
    config_dict = sanitize_for_hdf5(config.to_dict())
    readiness = enriched_summary.get("analysis_readiness", {})
    review_reasons = (
        readiness.get("review_reasons", [])
        if isinstance(readiness, dict)
        else []
    )
    review_summary = normalize_review_summary(
        enriched_summary,
        module="analysis",
        workflow_name=key_added,
        adata=adata,
        steps_executed=successful_steps,
        config=config_dict,
        warnings=review_reasons,
    )
    validate_review_summary_schema(review_summary, module="analysis", raise_on_error=True)
    validate_analysis_review_summary(review_summary, raise_on_error=True)
    return sanitize_for_hdf5(review_summary)


def _run_pseudobulk_first_de(
    adata: AnnData,
    condition_col: str,
    sample_col: str,
    cell_types: Optional[List[str]],
    config: AnalysisWorkflowConfig,
    *,
    contrasts: Optional[List[Tuple[str, str]]] = None,
    layer: Optional[str] = None,
    method: str = "auto",
    min_cells_per_sample: int = 5,
    min_samples_per_condition: int = 1,
    experimental_unit_col: Optional[str] = None,
    block_col: Optional[str] = None,
    design_covariates: Optional[List[str]] = None,
    key_added: str = "pseudobulk_first",
) -> Tuple[AnnData, Dict[str, Any]]:
    """Run sample-level pseudobulk DE per cell type after clustering/annotation."""
    if condition_col not in adata.obs.columns:
        raise KeyError(f"Condition column '{condition_col}' not found in adata.obs")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' not found in adata.obs")

    cell_type_key = _resolve_cell_type_key(adata, config)
    condition_values = sorted(pd.unique(adata.obs[condition_col].astype(str)).tolist())
    if len(condition_values) < 2:
        raise ValueError(
            f"Pseudobulk DE requires at least two conditions; found {condition_values}"
        )

    if contrasts is None:
        contrasts = list(itertools.combinations(condition_values, 2))

    group_names = cell_types
    if group_names is None:
        group_names = sorted(pd.unique(adata.obs[cell_type_key].astype(str)).tolist())
    if not group_names:
        raise ValueError("No cell types selected for pseudobulk-first DE.")

    pb_config = PseudobulkDEConfig(
        sample_col=sample_col,
        condition_key=condition_col,
        contrasts=contrasts,
        groupby=cell_type_key,
        group_names=group_names,
        layer=layer,
        use_raw=False,
        method=method,
        min_cells_per_sample=min_cells_per_sample,
        min_samples_per_condition=min_samples_per_condition,
        experimental_unit_col=experimental_unit_col,
        block_col=block_col,
        design_covariates=list(design_covariates or []),
        key_added=f"{key_added}_de",
        fallback_to_cell_level=False,
        single_sample_mode="descriptive",
    )
    results_df = run_pseudobulk_de(adata, config=pb_config)
    de_design = (
        adata.uns.get("sclucid", {})
        .get("analysis", {})
        .get("de", {})
        .get(f"{key_added}_de_design", {})
    )

    per_cell_type_results: Dict[str, pd.DataFrame] = {}
    if not results_df.empty and "group" in results_df.columns:
        for cell_type, sub_df in results_df.groupby(results_df["group"].astype(str)):
            per_cell_type_results[cell_type] = sub_df.reset_index(drop=True)

    # Build contrast records with per-cell-type independent-unit information,
    # including blocked groups that produced no gene-level rows.
    raw_contrast_records = de_design.get("contrasts", [])
    if isinstance(raw_contrast_records, dict) and all(
        str(key).isdigit() for key in raw_contrast_records
    ):
        raw_contrast_records = [
            raw_contrast_records[key]
            for key in sorted(raw_contrast_records, key=lambda key: int(str(key)))
        ]
    contrast_records: List[Dict[str, Any]] = [
        dict(record) for record in raw_contrast_records if isinstance(record, dict)
    ]
    under_replicated = False
    warning_message = ""
    for record in contrast_records:
        condition1 = str(record.get("condition1", ""))
        condition2 = str(record.get("condition2", ""))
        unit_counts = record.get("experimental_units_per_condition", {})
        n1 = int(unit_counts.get(condition1, 0)) if isinstance(unit_counts, dict) else 0
        n2 = int(unit_counts.get(condition2, 0)) if isinstance(unit_counts, dict) else 0
        if min(n1, n2) < 3:
            under_replicated = True
            if min(n1, n2) < 2:
                interpretation = (
                    "The result is descriptive only and cannot support formal inference."
                )
            else:
                interpretation = (
                    "The sample-level model is estimable but low-powered; effect uncertainty "
                    "and sensitivity require review."
                )
            warning_message = (
                f"Contrast '{condition2}_vs_{condition1}' in "
                f"'{record.get('group', 'all')}' has fewer than 3 independent units per "
                f"condition ({n1} vs {n2}). {interpretation}"
            )

    overall_valid = (
        bool(per_cell_type_results)
        and bool(contrast_records)
        and all(record.get("status") != "BLOCKED" for record in contrast_records)
        and all(
            bool(df["valid_for_publication_inference"].all())
            for df in per_cell_type_results.values()
        )
    )

    analysis_ns = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})
    analysis_ns[key_added] = {
        "contrasts": sanitize_for_hdf5(contrast_records),
        "per_cell_type_results": {
            ct: sanitize_for_hdf5(df) for ct, df in per_cell_type_results.items()
        },
        "inference_level": "sample_level",
        "valid_for_publication_inference": overall_valid,
        "design": sanitize_for_hdf5(de_design),
        "params": sanitize_for_hdf5(pb_config.to_dict()),
    }

    meta = {
        "under_replicated": under_replicated,
        "warning_message": warning_message,
        "valid_for_publication_inference": overall_valid,
    }
    return adata, meta


def run_pseudobulk_first_analysis(
    adata: AnnData,
    condition_col: str,
    sample_col: str,
    cell_types: Optional[List[str]] = None,
    *,
    config: Optional[AnalysisWorkflowConfig] = None,
    contrasts: Optional[List[Tuple[str, str]]] = None,
    layer: Optional[str] = None,
    method: str = "auto",
    min_cells_per_sample: int = 5,
    min_samples_per_condition: int = 1,
    experimental_unit_col: Optional[str] = None,
    block_col: Optional[str] = None,
    design_covariates: Optional[List[str]] = None,
    key_added: str = "pseudobulk_first",
    show_progress: bool = True,
    **kwargs,
) -> AnnData:
    """
    Run a pseudobulk-first analysis workflow.

    This workflow first runs standard clustering and annotation, then immediately
    aggregates cells to sample-level pseudobulk per cell type for condition DE.
    Cell-level condition DE is avoided as the primary inference output.

    Parameters
    ----------
    adata : AnnData
        Preprocessed single-cell data.
    condition_col : str
        Column in ``adata.obs`` with condition labels.
    sample_col : str
        Column in ``adata.obs`` with biological sample identifiers.
    cell_types : list of str, optional
        Subset of cell types to test. If None, all observed cell types are used.
    config : AnalysisWorkflowConfig, optional
        Workflow configuration. If None, defaults are used.
    contrasts : list of tuple[str, str], optional
        Condition contrasts to test. If None, all pairwise contrasts are used.
    layer : str, optional
        Layer containing count values for pseudobulk aggregation.
    method : str, default="auto"
        Pseudobulk DE method passed to ``run_pseudobulk_de``.
    min_cells_per_sample : int, default=5
        Minimum cells per sample to include the sample.
    min_samples_per_condition : int, default=1
        Minimum samples per condition required to test a contrast.
    experimental_unit_col : str, optional
        Independent biological unit used to count replicates.
    block_col : str, optional
        Paired/repeated-measures unit included in the sample-level design.
    design_covariates : list of str, optional
        Explicit sample-level adjustment covariates.
    key_added : str, default="pseudobulk_first"
        Namespace for stored results.
    show_progress : bool, default=True
        Show progress bars.
    **kwargs
        Additional overrides for ``AnalysisWorkflowConfig``.

    Returns:
    -------
    AnnData
        Annotated data with results in
        ``adata.uns["sclucid"]["analysis"][key_added]``.
    """
    if config is None:
        from .config import AnalysisWorkflowConfig as DefaultConfig

        config = apply_config_overrides(DefaultConfig(), **kwargs)
    else:
        config = apply_config_overrides(config, **kwargs)

    # Prevent recursion when the function is called from run_standard_analysis.
    if getattr(config, "pseudobulk_first", False):
        config = config.model_copy(update={"pseudobulk_first": False})

    adata = run_standard_analysis(adata, config=config, show_progress=show_progress)

    adata, meta = _run_pseudobulk_first_de(
        adata,
        condition_col=condition_col,
        sample_col=sample_col,
        cell_types=cell_types,
        config=config,
        contrasts=contrasts,
        layer=layer,
        method=method,
        min_cells_per_sample=min_cells_per_sample,
        min_samples_per_condition=min_samples_per_condition,
        experimental_unit_col=experimental_unit_col,
        block_col=block_col,
        design_covariates=design_covariates,
        key_added=key_added,
    )

    if meta["under_replicated"]:
        warnings.warn(meta["warning_message"], UserWarning, stacklevel=2)

    cluster_key = (
        config.clustering.key_added
        if config.clustering and config.clustering.key_added in adata.obs.columns
        else _default_groupby_key(adata)
    )
    successful_steps = list(
        adata.uns.get("sclucid", {}).get("analysis", {}).get("steps_executed", [])
    )
    if f"{key_added}_analysis" not in successful_steps:
        successful_steps.append(f"{key_added}_analysis")

    review_summary = _build_pseudobulk_first_review_summary(
        adata,
        config=config,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
        key_added=key_added,
        under_replicated=meta["under_replicated"],
        warning_message=meta["warning_message"],
    )
    adata.uns["sclucid"]["analysis"][f"{key_added}_review_summary"] = review_summary
    return adata


def run_standard_analysis(
    adata: AnnData,
    config: Optional[AnalysisWorkflowConfig] = None,
    *,
    context: Optional[Union[AnalysisContext, Dict[str, Any]]] = None,
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
    context : AnalysisContext or dict, optional
        Resolved project metadata used for sample-aware analysis defaults.
        Explicit analysis and proportion configuration remains authoritative.
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
    step_results: List[StepResult] = []
    markers_df: Optional[pd.DataFrame] = None
    annotation_review_table: Optional[pd.DataFrame] = None

    # Determine cluster key from config or existing analysis state. Downstream
    # marker/characterization steps require this key if clustering is skipped.
    from .config import ClusteringConfig

    if isinstance(config.clustering, ClusteringConfig):
        cluster_config = config.clustering
    elif isinstance(config.clustering, dict):
        cluster_config = ClusteringConfig(**config.clustering)
        config.clustering = cluster_config
    else:
        cluster_config = ClusteringConfig()
        config.clustering = cluster_config
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
                    recommended is not None
                    and cluster_config.method in {"leiden", "louvain"}
                ):
                    if getattr(config, "use_recommended_resolution", False):
                        cluster_config = cluster_config.model_copy(
                            update={"resolution": float(recommended)}
                        )
                        config.clustering = cluster_config
                        log.info(
                            "  Applying recommended resolution: "
                            f"{cluster_config.resolution:g}"
                        )
                    else:
                        log.info(
                            f"Review recommends resolution {recommended:g}; "
                            "set use_recommended_resolution=True to apply automatically."
                        )
                log.info(f"  Reviewed {len(review_df)} clustering resolution candidate(s)")
                successful_steps.append(step_name)
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="heuristic",
                        outputs={"n_candidates": int(len(review_df))},
                    )
                )

            # Step 1: Clustering
            elif step_name == "clustering":
                log.info("Step: Clustering")
                adata = cluster_cells(adata, cluster_config)
                log.info(f"  Clustering complete: {adata.obs[cluster_key].nunique()} clusters")
                successful_steps.append(step_name)
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="validated_core",
                        outputs={
                            "n_clusters": int(adata.obs[cluster_key].nunique()),
                            "cluster_key": cluster_key,
                        },
                    )
                )

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
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="validated_core",
                        outputs={"n_marker_rows": int(len(markers_df))},
                    )
                )

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
                annotation_config = (
                    config.annotation
                    if isinstance(config.annotation, AnnotationConfig)
                    else AnnotationConfig(**config.annotation)
                    if isinstance(config.annotation, dict)
                    else AnnotationConfig()
                )
                _sync_default_annotation_aliases(
                    adata,
                    annotation_key=annotation_config.key_added,
                    lineage_key=annotation_config.lineage_key,
                )
                n_annotated = (
                    adata.obs["cell_type"].notna().sum() if "cell_type" in adata.obs else n_annotated
                )
                log.info(f"  Annotated {n_annotated}/{len(adata)} cells")
                successful_steps.append(step_name)
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="heuristic",
                        outputs={
                            "n_annotated": int(n_annotated),
                            "n_total": int(len(adata)),
                        },
                    )
                )

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
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="heuristic",
                        outputs={"review_table_rows": int(annotation_review_table.shape[0])},
                    )
                )

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
                _sync_default_annotation_aliases(
                    adata,
                    annotation_key=annotation_config.key_added,
                    lineage_key=annotation_config.lineage_key,
                )
                log.info(f"  Applied consensus labels to obs['{annotation_config.key_added}']")
                successful_steps.append(step_name)
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="heuristic",
                        outputs={"final_key": annotation_config.key_added},
                    )
                )

            # Step 6: Malignancy interpretation (deprecated in analysis workflow;
            # prefer ``run_tumor_analysis`` or ``post_analysis_hooks``).
            elif step_name == "malignancy_interpretation":
                log.info("Step: Malignancy interpretation")
                warnings.warn(
                    "Running malignancy_interpretation inside run_standard_analysis is "
                    "deprecated. Use run_tumor_analysis() or config.post_analysis_hooks instead.",
                    FutureWarning,
                    stacklevel=2,
                )
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
                from ..tumor.malignancy import run_malignancy_interpretation

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
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="heuristic",
                        outputs={"review_table_rows": int(malignancy_table.shape[0])},
                        warnings=["deprecated: run via tumor workflow instead"],
                    )
                )

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
                    step_results.append(
                        StepResult(
                            name=step_name,
                            status="completed",
                            evidence_level="validated_core",
                            outputs={"groupby": active_cluster_key},
                        )
                    )
                except Exception as e:
                    if on_error == "skip":
                        log.warning(f"  Characterization failed: {e}. Skipping...")
                        step_results.append(
                            StepResult.from_exception(
                                name=step_name,
                                exc=e,
                                degraded=True,
                                evidence_level="unavailable",
                            )
                        )
                    else:
                        raise

            # Step 8: Cell type proportion analysis
            elif step_name == "proportion":
                log.info("Step: Cell type proportion analysis")
                proportion_config = config.proportion
                if isinstance(proportion_config, dict):
                    proportion_config = ProportionConfig(**proportion_config)
                elif proportion_config is None:
                    proportion_config = _build_context_aware_proportion_config(
                        adata,
                        config,
                        context=context,
                        **kwargs,
                    )
                config.proportion = proportion_config
                result = analyze_celltype_proportion(
                    adata,
                    method=getattr(config, "proportion_method", None),
                    config=proportion_config,
                    sample_col=proportion_config.sample_col,
                    condition_col=proportion_config.condition_col,
                    celltype_col=proportion_config.celltype_col,
                    out_dir=proportion_config.out_dir,
                    return_type="anndata",
                )
                if isinstance(result, AnnData):
                    adata = result
                prop_ns = adata.uns.get("sclucid", {}).get("proportion", {})
                stat_df = prop_ns.get("stat_df")
                successful_steps.append(step_name)
                step_results.append(
                    StepResult(
                        name=step_name,
                        status="completed",
                        evidence_level="validated_core"
                        if proportion_config.require_biological_replicates
                        else "exploratory",
                        outputs={
                            "method": prop_ns.get("method", getattr(config, "proportion_method", "auto")),
                            "celltype_col": proportion_config.celltype_col,
                            "sample_col": proportion_config.sample_col,
                            "condition_col": proportion_config.condition_col,
                            "n_result_rows": int(stat_df.shape[0])
                            if hasattr(stat_df, "shape")
                            else None,
                        },
                    )
                )

            # Step 9: Sample-level pseudobulk DE per cell type
            elif step_name == "pseudobulk_first":
                log.info("Step: Pseudobulk-first differential expression")
                resolved_context = infer_analysis_context(adata, context=context)
                pb_condition_col = kwargs.get("condition_col") or resolved_context.condition_key
                pb_sample_col = (
                    kwargs.get("sample_col")
                    or resolved_context.sample_key
                    or resolved_context.experimental_unit_key
                )
                pb_experimental_unit_col = (
                    kwargs.get("experimental_unit_col")
                    or resolved_context.experimental_unit_key
                    or pb_sample_col
                )
                pb_block_col = kwargs.get("block_col") or resolved_context.paired_key
                if "design_covariates" in kwargs:
                    pb_design_covariates = list(kwargs.get("design_covariates") or [])
                else:
                    # Batch metadata is a design candidate, not an automatic
                    # adjustment. Applying it without a rank/confounding review
                    # can remove the condition effect the user intends to test.
                    pb_design_covariates = []
                if not pb_condition_col or not pb_sample_col:
                    raise ValueError(
                        "Step 'pseudobulk_first' requires condition_col and sample_col arguments."
                    )
                adata, pb_meta = _run_pseudobulk_first_de(
                    adata,
                    condition_col=pb_condition_col,
                    sample_col=pb_sample_col,
                    cell_types=kwargs.get("cell_types"),
                    config=config,
                    contrasts=kwargs.get("contrasts"),
                    layer=kwargs.get("layer"),
                    method=kwargs.get("method", "auto"),
                    min_cells_per_sample=kwargs.get("min_cells_per_sample", 5),
                    min_samples_per_condition=kwargs.get("min_samples_per_condition", 1),
                    experimental_unit_col=pb_experimental_unit_col,
                    block_col=pb_block_col,
                    design_covariates=pb_design_covariates,
                    key_added=kwargs.get("key_added", "pseudobulk_first"),
                )
                successful_steps.append("pseudobulk_first")
                step_results.append(
                    StepResult(
                        name="pseudobulk_first",
                        status="completed",
                        evidence_level="validated_core",
                        outputs={
                            "n_cell_types_tested": len(
                                adata.uns.get("sclucid", {})
                                .get("analysis", {})
                                .get("pseudobulk_first", {})
                                .get("per_cell_type_results", {})
                            ),
                            "valid_for_publication_inference": pb_meta[
                                "valid_for_publication_inference"
                            ],
                        },
                        warnings=[pb_meta["warning_message"]]
                        if pb_meta["under_replicated"]
                        else [],
                    )
                )
                if pb_meta["under_replicated"]:
                    warnings.warn(pb_meta["warning_message"], UserWarning, stacklevel=2)

    except Exception as e:
        error_msg = f"Workflow failed at step '{current_step}': {str(e)}"
        log.error(error_msg)
        import traceback

        log.error(traceback.format_exc())

        if current_step is not None:
            step_results.append(
                StepResult.from_exception(
                    name=current_step,
                    exc=e,
                    degraded=False,
                    evidence_level="unavailable",
                )
            )

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
    config_dict = sanitize_for_hdf5(config.to_dict())
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})[
        UnsKeys.WORKFLOW_CONFIG
    ] = config_dict
    adata.uns["sclucid"]["analysis"][UnsKeys.STEPS_EXECUTED] = successful_steps
    adata.uns["sclucid"]["analysis"]["step_results"] = step_results_to_storage(step_results)

    # Run optional post-analysis hooks (e.g., tumor interpretation without
    # hard-coding a tumor import inside the analysis module).
    hooks = getattr(config, "post_analysis_hooks", None) or []
    for hook in hooks:
        try:
            if callable(hook):
                adata = hook(adata, config)
        except Exception as hook_exc:
            log.warning(f"Post-analysis hook failed: {hook_exc}")
            step_results.append(
                StepResult.from_exception(
                    name="post_analysis_hook",
                    exc=hook_exc,
                    degraded=True,
                    evidence_level="unavailable",
                )
            )
    # Re-save in case hooks added more results
    adata.uns["sclucid"]["analysis"]["step_results"] = step_results_to_storage(step_results)

    # Build and store review summary
    enriched_summary = enrich_analysis_review_summary(
        _build_analysis_review_summary(
            adata, config, successful_steps, cluster_key, step_results=step_results
        ),
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
        config=config_dict,
        warnings=(
            enriched_summary.get("analysis_readiness", {}).get("review_reasons", [])
            if isinstance(enriched_summary.get("analysis_readiness"), dict)
            else []
        ),
    )
    validate_review_summary_schema(review_summary, module="analysis", raise_on_error=True)
    validate_analysis_review_summary(review_summary, raise_on_error=True)
    review_summary = sanitize_for_hdf5(review_summary)
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
    step_results: Optional[List[StepResult]] = None,
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
    if step_results:
        from ..utils import summarize_step_results

        summary["step_results"] = summarize_step_results(step_results)
        summary["deprecated_steps_used"] = [
            r.name for r in step_results if "deprecated" in " ".join(r.warnings)
        ]

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

    # Pseudobulk-first summary
    if "pseudobulk_first" in successful_steps:
        pb_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("pseudobulk_first", {})
        summary["pseudobulk_first"] = {
            "n_cell_types_tested": len(pb_ns.get("per_cell_type_results", {})),
            "n_contrasts": len(pb_ns.get("contrasts", [])),
            "inference_level": pb_ns.get("inference_level"),
            "decision": "sample_level_de_primary_output",
            "valid_for_publication_inference": pb_ns.get("valid_for_publication_inference", False),
        }
        summary["artifacts"][
            "pseudobulk_first"
        ] = 'adata.uns["sclucid"]["analysis"]["pseudobulk_first"]'
        if not pb_ns.get("valid_for_publication_inference", False):
            summary["warnings"].append("pseudobulk_first_results_not_valid_for_publication_inference")

    # Cell type proportion summary
    if "proportion" in successful_steps:
        prop_ns = adata.uns.get("sclucid", {}).get("proportion", {})
        prop_config = config.proportion
        if isinstance(prop_config, dict):
            prop_config = ProportionConfig(**prop_config)
        stat_df = prop_ns.get("stat_df")
        summary["proportion"] = {
            "method": prop_ns.get("method", getattr(config, "proportion_method", None) or "auto"),
            "decision": "sample_level_celltype_composition_review",
            "celltype_col": getattr(prop_config, "celltype_col", None),
            "sample_col": getattr(prop_config, "sample_col", None),
            "condition_col": getattr(prop_config, "condition_col", None),
            "n_result_rows": int(stat_df.shape[0]) if hasattr(stat_df, "shape") else None,
        }
        summary["artifacts"][
            "proportion"
        ] = 'adata.uns["sclucid"]["proportion"]'

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

    steps = [_CUSTOM_STEP_ALIASES.get(step, step) for step in steps]
    invalid = set(steps) - set(ANALYSIS_WORKFLOW_STEPS) - {"scoring"}
    if invalid:
        raise ValueError(
            f"Invalid step names: {sorted(invalid)}. Valid steps are: "
            f"{ANALYSIS_WORKFLOW_STEPS + ['scoring']}"
        )

    log.info(f"Running custom analysis with {len(steps)} steps...")

    # Initialize progress bar
    step_iterator = get_progress_bar(
        steps, desc="Custom Analysis", enabled=show_progress, total=len(steps), unit="step"
    )

    for i, step in enumerate(step_iterator, 1):
        log.info(f"Step {i}/{len(steps)}: {step}")

        if step == "clustering_review":
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

        elif step == "annotation_evidence":
            config = step_configs.get(step, {})
            groupby = config.pop("groupby", _default_groupby_key(adata))
            run_annotation_evidence(adata, groupby, **config)
            log.info("  Annotation evidence complete")

        elif step == "annotation_consensus":
            config = step_configs.get(step, {})
            groupby = config.pop("groupby", _default_groupby_key(adata))
            review_table = config.pop(
                "annotation_review_table",
                adata.uns.get("sclucid", {})
                .get("analysis", {})
                .get("annotation", {})
                .get("annotation_review_table"),
            )
            build_annotation_consensus(adata, groupby, review_table, **config)
            log.info("  Annotation consensus complete")

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

        elif step == "proportion":
            config = step_configs.get(step, {})
            adata_result = analyze_celltype_proportion(adata, return_type="anndata", **config)
            if isinstance(adata_result, AnnData):
                adata = adata_result
            log.info("  Proportion analysis complete")

        elif step == "pseudobulk_first":
            config = step_configs.get(step, {})
            required = {"condition_col", "sample_col"}
            missing = sorted(required - set(config))
            if missing:
                raise ValueError(
                    "Custom step 'pseudobulk_first' requires config values: "
                    f"{missing}"
                )
            adata, _ = _run_pseudobulk_first_de(adata, **config)
            log.info("  Pseudobulk-first DE complete")

        elif step == "malignancy_interpretation":
            raise ValueError(
                "Custom step 'malignancy_interpretation' is deprecated in analysis. "
                "Use scLucid.tumor.run_malignancy_interpretation or post_analysis_hooks."
            )

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
                labels = adata_temp.obs[cluster_key]
                if isinstance(labels.dtype, pd.CategoricalDtype):
                    labels = labels.cat.codes
                score = silhouette_score(
                    adata_temp.obsm["X_pca"], labels.astype(int)
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
