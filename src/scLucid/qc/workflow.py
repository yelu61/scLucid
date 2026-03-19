"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, List, Any, Iterable, TypeVar, Literal

from anndata import AnnData
import numpy as np

from .config import QCWorkflowConfig
from .doublet import predict_doublets
from .filtering import (
    filter_cells,
    generate_qc_report,
    mark_low_quality_cell,
)
from .metrics import calculate_qc_metric
from ..utils import (
    get_progress_bar,
    save_workflow_result,
    WorkflowError,
    PartialResultManager,
    WorkflowCheckpoint,
)

log = logging.getLogger(__name__)

T = TypeVar('T')

# Keep for backward compatibility
QCWorkflowError = WorkflowError
PartialQCResult = PartialResultManager


def _progress_bar(
    iterable: Iterable[T],
    desc: str = "Processing",
    enabled: bool = True,
    total: Optional[int] = None,
) -> Iterable[T]:
    """
    Wrap an iterable with a tqdm progress bar if enabled.

    Uses utils.get_progress_bar for consistency with other modules.

    Args:
        iterable: The iterable to wrap.
        desc: Description for the progress bar.
        enabled: Whether to show the progress bar.
        total: Total number of items (optional).

    Returns:
        The original iterable or a tqdm-wrapped iterable.
    """
    return get_progress_bar(
        iterable, desc=desc, enabled=enabled, total=total, unit="sample"
    )


def _setup_workflow(
    adata_in: AnnData, save_dir: Optional[str], overwrite: bool
) -> Tuple[AnnData, Optional[Path]]:
    """
    Prepares the AnnData object and results directory for a workflow.

    Args:
        adata_in: Input AnnData object
        save_dir: Directory to save results. If None, no files will be saved.
        overwrite: If True, overwrite existing results directory

    Returns:
        Tuple of (adata_copy, results_path). results_path is None if save_dir is None.
    """
    adata = adata_in.copy()

    if save_dir is None:
        log.info("No save_dir specified. Running without file output.")
        return adata, None

    results_path = Path(save_dir)
    if results_path.exists() and not overwrite:
        log.warning(
            f"{save_dir} already exists. Old results may be overwritten. "
            "Consider setting overwrite=True in the config."
        )
    results_path.mkdir(parents=True, exist_ok=True)
    return adata, results_path


def _process_sample_qc(
    sample_adata: AnnData,
    config: QCWorkflowConfig,
    sample_name: str,
) -> AnnData:
    """
    Process QC for a single sample.

    This function is designed to be called in parallel for multiple samples.

    Args:
        sample_adata: AnnData object for a single sample
        config: QC workflow configuration
        sample_name: Name of the sample

    Returns:
        AnnData object with QC metrics computed
    """
    # Compute QC metrics. Fall back gracefully when cell cycle scoring is not feasible
    # (e.g., gene identifiers do not match marker lists).
    try:
        sample_adata = calculate_qc_metric(
            sample_adata,
            sample_key=None,  # Single sample, no sample key needed
            reporting_config=config.metrics_reporting_config,
            calculate_cell_cycle=True,
            cell_cycle_species=config.species,
        )
    except Exception as e:
        log.warning(
            f"Cell cycle scoring failed for sample '{sample_name}' ({e}). "
            "Retrying QC metrics without cell cycle scoring."
        )
        sample_adata = calculate_qc_metric(
            sample_adata,
            sample_key=None,
            reporting_config=config.metrics_reporting_config,
            calculate_cell_cycle=False,
        )

    return sample_adata


def _process_sample_doublet(
    sample_adata: AnnData,
    config: QCWorkflowConfig,
    save_dir: str,
    sample_name: str,
) -> AnnData:
    """
    Process doublet detection for a single sample.

    Args:
        sample_adata: AnnData object for a single sample
        config: QC workflow configuration
        save_dir: Directory to save doublet results
        sample_name: Name of the sample

    Returns:
        AnnData object with doublet predictions
    """
    # Update config save dir for this sample
    doublet_config = config.doublet_config
    doublet_config.save_dir = str(Path(save_dir) / sample_name)

    sample_adata = predict_doublets(
        sample_adata,
        config=doublet_config,
        sample_key=None,  # Single sample
    )

    return sample_adata


def _merge_sample_results(
    sample_results: List[Tuple[str, AnnData]],
    original_obs_names: List[str],
) -> AnnData:
    """
    Merge results from parallel sample processing.

    Args:
        sample_results: List of (sample_name, adata) tuples
        original_obs_names: Original observation names to preserve order

    Returns:
        Merged AnnData object
    """
    # Create a mapping from obs_name to row index
    obs_to_idx = {name: i for i, name in enumerate(original_obs_names)}

    # Initialize merged object
    merged_adata = sample_results[0][1].copy()
    merged_adata = merged_adata[original_obs_names].copy()

    # Merge QC metrics from each sample
    for sample_name, sample_adata in sample_results:
        # Get indices for this sample
        sample_idx = [i for i, name in enumerate(original_obs_names) if name in sample_adata.obs_names]

        # Merge each QC metric column
        for col in sample_adata.obs.columns:
            if col not in merged_adata.obs:
                # Initialize column
                merged_adata.obs[col] = np.nan

            # Copy values for this sample's cells
            for i, obs_name in enumerate(sample_adata.obs_names):
                if obs_name in obs_to_idx:
                    merged_adata.obs.loc[obs_name, col] = sample_adata.obs.loc[obs_name, col]

    return merged_adata


def _safe_parallel_process(
    process_func,
    samples: List[str],
    sample_data_func,
    config: QCWorkflowConfig,
    n_jobs: int,
    step_name: str = "processing",
    show_progress: bool = True,
) -> List[Tuple[str, Any]]:
    """
    Safely execute parallel processing with error handling and progress tracking.

    This wrapper ensures that:
    1. Individual sample failures don't crash the entire workflow
    2. Failed samples are logged and reported
    3. Successful samples are still merged
    4. Sequential fallback is available
    5. Progress bar shows processing status

    Args:
        process_func: Function to process each sample
        samples: List of sample names
        sample_data_func: Function to get sample data (adata[mask])
        config: QC workflow configuration
        n_jobs: Number of parallel jobs
        step_name: Name of the processing step for logging
        show_progress: Whether to show progress bar

    Returns:
        List of (sample_name, result) tuples. Failed samples have result=None.
    """
    from joblib import Parallel, delayed

    results = []
    failed_samples = []

    # Wrap samples with progress bar
    sample_iterator = _progress_bar(
        samples, desc=f"{step_name}", enabled=show_progress, total=len(samples)
    )

    # Process each sample with error handling
    for sample in sample_iterator:
        try:
            sample_data = sample_data_func(sample)
            result = process_func(sample_data, config, sample)
            results.append((sample, result))
        except Exception as e:
            log.error(f"Failed to process sample '{sample}' in {step_name}: {e}")
            failed_samples.append((sample, str(e)))
            results.append((sample, None))

    # Report summary if there were failures
    if failed_samples:
        log.warning("=" * 60)
        log.warning(f"PARALLEL PROCESSING SUMMARY - {step_name.upper()}")
        log.warning("=" * 60)
        log.warning(f"Total samples: {len(samples)}")
        log.warning(f"Successful: {len(samples) - len(failed_samples)}")
        log.warning(f"Failed: {len(failed_samples)}")
        for sample, error in failed_samples:
            log.warning(f"  - {sample}: {error}")
        log.warning("=" * 60)

    return results


def _run_qc_workflow(
    adata: AnnData,
    config: QCWorkflowConfig,
    results_path: Optional[Path],
    show_progress: bool = True,
) -> AnnData:
    """
    Run QC workflow.

    Args:
        adata: Input AnnData object
        config: QC workflow configuration
        results_path: Path to results directory. If None, no files will be saved.
        show_progress: Whether to show progress bars for multi-sample processing

    Returns:
        AnnData object with QC completed
    """
    # --- 1. QC metric calculation ---
    samples = adata.obs[config.sample_key].unique()

    if config.use_parallel and config.n_jobs != 1 and len(samples) > 1:
        log.info(f"Processing {len(samples)} samples in parallel with {config.n_jobs} jobs")

        results = _safe_parallel_process(
            process_func=_process_sample_qc,
            samples=list(samples),
            sample_data_func=lambda s: adata[adata.obs[config.sample_key] == s].copy(),
            config=config,
            n_jobs=config.n_jobs,
            step_name="QC metric calculation",
            show_progress=show_progress,
        )

        # Filter out failed samples (None results)
        successful_results = [(s, r) for s, r in results if r is not None]
        if not successful_results:
            raise RuntimeError("All samples failed QC metric calculation")

        if len(successful_results) < len(results):
            log.warning(f"Proceeding with {len(successful_results)}/{len(results)} successful samples")

        # Merge results
        adata = _merge_sample_results(
            successful_results,
            adata.obs_names.tolist(),
        )
    else:
        # Sequential processing
        if len(samples) == 1:
            adata = _process_sample_qc(adata, config, samples[0])
        else:
            adata = calculate_qc_metric(
                adata,
                sample_key=config.sample_key,
                reporting_config=config.metrics_reporting_config,
                calculate_cell_cycle=True,
                cell_cycle_species=config.species,
            )

    # --- 2. Doublet detection ---
    if results_path is not None:
        config.doublet_config.save_dir = str(results_path / "doublet")

    if config.use_parallel and config.n_jobs != 1 and len(samples) > 1:
        # Parallel doublet detection with error handling
        results = _safe_parallel_process(
            process_func=lambda data, cfg, name: _process_sample_doublet(
                data, cfg, config.doublet_config.save_dir, name
            ),
            samples=list(samples),
            sample_data_func=lambda s: adata[adata.obs[config.sample_key] == s].copy(),
            config=config,
            n_jobs=config.n_jobs,
            step_name="doublet detection",
            show_progress=show_progress,
        )

        # Filter out failed samples
        successful_results = [(s, r) for s, r in results if r is not None]
        if not successful_results:
            raise RuntimeError("All samples failed doublet detection")

        if len(successful_results) < len(results):
            log.warning(f"Proceeding with {len(successful_results)}/{len(results)} successful samples for doublet detection")

        # Merge doublet predictions
        adata = _merge_sample_results(
            successful_results,
            adata.obs_names.tolist(),
        )
    else:
        # Standard doublet detection
        adata = predict_doublets(
            adata,
            config=config.doublet_config,
            sample_key=config.sample_key,
        )

    # --- 3. Low-quality cell marking ---
    if results_path is not None:
        config.marking_config.save_dir = str(results_path / "low_quality")
    adata = mark_low_quality_cell(
        adata,
        config=config.marking_config,
        sample_key=config.sample_key,
    )

    return adata


def run_standard_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    overwrite: bool = False,
    *,
    show_progress: bool = True,
    # Error recovery
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: str = "raise",
    # Resume
    resume_from: Optional[str] = None,
) -> AnnData:
    """
    Run a standard single-cell RNA-seq QC workflow driven by a configuration object.

    If no config is provided, sensible defaults are used. This workflow includes:
    1. QC metric calculation (with cell cycle scoring).
    2. Doublet detection (Scrublet + basic heuristics).
    3. Low-quality cell marking (standard thresholds).
    4. Filtering of marked cells.
    5. Generation of a final report.

    New features in v0.3:
    - Error recovery with partial result saving
    - Resume from checkpoint
    - Consistent with preprocess/analysis workflows

    Args:
        adata_in: Input AnnData object (raw or pre-normalized).
        config: A QCWorkflowConfig object. If None, a default config is created.
        overwrite: If True, overwrite existing results directory.
        show_progress: If True, show progress bars for multi-sample processing.
        error_recovery: If True, enable error recovery mode.
        recovery_save_dir: Directory to save partial results on error.
        on_error: How to handle errors: "raise", "skip", or "save".
        resume_from: Path to checkpoint directory to resume from.

    Returns:
        Filtered AnnData object after QC.

    Examples:
        >>> # Basic usage with progress bar
        >>> adata_filtered = run_standard_qc(adata, show_progress=True)
        >>>
        >>> # With error recovery
        >>> adata_filtered = run_standard_qc(
        ...     adata,
        ...     error_recovery=True,
        ...     recovery_save_dir="./recovery",
        ...     on_error="save"
        ... )
        >>>
        >>> # Resume from checkpoint
        >>> adata_filtered = run_standard_qc(
        ...     adata,
        ...     resume_from="./recovery",
        ...     show_progress=True
        ... )
    """
    if config is None:
        log.info("No QCWorkflowConfig provided, using standard defaults.")
        config = QCWorkflowConfig()

    # Validate error recovery settings
    if error_recovery and on_error == "save" and not recovery_save_dir:
        raise ValueError("recovery_save_dir is required when error_recovery=True and on_error='save'")

    # Handle resume from checkpoint
    completed_steps: List[str] = []
    if resume_from:
        manager = PartialResultManager(resume_from)
        adata, checkpoint, _ = manager.load()
        completed_steps = checkpoint.completed_steps
        log.info(f"Resumed from checkpoint. Completed steps: {completed_steps}")
    else:
        adata = adata_in.copy()

    adata, results_path = _setup_workflow(adata, config.save_dir, overwrite)

    log.info("=" * 60)
    log.info("=== Starting Standard QC Workflow ===")
    log.info("=" * 60)
    log.info(f"Error recovery: {error_recovery}")
    log.info(f"Show progress: {show_progress}")

    # Track execution
    successful_steps: List[str] = []
    current_step = None

    try:
        # --- Step 1: QC Metrics ---
        if "qc_metrics" not in completed_steps:
            current_step = "qc_metrics"
            log.info("Step: QC Metrics Calculation")
            adata = _run_qc_workflow(adata, config, results_path, show_progress=show_progress)
            successful_steps.append("qc_metrics")
        else:
            log.info("Step: QC Metrics (already completed, skipping)")

        # --- Step 2: Filtering ---
        if "filtering" not in completed_steps:
            current_step = "filtering"
            log.info("Step: Cell Filtering")
            adata_filtered = filter_cells(adata, config=config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (already completed, using cached)")
            adata_filtered = adata

        # --- Step 3: Reporting ---
        if results_path is not None and "reporting" not in completed_steps:
            current_step = "reporting"
            log.info("Step: Report Generation")
            generate_qc_report(
                adata_filtered,
                save_dir=results_path / "report",
                sample_key=config.sample_key,
                adata_before=adata,
            )
            successful_steps.append("reporting")

    except Exception as e:
        error_msg = f"QC workflow failed at step '{current_step}': {str(e)}"
        log.error(error_msg)
        import traceback
        log.error(traceback.format_exc())

        if error_recovery and on_error in ["raise", "save"]:
            # Save partial results
            save_dir = recovery_save_dir or (str(results_path / "recovery") if results_path else "./recovery")
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(
                completed_steps=successful_steps,
                failed_step=current_step,
                error_message=str(e),
            )
            manager.save(adata, checkpoint, config)

            if on_error == "save":
                log.warning(f"QC failed but partial results saved to: {save_dir}")
                log.warning(f"To resume, use: run_standard_qc(adata, resume_from='{save_dir}')")
                return adata

        raise WorkflowError(error_msg, step_name=current_step or "unknown", original_error=e)

    # Save workflow result using standardized storage
    save_workflow_result(
        adata_filtered,
        module="qc",
        workflow_name="standard",
        steps=successful_steps,
        config=config.to_dict()
    )

    log.info("=" * 60)
    log.info("=== Standard QC Workflow Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)

    return adata_filtered


def run_advanced_qc(
    adata_in: AnnData,
    config: QCWorkflowConfig,
    overwrite: bool = False,
    *,
    show_progress: bool = True,
    # Error recovery
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: str = "raise",
    # Resume
    resume_from: Optional[str] = None,
) -> AnnData:
    """
    Run an advanced, fully configurable single-cell RNA-seq QC workflow.

    This workflow is entirely controlled by the provided QCWorkflowConfig object,
    allowing fine-grained control over every step.

    New features in v0.3:
    - Error recovery with partial result saving
    - Resume from checkpoint
    - Consistent with preprocess/analysis workflows

    Args:
        adata_in: Input AnnData object.
        config: A fully populated QCWorkflowConfig object.
        overwrite: If True, overwrite existing results directory.
        show_progress: If True, show progress bars for multi-sample processing.
        error_recovery: If True, enable error recovery mode.
        recovery_save_dir: Directory to save partial results on error.
        on_error: How to handle errors: "raise", "skip", or "save".
        resume_from: Path to checkpoint directory to resume from.

    Returns:
        Filtered AnnData object after QC.
    """
    # Validate error recovery settings
    if error_recovery and on_error == "save" and not recovery_save_dir:
        raise ValueError("recovery_save_dir is required when error_recovery=True and on_error='save'")

    # Handle resume from checkpoint
    completed_steps: List[str] = []
    if resume_from:
        manager = PartialResultManager(resume_from)
        adata, checkpoint, _ = manager.load()
        completed_steps = checkpoint.completed_steps
        log.info(f"Resumed from checkpoint. Completed steps: {completed_steps}")
    else:
        adata = adata_in.copy()

    adata, results_path = _setup_workflow(adata, config.save_dir, overwrite)

    log.info("=" * 60)
    log.info("=== Starting Advanced QC Workflow ===")
    log.info("=" * 60)
    if results_path is not None:
        log.info(f"Results will be saved in: {results_path}")
    else:
        log.info("Running without file output")
    log.info(f"Error recovery: {error_recovery}")

    # Track execution
    successful_steps: List[str] = []
    current_step = None

    try:
        # --- Step 1: QC Metrics ---
        if "qc_metrics" not in completed_steps:
            current_step = "qc_metrics"
            log.info("Step: QC Metrics Calculation")
            adata = _run_qc_workflow(adata, config, results_path, show_progress=show_progress)
            successful_steps.append("qc_metrics")
        else:
            log.info("Step: QC Metrics (already completed, skipping)")

        # --- Step 2: Filtering ---
        if "filtering" not in completed_steps:
            current_step = "filtering"
            log.info("Step: Cell Filtering")
            adata_filtered = filter_cells(adata, config=config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (already completed, using cached)")
            adata_filtered = adata

        # --- Step 3: Reporting ---
        if results_path is not None and "reporting" not in completed_steps:
            current_step = "reporting"
            log.info("Step: Report Generation")
            generate_qc_report(
                adata_filtered,
                save_dir=results_path / "report",
                sample_key=config.sample_key,
                adata_before=adata,
            )
            successful_steps.append("reporting")

    except Exception as e:
        error_msg = f"Advanced QC workflow failed at step '{current_step}': {str(e)}"
        log.error(error_msg)
        import traceback
        log.error(traceback.format_exc())

        if error_recovery and on_error in ["raise", "save"]:
            # Save partial results
            save_dir = recovery_save_dir or (str(results_path / "recovery") if results_path else "./recovery")
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(
                completed_steps=successful_steps,
                failed_step=current_step,
                error_message=str(e),
            )
            manager.save(adata, checkpoint, config)

            if on_error == "save":
                log.warning(f"QC failed but partial results saved to: {save_dir}")
                log.warning(f"To resume, use: run_advanced_qc(adata, config, resume_from='{save_dir}')")
                return adata

        raise WorkflowError(error_msg, step_name=current_step or "unknown", original_error=e)

    # Save workflow result using standardized storage
    save_workflow_result(
        adata_filtered,
        module="qc",
        workflow_name="advanced",
        steps=successful_steps,
        config=config.to_dict()
    )

    log.info("=" * 60)
    log.info("=== Advanced QC Workflow Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)

    return adata_filtered


__all__ = [
    "run_standard_qc",
    "run_advanced_qc",
    "QCWorkflowError",
    "PartialQCResult",
]
