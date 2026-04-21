"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, List, Any, Iterable, TypeVar, Literal

from anndata import AnnData
import numpy as np

from .config import QCWorkflowConfig
from .doublet import predict_doublets
from .filtering import (
    AdaptiveThresholdCalculator,
    filter_cells,
    generate_qc_report,
    mark_low_quality_cell,
)
from .metrics import calculate_qc_metric
from ..utils import (
    get_progress_bar,
    save_result,
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


def _prepare_runtime_qc_config(
    config: Optional[QCWorkflowConfig],
    tissue_type: str,
) -> QCWorkflowConfig:
    """
    Build a runtime config copy for workflow execution.

    Workflow entrypoints should never mutate the caller's config object. This helper
    centralizes the deep-copy behavior and fills in default tissue context only on
    the runtime copy.
    """
    runtime_config = config.model_copy(deep=True) if config is not None else QCWorkflowConfig()
    if runtime_config.tissue_type is None:
        runtime_config.tissue_type = tissue_type
    return runtime_config


def _ensure_sample_key(
    adata: AnnData,
    config: QCWorkflowConfig,
    warnings_list: Optional[List[str]] = None,
) -> None:
    """
    Ensure the configured sample key exists in adata.obs.

    Real single-sample datasets often arrive without a batch/sample column. In that
    case we create a synthetic single-sample label so the rest of the QC stack can
    continue to use a uniform multi-sample-aware code path.
    """
    if config.sample_key in adata.obs.columns:
        return

    candidate_keys = [
        "sampleID",
        "sample",
        "Sample",
        "orig.ident",
        "orig_ident",
        "patient",
        "patient_id",
        "donor",
        "donor_id",
        "batch",
        "Batch",
    ]
    for candidate in candidate_keys:
        if candidate in adata.obs.columns:
            original_key = config.sample_key
            config.sample_key = candidate
            msg = (
                f"Sample key '{original_key}' not found; using detected obs column "
                f"'{candidate}' for sample-aware QC."
            )
            log.info(msg)
            if warnings_list is not None:
                warnings_list.append(msg)
            return

    synthetic_sample = "sample_1"
    adata.obs[config.sample_key] = synthetic_sample
    msg = (
        f"Sample key '{config.sample_key}' not found in adata.obs; "
        f"created synthetic single-sample labels ('{synthetic_sample}')."
    )
    log.info(msg)
    if warnings_list is not None:
        warnings_list.append(msg)


def _add_tumor_aware_flags(
    adata: AnnData,
    tissue_type: str,
) -> None:
    """
    Store tumor-aware QC flags when tissue_type indicates tumor/cancer.

    Tumor tissues often have elevated mitochondrial content and other
    characteristics that should be flagged rather than aggressively filtered.
    """
    if not tissue_type or "tumor" not in tissue_type.lower() and "cancer" not in tissue_type.lower():
        return

    flags: Dict[str, Any] = {"tissue_type": tissue_type, "tumor_aware_enabled": True}

    if "pct_counts_mt" in adata.obs.columns:
        mt_values = adata.obs["pct_counts_mt"].values
        high_mt_frac = float(np.mean(mt_values > 10.0))
        flags["high_mt_population_flagged"] = high_mt_frac > 0.25
        flags["mean_pct_counts_mt"] = float(np.mean(mt_values))
        flags["fraction_mt_above_10pct"] = high_mt_frac

    if "pct_counts_ribo" in adata.obs.columns:
        ribo_values = adata.obs["pct_counts_ribo"].values
        flags["mean_pct_counts_ribo"] = float(np.mean(ribo_values))

    flags["note"] = (
        "Tumor-aware QC active: elevated mitochondrial content is flagged "
        "rather than hard-filtered. Review thresholds manually."
    )
    save_result(adata, "qc", "tumor_aware_flags", flags)
    log.info(f"Tumor-aware QC flags stored: {list(flags.keys())}")


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
    save_dir: Optional[str],
    sample_name: str,
) -> AnnData:
    """
    Process doublet detection for a single sample.

    Args:
        sample_adata: AnnData object for a single sample
        config: QC workflow configuration
        save_dir: Directory to save doublet results. If None, no files are saved.
        sample_name: Name of the sample

    Returns:
        AnnData object with doublet predictions
    """
    # Update config save dir for this sample
    doublet_config = config.doublet_config
    if save_dir is not None:
        doublet_config.save_dir = str(Path(save_dir) / sample_name)
    else:
        doublet_config.save_dir = None

    sample_adata = predict_doublets(
        sample_adata,
        config=doublet_config,
        sample_key=config.sample_key,
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
    import anndata as ad

    sample_adatas = [sample_adata for _, sample_adata in sample_results]
    merged_adata = ad.concat(sample_adatas, merge="same")
    merged_adata = merged_adata[original_obs_names].copy()
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


def _is_tumor_aware(tissue_type: Optional[str]) -> bool:
    if not tissue_type:
        return False
    return "tumor" in tissue_type.lower() or "cancer" in tissue_type.lower()


def _apply_qc_recommendations(
    config: QCWorkflowConfig,
    recommendation: Any,
) -> Tuple[QCWorkflowConfig, QCWorkflowConfig]:
    """Apply intelligent QC recommendation fields to a deep copy of the config.

    Returns:
        Tuple of (applied_config, original_config_snapshot).
        Only fields that were not explicitly set by the user are filled from recommendation.
    """
    # Deep copy to avoid mutating the caller's original config
    original = config.model_copy(deep=True)
    if recommendation is None:
        return original, original

    applied = config.model_copy(deep=True)
    rec_dict = recommendation.to_dict() if hasattr(recommendation, "to_dict") else {}

    def _is_user_set(obj, field_name: str) -> bool:
        return field_name in getattr(obj, "model_fields_set", set())

    # min_genes -> marking_config.thresholds.min_genes
    min_genes_rec = rec_dict.get("min_genes")
    if isinstance(min_genes_rec, dict) and min_genes_rec.get("threshold") is not None:
        if not (
            _is_user_set(applied, "marking_config")
            and _is_user_set(applied.marking_config, "thresholds")
            and _is_user_set(applied.marking_config.thresholds, "min_genes")
        ):
            applied.marking_config.thresholds.min_genes = int(min_genes_rec["threshold"])

    # max_mt_percent -> marking_config.thresholds.pc_mt
    mt_rec = rec_dict.get("max_mt_percent")
    if isinstance(mt_rec, dict) and mt_rec.get("threshold") is not None:
        if not (
            _is_user_set(applied, "marking_config")
            and _is_user_set(applied.marking_config, "thresholds")
            and _is_user_set(applied.marking_config.thresholds, "pc_mt")
        ):
            applied.marking_config.thresholds.pc_mt = float(mt_rec["threshold"])

    # n_counts -> marking_config.thresholds.min_counts
    counts_rec = rec_dict.get("n_counts")
    if isinstance(counts_rec, dict) and counts_rec.get("threshold") is not None:
        if not (
            _is_user_set(applied, "marking_config")
            and _is_user_set(applied.marking_config, "thresholds")
            and _is_user_set(applied.marking_config.thresholds, "min_counts")
        ):
            applied.marking_config.thresholds.min_counts = int(counts_rec["threshold"])

    # doublet_threshold -> doublet_config.score_threshold
    doublet_rec = rec_dict.get("doublet_threshold")
    if (
        isinstance(doublet_rec, dict)
        and doublet_rec.get("threshold") is not None
        and doublet_rec.get("confidence", 0) > 0
        and doublet_rec.get("method") != "no_doublet_scores"
    ):
        if not (
            _is_user_set(applied, "doublet_config")
            and _is_user_set(applied.doublet_config, "score_threshold")
        ):
            applied.doublet_config.score_threshold = float(doublet_rec["threshold"])

    return applied, original


def _compute_sample_thresholds(
    adata: AnnData,
    config: QCWorkflowConfig,
) -> Tuple[Dict[str, Any], List[str]]:
    """Compute per-sample adaptive thresholds when hierarchical/independent mode is active."""
    warnings: List[str] = []
    sample_thresholds: Dict[str, Any] = {}
    try:
        calculator = AdaptiveThresholdCalculator(adata, config.sample_key)
        metrics_to_compute = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        for metric in metrics_to_compute:
            try:
                thresholds = calculator._suggest_adaptive_thresholds(
                    metric, method=config.threshold_mode, percentile=95.0
                )
                for sample, th in thresholds.items():
                    sample_thresholds.setdefault(sample, {})[metric] = th
            except Exception as e:
                warnings.append(f"Adaptive threshold calculation failed for {metric}: {e}")
    except Exception as e:
        warnings.append(f"AdaptiveThresholdCalculator initialization failed: {e}")
    return sample_thresholds, warnings


def _diff_qc_recommendations(
    recommendation: Any,
    original_config: QCWorkflowConfig,
) -> Dict[str, Any]:
    """Compare recommended values against the original user config.

    This captures genuine user-vs-recommendation divergence.
    """
    if recommendation is None:
        return {}
    diffs: Dict[str, Any] = {}
    rec_dict = recommendation.to_dict() if hasattr(recommendation, "to_dict") else {}
    cfg_dict = original_config.to_dict()

    mapping = {
        "min_genes": ("marking_config", "thresholds", "min_genes"),
        "max_mt_percent": ("marking_config", "thresholds", "pc_mt"),
        "n_counts": ("marking_config", "thresholds", "min_counts"),
        "doublet_threshold": ("doublet_config", "score_threshold"),
    }

    for param_name, path in mapping.items():
        rec_val = None
        param_rec = rec_dict.get(param_name)
        if isinstance(param_rec, dict):
            rec_val = param_rec.get("threshold")

        actual_val = cfg_dict
        for key in path:
            if isinstance(actual_val, dict):
                actual_val = actual_val.get(key)
            else:
                actual_val = None
                break

        if rec_val is not None and actual_val is not None and rec_val != actual_val:
            diffs[param_name] = {"recommended": rec_val, "actual": actual_val}

    return diffs


def _store_qc_trace(
    adata: AnnData,
    config: QCWorkflowConfig,
    original_config: QCWorkflowConfig,
    recommendation: Any,
    sample_thresholds: Dict[str, Any],
    filtering_summary: Dict[str, Any],
    warnings: List[str],
) -> None:
    """Store unified QC trace under adata.uns['sclucid']['qc']."""
    n_samples = int(adata.obs[config.sample_key].nunique()) if config.sample_key in adata.obs else 1
    save_result(adata, "qc", "context", {
        "sample_key": config.sample_key,
        "threshold_mode": config.threshold_mode,
        "n_samples": n_samples,
        "tissue_type": config.tissue_type,
        "use_recommendations": config.use_recommendations,
    })
    if recommendation is not None:
        save_result(adata, "qc", "recommendation", recommendation.to_dict())
    save_result(adata, "qc", "original_config", original_config.to_dict())
    save_result(adata, "qc", "applied_config", config.to_dict())
    save_result(adata, "qc", "user_overrides", _diff_qc_recommendations(recommendation, original_config))
    save_result(adata, "qc", "sample_thresholds", sample_thresholds)
    save_result(adata, "qc", "filtering_summary", filtering_summary)
    save_result(adata, "qc", "warnings", warnings)


def _run_qc_workflow(
    adata: AnnData,
    config: QCWorkflowConfig,
    results_path: Optional[Path],
    show_progress: bool = True,
) -> Tuple[AnnData, Any, Dict[str, Any], List[str]]:
    """
    Run QC workflow.

    Args:
        adata: Input AnnData object
        config: QC workflow configuration
        results_path: Path to results directory. If None, no files will be saved.
        show_progress: Whether to show progress bars for multi-sample processing

    Returns:
        Tuple of (AnnData object with QC completed, recommendation, sample_thresholds, warnings)
    """
    warnings_list: List[str] = []
    recommendation = None
    sample_thresholds: Dict[str, Any] = {}

    active_tissue_type = config.tissue_type or "auto"

    _ensure_sample_key(adata, config, warnings_list)

    # Snapshot the original user config before any modifications
    original_config = config.model_copy(deep=True)

    # --- 0. Recommendation & threshold policy ---
    if config.use_recommendations:
        try:
            from .intelligent_qc import recommend_intelligent_qc
            recommendation = recommend_intelligent_qc(adata, tissue_type=active_tissue_type)
            config, original_config = _apply_qc_recommendations(config, recommendation)
            log.info("Intelligent QC recommendations applied to config.")
        except Exception as e:
            msg = f"Intelligent QC recommendation failed: {e}"
            warnings_list.append(msg)
            log.warning(msg)

    n_samples = int(adata.obs[config.sample_key].nunique()) if config.sample_key in adata.obs else 1
    if n_samples > 1 and config.threshold_mode != "pooled":
        sample_thresholds, policy_warnings = _compute_sample_thresholds(adata, config)
        warnings_list.extend(policy_warnings)
        if sample_thresholds:
            log.info(f"Computed {config.threshold_mode} per-sample thresholds for {len(sample_thresholds)} samples.")

    # Tumor-aware filtering adjustment
    if _is_tumor_aware(active_tissue_type):
        if "outlier_mt" in config.filter_config.criteria_to_filter:
            config.filter_config.criteria_to_filter = [
                c for c in config.filter_config.criteria_to_filter if c != "outlier_mt"
            ]
            msg = "Tumor-aware QC: outlier_mt excluded from filtering criteria."
            warnings_list.append(msg)
            log.info(msg)

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
    if config.doublet_config.run_algorithm or config.doublet_config.use_heuristics:
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
        sample_thresholds=sample_thresholds,
    )

    return adata, config, recommendation, sample_thresholds, warnings_list, original_config


def run_standard_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    overwrite: bool = False,
    *,
    tissue_type: str = "unknown",
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
    runtime_config = _prepare_runtime_qc_config(config, tissue_type)

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

    adata, results_path = _setup_workflow(adata, runtime_config.save_dir, overwrite)

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
            adata, applied_config, recommendation, sample_thresholds, qc_warnings, original_config = _run_qc_workflow(
                adata, runtime_config, results_path, show_progress=show_progress
            )
            _add_tumor_aware_flags(adata, applied_config.tissue_type or tissue_type)
            successful_steps.append("qc_metrics")
        else:
            log.info("Step: QC Metrics (already completed, skipping)")
            applied_config = runtime_config
            recommendation = None
            sample_thresholds = {}
            qc_warnings = []
            original_config = runtime_config.model_copy(deep=True)

        # --- Step 2: Filtering ---
        if "filtering" not in completed_steps:
            current_step = "filtering"
            log.info("Step: Cell Filtering")
            adata_filtered = filter_cells(adata, config=applied_config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (already completed, using cached)")
            adata_filtered = adata

        filtering_summary = adata_filtered.uns.get("sclucid", {}).get("qc", {}).get("filtering_results", {})
        _store_qc_trace(
            adata_filtered,
            applied_config,
            original_config,
            recommendation,
            sample_thresholds,
            filtering_summary,
            qc_warnings,
        )

        # --- Step 3: Reporting ---
        if results_path is not None and "reporting" not in completed_steps:
            current_step = "reporting"
            log.info("Step: Report Generation")
            generate_qc_report(
                adata_filtered,
                save_dir=results_path / "report",
                sample_key=applied_config.sample_key,
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
            manager.save(adata, checkpoint, applied_config if 'applied_config' in locals() else runtime_config)

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
        config=applied_config.to_dict()
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
    tissue_type: str = "unknown",
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
    runtime_config = _prepare_runtime_qc_config(config, tissue_type)

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

    adata, results_path = _setup_workflow(adata, runtime_config.save_dir, overwrite)

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
            adata, applied_config, recommendation, sample_thresholds, qc_warnings, original_config = _run_qc_workflow(
                adata, runtime_config, results_path, show_progress=show_progress
            )
            _add_tumor_aware_flags(adata, applied_config.tissue_type or tissue_type)
            successful_steps.append("qc_metrics")
        else:
            log.info("Step: QC Metrics (already completed, skipping)")
            applied_config = runtime_config
            recommendation = None
            sample_thresholds = {}
            qc_warnings = []
            original_config = runtime_config.model_copy(deep=True)

        # --- Step 2: Filtering ---
        if "filtering" not in completed_steps:
            current_step = "filtering"
            log.info("Step: Cell Filtering")
            adata_filtered = filter_cells(adata, config=applied_config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (already completed, using cached)")
            adata_filtered = adata

        filtering_summary = adata_filtered.uns.get("sclucid", {}).get("qc", {}).get("filtering_results", {})
        _store_qc_trace(
            adata_filtered,
            applied_config,
            original_config,
            recommendation,
            sample_thresholds,
            filtering_summary,
            qc_warnings,
        )

        # --- Step 3: Reporting ---
        if results_path is not None and "reporting" not in completed_steps:
            current_step = "reporting"
            log.info("Step: Report Generation")
            generate_qc_report(
                adata_filtered,
                save_dir=results_path / "report",
                sample_key=applied_config.sample_key,
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
            manager.save(adata, checkpoint, applied_config if 'applied_config' in locals() else runtime_config)

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
        config=applied_config.to_dict()
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
