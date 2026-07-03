"""Runtime helpers for QC workflow orchestration."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TypeVar

import numpy as np
from anndata import AnnData

from ..utils import get_progress_bar, save_result
from .config import QCWorkflowConfig
from .doublet import predict_doublets
from .metrics import calculate_qc_metric

log = logging.getLogger(__name__)
T = TypeVar("T")


# Define workflow steps for flexible execution
QC_WORKFLOW_STEPS = [
    "qc_metrics",
    "filtering",
    "reporting",
]

def _restore_empty_config_values(value: Any) -> Any:
    """Restore storage-normalized empty strings to ``None`` for config parsing."""
    if value == "":
        return None
    if isinstance(value, dict):
        return {key: _restore_empty_config_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_empty_config_values(item) for item in value]
    return value


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
    return get_progress_bar(iterable, desc=desc, enabled=enabled, total=total, unit="sample")


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
    if (
        not tissue_type
        or "tumor" not in tissue_type.lower()
        and "cancer" not in tissue_type.lower()
    ):
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
            tissue_type=config.tissue_type,
            sample_context_key=config.sample_context_key,
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
            tissue_type=config.tissue_type,
            sample_context_key=config.sample_context_key,
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

    doublet_group_key = doublet_config.detection_group_key or config.sample_key

    sample_adata = predict_doublets(
        sample_adata,
        config=doublet_config,
        sample_key=doublet_group_key,
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




def _resolve_qc_steps(
    steps: Optional[List[str]],
    skip_steps: Optional[List[str]],
    completed_steps: Optional[List[str]] = None,
) -> List[str]:
    """Resolve which QC steps to run."""
    if steps is not None and skip_steps is not None:
        raise ValueError("Cannot specify both 'steps' and 'skip_steps'. Choose one.")

    if steps is not None:
        resolved = list(steps)
    elif skip_steps is not None:
        resolved = [s for s in QC_WORKFLOW_STEPS if s not in skip_steps]
    else:
        resolved = QC_WORKFLOW_STEPS.copy()

    invalid = set(resolved) - set(QC_WORKFLOW_STEPS)
    if invalid:
        raise ValueError(f"Invalid step names: {invalid}. Valid steps are: {QC_WORKFLOW_STEPS}")

    if completed_steps:
        resolved = [s for s in resolved if s not in completed_steps]

    return resolved
