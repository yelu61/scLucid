"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import json
import logging
import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, TypeVar

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from anndata import AnnData

from ..runtime import effective_n_jobs
from ..utils import (
    PartialResultManager,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    get_progress_bar,
    normalize_review_summary,
    record_artifact,
    save_result,
    save_workflow_result,
    validate_review_summary_schema,
)
from ..utils.context import is_tumor_context, resolve_cell_type_key
from .ambient import (
    AMBIENT_CORRECTED_COUNTS_LAYER,
    correct_ambient_rna_linear,
    diagnose_ambient_rna,
    diagnose_empty_droplets,
    infer_ambient_input_context,
    record_ambient_layer_contract,
)
from .ambient_backends import correct_ambient_rna as correct_ambient_rna_unified
from .benchmark import evaluate_qc_benchmark, export_qc_benchmark_report
from .config import FilterConfig, QCWorkflowConfig
from .decisions import build_qc_decisions, summarize_qc_decisions, score_qc_gene_panels
from .doublet import predict_doublets
from .filtering import (
    filter_cells,
    mark_low_quality_cell,
    resolve_qc_thresholds,
)
from .filtering.core import AdaptiveThresholdCalculator
from .metrics import calculate_qc_metric
from .reporting import generate_qc_report
from .trace import enrich_qc_review_summary, validate_qc_review_summary

log = logging.getLogger(__name__)

T = TypeVar("T")

# Define workflow steps for flexible execution
QC_WORKFLOW_STEPS = [
    "qc_metrics",
    "filtering",
    "reporting",
]

# Keep for backward compatibility
QCWorkflowError = WorkflowError
PartialQCResult = PartialResultManager


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


def _is_tumor_aware(tissue_type: Optional[str]) -> bool:
    return is_tumor_context(tissue_type)


def _user_explicitly_set_threshold(
    config: QCWorkflowConfig,
    threshold_name: str,
) -> bool:
    """Return whether a nested QC threshold was explicitly supplied by the user."""
    if "marking_config" not in getattr(config, "model_fields_set", set()):
        return False
    marking_config = getattr(config, "marking_config", None)
    if marking_config is None or "thresholds" not in getattr(marking_config, "model_fields_set", set()):
        return False
    thresholds = getattr(marking_config, "thresholds", None)
    return threshold_name in getattr(thresholds, "model_fields_set", set())


def _resolve_filter_config_for_decisions(config: QCWorkflowConfig) -> FilterConfig:
    """Return the filter config after applying optional qc_decision filtering."""
    mode = getattr(config, "qc_decision_filter_mode", "off")
    if mode == "off":
        return config.filter_config
    if mode == "replace":
        return FilterConfig(
            criteria_to_filter=["qc_remove"],
            combination_logic="any",
        )
    criteria = list(config.filter_config.criteria_to_filter or [])
    if "qc_remove" not in criteria:
        criteria.append("qc_remove")
    return config.filter_config.model_copy(update={"criteria_to_filter": criteria})


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

    intelligent_thresholds: Dict[str, Any] = {}
    threshold_map = {
        "min_genes": "min_genes",
        "max_mt_percent": "pc_mt",
        "n_counts": "min_counts",
    }
    for rec_key, threshold_key in threshold_map.items():
        rec_value = rec_dict.get(rec_key)
        if isinstance(rec_value, dict) and rec_value.get("threshold") is not None:
            intelligent_thresholds[threshold_key] = rec_value["threshold"]

    manual_thresholds: Dict[str, Any] = {}
    if _is_user_set(config, "marking_config") and _is_user_set(
        config.marking_config, "thresholds"
    ):
        for field_name in config.marking_config.thresholds.model_fields_set:
            manual_thresholds[field_name] = getattr(config.marking_config.thresholds, field_name)

    if intelligent_thresholds or manual_thresholds:
        # Thresholds the user set explicitly (tracked via model_fields_set) are
        # authoritative and must survive recommendations verbatim, so use the
        # manual_override policy rather than treating them as floor/ceiling bounds.
        applied.marking_config.thresholds = resolve_qc_thresholds(
            intelligent=intelligent_thresholds,
            manual=manual_thresholds,
            policy="manual_override",
        )

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
    explicit_field_checks = {
        "min_genes": [
            ("marking_config",),
            ("thresholds",),
            ("min_genes",),
        ],
        "max_mt_percent": [
            ("marking_config",),
            ("thresholds",),
            ("pc_mt",),
        ],
        "n_counts": [
            ("marking_config",),
            ("thresholds",),
            ("min_counts",),
        ],
        "doublet_threshold": [
            ("doublet_config",),
            ("score_threshold",),
        ],
    }

    def _is_explicit_user_path(config_obj: Any, fields: list[tuple[str, ...]]) -> bool:
        current = config_obj
        for field_path in fields:
            field_name = field_path[0]
            if current is None or field_name not in getattr(current, "model_fields_set", set()):
                return False
            current = getattr(current, field_name, None)
        return True

    for param_name, path in mapping.items():
        rec_val = None
        param_rec = rec_dict.get(param_name)
        if isinstance(param_rec, dict):
            rec_val = param_rec.get("threshold")

        if not _is_explicit_user_path(original_config, explicit_field_checks[param_name]):
            continue

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


def _build_qc_review_summary(
    config: QCWorkflowConfig,
    original_config: QCWorkflowConfig,
    recommendation: Any,
    sample_thresholds: Dict[str, Any],
    filtering_summary: Dict[str, Any],
    warnings: List[str],
    ambient_summary: Optional[Dict[str, Any]] = None,
    empty_droplet_summary: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a human-reviewable summary of the QC run.

    This distills the full trace into the artifacts a reviewer needs:
    what was recommended, what was actually applied, what the user
    overrode, per-sample thresholds, and any tumor-aware cautions.
    """
    summary: Dict[str, Any] = {}

    # --- Recommendation summary ---
    rec_summary: Dict[str, Any] = {"available": recommendation is not None}
    if recommendation is not None:
        rec_dict = recommendation.to_dict() if hasattr(recommendation, "to_dict") else {}
        rec_summary["overall_strategy"] = rec_dict.get("overall_strategy", "unknown")
        rec_summary["overall_confidence"] = rec_dict.get("overall_confidence")
        rec_summary["data_quality_score"] = rec_dict.get("data_quality_score")
        rec_summary["concerns"] = rec_dict.get("concerns", [])
        key_thresholds: Dict[str, Any] = {}
        for param, rec_key, path in [
            ("min_genes", "min_genes", ("marking_config", "thresholds", "min_genes")),
            ("max_mt_percent", "max_mt_percent", ("marking_config", "thresholds", "pc_mt")),
            ("n_counts", "n_counts", ("marking_config", "thresholds", "min_counts")),
            ("doublet_threshold", "doublet_threshold", ("doublet_config", "score_threshold")),
        ]:
            rec_val = (
                rec_dict.get(rec_key, {}).get("threshold")
                if isinstance(rec_dict.get(rec_key), dict)
                else None
            )
            cfg_val = original_config.to_dict()
            for key in path:
                cfg_val = cfg_val.get(key) if isinstance(cfg_val, dict) else None
            key_thresholds[param] = {
                "recommended": rec_val,
                "user_provided": cfg_val,
            }
        rec_summary["key_thresholds"] = key_thresholds
    summary["recommendation_summary"] = rec_summary

    # --- Applied threshold summary ---
    th = config.marking_config.thresholds
    summary["applied_threshold_summary"] = {
        "min_genes": th.min_genes,
        "max_genes": th.max_genes,
        "min_counts": th.min_counts,
        "max_counts": th.max_counts,
        "pc_mt": th.pc_mt,
        "pc_hb": th.pc_hb,
        "nmads": th.nmads,
    }

    # --- User override summary ---
    overrides = _diff_qc_recommendations(recommendation, original_config)
    summary["user_override_summary"] = {
        "overrides_detected": bool(overrides),
        "details": overrides,
        "note": (
            "User-specified thresholds take precedence over recommendations. "
            "Empty details means the user accepted all recommendations or no recommendation was generated."
        ),
    }

    # --- Sample-level threshold summary ---
    n_samples = len(sample_thresholds)
    summary["sample_threshold_summary"] = {
        "mode": config.threshold_mode,
        "n_samples_with_thresholds": n_samples,
        "per_sample": (
            {
                sample: {
                    metric: {
                        "lower": (
                            round(vals["lower"], 2)
                            if isinstance(vals.get("lower"), (int, float))
                            else vals.get("lower")
                        ),
                        "upper": (
                            round(vals["upper"], 2)
                            if isinstance(vals.get("upper"), (int, float))
                            else vals.get("upper")
                        ),
                    }
                    for metric, vals in thresholds.items()
                }
                for sample, thresholds in sample_thresholds.items()
            }
            if sample_thresholds
            else {}
        ),
        "note": (
            "Per-sample thresholds are only computed in hierarchical/independent mode with >1 sample. "
            "Pooled mode uses a single global threshold."
        ),
    }

    # --- Tumor-aware summary ---
    is_tumor = _is_tumor_aware(config.tissue_type)
    tumor_notes: List[str] = []
    if is_tumor:
        tumor_notes.append(
            "Tumor-aware QC is active: elevated mitochondrial content is flagged rather than hard-filtered."
        )
        if "outlier_mt" not in config.filter_config.criteria_to_filter:
            tumor_notes.append(
                "Mitochondrial outlier filtering was disabled for this tumor dataset."
            )
        if config.marking_config.thresholds.pc_mt is not None:
            tumor_notes.append(
                "The mitochondrial threshold is retained as a warning signal for review and reporting."
            )
    tumor_warnings = [note for note in tumor_notes if "disabled" in note or "warning" in note]
    summary["tumor_aware_summary"] = {
        "enabled": is_tumor,
        "tissue_type": config.tissue_type,
        "notes": tumor_notes,
        "warnings": tumor_warnings,
        "filtering_criteria": list(config.filter_config.criteria_to_filter),
        "mitochondrial_filtering_enabled": "outlier_mt" in config.filter_config.criteria_to_filter,
    }

    # --- Filtering summary ---
    fs = filtering_summary if isinstance(filtering_summary, dict) else {}
    summary["filtering_summary"] = {
        "initial_cells": fs.get("initial_cells"),
        "final_cells": fs.get("final_cells"),
        "removed_cells": fs.get("removed_cells"),
        "removed_fraction": fs.get("removed_fraction"),
        "criteria_used": fs.get("criteria_used", config.filter_config.criteria_to_filter),
        "criteria_counts": fs.get("criteria_counts", {}),
        "review_criteria_counts": fs.get("review_criteria_counts", {}),
    }

    # --- Warnings ---
    summary["warnings"] = warnings
    if ambient_summary is not None:
        summary["ambient_rna_summary"] = ambient_summary
    if empty_droplet_summary is not None:
        summary["empty_droplet_summary"] = empty_droplet_summary

    return summary


def _sync_ambient_corrected_layer_to_output(
    *,
    source: AnnData,
    target: AnnData,
    output_layer: str,
) -> bool:
    """Copy an ambient-corrected layer from the review source to the final output."""
    if output_layer not in source.layers:
        return False
    if not target.var_names.equals(source.var_names):
        return False
    try:
        corrected = source.layers[output_layer]
        target.layers[output_layer] = corrected[target.obs_names, :].copy()
        return True
    except Exception:
        try:
            row_index = source.obs_names.get_indexer(target.obs_names)
            if np.any(row_index < 0):
                return False
            target.layers[output_layer] = source.layers[output_layer][row_index, :].copy()
            return True
        except Exception:
            return False


def _store_qc_trace(
    adata: AnnData,
    config: QCWorkflowConfig,
    original_config: QCWorkflowConfig,
    recommendation: Any,
    sample_thresholds: Dict[str, Any],
    filtering_summary: Dict[str, Any],
    warnings: List[str],
    steps_executed: Optional[List[str]] = None,
    adata_before_filtering: Optional[AnnData] = None,
) -> None:
    """Store unified QC trace under adata.uns['sclucid']['qc']."""
    filtering_summary = dict(filtering_summary or {})
    review_input = adata_before_filtering if adata_before_filtering is not None else adata
    if (
        _is_tumor_aware(config.tissue_type)
        and config.marking_config.thresholds.pc_mt is not None
        and "pct_counts_mt" in review_input.obs
    ):
        mt_values = np.asarray(review_input.obs["pct_counts_mt"], dtype=float)
        mt_threshold = float(config.marking_config.thresholds.pc_mt)
        review_counts = dict(filtering_summary.get("review_criteria_counts", {}) or {})
        review_counts["outlier_mt"] = int(np.sum(mt_values >= mt_threshold))
        filtering_summary["review_criteria_counts"] = review_counts

    n_samples = int(adata.obs[config.sample_key].nunique()) if config.sample_key in adata.obs else 1
    context = {
        "sample_key": config.sample_key,
        "threshold_mode": config.threshold_mode,
        "n_samples": n_samples,
        "tissue_type": config.tissue_type,
        "use_recommendations": config.use_recommendations,
    }
    save_result(
        adata,
        "qc",
        "context",
        context,
    )
    if recommendation is not None:
        save_result(adata, "qc", "recommendation", recommendation.to_dict())
    save_result(adata, "qc", "original_config", original_config.to_dict())
    save_result(adata, "qc", "applied_config", config.to_dict())
    save_result(
        adata, "qc", "user_overrides", _diff_qc_recommendations(recommendation, original_config)
    )
    save_result(adata, "qc", "sample_thresholds", sample_thresholds)
    save_result(adata, "qc", "filtering_summary", filtering_summary)
    save_result(adata, "qc", "warnings", warnings)
    benchmark_summary = None
    if adata_before_filtering is not None:
        benchmark_summary = evaluate_qc_benchmark(
            adata_before_filtering,
            adata,
            tissue_type=config.tissue_type,
            tissue=config.tissue,
            sample_key=config.sample_key,
            cell_type_key=_detect_cell_type_key(adata_before_filtering),
        )
        save_result(adata, "qc", "benchmark_summary", benchmark_summary)

    ambient_input = adata_before_filtering if adata_before_filtering is not None else adata
    try:
        ambient_input_context = infer_ambient_input_context(ambient_input)
        ambient_summary = diagnose_ambient_rna(ambient_input)
        ambient_summary["input_context"] = ambient_input_context
        save_result(adata, "qc", "ambient_rna_summary", ambient_summary)
        empty_droplet_summary = diagnose_empty_droplets(ambient_input)
        save_result(adata, "qc", "empty_droplet_summary", empty_droplet_summary)
        if ambient_summary.get("risk_level") in {"moderate", "high"}:
            warnings.append(
                "Ambient RNA risk is "
                f"{ambient_summary.get('risk_level')}; inspect ambient_rna_summary "
                "and consider Python backends such as CellBender or scAR."
            )
            if config.ambient_correction in {"linear", "auto"}:
                try:
                    if config.ambient_correction == "auto":
                        correction_summary = correct_ambient_rna_unified(
                            ambient_input,
                            method="auto",
                            output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                            empty_droplet_key="likely_empty_droplet"
                            if "likely_empty_droplet" in ambient_input.obs.columns
                            else None,
                        )
                    else:
                        correction_summary = correct_ambient_rna_linear(
                            ambient_input,
                            output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                            empty_droplet_key="likely_empty_droplet"
                            if "likely_empty_droplet" in ambient_input.obs.columns
                            else None,
                        )
                    if correction_summary.get("corrected"):
                        correction_summary["output_layer_synced_to_filtered_adata"] = (
                            _sync_ambient_corrected_layer_to_output(
                                source=ambient_input,
                                target=adata,
                                output_layer=str(
                                    correction_summary.get(
                                        "output_layer",
                                        AMBIENT_CORRECTED_COUNTS_LAYER,
                                    )
                                ),
                            )
                        )
                        backend = correction_summary.get("backend", "linear")
                        warnings.append(
                            f"Applied {backend} ambient RNA correction to layer "
                            f"'{correction_summary.get('output_layer')}'. "
                            f"Removed {correction_summary.get('removed_counts', 0):.0f} "
                            f"counts (mean rho = {correction_summary.get('mean_rho', 0):.3f})."
                        )
                    save_result(adata, "qc", "ambient_correction_summary", correction_summary)
                    record_ambient_layer_contract(
                        adata,
                        input_context=ambient_input_context,
                        correction_summary=correction_summary,
                        output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                    )
                except Exception as corr_exc:
                    warnings.append(f"Ambient RNA correction failed: {corr_exc}")
                    record_ambient_layer_contract(
                        adata,
                        input_context=ambient_input_context,
                        correction_summary={
                            "corrected": False,
                            "reason": f"ambient correction failed: {corr_exc}",
                            "review_required": True,
                        },
                        output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                    )
            else:
                record_ambient_layer_contract(
                    adata,
                    input_context=ambient_input_context,
                    correction_summary={"corrected": False, "reason": "correction_not_requested"},
                    output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                )
        if "ambient_layer_contract" not in adata.uns.get("sclucid", {}).get("qc", {}):
            record_ambient_layer_contract(
                adata,
                input_context=ambient_input_context,
                correction_summary={
                    "corrected": False,
                    "reason": "ambient_risk_below_correction_threshold",
                },
                output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
            )
        save_result(adata, "qc", "warnings", warnings)
    except Exception as exc:
        warnings.append(f"Ambient RNA diagnostics failed: {exc}")
        ambient_summary = {
            "available": False,
            "risk_level": "unknown",
            "reason": f"ambient diagnostic failed: {exc}",
            "review_required": True,
        }
        empty_droplet_summary = {
            "available": False,
            "risk_level": "unknown",
            "reason": f"empty-droplet diagnostic failed: {exc}",
            "review_required": True,
        }
        save_result(adata, "qc", "ambient_rna_summary", ambient_summary)
        save_result(adata, "qc", "empty_droplet_summary", empty_droplet_summary)
        save_result(adata, "qc", "warnings", warnings)

    # Build and store the review-facing summary
    base_review_summary = _build_qc_review_summary(
        config,
        original_config,
        recommendation,
        sample_thresholds,
        filtering_summary,
        warnings,
        ambient_summary=ambient_summary,
        empty_droplet_summary=empty_droplet_summary,
    )
    if benchmark_summary is not None:
        base_review_summary["benchmark_summary"] = benchmark_summary

    review_summary = normalize_review_summary(
        enrich_qc_review_summary(
            base_review_summary,
            adata=adata,
            config=config,
            original_config=original_config,
            recommendation=recommendation,
            sample_thresholds=sample_thresholds,
            filtering_summary=filtering_summary,
            warnings=warnings,
            context=context,
            steps_executed=steps_executed,
            adata_before_filtering=adata_before_filtering,
        ),
        module="qc",
        workflow_name="standard",
        adata=adata,
        steps_executed=steps_executed or [],
        config=config.to_dict(),
        warnings=warnings,
    )
    validate_review_summary_schema(review_summary, module="qc", raise_on_error=True)
    validate_qc_review_summary(review_summary, raise_on_error=True)
    save_result(adata, "qc", UnsKeys.REVIEW_SUMMARY, review_summary)
    return review_summary


def _detect_cell_type_key(adata: AnnData) -> Optional[str]:
    """Detect a likely cell type annotation column for benchmark stratification."""
    return resolve_cell_type_key(adata)


def _export_qc_review_summary(
    review_summary: Dict[str, Any],
    save_dir: Path,
    adata: Optional[AnnData] = None,
) -> None:
    """Export review summary as JSON and Markdown sidecars."""
    save_dir.mkdir(parents=True, exist_ok=True)

    # JSON sidecar
    json_path = save_dir / "qc_review_summary.json"
    json_path.write_text(json.dumps(review_summary, indent=2, default=str), encoding="utf-8")

    # Markdown sidecar
    md_lines = [
        "# QC Review Summary",
        "",
        "## Recommendation Summary",
        "",
    ]
    rec = review_summary.get("recommendation_summary", {})
    if rec.get("available"):
        md_lines.append(f"- **Strategy**: {rec.get('overall_strategy', 'unknown')}")
        md_lines.append(f"- **Confidence**: {rec.get('overall_confidence')}")
        md_lines.append(f"- **Data Quality Score**: {rec.get('data_quality_score')}")
        if rec.get("concerns"):
            md_lines.append("- **Concerns**:")
            for c in rec["concerns"]:
                md_lines.append(f"  - {c}")
        md_lines.append("")
        md_lines.append("| Parameter | Recommended | User Provided |")
        md_lines.append("|-----------|-------------|---------------|")
        for param, vals in rec.get("key_thresholds", {}).items():
            md_lines.append(
                f"| {param} | {vals.get('recommended')} | {vals.get('user_provided')} |"
            )
    else:
        md_lines.append(
            "- No recommendation was generated (recommendations disabled or engine failed)."
        )
    md_lines.append("")

    readiness = review_summary.get("qc_readiness", {})
    md_lines.extend(
        [
            "## QC Readiness",
            "",
            f"- **Status**: {readiness.get('status')}",
            f"- **Score**: {readiness.get('score')}",
            f"- **Verdict**: {readiness.get('verdict')}",
            "",
        ]
    )
    if readiness.get("blockers"):
        md_lines.append("- **Blockers**:")
        for blocker in readiness.get("blockers", []):
            md_lines.append(f"  - {blocker}")
    if readiness.get("review_reasons"):
        md_lines.append("- **Review reasons**:")
        for reason in readiness.get("review_reasons", []):
            md_lines.append(f"  - {reason}")
    md_lines.append("")

    action_items = review_summary.get("review_action_items", [])
    if isinstance(action_items, dict):
        action_items = list(action_items.values())
    if action_items:
        md_lines.extend(
            [
                "## Review Action Items",
                "",
                "| Priority | Action | Rationale |",
                "|----------|--------|-----------|",
            ]
        )
        for item in action_items:
            md_lines.append(
                "| {priority} | {action} | {rationale} |".format(
                    priority=item.get("priority"),
                    action=item.get("action"),
                    rationale=item.get("rationale"),
                )
            )
        md_lines.append("")

    md_lines.extend(
        [
            "## Threshold Reviewer Table",
            "",
            "| Parameter | Recommended | Applied | Source | Confidence | Affected Cells | Risk Note | Review Required |",
            "|-----------|-------------|---------|--------|------------|----------------|-----------|-----------------|",
        ]
    )
    decision_table = review_summary.get(
        "threshold_reviewer_table",
        review_summary.get("decision_table", []),
    )
    if isinstance(decision_table, dict):
        decision_table = list(decision_table.values())
    for row in decision_table:
        md_lines.append(
            "| {parameter} | {recommended} | {applied} | {source} | {confidence} | {affected} | {risk_note} | {review_required} |".format(
                parameter=row.get("parameter"),
                recommended=row.get("recommended"),
                applied=row.get("applied"),
                source=row.get("source"),
                confidence=row.get("confidence"),
                affected=row.get("affected_cells"),
                risk_note=str(row.get("risk_note") or "").replace("|", "/"),
                review_required=row.get("review_required"),
            )
        )
    md_lines.append("")

    qc_decisions = review_summary.get("qc_decision_summary", {})
    if qc_decisions:
        decision_counts = qc_decisions.get("decision_counts", {})
        evidence_summary = qc_decisions.get("evidence_summary", {})
        md_lines.extend(
            [
                "## QC Decision Summary",
                "",
                f"- **Policy**: {qc_decisions.get('policy')}",
                f"- **Cells**: {qc_decisions.get('n_cells')}",
                f"- **Review required cells**: {qc_decisions.get('review_required_cells')}",
                f"- **Risk note**: {qc_decisions.get('risk_note')}",
                "",
                "| Decision | Cells |",
                "|----------|-------|",
            ]
        )
        if isinstance(decision_counts, dict):
            for decision, count in decision_counts.items():
                md_lines.append(f"| {decision} | {count} |")
        md_lines.extend(["", "| Evidence | Cells |", "|----------|-------|"])
        if isinstance(evidence_summary, dict):
            for evidence, count in evidence_summary.items():
                md_lines.append(f"| {evidence} | {count} |")
        md_lines.append("")

    health = review_summary.get("output_health", {})
    md_lines.extend(
        [
            "## Output Health",
            "",
            f"- **Status**: {health.get('status')}",
            f"- **Cells**: {health.get('n_cells')}",
            f"- **Genes**: {health.get('n_genes')}",
        ]
    )
    issues = health.get("issues", [])
    if isinstance(issues, dict):
        issues = list(issues.values())
    if issues:
        md_lines.append("- **Issues**:")
        for issue in issues:
            md_lines.append(f"  - {issue}")
    md_lines.append("")

    benchmark = review_summary.get("benchmark_summary", {})
    if benchmark:
        retention = benchmark.get("retention", {})
        marker = benchmark.get("marker_fidelity", {})
        md_lines.extend(
            [
                "## Benchmark Summary",
                "",
                f"- **Profile**: {benchmark.get('profile_label')} ({benchmark.get('profile')})",
                f"- **Status**: {benchmark.get('status')}",
                f"- **Retention rate**: {retention.get('retention_rate')}",
                f"- **Marker fidelity**: {marker.get('overall_marker_fidelity')}",
                "",
            ]
        )

    md_lines.extend(
        [
            "## Applied Thresholds",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
        ]
    )
    for param, val in review_summary.get("applied_threshold_summary", {}).items():
        md_lines.append(f"| {param} | {val} |")
    md_lines.append("")

    ov = review_summary.get("user_override_summary", {})
    md_lines.extend(
        [
            "## User Overrides",
            "",
            f"- **Overrides detected**: {ov.get('overrides_detected', False)}",
        ]
    )
    if ov.get("details"):
        md_lines.append("- **Details**:")
        for param, vals in ov["details"].items():
            md_lines.append(
                f"  - {param}: recommended={vals.get('recommended')}, user={vals.get('actual')}"
            )
    md_lines.append("")

    st = review_summary.get("sample_threshold_summary", {})
    md_lines.extend(
        [
            "## Sample-Level Thresholds",
            "",
            f"- **Mode**: {st.get('mode')}",
            f"- **Samples with thresholds**: {st.get('n_samples_with_thresholds', 0)}",
            "",
        ]
    )
    if st.get("per_sample"):
        md_lines.append("```json")
        md_lines.append(json.dumps(st["per_sample"], indent=2, default=str))
        md_lines.append("```")
    md_lines.append("")

    ta = review_summary.get("tumor_aware_summary", {})
    md_lines.extend(
        [
            "## Tumor-Aware QC",
            "",
            f"- **Enabled**: {ta.get('enabled', False)}",
        ]
    )
    if ta.get("notes"):
        for note in ta["notes"]:
            md_lines.append(f"- {note}")
    md_lines.append("")

    doublet_summary = review_summary.get("doublet_evidence_summary", {})
    if doublet_summary:
        final_doublets = (
            doublet_summary.get("predictions", {}).get("predicted_doublet", {})
            if isinstance(doublet_summary.get("predictions"), dict)
            else {}
        )
        predicted_fraction = final_doublets.get("fraction")
        predicted_fraction_str = (
            f"{predicted_fraction:.1%}" if isinstance(predicted_fraction, (int, float)) else "N/A"
        )
        md_lines.extend(
            [
                "## Doublet Evidence",
                "",
                f"- **Status**: {doublet_summary.get('status')}",
                f"- **Predicted doublets**: {final_doublets.get('count', 'N/A')} ({predicted_fraction_str})",
                f"- **Review required**: {doublet_summary.get('review_required', False)}",
            ]
        )
        notes = doublet_summary.get("notes", [])
        if isinstance(notes, dict):
            notes = list(notes.values())
        for note in notes:
            md_lines.append(f"- {note}")
        method_keys = doublet_summary.get("method_metadata_keys", [])
        if method_keys:
            md_lines.append(f"- **Method metadata**: {method_keys}")
        benchmark_decision = doublet_summary.get("benchmark_decision", {})
        if isinstance(benchmark_decision, dict) and benchmark_decision:
            md_lines.extend(
                [
                    "- **Benchmark decision**:",
                    f"  - recommended_default_mode: {benchmark_decision.get('recommended_default_mode', 'N/A')}",
                    f"  - recommended_primary_method: {benchmark_decision.get('recommended_primary_method', 'N/A')}",
                    f"  - recommended_algorithm_weight: {benchmark_decision.get('recommended_algorithm_weight', 'N/A')}",
                    f"  - review_required: {benchmark_decision.get('review_required', False)}",
                ]
            )
            if benchmark_decision.get("risk_note"):
                md_lines.append(f"  - risk_note: {benchmark_decision['risk_note']}")
        benchmark_evidence = doublet_summary.get("benchmark_evidence", {})
        if benchmark_evidence:
            md_lines.append("- **Benchmark evidence**:")
            if isinstance(benchmark_evidence, dict):
                for key, value in benchmark_evidence.items():
                    md_lines.append(f"  - {key}: {value}")
            else:
                md_lines.append(f"  - {benchmark_evidence}")
        md_lines.append("")

    downstream = review_summary.get("downstream_preprocess_recommendations", {})
    md_lines.extend(
        [
            "## Downstream Preprocess Recommendations",
            "",
            f"- **Status**: {downstream.get('status')}",
            f"- **Ready for preprocess**: {downstream.get('ready_for_preprocess')}",
            "",
        ]
    )
    _recs = downstream.get("recommendations", [])
    if isinstance(_recs, dict):
        _recs = list(_recs.values())
    for item in _recs:
        md_lines.append(
            "- **{target}** ({priority}): {recommendation}".format(
                target=item.get("target"),
                priority=item.get("priority"),
                recommendation=item.get("recommendation"),
            )
        )
    md_lines.append("")

    fs = review_summary.get("filtering_summary", {})
    _removed_frac = fs.get("removed_fraction")
    _removed_frac_str = f"{_removed_frac:.1%}" if isinstance(_removed_frac, (int, float)) else "N/A"
    md_lines.extend(
        [
            "## Filtering Results",
            "",
            f"- **Initial cells**: {fs.get('initial_cells')}",
            f"- **Final cells**: {fs.get('final_cells')}",
            f"- **Removed**: {fs.get('removed_cells')} ({_removed_frac_str})",
            f"- **Criteria used**: {fs.get('criteria_used', [])}",
            "",
        ]
    )

    warnings = review_summary.get("warnings", [])
    if isinstance(warnings, dict):
        warnings = list(warnings.values())
    if warnings:
        md_lines.extend(
            [
                "## Warnings",
                "",
            ]
        )
        for w in warnings:
            md_lines.append(f"- {w}")
        md_lines.append("")

    md_path = save_dir / "qc_review_summary.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    if adata is not None:
        record_artifact(
            adata,
            "qc",
            "qc_review_summary_json",
            str(json_path),
            kind="json",
            description="QC review summary JSON sidecar",
        )
        record_artifact(
            adata,
            "qc",
            "qc_review_summary_md",
            str(md_path),
            kind="md",
            description="QC review summary Markdown sidecar",
        )
    log.info(f"QC review summary exported to {json_path} and {md_path}")


def _run_qc_workflow(
    adata: AnnData,
    config: QCWorkflowConfig,
    results_path: Optional[Path],
    show_progress: bool = True,
) -> Tuple[AnnData, QCWorkflowConfig, Any, Dict[str, Any], List[str], QCWorkflowConfig]:
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
            log.info(
                f"Computed {config.threshold_mode} per-sample thresholds for {len(sample_thresholds)} samples."
            )

    # Tumor-aware filtering adjustment
    if _is_tumor_aware(active_tissue_type):
        # The validation benchmark found that tumor-aware MT review thresholds
        # should not fall below the conventional 20% guardrail.  In tumor data,
        # MT% is treated as a review signal rather than an automatic deletion
        # criterion unless the user explicitly opts into a stricter threshold.
        pc_mt = config.marking_config.thresholds.pc_mt
        if (
            pc_mt is not None
            and pc_mt < 20.0
            and not _user_explicitly_set_threshold(original_config, "pc_mt")
        ):
            config.marking_config.thresholds.pc_mt = 20.0
            msg = (
                "Tumor-aware QC: pc_mt guardrail raised to 20.0 for review-only "
                "mitochondrial signaling; validated tumor-aware policy avoids "
                "mechanical high-MT deletion."
            )
            warnings_list.append(msg)
            log.info(msg)
        elif (
            pc_mt is not None
            and pc_mt < 20.0
            and _user_explicitly_set_threshold(original_config, "pc_mt")
        ):
            msg = (
                "Tumor-aware QC: user-provided pc_mt below 20.0 retained; review "
                "high-MT malignant/stress/program preservation before filtering."
            )
            warnings_list.append(msg)
            log.warning(msg)

        # Relax MAD threshold for tumor data (cells more heterogeneous)
        if config.marking_config.thresholds.nmads < 4.5:
            config.marking_config.thresholds.nmads = 4.5
            msg = "Tumor-aware QC: nmads relaxed to 4.5 for heterogeneous tumor populations."
            warnings_list.append(msg)
            log.info(msg)

        if "outlier_mt" in config.filter_config.criteria_to_filter:
            config.filter_config.criteria_to_filter = [
                c for c in config.filter_config.criteria_to_filter if c != "outlier_mt"
            ]
            msg = "Tumor-aware QC: outlier_mt excluded from filtering criteria."
            warnings_list.append(msg)
            log.info(msg)

        # Write MT review band from recommendation into obs when available.
        if (
            recommendation is not None
            and hasattr(recommendation, "max_mt_percent")
            and recommendation.max_mt_percent is not None
        ):
            mt_pct = np.asarray(adata.obs.get("pct_counts_mt", np.zeros(adata.n_obs)))
            review_lower = recommendation.max_mt_percent.evidence.get("review_band_lower")
            hard_threshold = recommendation.max_mt_percent.threshold
            if review_lower is not None and hard_threshold is not None:
                adata.obs["mt_review_flag"] = (
                    (mt_pct >= review_lower) & (mt_pct < hard_threshold)
                ).astype(int)
                adata.obs["mt_hard_fail"] = (mt_pct >= hard_threshold).astype(int)
                warnings_list.append(
                    "Tumor-aware QC: added 'mt_review_flag' and 'mt_hard_fail' obs columns "
                    "for metabolic-state review."
                )

    # Blood tissue: ensure hemoglobin outlier detection is active
    if active_tissue_type and (
        "blood" in active_tissue_type.lower()
        or "pbmc" in active_tissue_type.lower()
        or "bone marrow" in active_tissue_type.lower()
    ):
        if "outlier_hb" not in config.filter_config.criteria_to_filter:
            config.filter_config.criteria_to_filter.append("outlier_hb")
            msg = "Blood tissue detected: outlier_hb added to filtering criteria."
            warnings_list.append(msg)
            log.info(msg)

    # --- 1. QC metric calculation ---
    samples = adata.obs[config.sample_key].unique()

    active_n_jobs = effective_n_jobs(config.n_jobs, max_jobs=len(samples))
    if config.use_parallel and active_n_jobs != 1 and len(samples) > 1:
        log.info(f"Processing {len(samples)} samples in parallel with {active_n_jobs} jobs")

        results = _safe_parallel_process(
            process_func=_process_sample_qc,
            samples=list(samples),
            sample_data_func=lambda s: adata[adata.obs[config.sample_key] == s].copy(),
            config=config,
            n_jobs=active_n_jobs,
            step_name="QC metric calculation",
            show_progress=show_progress,
        )

        # Filter out failed samples (None results)
        successful_results = [(s, r) for s, r in results if r is not None]
        if not successful_results:
            raise RuntimeError("All samples failed QC metric calculation")

        if len(successful_results) < len(results):
            log.warning(
                f"Proceeding with {len(successful_results)}/{len(results)} successful samples"
            )

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
            try:
                adata = calculate_qc_metric(
                    adata,
                    sample_key=config.sample_key,
                    reporting_config=config.metrics_reporting_config,
                    calculate_cell_cycle=True,
                    cell_cycle_species=config.species,
                )
            except Exception as e:
                log.warning(
                    f"Cell cycle scoring failed for multi-sample QC metrics ({e}). "
                    "Retrying without cell cycle scoring."
                )
                adata = calculate_qc_metric(
                    adata,
                    sample_key=config.sample_key,
                    reporting_config=config.metrics_reporting_config,
                    calculate_cell_cycle=False,
                )

    # --- 2. Doublet detection ---
    if config.doublet_config.run_algorithm or config.doublet_config.use_heuristics:
        if results_path is not None:
            config.doublet_config.save_dir = str(results_path / "doublet")

        doublet_group_key = config.doublet_config.detection_group_key or config.sample_key
        if doublet_group_key not in adata.obs:
            log.warning(
                "Doublet detection group key '%s' not found; falling back to sample_key '%s'.",
                doublet_group_key,
                config.sample_key,
            )
            doublet_group_key = config.sample_key
        doublet_groups = adata.obs[doublet_group_key].unique()
        doublet_n_jobs = effective_n_jobs(config.n_jobs, max_jobs=len(doublet_groups))

        if config.use_parallel and doublet_n_jobs != 1 and len(doublet_groups) > 1:
            # Parallel doublet detection with error handling
            results = _safe_parallel_process(
                process_func=lambda data, cfg, name: _process_sample_doublet(
                    data, cfg, config.doublet_config.save_dir, name
                ),
                samples=list(doublet_groups),
                sample_data_func=lambda s: adata[adata.obs[doublet_group_key] == s].copy(),
                config=config,
                n_jobs=doublet_n_jobs,
                step_name="doublet detection",
                show_progress=show_progress,
            )

            # Filter out failed samples
            successful_results = [(s, r) for s, r in results if r is not None]
            if not successful_results:
                raise RuntimeError("All samples failed doublet detection")

            if len(successful_results) < len(results):
                log.warning(
                    f"Proceeding with {len(successful_results)}/{len(results)} successful samples for doublet detection"
                )

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
                sample_key=doublet_group_key,
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


def run_standard_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    overwrite: bool = False,
    *,
    tissue_type: str = "unknown",
    show_progress: bool = True,
    # Step control
    steps: Optional[List[str]] = None,
    skip_steps: Optional[List[str]] = None,
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

    Default path semantics:
        - ``use_recommendations=True``: intelligent QC recommendations are applied
          to thresholds that the caller did not explicitly set.
        - ``threshold_mode="hierarchical"``: per-sample thresholds are computed when
          multiple samples are present.
        - Tumor-aware adjustment is active when ``tissue_type`` contains "tumor" or
          "cancer" (e.g. ``outlier_mt`` is excluded from filtering).
        - A reviewer-facing summary is stored in
          ``adata.uns['sclucid']['qc']['review_summary']`` and written to disk as
          ``qc_review_summary.json`` / ``qc_review_summary.md`` when ``save_dir`` is set.

    New features in v0.4:
    - Step control via ``steps`` or ``skip_steps`` (consistent with preprocess/analysis)
    - Error recovery with partial result saving
    - Resume from checkpoint

    Args:
        adata_in: Input AnnData object (raw or pre-normalized).
        config: A QCWorkflowConfig object. If None, a default config is created.
        overwrite: If True, overwrite existing results directory.
        show_progress: If True, show progress bars for multi-sample processing.
        steps: Specific steps to run. See ``QC_WORKFLOW_STEPS`` for valid names.
        skip_steps: Steps to skip (alternative to specifying ``steps``).
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
        >>> # Skip reporting step
        >>> adata_filtered = run_standard_qc(adata, skip_steps=["reporting"])
        >>>
        >>> # Run only metrics and filtering
        >>> adata_filtered = run_standard_qc(adata, steps=["qc_metrics", "filtering"])
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
        raise ValueError(
            "recovery_save_dir is required when error_recovery=True and on_error='save'"
        )

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

    # Resolve steps
    steps_to_run = _resolve_qc_steps(steps, skip_steps, completed_steps)
    log.info(f"Steps to run: {steps_to_run}")

    # Track execution
    successful_steps: List[str] = []
    current_step = None

    try:
        # --- Step 1: QC Metrics ---
        if "qc_metrics" in steps_to_run:
            current_step = "qc_metrics"
            log.info("Step: QC Metrics Calculation")
            (
                adata,
                applied_config,
                recommendation,
                sample_thresholds,
                qc_warnings,
                original_config,
            ) = _run_qc_workflow(adata, runtime_config, results_path, show_progress=show_progress)
            _add_tumor_aware_flags(adata, applied_config.tissue_type or tissue_type)
            successful_steps.append("qc_metrics")
        else:
            log.info("Step: QC Metrics (skipped)")
            applied_config = runtime_config
            recommendation = None
            sample_thresholds = {}
            qc_warnings = []
            original_config = runtime_config.model_copy(deep=True)

        # --- Step 2: Filtering ---
        if "filtering" in steps_to_run:
            current_step = "filtering"
            log.info("Step: Cell Filtering")
            if applied_config.run_decision_engine:
                log.info("Building evidence-based QC decisions before filtering")
                decision_summary = build_qc_decisions(
                    adata,
                    tissue_type=applied_config.tissue_type or tissue_type,
                    policy=applied_config.qc_decision_policy,
                    score_layer=applied_config.qc_decision_score_layer,
                )
                adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
                    "qc_decision_summary"
                ] = decision_summary
            active_filter_config = _resolve_filter_config_for_decisions(applied_config)
            adata_filtered = filter_cells(adata, config=active_filter_config, copy=True)
            adata_filtered.uns.setdefault("sclucid", {}).setdefault("qc", {})[
                "qc_decision_filter_mode"
            ] = applied_config.qc_decision_filter_mode
            if applied_config.run_decision_engine and "qc_decision" in adata_filtered.obs:
                adata_filtered.uns.setdefault("sclucid", {}).setdefault("qc", {})[
                    "qc_decision_summary"
                ] = summarize_qc_decisions(
                    adata_filtered,
                    tissue_type=applied_config.tissue_type or tissue_type,
                    policy=applied_config.qc_decision_policy,
                )
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (skipped)")
            if applied_config.run_decision_engine:
                build_qc_decisions(
                    adata,
                    tissue_type=applied_config.tissue_type or tissue_type,
                    policy=applied_config.qc_decision_policy,
                    score_layer=applied_config.qc_decision_score_layer,
                )
            adata_filtered = adata

        filtering_summary = (
            adata_filtered.uns.get("sclucid", {}).get("qc", {}).get("filtering_results", {})
        )
        review_summary = _store_qc_trace(
            adata_filtered,
            applied_config,
            original_config,
            recommendation,
            sample_thresholds,
            filtering_summary,
            qc_warnings,
            steps_executed=successful_steps,
            adata_before_filtering=adata,
        )
        if results_path is not None:
            _export_qc_review_summary(review_summary, results_path, adata_filtered)
            benchmark_summary = review_summary.get("benchmark_summary")
            if isinstance(benchmark_summary, dict):
                export_qc_benchmark_report(benchmark_summary, results_path)

        # --- Step 3: Reporting ---
        if results_path is not None and "reporting" in steps_to_run:
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
            save_dir = recovery_save_dir or (
                str(results_path / "recovery") if results_path else "./recovery"
            )
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(
                completed_steps=successful_steps,
                failed_step=current_step,
                error_message=str(e),
            )
            manager.save(
                adata,
                checkpoint,
                applied_config if "applied_config" in locals() else runtime_config,
            )

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
        config=applied_config.to_dict(),
    )

    log.info("=" * 60)
    log.info("=== Standard QC Workflow Complete! ===")
    log.info(f"Completed steps: {successful_steps}")
    log.info("=" * 60)

    return adata_filtered


def recommend_qc_policy(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    *,
    tissue_type: str = "unknown",
    show_progress: bool = False,
) -> Dict[str, Any]:
    """Recommend a QC policy without filtering or mutating the input object.

    The returned bundle follows the canonical QC decision flow:
    profile dataset, propose candidate thresholds, score biological risk,
    choose/recommend policy, and emit a reviewer table. Filtering is not
    applied; use :func:`apply_qc_policy` or :func:`run_qc` for execution.
    """
    diagnostic_config = (
        config.model_copy(deep=True)
        if config is not None and hasattr(config, "model_copy")
        else config
    )
    if diagnostic_config is not None:
        diagnostic_config.save_dir = None
    diagnostic = run_standard_qc(
        adata_in,
        config=diagnostic_config,
        tissue_type=tissue_type,
        show_progress=show_progress,
        steps=["qc_metrics"],
    )
    qc_ns = diagnostic.uns.get("sclucid", {}).get("qc", {})
    review = qc_ns.get("review_summary", {}).get("data", {})
    applied_config = _restore_empty_config_values(qc_ns.get("applied_config", {}).get("data", {}))
    policy_flow = review.get("policy_flow", [])
    if isinstance(policy_flow, dict):
        policy_flow = [
            policy_flow[key]
            for key in sorted(policy_flow, key=lambda item: int(item) if str(item).isdigit() else str(item))
        ]
    return {
        "schema_version": "qc_policy_bundle_v1",
        "mode": "recommend_only",
        "policy_flow": policy_flow,
        "review_summary": review,
        "decision_table": review.get("decision_table", []),
        "recommended_threshold_summary": review.get("recommended_threshold_summary", {}),
        "filtering_policy_summary": review.get("qc_filtering_policy_summary", {}),
        "recommended_execution": {
            "entrypoint": "scLucid.qc.run_iterative_qc",
            "final_filter_policy": "decision_remove",
            "rationale": (
                "Reviewer-first QC should make final exclusion decisions from "
                "qc_decision == 'remove' after recording ambiguous cells as review "
                "or sensitivity_only evidence."
            ),
            "compatibility_path": (
                "scLucid.qc.apply_qc_policy remains available for legacy threshold-based "
                "execution and records qc_filtering_policy_summary for audit."
            ),
        },
        "doublet_evidence_summary": review.get("doublet_evidence_summary", {}),
        "applied_config": applied_config,
    }


def apply_qc_policy(
    adata_in: AnnData,
    policy: Optional[Dict[str, Any]] = None,
    config: Optional[QCWorkflowConfig] = None,
    *,
    tissue_type: str = "unknown",
    show_progress: bool = True,
    overwrite: bool = False,
) -> AnnData:
    """Apply a recommended or user-provided QC policy to an AnnData object."""
    if config is None and policy is not None:
        policy_config = policy.get("applied_config") or policy.get("config")
        if isinstance(policy_config, dict):
            config = QCWorkflowConfig(**_restore_empty_config_values(policy_config))
    return run_standard_qc(
        adata_in,
        config=config,
        tissue_type=tissue_type,
        show_progress=show_progress,
        overwrite=overwrite,
    )


def _as_bool_obs(adata: AnnData, column: str) -> pd.Series:
    """Return a boolean obs column or an all-False series."""
    if column not in adata.obs:
        return pd.Series(False, index=adata.obs_names)
    return adata.obs[column].fillna(False).astype(bool)


def _as_float_obs(adata: AnnData, column: str) -> pd.Series:
    """Return a numeric obs column or an all-NaN series."""
    if column not in adata.obs:
        return pd.Series(np.nan, index=adata.obs_names, dtype=float)
    return pd.to_numeric(adata.obs[column], errors="coerce")


def _sample_for_quick_review(
    adata: AnnData,
    *,
    max_cells: int,
    random_state: int,
) -> AnnData:
    """Return a deterministic subset for quick biology review."""
    if max_cells <= 0 or adata.n_obs <= max_cells:
        return adata.copy()
    rng = np.random.default_rng(random_state)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _run_quick_biology_review(
    adata: AnnData,
    *,
    sample_key: str,
    max_cells: int = 2000,
    n_top_genes: int = 1000,
    n_pcs: int = 20,
    n_neighbors: int = 10,
    resolution: float = 0.5,
    random_state: int = 0,
) -> Dict[str, Any]:
    """Run a temporary quick embedding/cluster review for iterative QC."""
    if adata.n_obs < 20 or adata.n_vars < 20:
        return {
            "schema_version": "quick_biology_review_v1",
            "status": "skipped_too_few_cells",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "review_required": False,
            "review_findings": [],
            "cluster_qc_table": [],
        }

    try:
        tmp = _sample_for_quick_review(adata, max_cells=max_cells, random_state=random_state)
        tmp.X = tmp.X.copy()

        # Avoid double-normalising if the input matrix is already log1p-normalised.
        # A raw count matrix typically has a max value well above 100; log1p(CPM)
        # or log1p-normalised data is usually below 10-20.
        x_max = float(np.asarray(tmp.X).max()) if not scipy.sparse.issparse(tmp.X) else float(tmp.X.max())
        if x_max > 100:
            sc.pp.normalize_total(tmp, target_sum=1e4)
            sc.pp.log1p(tmp)
        else:
            log.debug(
                "Quick biology review: input max value %.2f suggests already-normalised data; "
                "skipping normalize_total/log1p.",
                x_max,
            )
        if tmp.n_vars > 50:
            sc.pp.highly_variable_genes(
                tmp,
                n_top_genes=min(n_top_genes, tmp.n_vars),
                flavor="seurat",
            )
            if "highly_variable" in tmp.var and bool(tmp.var["highly_variable"].any()):
                tmp = tmp[:, tmp.var["highly_variable"].to_numpy()].copy()
        sc.pp.scale(tmp, max_value=10)
        effective_pcs = min(n_pcs, tmp.n_obs - 1, tmp.n_vars - 1)
        if effective_pcs < 2:
            raise ValueError("Too few cells/genes for PCA in quick biology review.")
        sc.tl.pca(tmp, n_comps=effective_pcs, random_state=random_state)
        sc.pp.neighbors(
            tmp,
            n_neighbors=min(n_neighbors, max(2, tmp.n_obs - 1)),
            n_pcs=effective_pcs,
            random_state=random_state,
        )
        sc.tl.umap(tmp, random_state=random_state)
        sc.tl.leiden(
            tmp,
            resolution=resolution,
            key_added="quick_qc_leiden",
            random_state=random_state,
        )

        obs = tmp.obs
        cluster_rows: list[dict[str, Any]] = []
        review_findings: list[dict[str, Any]] = []
        for cluster, group in obs.groupby("quick_qc_leiden", observed=True):
            n_cells = int(group.shape[0])
            qc_remove_frac = float(_as_bool_obs(tmp, "qc_remove").loc[group.index].mean())
            review_frac = float(_as_bool_obs(tmp, "qc_review_required").loc[group.index].mean())
            doublet_frac = float(_as_bool_obs(tmp, "predicted_doublet").loc[group.index].mean())
            stress_frac = float(_as_bool_obs(tmp, "stress_high").loc[group.index].mean())
            ambient_frac = float(_as_bool_obs(tmp, "ambient_risk").loc[group.index].mean())
            mean_mt = _as_float_obs(tmp, "pct_counts_mt").loc[group.index].mean()
            dominant_sample_fraction = None
            dominant_sample = None
            if sample_key in group:
                sample_counts = group[sample_key].astype(str).value_counts(normalize=True)
                if not sample_counts.empty:
                    dominant_sample = str(sample_counts.index[0])
                    dominant_sample_fraction = float(sample_counts.iloc[0])
            flags: list[str] = []
            if qc_remove_frac >= 0.30 or review_frac >= 0.50:
                flags.append("qc_enriched_cluster")
            if stress_frac >= 0.50:
                flags.append("stress_enriched_cluster")
            if doublet_frac >= 0.20:
                flags.append("doublet_enriched_cluster")
            if ambient_frac >= 0.30:
                flags.append("ambient_enriched_cluster")
            if dominant_sample_fraction is not None and dominant_sample_fraction >= 0.80 and n_cells >= 10:
                flags.append("sample_dominated_cluster")

            row = {
                "cluster": str(cluster),
                "n_cells": n_cells,
                "qc_remove_fraction": qc_remove_frac,
                "qc_review_required_fraction": review_frac,
                "doublet_fraction": doublet_frac,
                "stress_high_fraction": stress_frac,
                "ambient_risk_fraction": ambient_frac,
                "mean_pct_counts_mt": None if pd.isna(mean_mt) else float(mean_mt),
                "dominant_sample": dominant_sample,
                "dominant_sample_fraction": dominant_sample_fraction,
                "review_flags": flags,
            }
            cluster_rows.append(row)
            for flag in flags:
                review_findings.append(
                    {
                        "scope": "cluster",
                        "cluster": str(cluster),
                        "finding": flag,
                        "review_required": True,
                    }
                )

        stress_mask = _as_bool_obs(tmp, "stress_high")
        if sample_key in obs and int(stress_mask.sum()) >= 10:
            stress_samples = obs.loc[stress_mask, sample_key].astype(str).value_counts(normalize=True)
            if not stress_samples.empty and float(stress_samples.iloc[0]) >= 0.60:
                review_findings.append(
                    {
                        "scope": "sample",
                        "finding": "stress_high_sample_bias",
                        "sample": str(stress_samples.index[0]),
                        "fraction_of_stress_high_cells": float(stress_samples.iloc[0]),
                        "review_required": True,
                    }
                )

        doublet_mask = _as_bool_obs(tmp, "predicted_doublet")
        if int(doublet_mask.sum()) >= 5:
            doublet_clusters = obs.loc[doublet_mask, "quick_qc_leiden"].astype(str).nunique()
            review_findings.append(
                {
                    "scope": "embedding",
                    "finding": "doublet_boundary_review",
                    "doublet_positive_clusters": int(doublet_clusters),
                    "review_required": bool(doublet_clusters > 1),
                }
            )

        ambient_mask = _as_bool_obs(tmp, "ambient_risk")
        if int(ambient_mask.sum()) >= 5:
            ambient_clusters = obs.loc[ambient_mask, "quick_qc_leiden"].astype(str).nunique()
            if ambient_clusters >= max(2, obs["quick_qc_leiden"].nunique() // 2):
                review_findings.append(
                    {
                        "scope": "embedding",
                        "finding": "ambient_marker_widespread_leakage",
                        "ambient_positive_clusters": int(ambient_clusters),
                        "review_required": True,
                    }
                )

        return {
            "schema_version": "quick_biology_review_v1",
            "status": "complete",
            "n_cells_reviewed": int(tmp.n_obs),
            "n_genes_reviewed": int(tmp.n_vars),
            "cluster_key": "quick_qc_leiden",
            "n_clusters": int(obs["quick_qc_leiden"].nunique()),
            "umap_computed": bool("X_umap" in tmp.obsm),
            "review_required": any(item.get("review_required") for item in review_findings),
            "review_findings": review_findings,
            "cluster_qc_table": cluster_rows,
            "note": (
                "Quick review was computed on a temporary normalized/log1p subset and "
                "does not modify formal preprocessing layers."
            ),
        }
    except Exception as exc:
        return {
            "schema_version": "quick_biology_review_v1",
            "status": "failed",
            "error": str(exc),
            "review_required": True,
            "review_findings": [
                {
                    "scope": "workflow",
                    "finding": "quick_biology_review_failed",
                    "review_required": True,
                }
            ],
            "cluster_qc_table": [],
        }


def _refine_qc_decisions_from_review(
    adata: AnnData,
    quick_biology_review: Dict[str, Any],
    sample_key: str,
) -> Dict[str, Any]:
    """Use quick-biology-review findings to refine qc_decision labels.

    The quick review runs on a temporary embedding and may flag:

    - ``qc_enriched_cluster``: clusters with >30% qc_remove or >50% review_required
    - ``stress_enriched_cluster`` / ``stress_high_sample_bias``: stress signals
    - ``doublet_enriched_cluster`` / ``doublet_boundary_review``: doublet signals
    - ``ambient_enriched_cluster`` / ``ambient_marker_widespread_leakage``: ambient leakage
    - ``sample_dominated_cluster``: technical batch/sample segregation

    This function does **not** delete cells. It upgrades the decision label of
    affected cells to ``review`` or ``sensitivity_only`` when the quick review
    suggests the original decision may be missing biology or over-filtering.
    """
    if quick_biology_review.get("status") != "complete":
        return {"refined": False, "reason": "review_not_complete", "changes": {}}

    findings = quick_biology_review.get("review_findings", [])
    if not findings:
        return {"refined": False, "reason": "no_findings", "changes": {}}

    original_decisions = adata.obs["qc_decision"].copy()
    changes: Dict[str, int] = {}

    # Identify affected samples flagged by review findings.
    affected_samples: set[str] = set()
    for finding in findings:
        if finding.get("finding") == "stress_high_sample_bias":
            affected_samples.add(str(finding.get("sample", "")))

    if affected_samples and sample_key in adata.obs:
        sample_mask = adata.obs[sample_key].astype(str).isin(affected_samples)
        # Stress-high cells in affected samples: keep under sensitivity review.
        stress_mask = _as_bool_obs(adata, "stress_high") & sample_mask
        upgrade_mask = stress_mask & adata.obs["qc_decision"].isin(["keep"])
        if upgrade_mask.any():
            adata.obs.loc[upgrade_mask, "qc_decision"] = "sensitivity_only"
            adata.obs.loc[upgrade_mask, "qc_review_required"] = True
            changes["stress_sample_to_sensitivity"] = int(upgrade_mask.sum())

    # Clusters enriched for QC issues: ensure cells are at least under review.
    cluster_table = quick_biology_review.get("cluster_qc_table", [])
    review_clusters: set[str] = set()
    for row in cluster_table:
        flags = row.get("review_flags", [])
        if any(f in flags for f in ("qc_enriched_cluster", "ambient_enriched_cluster")):
            review_clusters.add(str(row.get("cluster", "")))

    if review_clusters:
        # Map quick-review clusters back to full adata via the sampled subset.
        # The quick review stores its cluster assignments only in the temporary
        # object, so we approximate by upgrading cells that share the flagged
        # QC evidence across the whole dataset.
        for flag in ("qc_high_mt", "ambient_risk", "qc_low_complexity"):
            if flag in adata.obs.columns:
                upgrade_mask = (
                    adata.obs[flag].astype(bool)
                    & adata.obs["qc_decision"].isin(["keep"])
                )
                if upgrade_mask.any():
                    adata.obs.loc[upgrade_mask, "qc_decision"] = "review"
                    adata.obs.loc[upgrade_mask, "qc_review_required"] = True
                    changes.setdefault(f"{flag}_to_review", 0)
                    changes[f"{flag}_to_review"] += int(upgrade_mask.sum())

    # Doublet-enriched clusters: doublet-like cells already flagged by the
    # algorithm should remain under review if they were about to be kept.
    doublet_upgrade = (
        _as_bool_obs(adata, "predicted_doublet")
        & adata.obs["qc_decision"].isin(["keep"])
    )
    if doublet_upgrade.any():
        adata.obs.loc[doublet_upgrade, "qc_decision"] = "review"
        adata.obs.loc[doublet_upgrade, "qc_review_required"] = True
        changes["doublet_to_review"] = int(doublet_upgrade.sum())

    # Recompute qc_remove to respect refined decisions.
    adata.obs["qc_remove"] = adata.obs["qc_decision"].eq("remove").to_numpy(bool)

    n_changed = int((adata.obs["qc_decision"] != original_decisions).sum())
    return {
        "refined": n_changed > 0,
        "n_changed": n_changed,
        "changes": changes,
    }


def _build_iterative_qc_summary(
    adata: AnnData,
    *,
    final_filter_policy: str,
    run_quick_review: bool,
    quick_biology_review: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build the phase-level audit summary for iterative QC."""
    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    decision_summary = qc_ns.get("qc_decision_summary", {})
    filtering_summary = qc_ns.get("filtering_results", {})
    phases = [
        {
            "phase": "lenient_cell_screen",
            "status": "complete",
            "outputs": [
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "outlier_*",
            ],
            "note": "Extreme barcode/cell quality evidence was marked before final filtering.",
        },
        {
            "phase": "doublet_contamination_stress_evidence",
            "status": "complete",
            "outputs": [
                "predicted_doublet",
                "qc_decision",
                "qc_reason",
                "qc_confidence",
                "stress_score",
                "hemoglobin_score",
                "platelet_score",
            ],
            "note": "Evidence columns are review signals; ambiguous biology is not automatically deleted.",
        },
        {
            "phase": "quick_biology_review",
            "status": "not_run"
            if not run_quick_review
            else (quick_biology_review or {}).get("status", "not_run"),
            "outputs": ["quick_biology_review"] if run_quick_review else [],
            "note": (
                "Quick embedding/cluster review was not requested."
                if not run_quick_review
                else "Temporary embedding review summarized QC/stress/doublet/ambient patterns."
            ),
        },
        {
            "phase": "final_qc_decision",
            "status": "complete",
            "outputs": ["qc_decision_summary", "filtering_results"],
            "note": f"Final filter policy: {final_filter_policy}.",
        },
    ]
    return {
        "schema_version": "iterative_qc_summary_v1",
        "final_filter_policy": final_filter_policy,
        "run_quick_review": bool(run_quick_review),
        "phases": phases,
        "quick_biology_review": quick_biology_review or {},
        "qc_decision_summary": decision_summary,
        "filtering_summary": filtering_summary,
        "recommended_next_step": (
            "Run formal preprocessing on retained cells; inspect review/sensitivity_only cells "
            "before irreversible exclusion in tumor or fragile-cell contexts."
        ),
    }


def run_iterative_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    *,
    tissue_type: str = "unknown",
    final_filter_policy: Literal["decision_remove", "legacy", "none"] = "decision_remove",
    run_quick_review: bool = False,
    quick_review_max_cells: int = 2000,
    quick_review_n_top_genes: int = 1000,
    quick_review_n_pcs: int = 20,
    quick_review_n_neighbors: int = 10,
    quick_review_resolution: float = 0.5,
    show_progress: bool = True,
    overwrite: bool = False,
) -> AnnData:
    """Run the reviewer-first iterative QC entrypoint.

    This first implementation formalizes the phase contract and routes final
    filtering through the unified ``qc_decision`` schema when requested. It is
    deliberately conservative: ambiguous high-MT/stress/doublet-like cells are
    marked for review unless the joint decision engine labels them ``remove``.
    """
    active_config = (
        config.model_copy(deep=True)
        if config is not None and hasattr(config, "model_copy")
        else QCWorkflowConfig()
    )
    active_config.run_decision_engine = True
    if final_filter_policy == "decision_remove":
        active_config.qc_decision_filter_mode = "replace"
        steps = None
    elif final_filter_policy == "legacy":
        active_config.qc_decision_filter_mode = "off"
        steps = None
    elif final_filter_policy == "none":
        active_config.qc_decision_filter_mode = "off"
        steps = ["qc_metrics"]
    else:
        raise ValueError(
            "final_filter_policy must be one of: 'decision_remove', 'legacy', 'none'"
        )

    result = run_standard_qc(
        adata_in,
        config=active_config,
        tissue_type=tissue_type,
        show_progress=show_progress,
        overwrite=overwrite,
        steps=steps,
    )
    quick_biology_review = None
    refinement = {"refined": False}
    if run_quick_review:
        quick_biology_review = _run_quick_biology_review(
            result,
            sample_key=active_config.sample_key,
            max_cells=quick_review_max_cells,
            n_top_genes=quick_review_n_top_genes,
            n_pcs=quick_review_n_pcs,
            n_neighbors=quick_review_n_neighbors,
            resolution=quick_review_resolution,
            random_state=active_config.doublet_config.random_state,
        )
        refinement = _refine_qc_decisions_from_review(
            result,
            quick_biology_review,
            sample_key=active_config.sample_key,
        )
        # If decisions were refined and final filtering uses qc_decision,
        # re-apply filtering so the refined labels take effect.
        if refinement.get("refined") and final_filter_policy == "decision_remove":
            from .decisions import summarize_qc_decisions

            summarize_qc_decisions(result, tissue_type=tissue_type, policy=active_config.qc_decision_policy)
            keep_mask = ~result.obs["qc_remove"].fillna(False).to_numpy(bool)
            result = result[keep_mask, :].copy()
            log.info(
                "Iterative QC: re-applied filtering after quick review; %d cells retained.",
                result.n_obs,
            )
    summary = _build_iterative_qc_summary(
        result,
        final_filter_policy=final_filter_policy,
        run_quick_review=run_quick_review,
        quick_biology_review=quick_biology_review,
    )
    summary["decision_refinement"] = refinement
    result.uns.setdefault("sclucid", {}).setdefault("qc", {})[
        "iterative_qc_summary"
    ] = summary
    return result


def run_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    **kwargs: Any,
) -> AnnData:
    """Canonical user-facing QC entrypoint.

    This is a thin alias for :func:`run_standard_qc`; it exists to make the
    public API read as ``recommend_qc_policy`` -> ``apply_qc_policy`` or
    ``run_qc`` for the one-step workflow.
    """
    return run_standard_qc(adata_in, config=config, **kwargs)


def run_advanced_qc(
    adata_in: AnnData,
    config: QCWorkflowConfig,
    overwrite: bool = False,
    *,
    tissue_type: str = "unknown",
    show_progress: bool = True,
    # Step control
    steps: Optional[List[str]] = None,
    skip_steps: Optional[List[str]] = None,
    # Error recovery
    error_recovery: bool = False,
    recovery_save_dir: Optional[str] = None,
    on_error: str = "raise",
    # Resume
    resume_from: Optional[str] = None,
) -> AnnData:
    """
    Run an advanced, fully configurable single-cell RNA-seq QC workflow.

    .. deprecated::
        ``run_advanced_qc`` is a thin compatibility wrapper around
        :func:`run_standard_qc`. Pass a fully populated ``QCWorkflowConfig`` to
        ``run_standard_qc`` instead; this alias will be removed in a future release.

    This workflow is entirely controlled by the provided QCWorkflowConfig object,
    allowing fine-grained control over every step.

    Reviewer-facing outputs:
        - ``adata.uns['sclucid']['qc']['review_summary']`` contains a structured
          digest of recommendations, applied thresholds, user overrides,
          sample-level thresholds, tumor-aware flags, and filtering results.
        - When ``save_dir`` is set, ``qc_review_summary.json`` and
          ``qc_review_summary.md`` sidecars are written alongside the report.

    New features in v0.4:
    - Step control via ``steps`` or ``skip_steps`` (consistent with preprocess/analysis)
    - Error recovery with partial result saving
    - Resume from checkpoint

    Args:
        adata_in: Input AnnData object.
        config: A fully populated QCWorkflowConfig object.
        overwrite: If True, overwrite existing results directory.
        show_progress: If True, show progress bars for multi-sample processing.
        steps: Specific steps to run. See ``QC_WORKFLOW_STEPS`` for valid names.
        skip_steps: Steps to skip (alternative to specifying ``steps``).
        error_recovery: If True, enable error recovery mode.
        recovery_save_dir: Directory to save partial results on error.
        on_error: How to handle errors: "raise", "skip", or "save".
        resume_from: Path to checkpoint directory to resume from.

    Returns:
        Filtered AnnData object after QC.
    """
    warnings.warn(
        "run_advanced_qc is a compatibility wrapper; use run_standard_qc with "
        "QCWorkflowConfig for the maintained QC workflow.",
        FutureWarning,
        stacklevel=2,
    )
    return run_standard_qc(
        adata_in,
        config=config,
        overwrite=overwrite,
        tissue_type=tissue_type,
        show_progress=show_progress,
        steps=steps,
        skip_steps=skip_steps,
        error_recovery=error_recovery,
        recovery_save_dir=recovery_save_dir,
        on_error=on_error,
        resume_from=resume_from,
    )


__all__ = [
    "run_qc",
    "run_iterative_qc",
    "recommend_qc_policy",
    "apply_qc_policy",
    "run_standard_qc",
    "QCWorkflowError",
    "QC_WORKFLOW_STEPS",
]
