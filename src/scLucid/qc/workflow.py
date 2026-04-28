"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import json
import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TypeVar

import numpy as np
from anndata import AnnData

from ..runtime import effective_n_jobs
from ..utils import (
    PartialResultManager,
    UnsKeys,
    WorkflowCheckpoint,
    WorkflowError,
    get_progress_bar,
    normalize_review_summary,
    save_result,
    save_workflow_result,
    validate_review_summary_schema,
)
from .benchmark import evaluate_qc_benchmark, export_qc_benchmark_report
from .config import QCWorkflowConfig
from .doublet import predict_doublets
from .filtering import (
    AdaptiveThresholdCalculator,
    filter_cells,
    generate_qc_report,
    mark_low_quality_cell,
)
from .metrics import calculate_qc_metric
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
    }

    # --- Warnings ---
    summary["warnings"] = warnings

    return summary


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

    # Build and store the review-facing summary
    base_review_summary = _build_qc_review_summary(
        config,
        original_config,
        recommendation,
        sample_thresholds,
        filtering_summary,
        warnings,
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
    for key in ["cell_type", "celltype", "cell_type_major", "annotation", "cell_annotation"]:
        if key in adata.obs:
            return key
    return None


def _export_qc_review_summary(
    review_summary: Dict[str, Any],
    save_dir: Path,
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
            "## Decision Table",
            "",
            "| Parameter | Applied | Source | Filter Enabled | Method | Confidence |",
            "|-----------|---------|--------|----------------|--------|------------|",
        ]
    )
    for row in review_summary.get("decision_table", []):
        md_lines.append(
            "| {parameter} | {applied} | {source} | {enabled} | {method} | {confidence} |".format(
                parameter=row.get("parameter"),
                applied=row.get("applied"),
                source=row.get("source"),
                enabled=row.get("is_filtering_enabled"),
                method=row.get("recommendation_method"),
                confidence=row.get("confidence"),
            )
        )
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
    if health.get("issues"):
        md_lines.append("- **Issues**:")
        for issue in health.get("issues", []):
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
    for item in downstream.get("recommendations", []):
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

    if review_summary.get("warnings"):
        md_lines.extend(
            [
                "## Warnings",
                "",
            ]
        )
        for w in review_summary["warnings"]:
            md_lines.append(f"- {w}")
        md_lines.append("")

    md_path = save_dir / "qc_review_summary.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    log.info(f"QC review summary exported to {json_path} and {md_path}")


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
            log.info(
                f"Computed {config.threshold_mode} per-sample thresholds for {len(sample_thresholds)} samples."
            )

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

        if config.use_parallel and active_n_jobs != 1 and len(samples) > 1:
            # Parallel doublet detection with error handling
            results = _safe_parallel_process(
                process_func=lambda data, cfg, name: _process_sample_doublet(
                    data, cfg, config.doublet_config.save_dir, name
                ),
                samples=list(samples),
                sample_data_func=lambda s: adata[adata.obs[config.sample_key] == s].copy(),
                config=config,
                n_jobs=active_n_jobs,
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
            adata_filtered = filter_cells(adata, config=applied_config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (skipped)")
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
            _export_qc_review_summary(review_summary, results_path)
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
    log.info("=== Starting Advanced QC Workflow ===")
    log.info("=" * 60)
    if results_path is not None:
        log.info(f"Results will be saved in: {results_path}")
    else:
        log.info("Running without file output")
    log.info(f"Error recovery: {error_recovery}")

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
            adata_filtered = filter_cells(adata, config=applied_config.filter_config, copy=True)
            successful_steps.append("filtering")
        else:
            log.info("Step: Filtering (skipped)")
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
            _export_qc_review_summary(review_summary, results_path)
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
        error_msg = f"Advanced QC workflow failed at step '{current_step}': {str(e)}"
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
                log.warning(
                    f"To resume, use: run_advanced_qc(adata, config, resume_from='{save_dir}')"
                )
                return adata

        raise WorkflowError(error_msg, step_name=current_step or "unknown", original_error=e)

    # Save workflow result using standardized storage
    save_workflow_result(
        adata_filtered,
        module="qc",
        workflow_name="advanced",
        steps=successful_steps,
        config=applied_config.to_dict(),
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
    "QC_WORKFLOW_STEPS",
]
