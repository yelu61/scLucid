"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from ..runtime import effective_n_jobs
from ..utils import (
    PartialResultManager,
    WorkflowCheckpoint,
    WorkflowError,
    save_workflow_result,
)
from ..utils.context import is_tumor_context
from .artifacts import (
    record_qc_decision_artifact,
    record_threshold_decision,
    record_threshold_recommendation,
)
from .config import FilterConfig, QCWorkflowConfig
from .doublet import predict_doublets
from .filtering import filter_cells
from .metrics import calculate_qc_metric
from .policy.adaptive_threshold import AdaptiveThresholdCalculator
from .policy.decisions import build_qc_decisions, summarize_qc_decisions, score_qc_gene_panels
from .policy.marking import mark_low_quality_cell
from .policy.thresholds import _resolve_qc_thresholds
from .reporting import generate_qc_report
from .workflow_artifacts import (
    _export_qc_review_summary,
    _store_qc_trace,
)
from .workflow_runtime import (
    QC_WORKFLOW_STEPS,
    _add_tumor_aware_flags,
    _ensure_sample_key,
    _merge_sample_results,
    _prepare_runtime_qc_config,
    _process_sample_doublet,
    _process_sample_qc,
    _resolve_qc_steps,
    _restore_empty_config_values,
    _safe_parallel_process,
    _setup_workflow,
)
from .workflow_review import (
    _build_iterative_qc_summary,
    _refine_qc_decisions_from_review,
    _run_quick_biology_review,
)

log = logging.getLogger(__name__)


def _is_tumor_aware(tissue_type: Optional[str]) -> bool:
    """Return whether the workflow should apply tumor-aware QC safeguards."""
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


def _threshold_payload(value: Any) -> Dict[str, Any]:
    """Convert threshold-like objects to plain dictionaries for QC artifacts."""
    if value is None:
        return {}
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if hasattr(value, "model_dump"):
        return value.model_dump()
    if isinstance(value, dict):
        return dict(value)
    return {}


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
        applied.marking_config.thresholds = _resolve_qc_thresholds(
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
            from .policy.intelligent_qc import recommend_intelligent_qc

            recommendation = recommend_intelligent_qc(adata, tissue_type=active_tissue_type)
            config, original_config = _apply_qc_recommendations(config, recommendation)
            record_threshold_recommendation(
                adata,
                source="intelligent_qc",
                payload=recommendation.to_dict() if hasattr(recommendation, "to_dict") else {},
            )
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
            record_threshold_recommendation(
                adata,
                source=f"{config.threshold_mode}_sample_thresholds",
                payload=sample_thresholds,
            )
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
                    tissue_type=config.tissue_type,
                    sample_context_key=config.sample_context_key,
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
                    tissue_type=config.tissue_type,
                    sample_context_key=config.sample_context_key,
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
    record_threshold_decision(
        adata,
        resolved_thresholds=config.marking_config.thresholds.to_dict(),
        policy=(
            "intelligent_recommendation_with_user_overrides"
            if recommendation is not None and config.use_recommendations
            else "configured_thresholds"
        ),
        sources={
            "intelligent_qc": recommendation.to_dict()
            if recommendation is not None and hasattr(recommendation, "to_dict")
            else {},
            "manual_config": original_config.marking_config.thresholds.to_dict(),
            "applied_config": config.marking_config.thresholds.to_dict(),
        },
        sample_key=config.sample_key,
        sample_thresholds=sample_thresholds,
        notes=warnings_list,
    )
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
                filtered_decision_summary = summarize_qc_decisions(
                    adata_filtered,
                    tissue_type=applied_config.tissue_type or tissue_type,
                    policy=applied_config.qc_decision_policy,
                )
                adata_filtered.uns.setdefault("sclucid", {}).setdefault("qc", {})[
                    "qc_decision_summary"
                ] = filtered_decision_summary
                record_qc_decision_artifact(
                    adata_filtered,
                    summary=filtered_decision_summary,
                    evidence_columns=[
                        "qc_low_counts",
                        "qc_low_genes",
                        "qc_high_mt",
                        "qc_low_complexity",
                        "qc_high_hb",
                        "platelet_contamination",
                        "hemoglobin_contamination",
                        "ambient_risk",
                        "stress_high",
                        "apoptosis_high",
                        "predicted_doublet",
                    ],
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
            from .policy.decisions import summarize_qc_decisions

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
