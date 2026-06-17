"""Reusable QC threshold decision layer.

This module hosts the QC threshold decision orchestrator, which combines
distribution-based threshold suggestions, optional intelligent-QC thresholds,
optional adaptive marking, unified low-quality marking, and optional
filtering/audit. It deliberately does **not** compute QC metrics or detect
doublets — it operates on metrics already present in ``adata.obs``.

Kept separate from ``suggestions.py`` (pure threshold suggestion) so that the
orchestration logic does not bloat the suggestion helpers.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, Dict, List, Literal, Optional

from anndata import AnnData

from .suggestions import resolve_qc_thresholds, suggest_qc_thresholds

log = logging.getLogger(__name__)

__all__ = [
    "run_qc_threshold_decision",
    "run_qc_decision_workflow",
]


def run_qc_threshold_decision(
    adata: AnnData,
    *,
    intelligent_thresholds: Optional[Dict[str, Any]] = None,
    manual_thresholds: Optional[Dict[str, Any]] = None,
    threshold_method: Literal["mad", "iqr", "percentile"] = "mad",
    threshold_policy: Literal[
        "intelligent_then_mad", "mad_then_intelligent", "manual_override"
    ] = "intelligent_then_mad",
    plot_distributions: Optional[bool] = None,
    sample_key: Optional[str] = None,
    group_key: Optional[str] = None,
    adaptive: bool = False,
    adaptive_batch_key: Optional[str] = None,
    adaptive_metrics: Optional[List[str]] = None,
    criteria_to_filter: Optional[List[str]] = None,
    filter_cells_result: bool = False,
    save_dir: Optional[str] = None,
    **threshold_kwargs,
) -> Dict[str, Any]:
    """Run the reusable QC threshold decision layer.

    The function combines distribution-based threshold suggestions, optional
    intelligent QC thresholds, optional adaptive marking, unified marking and
    optional filtering/audit. It does not calculate QC metrics or doublets.
    """
    from ..config import FilterConfig, MarkingConfig
    from .core import audit_filtering, filter_cells, mark_low_quality_cell

    if sample_key is None:
        for candidate in ("sampleID", "sample_id", "sample", "batch"):
            if candidate in adata.obs.columns:
                sample_key = candidate
                break
    if sample_key is None:
        adata.obs["_sclucid_sample"] = "sample_0"
        sample_key = "_sclucid_sample"

    threshold_table, suggested = suggest_qc_thresholds(
        adata,
        method=threshold_method,
        plot_distributions=bool(save_dir) if plot_distributions is None else plot_distributions,
        save_dir=save_dir,
        **threshold_kwargs,
    )
    mad_thresholds = suggested.to_dict()
    resolved = resolve_qc_thresholds(
        intelligent=intelligent_thresholds,
        mad=mad_thresholds,
        manual=manual_thresholds,
        policy=threshold_policy,
    )

    if adaptive:
        from .core import mark_low_quality_cells_adaptive

        adata = mark_low_quality_cells_adaptive(
            adata,
            batch_key=adaptive_batch_key or sample_key,
            metrics=adaptive_metrics,
        )

    marking_cfg = MarkingConfig(thresholds=resolved, save_dir=save_dir)
    adata = mark_low_quality_cell(adata, sample_key=sample_key, config=marking_cfg)

    filtered = None
    filtering_audit = None
    if filter_cells_result:
        filter_cfg = FilterConfig(criteria_to_filter=criteria_to_filter or [])
        filtered = filter_cells(adata, config=filter_cfg, copy=True)
        filtering_audit = audit_filtering(
            adata,
            filtered,
            sample_key=sample_key,
            group_key=group_key,
        )

    decision = {
        "adata": adata,
        "filtered": filtered,
        "threshold_table": threshold_table,
        "suggested_thresholds": suggested,
        "resolved_thresholds": resolved,
        "filtering_audit": filtering_audit,
        "criteria_to_filter": criteria_to_filter or [],
        "threshold_policy": threshold_policy,
    }
    decision_record = {
        "resolved_thresholds": resolved.to_dict(),
        "criteria_to_filter": criteria_to_filter or [],
        "adaptive": bool(adaptive),
        "filter_cells_result": bool(filter_cells_result),
    }
    qc_uns = adata.uns.setdefault("sclucid", {}).setdefault("qc", {})
    qc_uns["qc_threshold_decision"] = decision_record
    qc_uns["qc_decision_workflow"] = decision_record
    return decision


def run_qc_decision_workflow(*args, **kwargs) -> Dict[str, Any]:
    """Compatibility alias for :func:`run_qc_threshold_decision`."""
    warnings.warn(
        "run_qc_decision_workflow is transitional; use run_qc_threshold_decision "
        "for threshold decision/marking helpers.",
        FutureWarning,
        stacklevel=2,
    )
    return run_qc_threshold_decision(*args, **kwargs)
