"""QC threshold recommendation and decision chain."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from .adaptive_threshold import (
    THRESHOLD_RESULT_SCHEMA_VERSION,
    build_threshold_result,
    compute_mad_bounds,
    infer_qc_metric_type,
)
from ..artifacts import record_threshold_decision, record_threshold_recommendation
from ..config import QCThresholds

log = logging.getLogger(__name__)

THRESHOLD_RECOMMENDATION_BUNDLE_SCHEMA_VERSION = "qc_threshold_recommendation_bundle_v1"

__all__ = [
    "recommend_qc_thresholds",
    "decide_qc_thresholds",
    "apply_qc_threshold_decision",
    "run_qc_threshold_decision",
]


# --- Helper Functions ---


def _threshold_metric_map(top_gene_cols: List[str]) -> Dict[str, Any]:
    """Map obs metric names to QCThresholds keys."""
    metric_map: Dict[str, Any] = {
        "n_genes_by_counts": ("min_genes", "max_genes"),
        "total_counts": ("min_counts", "max_counts"),
        "pct_counts_mt": "pc_mt",
        "pct_counts_hb": "pc_hb",
    }
    for col in top_gene_cols:
        metric_map[col] = f"pc_{col.split('pct_counts_in_')[-1]}"
    return metric_map


def _threshold_metrics(adata: AnnData) -> Dict[str, str]:
    """Return QC metrics available for threshold recommendation."""
    metrics = {
        "n_genes_by_counts": "Gene counts per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial percentage",
    }
    if "pct_counts_hb" in adata.obs.columns:
        metrics["pct_counts_hb"] = "Hemoglobin percentage"

    top_gene_cols = [
        col for col in adata.obs.columns if re.match(r"pct_counts_in_top_\d+_genes", col)
    ]
    for col in top_gene_cols:
        metrics[col] = (
            col.replace("_", " ").replace("pct counts in ", "").replace(" genes", "").title()
        )
    return metrics


def _is_count_metric(metric: str) -> bool:
    """Return True for QC metrics represented as count-like lower/upper bounds."""
    return metric in {"n_genes_by_counts", "total_counts"}


def _is_review_only_metric(metric: str) -> bool:
    """Return True when thresholds should not become default hard filters."""
    metric_type = infer_qc_metric_type(metric)
    return metric_type == "review_evidence" or metric.startswith("pct_counts_in_top_")


def _candidate_record(
    *,
    metric: str,
    level_name: str,
    direction: str,
    threshold: float,
    values: np.ndarray,
    threshold_key: Optional[str],
    method: str,
    method_label: str,
    hard_filter_candidate: bool,
    fallback_reason: Optional[str] = None,
) -> Dict[str, Any]:
    """Build one structured threshold candidate record."""
    result = build_threshold_result(
        metric_name=metric,
        direction=direction,
        method=method,
        threshold=threshold,
        values=values,
        metric_type=infer_qc_metric_type(metric),
        fallback_reason=fallback_reason,
    )
    result.update(
        {
            "threshold_result_schema_version": THRESHOLD_RESULT_SCHEMA_VERSION,
            "level": level_name,
            "method_label": method_label,
            "threshold_key": threshold_key,
            "hard_filter_candidate": bool(hard_filter_candidate),
            "review_only": bool(not hard_filter_candidate),
            "converted_to_qc_threshold": bool(threshold_key and hard_filter_candidate),
        }
    )
    if result["review_only"]:
        result["review_note"] = (
            result["review_note"]
            + " This candidate is kept for review and is not promoted to the default "
            "hard-filter threshold."
        )
    return result


def _add_candidate(
    *,
    candidates: Dict[str, List[Dict[str, Any]]],
    table: Dict[str, Dict[str, Any]],
    metric: str,
    level_name: str,
    direction: str,
    threshold: float,
    values: np.ndarray,
    threshold_key: Optional[str],
    method: str,
    method_label: str,
    hard_filter_candidate: bool,
) -> None:
    """Append a candidate and update the legacy threshold table when appropriate."""
    if threshold_key and hard_filter_candidate:
        table.setdefault(level_name, {})[threshold_key] = threshold
    candidates.setdefault(metric, []).append(
        _candidate_record(
            metric=metric,
            level_name=level_name,
            direction=direction,
            threshold=threshold,
            values=values,
            threshold_key=threshold_key,
            method=method,
            method_label=method_label,
            hard_filter_candidate=hard_filter_candidate,
        )
    )


def _build_threshold_dataframe(all_suggestions: Dict[str, Dict[str, Any]]) -> pd.DataFrame:
    """Create the legacy threshold table from structured suggestions."""
    suggested_thresholds_df = pd.DataFrame.from_dict(all_suggestions, orient="index")
    if suggested_thresholds_df.empty:
        return suggested_thresholds_df

    cols_order = [
        "min_genes",
        "max_genes",
        "min_counts",
        "max_counts",
        "pc_mt",
        "pc_hb",
    ]
    top_gene_cols_sorted = sorted(
        [c for c in suggested_thresholds_df.columns if c.startswith("pc_top_")]
    )
    final_cols = [
        c for c in cols_order if c in suggested_thresholds_df.columns
    ] + top_gene_cols_sorted
    return suggested_thresholds_df[final_cols]


def _thresholds_from_dataframe(suggested_thresholds_df: pd.DataFrame) -> QCThresholds:
    """Create a legacy QCThresholds object from the first recommendation row."""
    default_thresholds_obj = QCThresholds()
    if suggested_thresholds_df.empty:
        return default_thresholds_obj

    default_series = suggested_thresholds_df.iloc[0]
    pc_top_genes_dict = {k: v for k, v in default_series.items() if k.startswith("pc_top_")}
    final_kwargs = {
        k: v for k, v in default_series.items() if not k.startswith("pc_top_") and pd.notna(v)
    }
    final_kwargs["pc_top_genes"] = pc_top_genes_dict
    return QCThresholds(**final_kwargs)


def _plot_threshold_recommendations(
    *,
    metrics: Dict[str, str],
    candidates: Dict[str, List[Dict[str, Any]]],
    save_dir: Optional[str],
) -> None:
    """Plot threshold candidates without participating in threshold selection."""
    n_metrics = len(metrics)
    n_cols = min(2, n_metrics)
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows), constrained_layout=True
    )
    axes = np.array(axes).flatten()

    for i, (metric, title) in enumerate(metrics.items()):
        if metric not in candidates:
            continue
        values = np.asarray(candidates[metric][0].get("_plot_values", []), dtype=float)
        if values.size == 0:
            continue
        ax = axes[i]
        ax.hist(values, bins=50, alpha=0.75, edgecolor="black")
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xlabel(metric.replace("_", " ").title(), fontsize=12)
        ax.set_ylabel("Frequency", fontsize=12)

        for candidate in candidates[metric]:
            threshold = candidate["threshold"]
            if not np.isfinite(threshold):
                continue
            color = "red" if candidate["direction"] == "lower" else "orange"
            label = (
                f"{candidate['threshold_key'] or candidate['direction']} "
                f"({candidate['method_label']}): {threshold:.1f}"
            )
            ax.axvline(
                x=threshold,
                color=color,
                linestyle="--",
                alpha=0.8,
                linewidth=1.5,
                label=label,
            )
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="upper right")
        ax.grid(axis="y", linestyle="--", alpha=0.6)

    for j in range(len(metrics), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Suggested QC Thresholds from Data Distribution",
        fontsize=18,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(
            save_path / "qc_threshold_suggestions.png",
            dpi=300,
            bbox_inches="tight",
        )
    plt.show()


def _serializable_candidates(
    candidates: Dict[str, List[Dict[str, Any]]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Drop plotting-only arrays from recommendation candidates."""
    clean: Dict[str, List[Dict[str, Any]]] = {}
    for metric, items in candidates.items():
        clean[metric] = []
        for item in items:
            clean[metric].append({k: v for k, v in item.items() if k != "_plot_values"})
    return clean


# --- Main Functions ---


def recommend_qc_thresholds(
    adata: AnnData,
    method: Literal["mad", "iqr", "percentile"] = "mad",
    mad_multipliers: Union[float, List[float]] = [3.0, 4.0, 5.0],
    iqr_multiplier: float = 1.5,
    percentile_range: Tuple[float, float] = (2.5, 97.5),
    plot_distributions: bool = False,
    save_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Recommend QC thresholds and return a structured recommendation bundle.

    This is the preferred threshold recommendation API. It separates statistical
    recommendation from threshold resolution, evidence marking, and filtering.
    The bundle includes both structured candidate evidence and the legacy table
    view needed by downstream decision code.
    """
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required QC columns: {missing_cols}")

    if method not in {"mad", "iqr", "percentile"}:
        raise ValueError(f"Unknown threshold recommendation method: {method}")

    log.info(f"Recommending QC thresholds using '{method}' method...")

    if isinstance(mad_multipliers, (int, float)):
        mad_multipliers = [mad_multipliers]

    metrics = _threshold_metrics(adata)
    top_gene_cols = [
        col for col in adata.obs.columns if re.match(r"pct_counts_in_top_\d+_genes", col)
    ]
    metric_map = _threshold_metric_map(top_gene_cols)
    all_suggestions: Dict[str, Dict[str, Any]] = {}
    candidates: Dict[str, List[Dict[str, Any]]] = {}

    for metric in metrics:
        data = adata.obs[metric].dropna()
        values = data.values
        is_count_metric = _is_count_metric(metric)
        hard_filter_candidate = not _is_review_only_metric(metric)

        if method == "mad":
            for multiplier in mad_multipliers:
                level_name = f"mad_x{multiplier}"
                lower_bound, upper_bound = compute_mad_bounds(
                    values, nmads=multiplier, direction="both"
                )
                if lower_bound == upper_bound:
                    log.warning(
                        f"MAD for metric '{metric}' is zero at MAD x{multiplier}. "
                        f"Falling back to percentile-based bounds."
                    )
                    lower_pct, upper_pct = data.quantile([0.05, 0.95])
                    lower_bound = float(lower_pct)
                    upper_bound = float(upper_pct)

                if is_count_metric:
                    lower_bound = max(0.0, lower_bound)
                    min_key, max_key = metric_map[metric]
                    _add_candidate(
                        candidates=candidates,
                        table=all_suggestions,
                        metric=metric,
                        level_name=level_name,
                        direction="lower",
                        threshold=int(lower_bound),
                        values=values,
                        threshold_key=min_key,
                        method="mad",
                        method_label=f"MAD x{multiplier}",
                        hard_filter_candidate=hard_filter_candidate,
                    )
                    _add_candidate(
                        candidates=candidates,
                        table=all_suggestions,
                        metric=metric,
                        level_name=level_name,
                        direction="upper",
                        threshold=int(upper_bound),
                        values=values,
                        threshold_key=max_key,
                        method="mad",
                        method_label=f"MAD x{multiplier}",
                        hard_filter_candidate=hard_filter_candidate,
                    )
                else:  # Percentage metric
                    key = metric_map.get(metric)
                    _add_candidate(
                        candidates=candidates,
                        table=all_suggestions,
                        metric=metric,
                        level_name=level_name,
                        direction="upper",
                        threshold=min(100.0, float(upper_bound)),
                        values=values,
                        threshold_key=key,
                        method="mad",
                        method_label=f"MAD x{multiplier}",
                        hard_filter_candidate=hard_filter_candidate,
                    )

        elif method == "iqr":
            level_name = f"iqr_x{iqr_multiplier}"
            q25, q75 = data.quantile([0.25, 0.75])
            iqr = q75 - q25
            upper_bound = q75 + iqr_multiplier * iqr

            if is_count_metric:
                lower_bound = max(0, q25 - iqr_multiplier * iqr)
                min_key, max_key = metric_map[metric]
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="lower",
                    threshold=int(lower_bound),
                    values=values,
                    threshold_key=min_key,
                    method="iqr",
                    method_label=f"IQR x{iqr_multiplier}",
                    hard_filter_candidate=hard_filter_candidate,
                )
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="upper",
                    threshold=int(upper_bound),
                    values=values,
                    threshold_key=max_key,
                    method="iqr",
                    method_label=f"IQR x{iqr_multiplier}",
                    hard_filter_candidate=hard_filter_candidate,
                )
            else:  # Percentage metric
                key = metric_map.get(metric)
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="upper",
                    threshold=min(100.0, float(upper_bound)),
                    values=values,
                    threshold_key=key,
                    method="iqr",
                    method_label=f"IQR x{iqr_multiplier}",
                    hard_filter_candidate=hard_filter_candidate,
                )

        elif method == "percentile":
            level_name = f"percentile_{percentile_range[0]}-{percentile_range[1]}"
            upper_bound = data.quantile(percentile_range[1] / 100)

            if is_count_metric:
                lower_bound = data.quantile(percentile_range[0] / 100)
                min_key, max_key = metric_map[metric]
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="lower",
                    threshold=int(lower_bound),
                    values=values,
                    threshold_key=min_key,
                    method="percentile",
                    method_label=f"{percentile_range[0]}-{percentile_range[1]} percentile",
                    hard_filter_candidate=hard_filter_candidate,
                )
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="upper",
                    threshold=int(upper_bound),
                    values=values,
                    threshold_key=max_key,
                    method="percentile",
                    method_label=f"{percentile_range[0]}-{percentile_range[1]} percentile",
                    hard_filter_candidate=hard_filter_candidate,
                )
            else:  # Percentage metric
                key = metric_map.get(metric)
                _add_candidate(
                    candidates=candidates,
                    table=all_suggestions,
                    metric=metric,
                    level_name=level_name,
                    direction="upper",
                    threshold=min(100.0, float(upper_bound)),
                    values=values,
                    threshold_key=key,
                    method="percentile",
                    method_label=f"{percentile_range[0]}-{percentile_range[1]} percentile",
                    hard_filter_candidate=hard_filter_candidate,
                )

        for candidate in candidates.get(metric, []):
            candidate.setdefault("_plot_values", values)

    if plot_distributions and metrics:
        _plot_threshold_recommendations(
            metrics=metrics,
            candidates=candidates,
            save_dir=save_dir,
        )

    suggested_thresholds_df = _build_threshold_dataframe(all_suggestions)
    default_thresholds_obj = _thresholds_from_dataframe(suggested_thresholds_df)
    serializable_candidates = _serializable_candidates(candidates)
    hard_filter_exclusions = [
        {
            "metric": metric,
            "reason": "review_only_threshold_candidate",
            "candidate_count": len(items),
        }
        for metric, items in serializable_candidates.items()
        if items and all(item.get("review_only") for item in items)
    ]

    bundle = {
        "schema_version": THRESHOLD_RECOMMENDATION_BUNDLE_SCHEMA_VERSION,
        "threshold_result_schema_version": THRESHOLD_RESULT_SCHEMA_VERSION,
        "method": method,
        "method_parameters": {
            "mad_multipliers": list(mad_multipliers),
            "iqr_multiplier": iqr_multiplier,
            "percentile_range": tuple(percentile_range),
        },
        "metrics": list(metrics.keys()),
        "candidate_thresholds": serializable_candidates,
        "hard_filter_exclusions": hard_filter_exclusions,
        "threshold_table": suggested_thresholds_df,
        "threshold_table_records": suggested_thresholds_df.to_dict(orient="index"),
        "suggested_thresholds": default_thresholds_obj,
        "suggested_thresholds_dict": default_thresholds_obj.to_dict(),
        "score_semantics": (
            "Threshold recommendations are candidates. Filtering starts only after "
            "threshold decision and mark-evidence steps."
        ),
    }

    log.info("Comparison of recommended QC thresholds:")
    log.info("\n" + suggested_thresholds_df.to_string())
    return bundle


def _resolve_qc_thresholds(
    *,
    intelligent: Optional[Dict[str, Any]] = None,
    mad: Optional[Dict[str, Any]] = None,
    manual: Optional[Dict[str, Any]] = None,
    policy: Literal["intelligent_then_mad", "mad_then_intelligent", "manual_override"] = "intelligent_then_mad",
) -> QCThresholds:
    """Merge QC threshold sources into a single ``QCThresholds`` object.

    Internal helper for combining data-driven recommendations with manual
    project overrides. Manual overrides are treated as safety bounds unless
    ``policy='manual_override'`` is selected:

    - ``manual.min_genes`` / ``manual.min_counts`` act as floors.
    - ``manual.max_genes`` / ``manual.max_counts`` / ``manual.pc_mt`` act as
      ceilings.

    Parameters
    ----------
    intelligent
        Threshold dict from ``recommend_intelligent_qc`` (keys like ``min_genes``,
        ``min_counts``, ``pc_mt``).
    mad
        Threshold dict from ``recommend_qc_thresholds`` (a ``QCThresholds.to_dict()``).
    manual
        Manual overrides. Missing keys are ignored.
    policy
        Which data-driven source takes priority when both are provided.

    Returns:
    -------
    QCThresholds
        Resolved thresholds.

    Examples:
    --------
    >>> thresholds = _resolve_qc_thresholds(
    ...     intelligent={"min_genes": 300, "pc_mt": 8.0},
    ...     mad={"min_genes": 200, "max_genes": 8000, "pc_mt": 10.0},
    ...     manual={"min_genes": 500, "pc_mt": 12.0},
    ... )
    """
    intelligent = intelligent or {}
    mad = mad or {}
    manual = manual or {}

    def _pick_data(key: str) -> Any:
        i_val = intelligent.get(key)
        m_val = mad.get(key)
        if policy == "manual_override":
            return manual.get(key) if manual.get(key) is not None else (i_val if i_val is not None else m_val)
        if policy == "intelligent_then_mad":
            return i_val if i_val is not None else m_val
        # mad_then_intelligent
        return m_val if m_val is not None else i_val

    min_genes = _pick_data("min_genes")
    max_genes = _pick_data("max_genes")
    min_counts = _pick_data("min_counts")
    max_counts = _pick_data("max_counts")
    pc_mt = _pick_data("pc_mt")
    pc_hb = _pick_data("pc_hb")
    pc_top_genes = _pick_data("pc_top_genes") or {}

    # Manual overrides as floors/ceilings
    if manual.get("min_genes") is not None and min_genes is not None:
        min_genes = max(min_genes, manual["min_genes"])
    elif manual.get("min_genes") is not None:
        min_genes = manual["min_genes"]

    if manual.get("min_counts") is not None and min_counts is not None:
        min_counts = max(min_counts, manual["min_counts"])
    elif manual.get("min_counts") is not None:
        min_counts = manual["min_counts"]

    if manual.get("max_genes") is not None and max_genes is not None:
        max_genes = min(max_genes, manual["max_genes"])
    elif manual.get("max_genes") is not None:
        max_genes = manual["max_genes"]

    if manual.get("max_counts") is not None and max_counts is not None:
        max_counts = min(max_counts, manual["max_counts"])
    elif manual.get("max_counts") is not None:
        max_counts = manual["max_counts"]

    if manual.get("pc_mt") is not None and pc_mt is not None:
        pc_mt = min(pc_mt, manual["pc_mt"])
    elif manual.get("pc_mt") is not None:
        pc_mt = manual["pc_mt"]

    if manual.get("pc_hb") is not None and pc_hb is not None:
        pc_hb = min(pc_hb, manual["pc_hb"])
    elif manual.get("pc_hb") is not None:
        pc_hb = manual["pc_hb"]

    kwargs: Dict[str, Any] = {}
    if min_genes is not None:
        kwargs["min_genes"] = min_genes
    if max_genes is not None:
        kwargs["max_genes"] = max_genes
    if min_counts is not None:
        kwargs["min_counts"] = min_counts
    if max_counts is not None:
        kwargs["max_counts"] = max_counts
    if pc_mt is not None:
        kwargs["pc_mt"] = pc_mt
    if pc_hb is not None:
        kwargs["pc_hb"] = pc_hb
    if pc_top_genes:
        kwargs["pc_top_genes"] = pc_top_genes

    return QCThresholds(**kwargs)


def _thresholds_to_dict(value: Any) -> Dict[str, Any]:
    """Convert threshold objects or mappings into plain dictionaries."""
    if value is None:
        return {}
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if hasattr(value, "model_dump"):
        return value.model_dump()
    return dict(value)


def _recommendation_payload(bundle: Dict[str, Any]) -> Dict[str, Any]:
    """Return a serializable threshold recommendation payload."""
    return {
        "schema_version": bundle.get("schema_version"),
        "threshold_result_schema_version": bundle.get("threshold_result_schema_version"),
        "method": bundle.get("method"),
        "method_parameters": bundle.get("method_parameters", {}),
        "metrics": bundle.get("metrics", []),
        "candidate_thresholds": bundle.get("candidate_thresholds", {}),
        "hard_filter_exclusions": bundle.get("hard_filter_exclusions", []),
        "threshold_table": bundle.get("threshold_table_records", {}),
        "suggested_thresholds": bundle.get("suggested_thresholds_dict", {}),
        "score_semantics": bundle.get("score_semantics"),
    }


def decide_qc_thresholds(
    adata: AnnData,
    *,
    intelligent_thresholds: Optional[Dict[str, Any]] = None,
    manual_thresholds: Optional[Dict[str, Any]] = None,
    threshold_method: Literal["mad", "iqr", "percentile"] = "mad",
    threshold_policy: Literal[
        "intelligent_then_mad", "mad_then_intelligent", "manual_override"
    ] = "intelligent_then_mad",
    plot_distributions: Optional[bool] = None,
    save_dir: Optional[str] = None,
    **threshold_kwargs,
) -> Dict[str, Any]:
    """Generate and resolve QC thresholds without mutating ``adata``."""
    recommendation = recommend_qc_thresholds(
        adata,
        method=threshold_method,
        plot_distributions=bool(save_dir) if plot_distributions is None else plot_distributions,
        save_dir=save_dir,
        **threshold_kwargs,
    )
    threshold_table = recommendation["threshold_table"]
    suggested = recommendation["suggested_thresholds"]
    distribution_thresholds = suggested.to_dict()
    resolved = _resolve_qc_thresholds(
        intelligent=intelligent_thresholds,
        mad=distribution_thresholds,
        manual=manual_thresholds,
        policy=threshold_policy,
    )
    record_threshold_recommendation(
        adata,
        source=f"distribution_{threshold_method}",
        payload=_recommendation_payload(recommendation),
    )
    if intelligent_thresholds:
        record_threshold_recommendation(
            adata,
            source="intelligent_qc",
            payload=_thresholds_to_dict(intelligent_thresholds),
        )
    if manual_thresholds:
        record_threshold_recommendation(
            adata,
            source="manual",
            payload=_thresholds_to_dict(manual_thresholds),
        )
    record_threshold_decision(
        adata,
        resolved_thresholds=resolved.to_dict(),
        policy=threshold_policy,
        sources={
            "distribution": distribution_thresholds,
            "intelligent": _thresholds_to_dict(intelligent_thresholds),
            "manual": _thresholds_to_dict(manual_thresholds),
        },
        notes=["Resolved thresholds are used for QC evidence marking."],
    )
    return {
        "recommendation_bundle": recommendation,
        "threshold_table": threshold_table,
        "suggested_thresholds": suggested,
        "resolved_thresholds": resolved,
        "threshold_policy": threshold_policy,
    }


def apply_qc_threshold_decision(
    adata: AnnData,
    *,
    resolved_thresholds: Any,
    sample_key: Optional[str] = None,
    group_key: Optional[str] = None,
    criteria_to_filter: Optional[List[str]] = None,
    filter_cells_result: bool = False,
    save_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Apply a resolved QC threshold decision and record execution outputs."""
    from ..config import FilterConfig, MarkingConfig
    from ..filtering import filter_cells
    from .marking import audit_filtering, mark_low_quality_cell

    if sample_key is None:
        for candidate in ("sampleID", "sample_id", "sample", "batch"):
            if candidate in adata.obs.columns:
                sample_key = candidate
                break
    if sample_key is None:
        adata.obs["_sclucid_sample"] = "sample_0"
        sample_key = "_sclucid_sample"

    marking_cfg = MarkingConfig(thresholds=resolved_thresholds, save_dir=save_dir)
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

    decision_record = {
        "schema_version": "qc_threshold_execution_v1",
        "resolved_thresholds": resolved_thresholds.to_dict()
        if hasattr(resolved_thresholds, "to_dict")
        else dict(resolved_thresholds),
        "criteria_to_filter": criteria_to_filter or [],
        "filter_cells_result": bool(filter_cells_result),
    }
    qc_uns = adata.uns.setdefault("sclucid", {}).setdefault("qc", {})
    qc_uns["qc_threshold_decision"] = decision_record
    record_threshold_decision(
        adata,
        resolved_thresholds=decision_record["resolved_thresholds"],
        policy="execution_applied",
        sources={"resolved": decision_record["resolved_thresholds"]},
        sample_key=sample_key,
        notes=[
            "Threshold decision was applied to mark evidence columns.",
            "Filtering was executed in this helper."
            if filter_cells_result
            else "Filtering was not executed.",
        ],
    )
    return {
        "adata": adata,
        "filtered": filtered,
        "filtering_audit": filtering_audit,
        "criteria_to_filter": criteria_to_filter or [],
        "sample_key": sample_key,
        "decision_record": decision_record,
    }


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
    criteria_to_filter: Optional[List[str]] = None,
    filter_cells_result: bool = False,
    save_dir: Optional[str] = None,
    **threshold_kwargs,
) -> Dict[str, Any]:
    """Run threshold recommendation, decision, marking, and optional filtering."""
    decision = decide_qc_thresholds(
        adata,
        intelligent_thresholds=intelligent_thresholds,
        manual_thresholds=manual_thresholds,
        threshold_method=threshold_method,
        threshold_policy=threshold_policy,
        plot_distributions=plot_distributions,
        save_dir=save_dir,
        **threshold_kwargs,
    )
    execution = apply_qc_threshold_decision(
        adata,
        resolved_thresholds=decision["resolved_thresholds"],
        sample_key=sample_key,
        group_key=group_key,
        criteria_to_filter=criteria_to_filter,
        filter_cells_result=filter_cells_result,
        save_dir=save_dir,
    )
    decision.update(execution)
    return decision
