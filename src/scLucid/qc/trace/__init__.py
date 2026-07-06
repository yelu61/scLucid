"""QC trace and review-summary schema helpers.

This module keeps the workflow-facing QC audit contract independent from the
execution code. It turns recommendations, applied config, and filtering output
into a compact machine-readable decision table and evidence chain.
"""

from __future__ import annotations

from collections.abc import Mapping
from enum import Enum
from typing import Any

import pandas as pd
from anndata import AnnData

from scLucid.utils.contracts import _review_payload

from ...utils.context import is_tumor_context as _shared_is_tumor_context

QC_TRACE_SCHEMA_VERSION = "1.0"
QC_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

QC_REQUIRED_REVIEW_SECTIONS = {
    "recommendation_summary",
    "recommended_threshold_summary",
    "applied_threshold_summary",
    "user_override_summary",
    "sample_threshold_summary",
    "tumor_aware_summary",
    "filtering_summary",
    "qc_filtering_policy_summary",
    "qc_retention_audit_summary",
    "doublet_evidence_summary",
    "ambient_evidence_summary",
    "post_annotation_qc_review",
    "qc_benchmark_scorecard",
    "warnings",
    "decision_table",
    "qc_reviewer_table",
    "evidence_chain",
    "execution_trace",
    "output_health",
    "downstream_preprocess_recommendations",
    "qc_readiness",
    "review_action_items",
    "reproducibility_manifest",
    "evidence_bundle",
    "module_maturity",
}

QC_REQUIRED_OBS_METRICS = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
]

QC_STABLE_ENTRYPOINTS = (
    "scLucid.qc.run_qc",
    "scLucid.qc.recommend_qc_policy",
    "scLucid.qc.apply_qc_policy",
    "scLucid.qc.run_standard_qc",
    "scLucid.qc.calculate_qc_metric",
    "scLucid.qc.recommend_intelligent_qc",
    "scLucid.qc.run_qc_threshold_decision",
    "scLucid.qc.build_qc_decisions",
    "scLucid.qc.build_post_annotation_qc_review",
    "scLucid.qc.filter_cells",
)

QC_EXPECTED_ARTIFACTS = (
    "qc_review_summary.json",
    "qc_review_summary.md",
    "qc_benchmark.json",
    "qc_benchmark.md",
)

_PARAMETER_SPECS = [
    {
        "parameter": "min_genes",
        "recommended_key": "min_genes",
        "applied_path": ("marking_config", "thresholds", "min_genes"),
        "obs_metric": "n_genes_by_counts",
        "filtering_flag": "outlier_min_genes",
        "direction": "lower_bound",
    },
    {
        "parameter": "max_genes",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "max_genes"),
        "obs_metric": "n_genes_by_counts",
        "filtering_flag": "outlier_max_genes",
        "direction": "upper_bound",
    },
    {
        "parameter": "n_counts",
        "recommended_key": "n_counts",
        "applied_path": ("marking_config", "thresholds", "min_counts"),
        "obs_metric": "total_counts",
        "filtering_flag": "outlier_min_counts",
        "direction": "lower_bound",
    },
    {
        "parameter": "max_counts",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "max_counts"),
        "obs_metric": "total_counts",
        "filtering_flag": "outlier_max_counts",
        "direction": "upper_bound",
    },
    {
        "parameter": "max_mt_percent",
        "recommended_key": "max_mt_percent",
        "applied_path": ("marking_config", "thresholds", "pc_mt"),
        "obs_metric": "pct_counts_mt",
        "filtering_flag": "outlier_mt",
        "direction": "upper_bound",
    },
    {
        "parameter": "max_hb_percent",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "pc_hb"),
        "obs_metric": "pct_counts_hb",
        "filtering_flag": "outlier_hb",
        "direction": "upper_bound",
    },
    {
        "parameter": "doublet_threshold",
        "recommended_key": "doublet_threshold",
        "applied_path": ("doublet_config", "score_threshold"),
        "obs_metric": "doublet_score",
        "filtering_flag": "predicted_doublet",
        "direction": "upper_bound",
    },
    {
        "parameter": "nmads",
        "recommended_key": None,
        "applied_path": ("marking_config", "thresholds", "nmads"),
        "obs_metric": "qc_metric_distribution",
        "filtering_flag": "outlier_qc_metrics",
        "direction": "mad_outlier",
    },
]


def _json_safe(value: Any) -> Any:
    """Convert common scientific/Python objects into JSON-safe values."""
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def _to_dict(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, dict):
        return _json_safe(value)
    if hasattr(value, "to_dict"):
        return _json_safe(value.to_dict())
    if hasattr(value, "model_dump"):
        return _json_safe(value.model_dump())
    return {}


def _get_nested(data: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = data
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _recommendation_value(rec_dict: Mapping[str, Any], key: str | None, field: str) -> Any:
    if key is None:
        return None
    rec = rec_dict.get(key)
    if not isinstance(rec, Mapping):
        return None
    return rec.get(field)


def _decision_source(
    *,
    parameter: str,
    applied_value: Any,
    recommended_value: Any,
    user_overrides: Mapping[str, Any],
) -> str:
    if parameter in user_overrides:
        return "user_override"
    if recommended_value is not None and applied_value == recommended_value:
        return "recommendation"
    if applied_value is None:
        return "disabled_or_not_available"
    return "default_or_config"

from .review import (
    build_ambient_evidence_summary,
    build_post_annotation_qc_review,
    build_qc_benchmark_scorecard,
    build_qc_decision_table,
    build_qc_reviewer_table,
    enrich_qc_decision_table_for_review,
    get_qc_module_contract,
)
from .summary import (
    build_qc_module_maturity_assessment,
    enrich_qc_review_summary,
    summarize_qc_review_summary,
    validate_qc_module_completeness,
    validate_qc_review_summary,
)

__all__ = [
    "QC_TRACE_SCHEMA_VERSION",
    "QC_MODULE_MATURITY_SCHEMA_VERSION",
    "QC_REQUIRED_REVIEW_SECTIONS",
    "QC_REQUIRED_OBS_METRICS",
    "QC_STABLE_ENTRYPOINTS",
    "build_qc_decision_table",
    "build_qc_reviewer_table",
    "build_ambient_evidence_summary",
    "build_post_annotation_qc_review",
    "build_qc_benchmark_scorecard",
    "build_qc_module_maturity_assessment",
    "enrich_qc_decision_table_for_review",
    "enrich_qc_review_summary",
    "get_qc_module_contract",
    "summarize_qc_review_summary",
    "validate_qc_module_completeness",
    "validate_qc_review_summary",
]
