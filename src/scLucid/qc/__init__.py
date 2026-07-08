"""Quality control public API for scLucid."""

import warnings
from collections.abc import Iterable
from importlib import import_module

__all__ = []


def _export(module: str, names: Iterable[str], *, optional: bool = False) -> bool:
    """Import names from a submodule without breaking package import."""
    try:
        loaded = import_module(f"{__name__}.{module}")
    except Exception as exc:
        level = "optional" if optional else "required"
        warnings.warn(
            f"Could not import {level} QC module '{module}': {exc}",
            ImportWarning,
        )
        return False

    found = False
    for name in names:
        if hasattr(loaded, name):
            globals()[name] = getattr(loaded, name)
            __all__.append(name)
            found = True
    return found


# Configuration
_export(
    "config",
    [
        "MetricsReportingConfig",
        "QCThresholds",
        "MarkerConfig",
        "DoubletConfig",
        "MarkingConfig",
        "FilterConfig",
        "QCWorkflowConfig",
    ],
)

# Core
_export("metrics", ["calculate_qc_metric"])
_export(
    "ambient",
    [
        "AMBIENT_CORRECTED_COUNTS_LAYER",
        "diagnose_ambient_rna",
        "register_external_ambient_result",
    ],
)
_export("ambient_backends", ["correct_ambient_rna", "cellbender_available"], optional=True)
_export("cycle", ["score_cell_cycle"])
_export(
    "artifacts",
    [
        "QC_ARTIFACT_CONTRACT_SCHEMA_VERSION",
        "QC_ARTIFACT_CONTRACT",
        "get_qc_artifact_contract",
        "record_qc_artifact_contract",
        "record_threshold_recommendation",
        "record_threshold_decision",
        "record_mark_evidence",
        "record_qc_decision_artifact",
        "record_filter_result",
        "record_benchmark_review",
        "cleanup_qc_intermediates",
    ],
)
_export(
    "policy.decisions",
    [
        "QC_DECISION_SCHEMA_VERSION",
        "QC_DECISION_VALUES",
        "build_qc_decisions",
        "score_qc_gene_panels",
        "summarize_qc_decisions",
    ],
)
_export(
    "doublet",
    [
        "generate_doublet_rates",
        "create_custom_marker_dict",
        "predict_doublets",
        "audit_doublets",
        "DOUBLET_OBS_COLUMNS",
    ],
)
_export(
    "filtering",
    [
        "filter_cells",
    ],
)
_export(
    "policy.thresholds",
    [
        "recommend_qc_thresholds",
        "decide_qc_thresholds",
        "apply_qc_threshold_decision",
        "run_qc_threshold_decision",
    ],
)

# Extended QC
_export(
    "policy.adaptive_threshold",
    [
        "AdaptiveThresholdLearner",
        "MultiMetricAdaptiveLearner",
        "THRESHOLD_RESULT_SCHEMA_VERSION",
        "build_threshold_result",
        "infer_qc_metric_type",
        "recommended_threshold_methods",
    ],
    optional=True,
)
_export(
    "reporting",
    [
        "generate_qc_report",
        "generate_qc_html_report",
    ],
    optional=True,
)
_export(
    "workflow",
    [
        "run_qc",
        "run_iterative_qc",
        "recommend_qc_policy",
        "apply_qc_policy",
        "run_standard_qc",
        "QC_WORKFLOW_STEPS",
        "QCWorkflowError",
    ],
)
_export(
    "policy.benchmark",
    [
        "QC_BENCHMARK_SCHEMA_VERSION",
        "BENCHMARK_PROFILES",
        "build_qc_benchmark_assessment",
        "compute_marker_fidelity",
        "compute_retention_metrics",
        "evaluate_qc_benchmark",
        "export_qc_benchmark_report",
        "infer_qc_benchmark_profile",
        "render_qc_benchmark_compact_markdown",
        "render_qc_benchmark_markdown",
    ],
)
_export(
    "trace",
    [
        "QC_TRACE_SCHEMA_VERSION",
        "QC_MODULE_MATURITY_SCHEMA_VERSION",
        "QC_REQUIRED_REVIEW_SECTIONS",
        "QC_REQUIRED_OBS_METRICS",
        "QC_STABLE_ENTRYPOINTS",
        "build_qc_decision_table",
        "build_qc_handoff_readiness",
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
    ],
)

# Intelligent QC
_export(
    "policy.intelligent_qc",
    [
        "IntelligentQCRecommender",
        "recommend_intelligent_qc",
        "QCRecommendation",
        "ThresholdRecommendation",
        "StrategyType",
    ],
    optional=True,
)
