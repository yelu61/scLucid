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
        "build_ambient_layer_contract",
        "diagnose_ambient_rna",
        "diagnose_empty_droplets",
        "infer_ambient_input_context",
        "record_ambient_correction_status",
        "record_ambient_layer_contract",
        "register_external_ambient_result",
        "correct_ambient_rna_linear",
    ],
)
_export("ambient_backends", ["correct_ambient_rna", "cellbender_available"], optional=True)
_export("cycle", ["score_cell_cycle"])
_export(
    "decisions",
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
        "predict_doublets_with_profiling",
        "audit_doublets",
    ],
)
_export(
    "filtering",
    [
        "suggest_qc_thresholds",
        "identify_outliers",
        "mark_low_quality_cell",
        "mark_low_quality_cells_adaptive",
        "filter_cells",
        "audit_filtering",
        "resolve_qc_thresholds",
        "decide_qc_thresholds",
        "apply_qc_threshold_decision",
        "run_qc_threshold_decision",
    ],
)

# Extended QC
_export(
    "adaptive_threshold",
    ["AdaptiveThresholdLearner", "MultiMetricAdaptiveLearner"],
    optional=True,
)
_export(
    "reporting",
    ["EnhancedQCReport", "generate_qc_report", "generate_qc_html_report", "InteractiveReportGenerator"],
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
    "benchmark",
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
    "intelligent_qc",
    [
        "IntelligentQCRecommender",
        "recommend_intelligent_qc",
        "QCRecommendation",
        "ThresholdRecommendation",
        "StrategyType",
    ],
    optional=True,
)

# Transitional aliases stay importable but are intentionally omitted from __all__.
try:
    _workflow_decision = import_module(f"{__name__}.filtering.workflow_decision")
    if hasattr(_workflow_decision, "run_qc_decision_workflow"):
        run_qc_decision_workflow = getattr(_workflow_decision, "run_qc_decision_workflow")
except Exception:
    pass

try:
    _workflow = import_module(f"{__name__}.workflow")
    if hasattr(_workflow, "run_advanced_qc"):
        run_advanced_qc = getattr(_workflow, "run_advanced_qc")
except Exception:
    pass
