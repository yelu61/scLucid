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
_export("cycle", ["score_cell_cycle"])
_export(
    "doublet",
    [
        "generate_doublet_rates",
        "create_custom_marker_dict",
        "predict_doublets",
        "predict_doublets_with_profiling",
    ],
)
_export(
    "filtering",
    [
        "suggest_qc_thresholds",
        "identify_outliers",
        "generate_qc_report",
        "mark_low_quality_cell",
        "mark_low_quality_cells_adaptive",
        "filter_cells",
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
    ["EnhancedQCReport", "generate_qc_html_report", "InteractiveReportGenerator"],
    optional=True,
)
_export("workflow", ["run_advanced_qc", "run_standard_qc"])

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
