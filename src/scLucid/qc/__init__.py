"""
Quality Control (QC) module for scLucid.

This module provides a complete workflow for single-cell data quality control,
from metric calculation to filtering and reporting.
"""

# --- Import Configuration Classes ---
# Make the config objects directly accessible to the user, e.g., sclucid.qc.QCThresholds
from .config import MetricsReportingConfig, QCThresholds, MarkerConfig, DoubletConfig, MarkingConfig, FilterConfig, QCWorkflowConfig

# --- Import Core Functions ---
from .metrics import calculate_qc_metric
from .cycle import score_cell_cycle
from .doublet import (
    generate_doublet_rates,
    create_custom_marker_dict,
    predict_doublets,
    predict_doublets_with_profiling
)
from .filtering import (
    suggest_qc_thresholds,
    identify_outliers,
    generate_qc_report,
    mark_low_quality_cell,
    mark_low_quality_cells_adaptive,
    filter_cells,
)

# --- Import High-Level Workflow Functions ---
from .workflow import run_advanced_qc, run_standard_qc

# --- Define Public API for the Module ---
__all__ = [
    # Configuration Classes
    "MetricsReportingConfig",
    "QCThresholds",
    "MarkerConfig",
    "DoubletConfig",
    "MarkingConfig",
    "FilterConfig",
    "QCWorkflowConfig",
    # Core Functions
    "calculate_qc_metric",
    "score_cell_cycle",
    "generate_doublet_rates",
    "create_custom_marker_dict",
    "predict_doublets",
    "predict_doublets_with_profiling",
    "suggest_qc_thresholds",
    "identify_outliers",
    "mark_low_quality_cell",
    "mark_low_quality_cells_adaptive",
    "filter_cells",
    "generate_qc_report",
    # Workflow Functions
    "run_standard_qc",
    "run_advanced_qc",
]
