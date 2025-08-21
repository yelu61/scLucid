"""
Quality Control (QC) module for scLucid.

This module provides a complete workflow for single-cell data quality control,
from metric calculation to filtering and reporting.
"""

# --- Import Configuration Classes ---
# Make the config objects directly accessible to the user, e.g., sclucid.qc.QCThresholds
from .config import QCThresholds, FilterConfig, MarkerConfig, DoubletConfig

# --- Import Core Functions ---
from .metrics import calculate_qc_metric
from .cycle import score_cell_cycle
from .doublet import (
    generate_doublet_rates, 
    predict_doublets, 
    export_doublet_stats,
    create_custom_marker_dict,
    create_doublet_marker_config_from_manager,
)
from .filtering import (
    mark_low_quality_cell,
    filter_cells,
    suggest_qc_thresholds,
    generate_qc_report,
)

# --- Import High-Level Workflow Functions ---
from .workflow import run_standard_qc, run_advanced_qc

# --- Define Public API for the Module ---
__all__ = [
    # Configuration Classes
    "QCThresholds",
    "FilterConfig",
    "MarkerConfig",
    "DoubletConfig",
    
    # Core Functions
    "calculate_qc_metric",
    "score_cell_cycle",
    "generate_doublet_rates",
    "predict_doublets",
    "export_doublet_stats",
    "create_custom_marker_dict",
    "create_doublet_marker_config_from_manager",
    "mark_low_quality_cell",
    "suggest_qc_thresholds",
    "filter_cells",
    "generate_qc_report",

    # Workflow Functions
    "run_standard_qc",
    "run_advanced_qc",
]