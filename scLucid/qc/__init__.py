"""
Quality control module for single-cell RNA-seq data.

This module provides functions for calculating QC metrics, identifying low-quality
cells and doublets, and filtering cells based on QC parameters.
"""

# Import and expose key functions from submodules
from .cycle import get_cell_cycle_genes, score_cell_cycle
from .doublet import (
    generate_doublet_rates,
    predict_doublets,
)
from .filtering import (
    filter_cells,
    generate_qc_report,
    mark_low_quality_cell,
    suggest_qc_thresholds,
)
from .metrics import calculate_qc_metric
from .workflow import run_advanced_qc, run_standard_qc

# Define what should be accessible when importing from this module
__all__ = [
    # calculating metrics
    "calculate_qc_metric",
    "get_cell_cycle_genes",
    "score_cell_cycle",
    # marking cells
    "suggest_qc_thresholds",
    "mark_low_quality_cell",
    "generate_doublet_rates",
    "predict_doublets",
    # filtering cells
    "filter_cells",
    # report
    "generate_qc_report",
    # workflow
    "run_standard_qc",
    "run_advanced_qc",
]
