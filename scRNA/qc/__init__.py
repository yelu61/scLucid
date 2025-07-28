"""
Quality control module for single-cell RNA-seq data.

This module provides functions for calculating QC metrics, identifying low-quality
cells and doublets, and filtering cells based on QC parameters.
"""

# Import and expose key functions from submodules
from .doublet import generate_doublet_rates, is_doublet
from .filtering import _identify_outliers, filter_cells, is_low_quality_cell
from .metrics import calculate_qc_metric

# Define what should be accessible when importing from this module
__all__ = [
    "calculate_qc_metric",
    "_identify_outliers",
    "is_low_quality_cell",
    "generate_doublet_rates",
    "is_doublet",
    "filter_cells",
]
