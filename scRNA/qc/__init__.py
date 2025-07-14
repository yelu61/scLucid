"""
Quality control module for single-cell RNA-seq data.

This module provides functions for calculating QC metrics, identifying low-quality
cells and doublets, and filtering cells based on QC parameters.
"""

# Import and expose key functions from submodules
from .doublet import is_doublet
from .filtering import filter_cells, identify_outliers, is_low_quality_cell
from .metrics import calculate_qc_metric

# Define what should be accessible when importing from this module
__all__ = [
    "calculate_qc_metric",
    "identify_outliers",
    "is_low_quality_cell",
    "is_doublet",
    "filter_cells",
]
