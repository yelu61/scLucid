"""
Quality control module for single-cell RNA-seq data.

This module provides functions for calculating QC metrics, identifying low-quality
cells and doublets, and filtering cells based on QC parameters.
"""

# Import and expose key functions from submodules
from .metrics import calculate_qc_metric
from .filtering import is_low_quality_cell, identify_outliers, filter_cells
from .doublet import is_doublet

# Define what should be accessible when importing from this module
__all__ = [
    "calculate_qc_metric",
    "identify_outliers",
    "is_low_quality_cell", 
    "is_doublet",
    "filter_cells"
]