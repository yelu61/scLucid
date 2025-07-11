"""
Quality control module for single-cell RNA-seq data.

This module provides functions for calculating QC metrics, identifying low-quality
cells and doublets, and filtering cells based on QC parameters.
"""

# Version tracking
__version__ = "0.1.0"

# Import and expose key functions from submodules
from .metrics import calculate_qc_metric
from .filtering import is_low_quality_cell, filter_low_quality_cells, identify_outliers
from .doublet import is_doublet, filter_doublets

# Define what should be accessible when importing from this module
__all__ = [
    "calculate_qc_metric",
    "identify_outliers",
    "is_low_quality_cell", 
    "filter_low_quality_cells",
    "is_doublet",
    "filter_doublets"
]