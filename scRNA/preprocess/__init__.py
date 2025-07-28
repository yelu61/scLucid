"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
batch correction, scaling, and other preprocessing steps.
"""

__version__ = "0.1.0"

# Import key functions
from .cycle import score_cell_cycle
from .hvg import annotate_hvg, select_hvg
from .integrate import batch_correction
from .normalize import normalize_data, regress_out
from .scale import scale_data

__all__ = [
    # Normalization
    "normalize_data",
    "score_cell_cycle",
    "regress_out",
    # Feature selection
    "annotate_hvg",
    "select_hvg",
    # Integration
    "batch_correction",
    # Scaling
    "scale_data",
]
