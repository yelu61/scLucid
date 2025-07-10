"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
batch correction, scaling, and other preprocessing steps.
"""

__version__ = "0.1.0"

# Import key functions
from .normalize import normalize_data, regress_out
from .hvg import annotate_hvg, select_hvg
from .scale import scale_data
from .cycle import score_cell_cycle
from .integrate import integrate_scanorama, integrate_scvi, integrate_harmony, batch_correction

__all__ = [
    # Normalization
    "normalize_data",
    "score_cell_cycle",
    "regress_out",
    # Feature selection
    "annotate_hvg",
    "select_hvg",
    # Integration
    "integrate_scanorama",
    "integrate_scvi",
    "integrate_harmony",
    "batch_correction",
    # Scaling
    "scale_data",
]