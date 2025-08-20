"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
batch correction, scaling, and other preprocessing steps.
"""

__version__ = "0.1.0"

# Import key functions
from .cycle import score_cell_cycle, get_cell_cycle_genes
from .hvg import annotate_hvg, select_hvg, evaluate_hvg_stability, plot_hvg_metrics
from .integrate import batch_correction
from .normalize import normalize_data, regress_out, plot_normalization_comparison
from .scale import scale_data, plot_scaling_effect

__all__ = [
    # Normalization
    "normalize_data",
    "plot_normalization_comparison",
    "score_cell_cycle",
    "get_cell_cycle_genes",
    "regress_out",
    # Feature selection
    "annotate_hvg",
    "select_hvg",
    "evaluate_hvg_stability",
    "plot_hvg_metrics",
    # Integration
    "batch_correction",
    # Scaling
    "scale_data",
    "plot_scaling_effect",
]
