"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
"""

# Import and expose key functions from submodules
from .utils import identify_outliers, use_layer_as_X
from .plotting import (
    plot_composition,
    plot_embedding,
    plot_enrichment,
    plot_marker_expression,
    plot_marker_heatmap,
    plot_feature_correlation,
    plot_multi_modality,
    plot_pseudotime,
    plot_spatial,
    plot_volcano
)


# Define what should be accessible when importing from this module
__all__ = [
    "use_layer_as_X",
    "identify_outliers",
    "plot_composition",
    "plot_embedding",
    "plot_marker_heatmap",
    "plot_enrichment",
    "plot_marker_expression",
    "plot_feature_correlation",
    "plot_multi_modality",
    "plot_pseudotime",
    "plot_spatial",
    "plot_volcano"
]
