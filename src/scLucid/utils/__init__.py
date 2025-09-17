"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
"""

# Import and expose key functions from submodules
from .marker_manager import (
    Manager,
    CellType,
    get_marker_manager,
    KNOWN_SPECIES,
    MARKER_FORMATS,
)
from .utils import load_10x_data,identify_outliers, use_layer_as_X, sanitize_for_hdf5, subset_adata, subset_from_annotations
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
    plot_volcano,
    plot_ridge,
    plot_coexpression,
)


# Define what should be accessible when importing from this module
__all__ = [
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "identify_outliers",
    "subset_adata",
    "subset_from_annotations",
    "Manager",
    "CellType",
    "get_marker_manager",
    "KNOWN_SPECIES",
    "MARKER_FORMATS",
    "plot_composition",
    "plot_embedding",
    "plot_marker_heatmap",
    "plot_enrichment",
    "plot_marker_expression",
    "plot_feature_correlation",
    "plot_multi_modality",
    "plot_pseudotime",
    "plot_spatial",
    "plot_volcano",
    "plot_ridge",
    "plot_coexpression"
]
