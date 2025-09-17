"""
General-Purpose Visualization Library for scLucid.

This module provides a comprehensive suite of plotting functions to generate
publication-quality visualizations for all stages of single-cell analysis,
from QC and preprocessing to final annotation and downstream analysis.
"""

# Import and expose all user-facing plotting functions from the submodule
from .main import (
    plot_composition,
    plot_coexpression,
    plot_differential_abundance,
    plot_embedding,
    plot_enrichment,
    plot_feature_correlation,
    plot_marker_expression,
    plot_marker_heatmap,
    plot_multi_modality,
    plot_pseudotime,
    plot_ridge,
    plot_spatial,
    plot_volcano,
)

# Define what should be accessible when a user does `from scLucid.plotting import *`
__all__ = [
    "plot_composition",
    "plot_coexpression",
    "plot_differential_abundance",
    "plot_embedding",
    "plot_enrichment",
    "plot_feature_correlation",
    "plot_marker_expression",
    "plot_marker_heatmap",
    "plot_multi_modality",
    "plot_pseudotime",
    "plot_ridge",
    "plot_spatial",
    "plot_volcano",
]