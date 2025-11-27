"""
General-Purpose Visualization Library for scLucid.

This module provides a comprehensive suite of plotting functions to generate
publication-quality visualizations for all stages of single-cell analysis,
from QC and preprocessing to final annotation and downstream analysis.
"""

# Import and expose all user-facing plotting functions from the submodule
from .main import (
    plot_embedding,
    plot_faceted_embedding,
    plot_dotplot,
    plot_stacked_violin,
    plot_split_violin_with_stats,
    plot_marker_expression,
    plot_faceted_feature,
    plot_marker_heatmap,
    plot_volcano,
    plot_ridge,
    plot_feature_correlation,
    plot_coexpression,
    plot_differential_abundance,
)

# Define what should be accessible when a user does `from scLucid.plotting import *`
__all__ = [
    "plot_embedding",
    "plot_faceted_embedding",
    "plot_dotplot",
    "plot_stacked_violin",
    "plot_split_violin_with_stats",
    "plot_marker_expression",
    "plot_faceted_feature",
    "plot_marker_heatmap",
    "plot_volcano",
    "plot_ridge",
    "plot_feature_correlation",
    "plot_coexpression",
    "plot_differential_abundance",
]