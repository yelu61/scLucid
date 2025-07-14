"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, and cell type annotation.
"""

# Import and expose key functions from submodules
from .annotation import annotate_clusters, score_cell_types
from .cluster import find_resolution, merge_clusters
from .differential import filter_markers, find_markers, get_conserved_markers
from .dimension import plot_composition, plot_embedding, plot_marker_heatmap
from .manager import Manager

# Define what should be accessible when importing from this module
__all__ = [
    "Manager",
    "annotate_clusters",
    "score_cell_types",
    "find_resolution",
    "merge_clusters",
    "find_markers",
    "filter_markers",
    "get_conserved_markers",
    "plot_embedding",
    "plot_marker_heatmap",
    "plot_composition",
]
