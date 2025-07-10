"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, and cell type annotation.
"""

# Import and expose key functions from submodules
from .manager import Manager
from .annotation import annotate_clusters, score_cell_types
from .cluster import marker_guided_clustering, evaluate_resolution
from .differential import find_markers, marker_enrichment_analysis
from .dimension import plot_marker_expression, plot_cell_type_composition

# Define what should be accessible when importing from this module
__all__ = [
    "Manager",
    "annotate_clusters",
    "score_cell_types",
    "marker_guided_clustering",
    "evaluate_resolution",
    "find_markers",
    "marker_enrichment_analysis",
    "plot_marker_expression",
    "plot_cell_type_composition"
]