"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, and cell type annotation.
"""

# Import and expose key functions from submodules
from .annotation import (
    annotate_clusters,
    evaluate_annotation,
    score_cell_types,
    transfer_labels,
)
from .cluster import (
    cluster_cells,
    find_resolution,
    merge_clusters,
    optimize_neighbors_pcs,
)
from .differential import (
    compare_groups,
    filter_markers,
    find_markers,
    get_conserved_markers,
    run_enrichment,
    visualize_markers,
)

# Define what should be accessible when importing from this module
__all__ = [
    # Clustering and dimensionality reduction
    "optimize_neighbors_pcs",
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
    # Annotation
    "annotate_clusters",
    "score_cell_types",
    "transfer_labels",
    "evaluate_annotation",
    # Differential expression
    "find_markers",
    "filter_markers",
    "get_conserved_markers",
    "run_enrichment",
    "compare_groups",
    "visualize_markers",
]
