"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, and cell type annotation.
"""

# Import and expose key functions from submodules
from .dimension import run_pca, run_umap, run_tsne
from .cluster import run_clustering, find_markers
from .differential import find_markers_between_groups, run_gsea
from .annotation import annotate_cell_types, score_gene_sets

# Define what should be accessible when importing from this module
__all__ = [
    "run_pca",
    "run_umap",
    "run_tsne",
    "run_clustering",
    "find_markers",
    "find_markers_between_groups",
    "run_gsea",
    "annotate_cell_types",
    "score_gene_sets"
]