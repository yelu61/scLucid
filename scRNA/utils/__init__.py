"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
"""

# Import and expose key functions from submodules
# from .plotting import plot_qc, plot_clusters, plot_gene_expression
# from .io import read_data, write_data, merge_datasets
# from .stats import calculate_statistics, perform_test
from .anndata_helpers import use_layer_as_X
from .plotting import (
    plot_composition,
    plot_embedding,
    plot_marker_heatmap,
    plot_enrichment,
)
# Define what should be accessible when importing from this module
__all__ = [
    "use_layer_as_X",
    "plot_composition",
    "plot_embedding",
    "plot_marker_heatmap",
    "plot_enrichment",
    ]
