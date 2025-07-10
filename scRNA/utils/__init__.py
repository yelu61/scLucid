"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
"""

# Import and expose key functions from submodules
from .plotting import plot_qc, plot_clusters, plot_gene_expression
from .io import read_data, write_data, merge_datasets
from .stats import calculate_statistics, perform_test

# Define what should be accessible when importing from this module
__all__ = [
    "plot_qc",
    "plot_clusters",
    "plot_gene_expression",
    "read_data",
    "write_data",
    "merge_datasets",
    "calculate_statistics",
    "perform_test"
]