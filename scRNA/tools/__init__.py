"""
Specialized tools for single-cell RNA-seq data analysis.

This module provides specialized analysis tools for tasks like
copy number variation analysis and trajectory analysis.
"""

# Import and expose key functions from submodules
from .infercnv import run_cnv_analysis, find_tumor
from .trajectory import run_trajectory_analysis, order_cells_pseudotime

# Define what should be accessible when importing from this module
__all__ = [
    "run_cnv_analysis",
    "find_tumor",
    "run_trajectory_analysis",
    "order_cells_pseudotime"
]