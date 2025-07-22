"""
scRNA: A comprehensive toolkit for single-cell RNA sequencing data analysis.

This package provides a modular framework for analyzing single-cell RNA-seq data,
including quality control, preprocessing, analysis, specialized tools, and workflows.
"""

__version__ = "0.1.0"
__author__ = "Ye Lu"

# Import main functionality for easy access
from . import analysis, preprocess, qc, tools, utils
# from . import workflows
from .config import settings, set_figure_params

# from . import workflows

# Optional
# from .workflows import run_standard_workflow

# Define what should be accessible when using `from scRNA import *`
__all__ = [
    # Configuration settings
    "settings",
    "set_figure_params",
    # Modules
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "utils",
    # "workflows",
    # "run_standard_workflow"
]
