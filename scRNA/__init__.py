"""
scRNA: A comprehensive toolkit for single-cell RNA sequencing data analysis.

This package provides a modular framework for analyzing single-cell RNA-seq data,
including quality control, preprocessing, analysis, specialized tools, and workflows.
"""

__version__ = "0.1.0"
__author__ = "Ye Lu"

# Import main functionality for easy access
from . import qc
from . import preprocess
from . import analysis
from . import tools
from . import utils
from . import workflows

# Import key functions for direct access from the package
from .qc import calculate_qc_metric, filter_low_quality_cells
from .preprocess import normalize_data, annotate_hvg, batch_correction
from .analysis import annotate_clusters
from .tools import run_cnv_analysis
#from .workflows import run_standard_workflow

# Define what should be accessible when using `from scRNA import *`
__all__ = [
    # Modules
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "utils",
    "workflows",
    
    # Key functions
    "calculate_qc_metric",
    "filter_low_quality_cells",
    "normalize_data",
    "annotate_hvg",
    "batch_correction",
    "annotate_clusters",
    "run_cnv_analysis",
    #"run_standard_workflow"
]