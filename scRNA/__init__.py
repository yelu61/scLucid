"""
scRNA: A comprehensive toolkit for single-cell RNA sequencing data analysis.

This package provides a modular framework for analyzing single-cell RNA-seq data,
including quality control, preprocessing, analysis, specialized tools, and workflows.
"""

__version__ = "0.1.0"
__author__ = "Ye Lu"

from .config import settings

# Import main functionality for easy access
from . import qc
from . import preprocess
from . import analysis
from . import tools
from . import utils
#from . import workflows

# Optional
#from .workflows import run_standard_workflow

# Define what should be accessible when using `from scRNA import *`
__all__ = [
    # Configuration settings
    "settings",  
    # Modules
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "utils",
    #"workflows",
    #"run_standard_workflow"
]