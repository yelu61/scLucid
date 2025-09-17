"""
scLucid: An integrated, flexible, and exploratory analysis system for single-cell genomics.
"""

__version__ = "0.1.0"
__author__ = "Ye LU"

# --- Core Modules ---
# Import the full modules for explicit access
from . import qc
from . import preprocess
from . import analysis
from . import tools
from . import markers
from . import plotting
from . import utils

# --- Import with convenient aliases for interactive use ---
from . import preprocess as pp
from . import analysis as al
from . import tools as tl
from . import utils as ut
from . import plotting as pl

# --- Configuration and Settings ---
from .settings import setup_logging, set_figure_params, reset_figure_params

# --- High-Level Workflows ---
# Expose the main workflow functions at the top level for easy access
from .qc.workflow import run_standard_qc, run_advanced_qc
from .preprocess.workflow import run_preprocessing
from .analysis import run_annotation, characterize_clusters

# Define the public API using __all__
__all__ = [
    # Convenient Aliases
    "pp",
    "al",
    "tl",
    "ut",
    "pl",
    
    # Full Modules
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "markers",
    "plotting",
    "utils",
    
    # Configuration
    "setup_logging",
    "set_figure_params",
    "reset_figure_params",
    
    # High-Level Workflows
    "run_standard_qc",
    "run_advanced_qc",
    "run_preprocessing",
    "run_annotation",
    "characterize_clusters",
]