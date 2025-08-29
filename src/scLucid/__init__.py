"""
scLucid: A comprehensive and flexible system for single-cell genomics analysis.
"""

__version__ = "0.1.0"
__author__ = "Ye LU"

# --- Core Modules ---
# Import the full modules for explicit access
from . import qc
from . import preprocess
from . import analysis
from . import tools
from . import utils
#from . import datasets

# --- Import with convenient aliases for interactive use ---
from . import preprocess as pp
from . import analysis as al
from . import tools as tl
from . import utils as ut

# --- Configuration and Settings ---
from .settings import setup_logging, set_figure_params, reset_figure_params

# --- High-Level Workflows ---
# Expose the main workflow functions at the top level for easy access
from .qc.workflow import run_standard_qc, run_advanced_qc
from .preprocess.workflow import run_preprocessing

# Define the public API using __all__
__all__ = [
    # Convenient Aliases
    "pp",
    "al",
    "tl",
    "ut",
    
    # Full Modules
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "utils",
    #"datasets",
    
    # Configuration
    "setup_logging",
    "set_figure_params",
    "reset_figure_params",
    
    # High-Level Workflows
    "run_standard_qc",
    "run_advanced_qc",
    "run_preprocessing",
]