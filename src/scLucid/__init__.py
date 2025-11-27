"""
scLucid: A Comprehensive System for Single-Cell Analysis
=========================================================

scLucid is a powerful and flexible Python toolkit for the analysis of
single-cell RNA-sequencing data.

Main modules:
- qc: Quality control
- preprocess: Data preprocessing  
- analysis: Clustering, annotation, differential expression
- plotting: Visualization
- marker: Marker gene management
- tools: Advanced analyses (velocity, CNV, trajectory, etc.)
- utils: Utility functions
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("sclucid")
except PackageNotFoundError:
    __version__ = "0.1"
    
# --- Core Modules ---
# Import the full modules for explicit access
from . import qc
from . import preprocess
from . import analysis
from . import markers
from . import plotting
from . import utils
# Make tools import optional
try:
    from . import tools
except ImportError as e:
    import warnings
    warnings.warn(
        f"Could not import tools module: {e}. "
        "Install with 'pip install sclucid[tools]' to use advanced features.",
        ImportWarning
    )
    tools = None

# --- Import with convenient aliases for interactive use ---
from . import preprocess as pp
from . import analysis as al
from . import tools as tl
from . import utils as ut
from . import plotting as pl

# --- Configuration and Settings ---
from .settings import setup_logging, set_figure_params, reset_figure_params
from .config import get_config, set_config, reset_config

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
    
    # Utility Functions
    'get_config',
    'set_config',
    'reset_config',
]