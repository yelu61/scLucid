"""
Specialized tools for single-cell RNA-seq data analysis.

This module provides specialized analysis tools for advanced single-cell analytics,
including CNV analysis, trajectory inference, RNA velocity, and gene regulatory
network inference. It also provides a bridge for running R-based tools.
"""

import logging

log = logging.getLogger(__name__)

# Import and expose key functions from submodules

# --- Python-native Tools ---
from .infercnv import run_cnv_analysis, find_tumor
from .scvelo import run_velocity_analysis, plot_velocity_results
from .trajectory import run_trajectory_analysis, plot_trajectory
from .scenic import run_scenic, analyze_scenic_results
from .cellphonedb import run_cellphonedb
from .sccoda import run_sccoda

# --- R Bridge ---
from .rtools import RTools

# Define what should be accessible when importing the 'tools' module
__all__ = [
    # CNV analysis
    "run_cnv_analysis",
    "find_tumor",
    # Trajectory analysis (unified interface)
    "run_trajectory_analysis",
    "plot_trajectory",
    # RNA velocity (can be called by trajectory or directly)
    "run_velocity_analysis",
    "plot_velocity_results",
    # Gene Regulatory Network analysis (pySCENIC)
    "run_scenic",
    "analyze_scenic_results",
    # Cell-cell communication analysis (CellPhoneDB)
    "run_cellphonedb",
    # Compositional analysis
    "run_sccoda",
    # R Tools Bridge (access R functions via this class)
    "RTools",
]
