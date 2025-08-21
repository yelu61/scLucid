"""
Specialized tools for advanced single-cell RNA-seq data analysis.

This module provides a suite of high-level functions for specialized analyses,
including CNV inference, trajectory and dynamics analysis, gene regulatory
network inference, compositional analysis, and cell-cell communication.
It also features a robust bridge to the R ecosystem.
"""

import logging

log = logging.getLogger(__name__)

# --- Python-native Tools (High-Level API) ---
from .infercnv import run_cnv_analysis, find_tumor
from .scvelo import run_velocity_analysis, plot_velocity_results
from .trajectory import run_trajectory_analysis, plot_trajectory
from .scenic import run_scenic, analyze_scenic_results
from .cellphonedb import run_cellphonedb
from .sccoda import run_sccoda

# --- R Bridge (Class-based access to R tools) ---
from .rtools import RTools

# Define the public API for the 'tools' module
__all__ = [
    # CNV Analysis
    "run_cnv_analysis",
    "find_tumor",
    
    # Trajectory & Dynamics (Unified Interface)
    "run_trajectory_analysis",
    "plot_trajectory",
    
    # RNA Velocity (Can be called by trajectory module or directly)
    "run_velocity_analysis",
    "plot_velocity_results",
    
    # Gene Regulatory Network Analysis (pySCENIC)
    "run_scenic",
    "analyze_scenic_results",
    
    # Cell-Cell Communication (CellPhoneDB)
    "run_cellphonedb",
    
    # Compositional Analysis
    "run_sccoda",
    
    # R Tools Bridge (Access R functions like CellChat, Monocle3 via this class)
    "RTools",
]
