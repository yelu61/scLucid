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
from .cellphonedb import (
    run_cellphonedb,
    run_cellphonedb_batch,
    run_cellphonedb_by_group,
    summarize_cellphonedb,
)
from .infercnv import find_tumor, run_cnv_analysis
from .sccoda import (
    plot_sccoda_proportion_with_significance,
    run_sccoda,
    run_sccoda_batch,
    summarize_sccoda,
)
from .pyscenic import (
    analyze_scenic_results,
    export_scenic_report,
    run_scenic,
    run_scenic_batch,
    run_scenic_by_group,
)
#from .trajectory import run_trajectory_analysis, plot_trajectory

# --- R Bridge (Class-based access to R tools) ---
# from .rtools import RTools

# Define the public API for the 'tools' module
__all__ = [
    # CNV Analysis
    "run_cnv_analysis",
    "find_tumor",
    # Trajectory & Dynamics (Unified Interface)
    #"run_trajectory_analysis",
    #"plot_trajectory",
    # Gene Regulatory Network Analysis (pySCENIC)
    "run_scenic",
    "run_scenic_batch",
    "run_scenic_by_group",
    "analyze_scenic_results",
    "export_scenic_report",
    # Cell-Cell Communication (CellPhoneDB)
    "run_cellphonedb",
    "run_cellphonedb_batch",
    "run_cellphonedb_by_group",
    "summarize_cellphonedb",
    # Compositional Analysis
    "run_sccoda",
    "run_sccoda_batch",
    "summarize_sccoda",
    "plot_sccoda_proportion_with_significance",
    
    # R Tools Bridge (Access R functions like CellChat, Monocle3 via this class)
    # "RTools",
]
