"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
scaling, batch correction, and other essential preprocessing steps.
"""

__version__ = "0.1.0"

# --- Configuration Objects ---
from .config import (
    NormalizationConfig,
    HVGConfig,
    IntegrationConfig,
    GraphConfig,
    PreprocessingWorkflowConfig,
)

# --- Core Functions ---
from .normalize import normalize_data, regress_out
from .hvg import find_hvgs, select_hvg_sets, suggest_hvg_choice
from .scale import scale_data
from .integrate import batch_correction
from .neighbors import optimize_neighbors_pcs

# --- Plotting & Evaluation Functions ---
from .normalize import plot_normalization_effect
from .hvg import plot_hvg_metrics, evaluate_hvg_stability
from .scale import plot_scaling_effect
from .integrate import evaluate_integration

# --- High-Level Workflow ---
from .workflow import run_preprocessing

# --- Public API Definition ---
__all__ = [
    # Configuration
    "PreprocessingWorkflowConfig",
    "NormalizationConfig",
    "HVGConfig",
    "IntegrationConfig",
    "GraphConfig",
    
    # Workflow
    "run_preprocessing",
    
    # Core Functions
    "normalize_data",
    "regress_out",
    "find_hvgs",
    "select_hvg_sets", 
    
    "scale_data",
    "batch_correction",
    
    # Plotting & Evaluation
    "plot_normalization_effect",
    "plot_hvg_metrics",
    "plot_scaling_effect",
    "suggest_hvg_choice",
    "evaluate_hvg_stability",
    "evaluate_integration",
    "optimize_neighbors_pcs",
]