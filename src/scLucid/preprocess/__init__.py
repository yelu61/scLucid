"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
scaling, batch correction, and other essential preprocessing steps.
"""

# --- Configuration Objects ---
from .config import (
    GraphConfig,
    HVGConfig,
    IntegrationConfig,
    NeighborsConfig,
    NormalizationConfig,
    ScalingConfig,
    WorkflowConfig,
)
from .hvg import (
    evaluate_hvg_stability,
    find_hvgs,
    plot_hvg_metrics,
    select_hvg_sets,
    suggest_hvg_choice,
)
from .integrate import batch_correction, evaluate_integration
from .neighbors import optimize_neighbors_pcs

# --- Core Functions ---
# --- Plotting & Evaluation Functions ---
from .normalize import normalize_data, plot_normalization_effect
from .scale import plot_scaling_effect, regress_out, scale_data

# --- High-Level Workflow ---
from .workflow import run_preprocessing

# --- Public API Definition ---
__all__ = [
    # Configuration
    "WorkflowConfig",
    "NormalizationConfig",
    "HVGConfig",
    "ScalingConfig",
    "IntegrationConfig",
    "NeighborsConfig",
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
