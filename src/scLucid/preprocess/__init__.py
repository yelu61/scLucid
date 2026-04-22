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
    PreprocessingWorkflowConfig,
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
from .workflow import (
    WORKFLOW_STEPS,
    PartialWorkflowResult,
    WorkflowError,
    run_preprocessing,
)

# --- Backend Abstraction ---
from .backend import (
    PreprocessingBackend,
    ScanpyBackend,
    RapidsBackend,
    get_backend,
    set_backend,
    list_available_backends,
)
from .gene_biotype import (
    annotate_gene_biotypes,
    filter_genes_by_biotype,
    get_biotype_statistics,
    recommend_biotype_strategy,
)
from .adaptive_normalize import (
    AdaptiveNormalizationConfig,
    adaptive_normalize,
    estimate_cell_size_factors,
    quality_aware_normalize,
)

# --- Intelligent Preprocessing ---
from .intelligent import (
    BatchCorrectionRecommendation,
    DataProfile,
    HVGRecommendation,
    IntelligentPreprocessConfig,
    IntelligentPreprocessRecommender,
    NeighborsRecommendation,
    PCARecommendation,
    PreprocessingStrategy,
    ResolutionRecommendation,
    recommend_intelligent_preprocessing,
    run_intelligent_preprocessing,
)

# --- Public API Definition ---
__all__ = [
    # Configuration
    "WorkflowConfig",
    "PreprocessingWorkflowConfig",
    "NormalizationConfig",
    "HVGConfig",
    "ScalingConfig",
    "IntegrationConfig",
    "NeighborsConfig",
    "GraphConfig",
    # Workflow
    "run_preprocessing",
    "WORKFLOW_STEPS",
    "WorkflowError",
    "PartialWorkflowResult",
    # Core Functions
    "normalize_data",
    "regress_out",
    "find_hvgs",
    "select_hvg_sets",
    "scale_data",
    "batch_correction",
    "annotate_gene_biotypes",
    "filter_genes_by_biotype",
    "get_biotype_statistics",
    "recommend_biotype_strategy",
    # Backend
    "PreprocessingBackend",
    "ScanpyBackend",
    "RapidsBackend",
    "get_backend",
    "set_backend",
    "list_available_backends",
    # Plotting & Evaluation
    "plot_normalization_effect",
    "plot_hvg_metrics",
    "plot_scaling_effect",
    "suggest_hvg_choice",
    "evaluate_hvg_stability",
    "evaluate_integration",
    "optimize_neighbors_pcs",
    # Intelligent Preprocessing
    "IntelligentPreprocessConfig",
    "IntelligentPreprocessRecommender",
    "PreprocessingStrategy",
    "DataProfile",
    "HVGRecommendation",
    "PCARecommendation",
    "NeighborsRecommendation",
    "ResolutionRecommendation",
    "BatchCorrectionRecommendation",
    "recommend_intelligent_preprocessing",
    "run_intelligent_preprocessing",
    # Adaptive Normalization
    "AdaptiveNormalizationConfig",
    "adaptive_normalize",
    "estimate_cell_size_factors",
    "quality_aware_normalize",
]
