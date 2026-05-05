"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
scaling, batch correction, and other essential preprocessing steps.
"""

# --- Configuration Objects ---
from .adaptive_normalize import (
    AdaptiveNormalizationConfig,
    adaptive_normalize,
    estimate_cell_size_factors,
    quality_aware_normalize,
)

# --- Backend Abstraction ---
from .backend import (
    PreprocessingBackend,
    RapidsBackend,
    ScanpyBackend,
    get_backend,
    list_available_backends,
    set_backend,
)
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
from .gene_biotype import (
    annotate_gene_biotypes,
    apply_gene_biotype_strategy,
    filter_genes_by_biotype,
    get_biotype_statistics,
    get_gene_biotype_cache_dir,
    list_gene_biotype_resources,
    load_gene_biotypes,
    recommend_biotype_strategy,
)
from .hvg import (
    evaluate_hvg_stability,
    find_hvgs,
    plot_hvg_metrics,
    select_hvg_sets,
    suggest_hvg_choice,
)
from .integrate import batch_correction, evaluate_integration

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
from .neighbors import optimize_neighbors_pcs

# --- Core Functions ---
# --- Plotting & Evaluation Functions ---
from .normalize import normalize_data, plot_normalization_effect
from .scale import plot_scaling_effect, regress_out, scale_data
from .trace import (
    PREPROCESS_REQUIRED_REVIEW_SECTIONS,
    PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION,
    PREPROCESS_STABLE_ENTRYPOINTS,
    PREPROCESS_TRACE_SCHEMA_VERSION,
    build_preprocess_module_maturity_assessment,
    build_qc_input_context,
    build_step_evidence_summary,
    enrich_preprocessing_review_summary,
    get_preprocess_module_contract,
    summarize_preprocess_review_summary,
    validate_preprocess_module_completeness,
    validate_preprocessing_review_summary,
)

# --- High-Level Workflow ---
from .workflow import (
    WORKFLOW_STEPS,
    PartialWorkflowResult,
    WorkflowError,
    run_preprocessing,
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
    "PREPROCESS_REQUIRED_REVIEW_SECTIONS",
    "PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION",
    "PREPROCESS_STABLE_ENTRYPOINTS",
    "PREPROCESS_TRACE_SCHEMA_VERSION",
    "build_preprocess_module_maturity_assessment",
    "build_qc_input_context",
    "build_step_evidence_summary",
    "enrich_preprocessing_review_summary",
    "get_preprocess_module_contract",
    "summarize_preprocess_review_summary",
    "validate_preprocess_module_completeness",
    "validate_preprocessing_review_summary",
    # Core Functions
    "apply_gene_biotype_strategy",
    "normalize_data",
    "regress_out",
    "find_hvgs",
    "select_hvg_sets",
    "scale_data",
    "batch_correction",
    "annotate_gene_biotypes",
    "filter_genes_by_biotype",
    "get_biotype_statistics",
    "load_gene_biotypes",
    "list_gene_biotype_resources",
    "get_gene_biotype_cache_dir",
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
