"""
Preprocessing module for single-cell RNA-seq data.

This module provides functions for normalization, feature selection,
scaling, batch correction, and other essential preprocessing steps.
"""

# --- Configuration Objects ---
from .adaptive_normalize import (
    AdaptiveNormalizationConfig,
    estimate_cell_size_factors,
)

from .config import (
    GeneBiotypeConfig,
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
    filter_genes_by_biotype,
    get_biotype_statistics,
    load_gene_biotypes,
    recommend_biotype_strategy,
)
from .hvg import (
    evaluate_hvg_stability,
    find_hvgs,
    plot_hvg_metrics,
    select_and_audit_hvgs,
    select_hvg_sets,
    suggest_hvg_choice,
)
from .integrate import (
    batch_correction,
    decide_integration,
    detect_integration_confounding,
    diagnose_integration_risk,
    evaluate_integration,
)

from .neighbors import optimize_neighbors_pcs, run_embedding_pipeline
from .policy import apply_preprocess_policy, recommend_preprocess_policy

# --- Core Functions ---
# --- Plotting & Evaluation Functions ---
from .normalize import normalize_data, plot_normalization_effect
from .scale import diagnose_cell_cycle_regression, plot_scaling_effect, regress_out, scale_data
from .trace import (
    PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION,
    PREPROCESS_REQUIRED_REVIEW_SECTIONS,
    PREPROCESS_STABLE_ENTRYPOINTS,
    PREPROCESS_TRACE_SCHEMA_VERSION,
    build_layer_transition_table,
    build_normalization_decision_policy,
    build_preprocess_decision_summary,
    build_preprocess_layer_contract,
    build_preprocess_method_semantics,
    build_preprocess_module_maturity_assessment,
    build_preprocess_reviewer_table,
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
    run_iterative_preprocessing,
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
    "GeneBiotypeConfig",
    # Workflow
    "run_preprocessing",
    "run_iterative_preprocessing",
    "recommend_preprocess_policy",
    "apply_preprocess_policy",
    "WORKFLOW_STEPS",
    "PartialWorkflowResult",
    "WorkflowError",
    "PREPROCESS_REQUIRED_REVIEW_SECTIONS",
    "PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION",
    "PREPROCESS_STABLE_ENTRYPOINTS",
    "PREPROCESS_TRACE_SCHEMA_VERSION",
    "build_layer_transition_table",
    "build_normalization_decision_policy",
    "build_preprocess_decision_summary",
    "build_preprocess_layer_contract",
    "build_preprocess_method_semantics",
    "build_preprocess_module_maturity_assessment",
    "build_preprocess_reviewer_table",
    "build_qc_input_context",
    "build_step_evidence_summary",
    "enrich_preprocessing_review_summary",
    "get_preprocess_module_contract",
    "summarize_preprocess_review_summary",
    "validate_preprocess_module_completeness",
    "validate_preprocessing_review_summary",
    # Core Functions
    "normalize_data",
    "regress_out",
    "find_hvgs",
    "select_and_audit_hvgs",
    "select_hvg_sets",
    "scale_data",
    "diagnose_cell_cycle_regression",
    "batch_correction",
    "diagnose_integration_risk",
    "detect_integration_confounding",
    "decide_integration",
    "annotate_gene_biotypes",
    "filter_genes_by_biotype",
    "get_biotype_statistics",
    "load_gene_biotypes",
    "recommend_biotype_strategy",
    # Plotting & Evaluation
    "plot_normalization_effect",
    "plot_hvg_metrics",
    "plot_scaling_effect",
    "suggest_hvg_choice",
    "evaluate_hvg_stability",
    "evaluate_integration",
    "optimize_neighbors_pcs",
    "run_embedding_pipeline",
    # Adaptive Normalization
    "AdaptiveNormalizationConfig",
    "estimate_cell_size_factors",
]

# Importing selected objects from a submodule can leave the submodule object on
# this namespace. Remove it so ``scl.pp.adaptive_normalize`` is not mistaken for
# the removed top-level function alias.
globals().pop("adaptive_normalize", None)
