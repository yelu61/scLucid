"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting, data input/output,
statistical calculations, workflow management, marker gene management,
data validation, and performance profiling.

Key Components:
- Workflow utilities: Progress tracking, error recovery, checkpoint management
- Marker management: Hierarchical cell type markers from TOML/JSON files
- Data loading: Built-in dataset loaders
- Storage management: Standardized adata.uns['sclucid'] storage (storage.py)
- Validation utilities: AnnData structure validation and analysis readiness checks
- Profiling utilities: Performance monitoring with time and memory tracking

Key exposed components:
- Manager: The main class for handling marker hierarchies.
- CellType: A dataclass representing a single cell type entry.
- get_marker_manager: A factory function to easily build a combined manager.
- KNOWN_SPECIES: A list of built-in species supported.
- MARKER_FORMATS: A list of supported marker file formats.
- Workflow utilities: get_progress_bar, PartialResultManager, WorkflowError
- Storage utilities: get_storage, save_result, load_result
- Validation utilities: validate_adata, ValidationError, assert_qc_ready
- Profiling utilities: PerformanceProfiler, profile_performance, memory_tracker
"""

# Import runtime safeguards
from ..runtime import (
    effective_n_jobs,
    is_ci_environment,
    run_joblib_or_sequential,
    setup_runtime_environment,
)
from .audit_report import export_audit_report

# Import and expose key functions from submodules
from .context import (
    AnalysisContext,
    DatasetProfile,
    DatasetType,
    ProjectContext,
    infer_analysis_context,
    infer_dataset_profile,
    is_multi_sample_hint,
    normalize_dataset_type,
)
from .contracts import (
    API_LAYER_CONTRACTS,
    API_LAYER_ORDER,
    MINIMAL_WORKFLOW_CONTRACT,
    REVIEW_SUMMARY_RECOMMENDED_KEYS,
    REVIEW_SUMMARY_REQUIRED_KEYS,
    SCHEMA_VERSION,
    SCLUCID_ROOT,
    STAGE_CONTRACTS,
    STAGE_ORDER,
    APILayerContract,
    AssayKeys,
    ContractError,
    ContractValidationResult,
    LayerKeys,
    LayerSemanticKeys,
    ModalityContractResult,
    ModalityKeys,
    Modules,
    ObsKeys,
    ObsmKeys,
    StageContract,
    UnsKeys,
    VarKeys,
    api_layer_contract_to_dict,
    build_config_lineage,
    ensure_sclucid_namespace,
    format_contract_error,
    get_api_layer_spec,
    get_contract_spec,
    get_minimal_workflow_contract,
    get_stage_contract,
    infer_anndata_semantics,
    module_namespace,
    normalize_review_summary,
    record_artifact,
    record_config_lineage,
    record_contract_result,
    record_error,
    register_anndata_semantics,
    stage_contract_to_dict,
    validate_all_stage_contracts,
    validate_modality_contract,
    validate_review_summary_schema,
    validate_stage_contract,
)
from .evidence import (
    EVIDENCE_SCHEMA_VERSION,
    DecisionRecord,
    EvidenceBundle,
    EvidenceItem,
    ReviewAction,
    model_to_dict,
)
from .helpers import (
    assess_matrix_semantics,
    build_metadata_dicts,
    merge_obs_metadata,
    print_sample_crosstab,
    sanitize_for_hdf5,
    subset_adata,
    subset_from_annotations,
    use_layer_as_X,
)
from .io import load_10x_data, read_10x, read_h5ad

# Import and expose key functions and classes from the submodule
from .manager import (
    KNOWN_SPECIES,
    MARKER_FORMATS,
    CellType,
    Manager,
    _get_cancer_markers,
    canonicalize_marker_label,
    get_gene_display_aliases,
    get_marker_aliases,
    get_marker_manager,
    load_gene_set_manager,
    load_gene_sets,
    load_marker_aliases,
)
from .manual_review import finalize_manual_review_summary
from .marker_sets import filter_marker_dict, flatten_marker_dict

# Import profiling utilities
from .profiling import (
    BenchmarkRunner,
    PerformanceProfiler,
    PerformanceStats,
    estimate_adata_memory,
    get_memory_usage,
    memory_tracker,
    profile_function,
    profile_performance,
)
from .resource_audit import (
    assert_trusted_resources,
    audit_curation_index,
    audit_geneset_resources,
    audit_marker_entry_quality,
    audit_marker_resources,
    audit_resource_manifest,
    build_resource_trust_report,
    classify_literature_resource_utility,
    load_marker_curation_literature_index,
    load_reference_index,
    load_resource_manifest,
)

# Import result cleanup utilities
from .result_cleanup import (
    clear_sclucid_results,
    compact_sclucid_uns,
    list_sclucid_modules,
)
from .step_result import (
    EvidenceLevel,
    StepResult,
    StepStatus,
    rollup_step_status,
    step_results_from_storage,
    step_results_to_storage,
    summarize_step_results,
)

# Import storage utilities (new simplified interface)
from .storage import (
    STORAGE_ROOT,
    VALID_MODULES,
    clear_storage,
    export_review_summary,
    get_storage,
    has_result,
    list_results,
    load_config,
    load_result,
    load_workflow_result,
    migrate_legacy_storage,
    save_result,
    save_workflow_result,
    write_h5ad_safe,
)

# Import validation utilities
from .validation import (
    ValidationError,
    assert_analysis_ready,
    assert_preprocessing_ready,
    assert_qc_ready,
    check_layer_consistency,
    validate_adata,
    validate_analysis_results,
    validate_config,
    validate_workflow_contract,
)
from .validation_scaffold import (
    COMPARATIVE_READINESS_LABEL,
    VALIDATION_SCAFFOLD_SCHEMA_VERSION,
    VALIDATION_SCOPE,
    build_qc_preprocess_validation,
    validation_table_to_dataframe,
    write_validation_outputs,
)

# Import workflow utilities
from .workflow_utils import (
    BaseWorkflow,
    PartialResultManager,
    RecoveryError,
    StepError,
    WorkflowCheckpoint,
    WorkflowError,
    WorkflowStepIterator,
    get_progress_bar,
    merge_partial_results,
    progress_decorator,
    with_error_recovery,
)

# Import data loading utilities
try:
    from .data_loader import (
        filter_by_species,  # noqa: F401
        filter_by_tissue_type,  # noqa: F401
        get_dataset_info,  # noqa: F401
        load_all_datasets,  # noqa: F401
        load_luad,  # noqa: F401
        load_melanoma,  # noqa: F401
        load_pbmc3k,  # noqa: F401
        print_dataset_summary,  # noqa: F401
    )

    _data_loader_available = True
except ImportError:
    _data_loader_available = False

# Define what should be accessible when importing from this module
__all__ = [
    # Validation
    "ValidationError",
    "validate_adata",
    "validate_config",
    "validate_analysis_results",
    "check_layer_consistency",
    "assess_matrix_semantics",
    "assert_preprocessing_ready",
    "assert_analysis_ready",
    "validate_workflow_contract",
    "write_h5ad_safe",
    # Helper functions
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
    "build_metadata_dicts",
    "AnalysisContext",
    "DatasetProfile",
    "DatasetType",
    "ProjectContext",
    "infer_analysis_context",
    "infer_dataset_profile",
    "is_multi_sample_hint",
    "normalize_dataset_type",
    "API_LAYER_CONTRACTS",
    "API_LAYER_ORDER",
    "APILayerContract",
    "AssayKeys",
    "ContractError",
    "ContractValidationResult",
    "LayerKeys",
    "LayerSemanticKeys",
    "MINIMAL_WORKFLOW_CONTRACT",
    "ModalityContractResult",
    "ModalityKeys",
    "Modules",
    "ObsmKeys",
    "ObsKeys",
    "SCHEMA_VERSION",
    "SCLUCID_ROOT",
    "STAGE_CONTRACTS",
    "STAGE_ORDER",
    "StageContract",
    "REVIEW_SUMMARY_RECOMMENDED_KEYS",
    "REVIEW_SUMMARY_REQUIRED_KEYS",
    "UnsKeys",
    "VarKeys",
    "api_layer_contract_to_dict",
    "build_config_lineage",
    "ensure_sclucid_namespace",
    "format_contract_error",
    "get_api_layer_spec",
    "get_contract_spec",
    "get_minimal_workflow_contract",
    "get_stage_contract",
    "infer_anndata_semantics",
    "module_namespace",
    "normalize_review_summary",
    "record_artifact",
    "record_contract_result",
    "record_config_lineage",
    "record_error",
    "register_anndata_semantics",
    "stage_contract_to_dict",
    "validate_all_stage_contracts",
    "validate_modality_contract",
    "validate_review_summary_schema",
    "validate_stage_contract",
    "effective_n_jobs",
    "is_ci_environment",
    "run_joblib_or_sequential",
    "setup_runtime_environment",
    "flatten_marker_dict",
    "filter_marker_dict",
    "finalize_manual_review_summary",
    "DecisionRecord",
    "EVIDENCE_SCHEMA_VERSION",
    "EvidenceBundle",
    "EvidenceItem",
    "ReviewAction",
    "model_to_dict",
    "EvidenceLevel",
    "StepResult",
    "StepStatus",
    "rollup_step_status",
    "step_results_from_storage",
    "step_results_to_storage",
    "summarize_step_results",
    # Result cleanup
    "clear_sclucid_results",
    "compact_sclucid_uns",
    "list_sclucid_modules",
    # Workflow utilities
    "get_progress_bar",
    "progress_decorator",
    "WorkflowError",
    "StepError",
    "RecoveryError",
    "WorkflowCheckpoint",
    "PartialResultManager",
    "WorkflowStepIterator",
    "BaseWorkflow",
    "with_error_recovery",
    "merge_partial_results",
    # Marker management
    "CellType",
    "Manager",
    "canonicalize_marker_label",
    "get_gene_display_aliases",
    "get_marker_aliases",
    "get_marker_manager",
    "load_gene_set_manager",
    "load_gene_sets",
    "load_marker_aliases",
    "_get_cancer_markers",
    "KNOWN_SPECIES",
    "MARKER_FORMATS",
    "assert_trusted_resources",
    "audit_curation_index",
    "audit_geneset_resources",
    "audit_marker_entry_quality",
    "audit_marker_resources",
    "audit_resource_manifest",
    "build_resource_trust_report",
    "classify_literature_resource_utility",
    "load_marker_curation_literature_index",
    "load_reference_index",
    "load_resource_manifest",
    # Storage management (new simplified interface)
    "get_storage",
    "save_result",
    "load_result",
    "load_config",
    "has_result",
    "list_results",
    "clear_storage",
    "migrate_legacy_storage",
    "save_workflow_result",
    "load_workflow_result",
    "export_review_summary",
    "STORAGE_ROOT",
    "VALID_MODULES",
    # Validation utilities
    "ValidationError",
    "validate_adata",
    "validate_config",
    "validate_analysis_results",
    "check_layer_consistency",
    "assert_qc_ready",
    "assert_preprocessing_ready",
    "assert_analysis_ready",
    "COMPARATIVE_READINESS_LABEL",
    "VALIDATION_SCAFFOLD_SCHEMA_VERSION",
    "VALIDATION_SCOPE",
    "build_qc_preprocess_validation",
    "validation_table_to_dataframe",
    "write_validation_outputs",
    # Audit reporting
    "export_audit_report",
    # Data loading
    "read_10x",
    "read_h5ad",
    "load_10x_data",
    "build_metadata_dicts",
    "print_sample_crosstab",
    "assess_matrix_semantics",
    # Profiling utilities
    "PerformanceStats",
    "PerformanceProfiler",
    "BenchmarkRunner",
    "profile_performance",
    "profile_function",
    "memory_tracker",
    "get_memory_usage",
    "estimate_adata_memory",
]

# Add data loader functions if available
if _data_loader_available:
    __all__.extend(
        [
            "load_pbmc3k",
            "load_luad",
            "load_melanoma",
            "load_all_datasets",
            "get_dataset_info",
            "print_dataset_summary",
            "filter_by_species",
            "filter_by_tissue_type",
        ]
    )
