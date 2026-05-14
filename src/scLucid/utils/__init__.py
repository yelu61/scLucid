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

# Import and expose key functions from submodules
from .context import (
    AnalysisContext,
    DatasetProfile,
    DatasetType,
    infer_analysis_context,
    infer_dataset_profile,
    is_multi_sample_hint,
    normalize_dataset_type,
)
from .contracts import (
    API_LAYER_CONTRACTS,
    API_LAYER_ORDER,
    SCHEMA_VERSION,
    SCLUCID_ROOT,
    STAGE_CONTRACTS,
    STAGE_ORDER,
    APILayerContract,
    ContractError,
    ContractValidationResult,
    LayerKeys,
    MINIMAL_WORKFLOW_CONTRACT,
    Modules,
    ObsKeys,
    ObsmKeys,
    REVIEW_SUMMARY_RECOMMENDED_KEYS,
    REVIEW_SUMMARY_REQUIRED_KEYS,
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
    module_namespace,
    normalize_review_summary,
    record_artifact,
    record_contract_result,
    record_config_lineage,
    record_error,
    stage_contract_to_dict,
    validate_all_stage_contracts,
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
    load_10x_data,
    merge_obs_metadata,
    sanitize_for_hdf5,
    subset_adata,
    subset_from_annotations,
    use_layer_as_X,
)

# Import and expose key functions and classes from the submodule
from .manager import (
    KNOWN_SPECIES,
    MARKER_FORMATS,
    CellType,
    Manager,
    _get_cancer_markers,
    get_marker_manager,
    load_gene_set_manager,
    load_gene_sets,
)
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

# Import result cleanup utilities
from .result_cleanup import (
    clear_sclucid_results,
    list_sclucid_modules,
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
)
from .audit_report import export_audit_report
from .helpers import read_10x
from .io import read_h5ad
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
    # Helper functions
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
    "AnalysisContext",
    "DatasetProfile",
    "DatasetType",
    "infer_analysis_context",
    "infer_dataset_profile",
    "is_multi_sample_hint",
    "normalize_dataset_type",
    "API_LAYER_CONTRACTS",
    "API_LAYER_ORDER",
    "APILayerContract",
    "ContractError",
    "ContractValidationResult",
    "LayerKeys",
    "MINIMAL_WORKFLOW_CONTRACT",
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
    "module_namespace",
    "normalize_review_summary",
    "record_artifact",
    "record_contract_result",
    "record_config_lineage",
    "record_error",
    "stage_contract_to_dict",
    "validate_all_stage_contracts",
    "validate_review_summary_schema",
    "validate_stage_contract",
    "effective_n_jobs",
    "is_ci_environment",
    "run_joblib_or_sequential",
    "setup_runtime_environment",
    "flatten_marker_dict",
    "filter_marker_dict",
    "DecisionRecord",
    "EVIDENCE_SCHEMA_VERSION",
    "EvidenceBundle",
    "EvidenceItem",
    "ReviewAction",
    "model_to_dict",
    # Result cleanup
    "clear_sclucid_results",
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
    "get_marker_manager",
    "load_gene_set_manager",
    "load_gene_sets",
    "_get_cancer_markers",
    "KNOWN_SPECIES",
    "MARKER_FORMATS",
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
