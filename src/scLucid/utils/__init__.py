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

# Import and expose key functions from submodules
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
    GeneSetManager,
    Manager,
    _get_cancer_markers,
    get_geneset_manager,
    get_marker_manager,
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
        filter_by_species,
        filter_by_tissue_type,
        get_dataset_info,
        load_all_datasets,
        load_luad,
        load_melanoma,
        load_pbmc3k,
        print_dataset_summary,
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
    "flatten_marker_dict",
    "filter_marker_dict",
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
    "GeneSetManager",
    "get_marker_manager",
    "get_geneset_manager",
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
