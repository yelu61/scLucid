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
from .marker_sets import flatten_marker_dict, filter_marker_dict

# Import workflow utilities
from .workflow_utils import (
    get_progress_bar,
    progress_decorator,
    WorkflowError,
    StepError,
    RecoveryError,
    WorkflowCheckpoint,
    PartialResultManager,
    WorkflowStepIterator,
    BaseWorkflow,
    with_error_recovery,
    merge_partial_results,
)

# Import validation utilities
from .validation import (
    ValidationError,
    validate_adata,
    validate_config,
    validate_analysis_results,
    check_layer_consistency,
    assert_qc_ready,
    assert_preprocessing_ready,
    assert_analysis_ready,
)

# Import profiling utilities
from .profiling import (
    PerformanceStats,
    PerformanceProfiler,
    BenchmarkRunner,
    profile_performance,
    profile_function,
    memory_tracker,
    get_memory_usage,
    estimate_adata_memory,
)

# Import storage utilities (new simplified interface)
from .storage import (
    get_storage,
    save_result,
    load_result,
    load_config,
    has_result,
    list_results,
    clear_storage,
    migrate_legacy_storage,
    save_workflow_result,
    load_workflow_result,
    STORAGE_ROOT,
    VALID_MODULES,
)

# Import and expose key functions and classes from the submodule
from .manager import (
    CellType,
    Manager,
    GeneSetManager,
    get_marker_manager,
    get_geneset_manager,
    _get_cancer_markers,
    KNOWN_SPECIES,
    MARKER_FORMATS,
)


# Import data loading utilities
try:
    from .data_loader import (
        load_pbmc3k,
        load_luad,
        load_melanoma,
        load_all_datasets,
        get_dataset_info,
        print_dataset_summary,
        filter_by_species,
        filter_by_tissue_type,
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
    __all__.extend([
        "load_pbmc3k",
        "load_luad",
        "load_melanoma",
        "load_all_datasets",
        "get_dataset_info",
        "print_dataset_summary",
        "filter_by_species",
        "filter_by_tissue_type",
    ])
