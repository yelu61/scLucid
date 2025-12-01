"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
This module provides a robust system for loading, managing, and querying hierarchical
cell type markers from TOML or JSON files. It is central to the biology-aware
features of the scLucid toolkit.

Key exposed components:
- Manager: The main class for handling marker hierarchies.
- CellType: A dataclass representing a single cell type entry.
- get_marker_manager: A factory function to easily build a combined manager.
- KNOWN_SPECIES: A list of built-in species supported.
- MARKER_FORMATS: A list of supported marker file formats.
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
    CellType,
    Manager,
    get_marker_manager,
    KNOWN_SPECIES,
    MARKER_FORMATS,
)

from .storage_manager import (
    StorageManager,
    check_storage_status,
    cleanup_storage,
    optimize_storage,
)

# Define what should be accessible when importing from this module
__all__ = [
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
    "CellType",
    "Manager",
    "get_marker_manager",
    "KNOWN_SPECIES",
    "MARKER_FORMATS",
    "StorageManager",
    "check_storage_status",
    "cleanup_storage",
    "optimize_storage",
]
