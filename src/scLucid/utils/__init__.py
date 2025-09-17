"""
Utility functions for single-cell RNA-seq data analysis.

This module provides general utility functions for plotting,
data input/output, and statistical calculations.
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

# Define what should be accessible when importing from this module
__all__ = [
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
]
