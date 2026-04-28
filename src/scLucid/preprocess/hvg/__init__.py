"""Highly variable gene selection subpackage."""

from .core import *  # noqa: F401, F403

# Re-export symbols that moved to focused submodules.
from .core import (  # noqa: F401
    _exclude_genes,
    _gene_type_detection,
    _get_hvg_input_matrix,
    _infer_species_from_gene_names,
    _validate_hvg_input_matrix,
    find_hvgs,
)
from .plotting import plot_hvg_metrics  # noqa: F401
from .selection import select_hvg_sets, suggest_hvg_choice  # noqa: F401
from .stability import evaluate_hvg_stability  # noqa: F401
