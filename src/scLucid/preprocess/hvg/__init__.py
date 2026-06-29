"""Highly variable gene selection subpackage."""

from .core import PROTECTED_GENE_PRESETS, find_hvgs  # noqa: F401
from .plotting import plot_hvg_metrics  # noqa: F401
from .selection import select_and_audit_hvgs, select_hvg_sets, suggest_hvg_choice  # noqa: F401
from .stability import evaluate_hvg_stability  # noqa: F401

__all__ = [
    "PROTECTED_GENE_PRESETS",
    "find_hvgs",
    "plot_hvg_metrics",
    "suggest_hvg_choice",
    "select_hvg_sets",
    "select_and_audit_hvgs",
    "evaluate_hvg_stability",
]
