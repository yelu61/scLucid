"""
Copy Number Variation (CNV) analysis for tumor single-cell data.

This module provides tools for:
- Inferring CNV from scRNA-seq data
- Identifying tumor cells based on CNV patterns
- Clonal evolution analysis
- CNV signature extraction
"""

from importlib import import_module
from typing import Any

from .clone_analysis import (
    CloneAnalyzer,
    calculate_clonal_diversity,
    infer_clonal_phylogeny,
)
from .clone_analysis import (
    identify_clones as identify_clones_from_cnv,
)
from .cnv_signature import (
    CNVSigExtractor,
    assign_cnv_signature,
    extract_cnv_signatures,
)
from .infercnv import (
    CNVAnalyzer,
    calculate_cnv_score,
    find_tumor_cells,
    identify_clones,
    infer_cnv,
    plot_aneuploid_proportion,
    plot_cnv_distribution,
    plot_cnv_heatmap,
    plot_per_chromosome_scores,
)


def _load_infercnvpy_backend():
    """Load the optional infercnvpy adapter only when it is called."""
    try:
        return import_module(f"{__name__}.infercnvpy")
    except ModuleNotFoundError as exc:
        if exc.name == "infercnvpy":
            raise ImportError(
                "run_cnv_analysis requires the optional 'infercnvpy' package. "
                "Install scLucid with the CNV extra before calling this backend."
            ) from exc
        raise


def run_cnv_analysis(*args: Any, **kwargs: Any):
    """Run the optional infercnvpy backend without loading it at core import time."""
    return _load_infercnvpy_backend().run_cnv_analysis(*args, **kwargs)


def find_tumor(*args: Any, **kwargs: Any):
    """Call the optional infercnvpy tumor classifier lazily."""
    return _load_infercnvpy_backend().find_tumor(*args, **kwargs)

__all__ = [
    "infer_cnv",
    "find_tumor_cells",
    "identify_clones",
    "calculate_cnv_score",
    "CNVAnalyzer",
    "infer_clonal_phylogeny",
    "calculate_clonal_diversity",
    "identify_clones_from_cnv",
    "CloneAnalyzer",
    "extract_cnv_signatures",
    "assign_cnv_signature",
    "CNVSigExtractor",
    "plot_cnv_distribution",
    "plot_cnv_heatmap",
    "plot_per_chromosome_scores",
    "plot_aneuploid_proportion",
    "run_cnv_analysis",
    "find_tumor",
]
