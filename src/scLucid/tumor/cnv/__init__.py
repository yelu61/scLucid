"""
Copy Number Variation (CNV) analysis for tumor single-cell data.

This module provides tools for:
- Inferring CNV from scRNA-seq data
- Identifying tumor cells based on CNV patterns
- Clonal evolution analysis
- CNV signature extraction
"""

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

try:
    from .infercnvpy import find_tumor, run_cnv_analysis
except ImportError:
    find_tumor = None
    run_cnv_analysis = None

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
]

if run_cnv_analysis is not None:
    __all__.extend(["run_cnv_analysis", "find_tumor"])
