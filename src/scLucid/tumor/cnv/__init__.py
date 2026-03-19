"""
Copy Number Variation (CNV) analysis for tumor single-cell data.

This module provides tools for:
- Inferring CNV from scRNA-seq data
- Identifying tumor cells based on CNV patterns
- Clonal evolution analysis
- CNV signature extraction
"""

from .infercnv import (
    infer_cnv,
    find_tumor_cells,
    identify_clones,
    calculate_cnv_score,
    CNVAnalyzer,
)

from .clone_analysis import (
    construct_phylogeny,
    calculate_clonal_diversity,
    track_clonal_evolution,
    CloneTracker,
)

from .cnv_signature import (
    extract_cnv_signatures,
    classify_cnv_pattern,
    CNVDictionary,
)

__all__ = [
    "infer_cnv",
    "find_tumor_cells",
    "identify_clones",
    "calculate_cnv_score",
    "CNVAnalyzer",
    "construct_phylogeny",
    "calculate_clonal_diversity",
    "track_clonal_evolution",
    "CloneTracker",
    "extract_cnv_signatures",
    "classify_cnv_pattern",
    "CNVDictionary",
]
