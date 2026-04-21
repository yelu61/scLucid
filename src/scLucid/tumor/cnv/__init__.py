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
    infer_clonal_phylogeny,
    calculate_clonal_diversity,
    identify_clones as identify_clones_from_cnv,
    CloneAnalyzer,
)

from .cnv_signature import (
    extract_cnv_signatures,
    assign_cnv_signature,
    CNVSigExtractor,
)

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
]
