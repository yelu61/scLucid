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
