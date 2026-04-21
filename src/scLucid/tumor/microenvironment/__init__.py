"""
Tumor Microenvironment (TME) analysis module.

Provides tools for:
- TME cell type deconvolution
- Cell-cell interaction analysis
- Ecosystem characterization
- Immune infiltration assessment
"""

from .deconvolution import (
    deconvolve_tme,
    estimate_stromal_content,
    analyze_immune_infiltration,
    TMEProfiler,
)

from .interaction import (
    InteractionAnalyzer,
    analyze_cell_interactions,
    find_dominant_interactions,
    score_immune_interactions,
)

from .ecosystem import (
    EcosystemAnalyzer,
    analyze_ecosystem_composition,
    calculate_tumor_microenvironment_score,
    compare_ecosystems,
)

__all__ = [
    "deconvolve_tme",
    "estimate_stromal_content",
    "analyze_immune_infiltration",
    "TMEProfiler",
    "InteractionAnalyzer",
    "analyze_cell_interactions",
    "find_dominant_interactions",
    "score_immune_interactions",
    "EcosystemAnalyzer",
    "analyze_ecosystem_composition",
    "calculate_tumor_microenvironment_score",
    "compare_ecosystems",
]
