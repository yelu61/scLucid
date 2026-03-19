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
    analyze_cell_interactions,
    calculate_communication_strength,
    InteractionAnalyzer,
)

from .ecosystem import (
    score_tme_ecosystem,
    identify_ecosystem_types,
    EcosystemClassifier,
)

__all__ = [
    "deconvolve_tme",
    "estimate_stromal_content",
    "analyze_immune_infiltration",
    "TMEProfiler",
    "analyze_cell_interactions",
    "calculate_communication_strength",
    "InteractionAnalyzer",
    "score_tme_ecosystem",
    "identify_ecosystem_types",
    "EcosystemClassifier",
]
