"""
Tumor evolution and phylogenetic analysis module.

Provides tools for:
- Phylogenetic tree construction
- Tumor progression trajectory analysis
- Metastasis tracking
- Clonal evolution visualization
"""

from .metastasis import (
    MetastasisTracker,
    analyze_dissemination,
    predict_metastasis_risk,
)
from .phylogeny import (
    PhylogenyBuilder,
    build_phylogenetic_tree,
    root_tree,
)
from .trajectory import (
    ProgressionAnalyzer,
    analyze_tumor_progression,
    identify_transition_states,
)

__all__ = [
    "build_phylogenetic_tree",
    "root_tree",
    "PhylogenyBuilder",
    "analyze_tumor_progression",
    "identify_transition_states",
    "ProgressionAnalyzer",
    "predict_metastasis_risk",
    "analyze_dissemination",
    "MetastasisTracker",
]
