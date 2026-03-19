"""
Tumor evolution and phylogenetic analysis module.

Provides tools for:
- Phylogenetic tree construction
- Tumor progression trajectory analysis
- Metastasis tracking
- Clonal evolution visualization
"""

from .phylogeny import (
    build_phylogenetic_tree,
    root_tree,
    PhylogenyBuilder,
)

from .trajectory import (
    analyze_tumor_progression,
    identify_transition_states,
    ProgressionAnalyzer,
)

from .metastasis import (
    predict_metastasis_risk,
    analyze_dissemination,
    MetastasisTracker,
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
