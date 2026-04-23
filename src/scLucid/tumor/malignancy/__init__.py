"""
Malignancy analysis module for identifying and characterizing malignant cells.

Provides tools for:
- Malignant cell scoring and classification
- Proliferation index calculation
- Cancer stem cell identification
- Metastatic potential estimation
"""

from .classification import (
    MalignancyClassifier,
    classify_malignant_cells,
    score_malignancy_potential,
)
from .scoring import (
    MalignancyScorer,
    calculate_proliferation_index,
    estimate_metastatic_potential,
    score_malignancy,
)
from .stemness import (
    StemnessAnalyzer,
    calculate_stemness_score,
    compare_stemness_between_groups,
    identify_cancer_stem_cells,
)

__all__ = [
    "score_malignancy",
    "calculate_proliferation_index",
    "estimate_metastatic_potential",
    "MalignancyScorer",
    "classify_malignant_cells",
    "score_malignancy_potential",
    "MalignancyClassifier",
    "calculate_stemness_score",
    "identify_cancer_stem_cells",
    "compare_stemness_between_groups",
    "StemnessAnalyzer",
]
