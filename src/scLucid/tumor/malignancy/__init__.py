"""
Malignancy analysis module for identifying and characterizing malignant cells.

Provides tools for:
- Malignant cell scoring and classification
- Proliferation index calculation
- Cancer stem cell identification
- Metastatic potential estimation
"""

from .scoring import (
    score_malignancy,
    calculate_proliferation_index,
    estimate_metastatic_potential,
    MalignancyScorer,
)

from .classification import (
    classify_malignant_cells,
    score_malignancy_potential,
    MalignancyClassifier,
)

from .stemness import (
    calculate_stemness_score,
    identify_cancer_stem_cells,
    compare_stemness_between_groups,
    StemnessAnalyzer,
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
