"""
Therapy response and resistance analysis module.

Provides tools for:
- Drug resistance mechanism identification
- Therapy response prediction
- Therapeutic target discovery
- Patient stratification
"""

from .resistance import (
    identify_resistance_mechanisms,
    score_drug_resistance,
    ResistanceAnalyzer,
)

from .prediction import (
    predict_therapy_response,
    stratify_patients,
    ResponsePredictor,
)

from .target import (
    discover_therapeutic_targets,
    prioritize_druggable_genes,
    TargetDiscovery,
)

__all__ = [
    "identify_resistance_mechanisms",
    "score_drug_resistance",
    "ResistanceAnalyzer",
    "predict_therapy_response",
    "stratify_patients",
    "ResponsePredictor",
    "discover_therapeutic_targets",
    "prioritize_druggable_genes",
    "TargetDiscovery",
]
