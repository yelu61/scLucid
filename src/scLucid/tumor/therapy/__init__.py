"""
Therapy response and resistance analysis module.

Provides tools for:
- Drug resistance mechanism identification
- Therapy response prediction
- Therapeutic target discovery
- Patient stratification
"""

from .prediction import (
    ResponsePredictor,
    predict_therapy_response,
    stratify_patients,
)
from .resistance import (
    ResistanceAnalyzer,
    identify_resistance_mechanisms,
    score_drug_resistance,
)
from .target import (
    TargetDiscovery,
    discover_therapeutic_targets,
    prioritize_druggable_genes,
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
