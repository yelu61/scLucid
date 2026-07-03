"""Doublet detection."""

from .core import (
    ALGORITHM_PRED_COL,
    ALGORITHM_SCORE_COL,
    COMBINED_SCORE_COL,
    EXPECTED_HETEROTYPIC_RATE_COL,
    EXPECTED_HOMOTYPIC_RATE_COL,
    EXPECTED_TOTAL_RATE_COL,
    FINAL_PRED_COL,
    HETEROTYPIC_RISK_COL,
    HEURISTIC_PRED_COL,
    HEURISTIC_SCORE_COL,
    HOMOTYPIC_RISK_COL,
    LINEAGE_SCORES_KEY,
    DOUBLET_OBS_COLUMNS,
    audit_doublets,
    create_custom_marker_dict,
    generate_doublet_rates,
)
from .ensemble import predict_doublets
from .profiler import DoubletEvidenceProfiler, predict_doublets_with_profiling

__all__ = [
    "ALGORITHM_PRED_COL",
    "ALGORITHM_SCORE_COL",
    "COMBINED_SCORE_COL",
    "EXPECTED_HETEROTYPIC_RATE_COL",
    "EXPECTED_HOMOTYPIC_RATE_COL",
    "EXPECTED_TOTAL_RATE_COL",
    "FINAL_PRED_COL",
    "HETEROTYPIC_RISK_COL",
    "HEURISTIC_PRED_COL",
    "HEURISTIC_SCORE_COL",
    "HOMOTYPIC_RISK_COL",
    "LINEAGE_SCORES_KEY",
    "DOUBLET_OBS_COLUMNS",
    "generate_doublet_rates",
    "create_custom_marker_dict",
    "audit_doublets",
    "predict_doublets",
    "DoubletEvidenceProfiler",
    "predict_doublets_with_profiling",
]
