"""Unified recommendation engine (split into submodules for maintainability)."""

from .core import RecommendationEngine, recommend_analysis_parameters

__all__ = [
    "RecommendationEngine",
    "recommend_analysis_parameters",
]
