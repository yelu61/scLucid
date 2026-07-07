"""Unified recommendation API for cross-stage parameter selection."""

from .config import RecommendationConfig
from .engine import RecommendationEngine, recommend_analysis_parameters
from .preprocess import (
    BatchCorrectionRecommendation,
    DataProfile,
    HVGRecommendation,
    IntelligentPreprocessConfig,
    IntelligentPreprocessRecommender,
    NeighborsRecommendation,
    PCARecommendation,
    PreprocessingStrategy,
    ResolutionRecommendation,
    recommend_intelligent_preprocessing,
    run_intelligent_preprocessing,
)
from .schema import (
    ParameterRecommendation,
    RecommendationSection,
    WorkflowRecommendations,
)

__all__ = [
    "RecommendationConfig",
    "RecommendationEngine",
    "recommend_analysis_parameters",
    "IntelligentPreprocessConfig",
    "IntelligentPreprocessRecommender",
    "PreprocessingStrategy",
    "DataProfile",
    "HVGRecommendation",
    "PCARecommendation",
    "NeighborsRecommendation",
    "ResolutionRecommendation",
    "BatchCorrectionRecommendation",
    "recommend_intelligent_preprocessing",
    "run_intelligent_preprocessing",
    "ParameterRecommendation",
    "RecommendationSection",
    "WorkflowRecommendations",
]
