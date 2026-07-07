"""
Intelligent preprocessing submodule for scLucid.

This submodule provides data-driven parameter selection and automated
preprocessing optimization capabilities.
"""

# Configuration
from .config import IntelligentPreprocessConfig

# Data classes
from .data_classes import (
    BatchCorrectionRecommendation,
    DataProfile,
    HVGRecommendation,
    NeighborsRecommendation,
    PCARecommendation,
    PreprocessingStrategy,
    ResolutionRecommendation,
)

# Main recommender
from .recommender import (
    IntelligentPreprocessRecommender,
    recommend_intelligent_preprocessing,
    run_intelligent_preprocessing,
)

__all__ = [
    # Configuration
    "IntelligentPreprocessConfig",
    # Data classes
    "DataProfile",
    "HVGRecommendation",
    "PCARecommendation",
    "NeighborsRecommendation",
    "ResolutionRecommendation",
    "BatchCorrectionRecommendation",
    "PreprocessingStrategy",
    # Functions
    "IntelligentPreprocessRecommender",
    "recommend_intelligent_preprocessing",
    "run_intelligent_preprocessing",
]
