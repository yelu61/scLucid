"""Configuration for the unified recommendation engine."""

from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import ConfigDict, Field, field_validator

from ..analysis.config import AnnotationConfig, ResolutionSearchConfig
from ..base_config import SclucidBaseConfig
from ..preprocess.intelligent.config import IntelligentPreprocessConfig
from ..qc.intelligent_qc import IntelligentQCConfig
from ..tumor.config import TumorAnalysisConfig


RecommendationModule = Literal["qc", "preprocess", "clustering", "annotation", "tumor"]


class RecommendationConfig(SclucidBaseConfig):
    """Controls which recommendation stages run and how they are configured."""

    model_config = ConfigDict(extra="ignore")

    modules: List[RecommendationModule] = Field(
        default_factory=lambda: ["qc", "preprocess", "clustering", "annotation"]
    )
    qc: Optional[IntelligentQCConfig] = Field(default=None)
    preprocess: Optional[IntelligentPreprocessConfig] = Field(default=None)
    annotation: Optional[AnnotationConfig] = Field(default=None)
    tumor: Optional[TumorAnalysisConfig] = Field(default=None)
    resolution_search: Optional[ResolutionSearchConfig] = Field(default=None)
    clustering_selection_strategy: Literal["elbow", "peak", "balanced"] = Field(
        default="balanced"
    )
    clustering_method: Literal["leiden", "louvain"] = Field(default="leiden")
    clustering_use_rep: str = Field(default="X_pca")
    prepare_clustering_rep_if_missing: bool = Field(default=True)

    @field_validator("modules")
    @classmethod
    def validate_modules(cls, value: List[RecommendationModule]) -> List[RecommendationModule]:
        """Ensure modules are unique and ordered."""
        if not value:
            raise ValueError("At least one recommendation module must be selected.")
        seen = set()
        ordered = []
        for item in value:
            if item not in seen:
                ordered.append(item)
                seen.add(item)
        return ordered
