"""Shared schema for unified parameter recommendations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from ..analysis.config import AnnotationConfig, ClusteringConfig


@dataclass
class ParameterRecommendation:
    """Standardized recommendation for a single tunable parameter."""

    name: str
    value: Any
    method: str
    confidence: float
    rationale: str
    ci_lower: Optional[Any] = None
    ci_upper: Optional[Any] = None
    evidence: Dict[str, Any] = field(default_factory=dict)
    alternatives: List[Any] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        return {
            "name": self.name,
            "value": self.value,
            "method": self.method,
            "confidence": self.confidence,
            "rationale": self.rationale,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "evidence": self.evidence,
            "alternatives": self.alternatives,
        }


@dataclass
class RecommendationSection:
    """Recommendation group for one analysis stage."""

    name: str
    summary: str
    confidence: float
    parameters: List[ParameterRecommendation] = field(default_factory=list)
    concerns: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    raw_result: Any = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def get_parameter(self, name: str) -> Optional[ParameterRecommendation]:
        """Return the first matching parameter recommendation."""
        return next((param for param in self.parameters if param.name == name), None)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        return {
            "name": self.name,
            "summary": self.summary,
            "confidence": self.confidence,
            "parameters": [param.to_dict() for param in self.parameters],
            "concerns": self.concerns,
            "notes": self.notes,
            "metadata": self.metadata,
        }


@dataclass
class WorkflowRecommendations:
    """Unified recommendation bundle spanning multiple pipeline stages."""

    sections: Dict[str, RecommendationSection]
    overall_confidence: float
    context: Dict[str, Any] = field(default_factory=dict)
    concerns: List[str] = field(default_factory=list)

    def get_section(self, name: str) -> Optional[RecommendationSection]:
        """Return a recommendation section by name."""
        return self.sections.get(name)

    def to_qc_thresholds(self) -> Dict[str, Any]:
        """Convert the QC section into threshold values."""
        section = self.get_section("qc")
        if section is None:
            return {}
        return {param.name: param.value for param in section.parameters}

    def to_preprocessing_config(self, base_config: Optional[Any] = None) -> Any:
        """Convert the preprocess section back to an executable config when possible."""
        section = self.get_section("preprocess")
        if section is None or section.raw_result is None:
            return None
        if hasattr(section.raw_result, "to_config"):
            return section.raw_result.to_config(base_config=base_config)
        raise TypeError("Preprocess recommendation does not expose to_config().")

    def to_clustering_config(
        self,
        *,
        method: str = "leiden",
        use_rep: str = "X_pca",
        key_added: Optional[str] = None,
        base_config: Optional[ClusteringConfig] = None,
    ) -> Optional[ClusteringConfig]:
        """Convert the clustering section into a clustering config."""
        section = self.get_section("clustering")
        if section is None:
            return None

        resolution_param = section.get_parameter("resolution")
        n_clusters_param = section.get_parameter("n_clusters")

        if base_config is None:
            config = ClusteringConfig(method=method, use_rep=use_rep, key_added=key_added)
        else:
            config = base_config.model_copy()
            config.method = method
            config.use_rep = use_rep
            config.key_added = key_added if key_added is not None else config.key_added

        if resolution_param is not None:
            config.resolution = float(resolution_param.value)

        if method == "kmeans" and n_clusters_param is not None:
            config.n_clusters = int(n_clusters_param.value)

        return config

    def to_annotation_config(
        self,
        *,
        base_config: Optional[AnnotationConfig] = None,
    ) -> Optional[AnnotationConfig]:
        """Convert the annotation section into an annotation config."""
        section = self.get_section("annotation")
        if section is None:
            return None

        if isinstance(section.raw_result, AnnotationConfig):
            return section.raw_result.model_copy()

        if base_config is None:
            config = AnnotationConfig()
        else:
            config = base_config.model_copy()

        for param in section.parameters:
            if hasattr(config, param.name):
                setattr(config, param.name, param.value)

        for key in ["cluster_key", "marker_species", "marker_tissue", "key_added"]:
            if key in section.metadata and hasattr(config, key):
                setattr(config, key, section.metadata[key])

        return config

    def to_tumor_config(
        self,
        *,
        base_config: Optional[Any] = None,
    ) -> Optional[Any]:
        """Convert the tumor section into a TumorAnalysisConfig."""
        from ..tumor.config import TumorAnalysisConfig

        section = self.get_section("tumor")
        if section is None:
            return None

        if isinstance(section.raw_result, TumorAnalysisConfig):
            return section.raw_result.model_copy()

        if base_config is None:
            config = TumorAnalysisConfig()
        else:
            config = base_config.model_copy()

        for param in section.parameters:
            if hasattr(config, param.name) and param.value is not None:
                setattr(config, param.name, param.value)

        return config

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        return {
            "sections": {name: section.to_dict() for name, section in self.sections.items()},
            "overall_confidence": self.overall_confidence,
            "context": self.context,
            "concerns": self.concerns,
        }
