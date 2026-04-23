"""
Data classes for intelligent preprocessing recommendations.
"""

import copy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional

import numpy as np
import pandas as pd


@dataclass
class DataProfile:
    """
    Characterization of input data for strategy selection.
    """

    # Basic characteristics
    n_cells: int
    n_genes: int
    sparsity: float
    median_counts_per_cell: float
    median_genes_per_cell: float

    # Classification flags
    is_sparse: bool
    is_small_dataset: bool
    is_medium_dataset: bool
    is_large_dataset: bool
    has_batch_info: bool = False
    n_batches: Optional[int] = None

    # Quality indicators
    data_quality_score: float = 0.0
    potential_issues: List[str] = field(default_factory=list)

    # Suggested strategy type
    strategy_type: Literal["minimal", "standard", "aggressive", "large_scale"] = "standard"

    @classmethod
    def from_adata(
        cls, adata: Any, batch_key: Optional[str] = None, config: Optional[Any] = None
    ) -> "DataProfile":
        """Analyze AnnData and create profile."""
        # Handle default config
        if config is None:
            sparse_threshold = 0.9
            small_threshold = 1000
            large_threshold = 50000
        else:
            sparse_threshold = config.sparse_data_threshold
            small_threshold = config.small_dataset_threshold
            large_threshold = config.large_dataset_threshold

        # Calculate sparsity
        if hasattr(adata.X, "nnz"):
            sparsity = 1 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
        else:
            sparsity = np.mean(adata.X == 0)

        # Count statistics
        total_counts = adata.obs.get("total_counts", np.array([0]))
        n_genes = adata.obs.get("n_genes_by_counts", np.array([0]))
        median_counts = np.median(total_counts) if len(total_counts) > 0 else 0
        median_genes = np.median(n_genes) if len(n_genes) > 0 else 0

        # Classify dataset size
        is_small = adata.n_obs < small_threshold
        is_large = adata.n_obs > large_threshold
        is_medium = not is_small and not is_large

        # Batch info
        has_batch = batch_key is not None and batch_key in adata.obs.columns
        n_batches = adata.obs[batch_key].nunique() if has_batch else None

        # Determine strategy
        if is_large:
            strategy = "large_scale"
        elif sparsity > sparse_threshold:
            strategy = "aggressive"
        elif is_small and not has_batch:
            strategy = "minimal"
        else:
            strategy = "standard"

        # Quality score
        quality_score = 100.0
        issues = []

        if median_genes < 200:
            quality_score -= 20
            issues.append("Low median genes per cell")
        if sparsity > 0.95:
            quality_score -= 10
            issues.append("Very sparse data")
        if adata.n_obs < 100:
            quality_score -= 30
            issues.append("Very small dataset")

        return cls(
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            sparsity=sparsity,
            median_counts_per_cell=median_counts,
            median_genes_per_cell=median_genes,
            is_sparse=sparsity > sparse_threshold,
            is_small_dataset=is_small,
            is_medium_dataset=is_medium,
            is_large_dataset=is_large,
            has_batch_info=has_batch,
            n_batches=n_batches,
            data_quality_score=max(0, quality_score),
            potential_issues=issues,
            strategy_type=strategy,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "sparsity": self.sparsity,
            "median_counts_per_cell": self.median_counts_per_cell,
            "median_genes_per_cell": self.median_genes_per_cell,
            "is_sparse": self.is_sparse,
            "is_small_dataset": self.is_small_dataset,
            "is_large_dataset": self.is_large_dataset,
            "has_batch_info": self.has_batch_info,
            "n_batches": self.n_batches,
            "data_quality_score": self.data_quality_score,
            "potential_issues": self.potential_issues,
            "strategy_type": self.strategy_type,
        }


@dataclass
class HVGRecommendation:
    """Data-driven HVG selection recommendation."""

    n_top_genes: int
    variance_explained: float
    ci_lower: int
    ci_upper: int
    method: str
    confidence: float
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_top_genes": self.n_top_genes,
            "variance_explained": self.variance_explained,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "method": self.method,
            "confidence": self.confidence,
            "evidence": self.evidence,
        }


@dataclass
class PCARecommendation:
    """Data-driven PCA dimensionality recommendation."""

    n_pcs: int
    variance_explained: float
    ci_lower: int
    ci_upper: int
    method: str
    confidence: float
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_pcs": self.n_pcs,
            "variance_explained": self.variance_explained,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "method": self.method,
            "confidence": self.confidence,
            "evidence": self.evidence,
        }


@dataclass
class NeighborsRecommendation:
    """Data-driven neighbors and PCs recommendation."""

    n_neighbors: int
    n_pcs: int
    silhouette_score: float
    ci_lower_neighbors: int
    ci_upper_neighbors: int
    ci_lower_pcs: int
    ci_upper_pcs: int
    method: str
    confidence: float
    search_results: Optional[pd.DataFrame] = None
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_neighbors": self.n_neighbors,
            "n_pcs": self.n_pcs,
            "silhouette_score": self.silhouette_score,
            "ci_lower_neighbors": self.ci_lower_neighbors,
            "ci_upper_neighbors": self.ci_upper_neighbors,
            "ci_lower_pcs": self.ci_lower_pcs,
            "ci_upper_pcs": self.ci_upper_pcs,
            "method": self.method,
            "confidence": self.confidence,
            "evidence": self.evidence,
        }


@dataclass
class ResolutionRecommendation:
    """Data-driven clustering resolution recommendation."""

    resolution: float
    n_clusters: int
    stability_score: float
    ci_lower: float
    ci_upper: float
    method: str
    confidence: float
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "resolution": self.resolution,
            "n_clusters": self.n_clusters,
            "stability_score": self.stability_score,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "method": self.method,
            "confidence": self.confidence,
            "evidence": self.evidence,
        }


@dataclass
class BatchCorrectionRecommendation:
    """Batch correction method recommendation."""

    needs_correction: bool
    severity_score: float
    recommended_method: Optional[str]
    alternative_methods: List[str] = field(default_factory=list)
    method_scores: Dict[str, float] = field(default_factory=dict)
    confidence: float = 0.0
    evidence: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "needs_correction": self.needs_correction,
            "severity_score": self.severity_score,
            "recommended_method": self.recommended_method,
            "alternative_methods": self.alternative_methods,
            "method_scores": self.method_scores,
            "confidence": self.confidence,
            "evidence": self.evidence,
        }


@dataclass
class PreprocessingStrategy:
    """Complete preprocessing strategy recommendation."""

    data_profile: DataProfile
    hvg: HVGRecommendation
    pca: PCARecommendation
    neighbors: NeighborsRecommendation
    resolution: ResolutionRecommendation
    batch_correction: Optional[BatchCorrectionRecommendation] = None
    overall_confidence: float = 0.0
    concerns: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "data_profile": self.data_profile.to_dict(),
            "hvg": self.hvg.to_dict(),
            "pca": self.pca.to_dict(),
            "neighbors": self.neighbors.to_dict(),
            "resolution": self.resolution.to_dict(),
            "batch_correction": self.batch_correction.to_dict() if self.batch_correction else None,
            "overall_confidence": self.overall_confidence,
            "concerns": self.concerns,
            "recommendations": self.recommendations,
        }

    def to_config(self, base_config: Optional[Any] = None) -> Any:
        """
        Convert recommendations to executable PreprocessingWorkflowConfig.

        Args:
            base_config: Optional base config to override

        Returns:
            PreprocessingWorkflowConfig with recommended parameters
        """
        from ..config import PreprocessingWorkflowConfig

        if base_config is None:
            config = PreprocessingWorkflowConfig()
        elif hasattr(base_config, "model_copy"):
            config = base_config.model_copy(deep=True)
        else:
            config = copy.deepcopy(base_config)

        # Apply HVG recommendation
        config.hvg = config.hvg.model_copy(update={"n_top_genes": self.hvg.n_top_genes})

        # Apply PCA/Graph recommendation
        config.graph = config.graph.model_copy(
            update={
                "n_pcs": self.neighbors.n_pcs,
                "n_neighbors": self.neighbors.n_neighbors,
            }
        )

        # Apply resolution recommendation
        config.scaling = config.scaling.model_copy()

        # Store resolution in uns for later use
        if not hasattr(config, "_resolution"):
            config._resolution = self.resolution.resolution

        # Apply batch correction recommendation
        if self.batch_correction and self.batch_correction.needs_correction:
            config.run_integration = True
            batch_key = self.batch_correction.evidence.get("batch_key")
            config.integration = config.integration.model_copy(
                update={
                    "method": self.batch_correction.recommended_method or "harmony",
                    "batch_key": batch_key or config.integration.batch_key,
                }
            )
        else:
            config.run_integration = False

        return config

    def to_review_summary(self) -> Dict[str, Any]:
        """
        Build a human-reviewable summary of the preprocessing recommendation.

        Returns:
            Structured dict with data profile, key recommendations,
            batch-correction decision, and any concerns.
        """
        summary: Dict[str, Any] = {
            "data_profile": self.data_profile.to_dict(),
            "overall_confidence": round(self.overall_confidence, 3),
            "concerns": self.concerns,
            "recommendations": self.recommendations,
        }

        summary["hvg"] = {
            "n_top_genes": self.hvg.n_top_genes,
            "variance_explained": round(self.hvg.variance_explained, 3),
            "confidence": round(self.hvg.confidence, 3),
            "ci": [self.hvg.ci_lower, self.hvg.ci_upper],
        }

        summary["pca"] = {
            "n_pcs": self.pca.n_pcs,
            "variance_explained": round(self.pca.variance_explained, 3),
            "confidence": round(self.pca.confidence, 3),
            "ci": [self.pca.ci_lower, self.pca.ci_upper],
        }

        summary["neighbors"] = {
            "n_neighbors": self.neighbors.n_neighbors,
            "n_pcs": self.neighbors.n_pcs,
            "silhouette_score": round(self.neighbors.silhouette_score, 3),
            "confidence": round(self.neighbors.confidence, 3),
        }

        summary["resolution"] = {
            "resolution": round(self.resolution.resolution, 2),
            "expected_clusters": self.resolution.n_clusters,
            "stability_score": round(self.resolution.stability_score, 3),
            "confidence": round(self.resolution.confidence, 3),
        }

        if self.batch_correction is not None:
            summary["batch_correction"] = {
                "needs_correction": self.batch_correction.needs_correction,
                "severity_score": round(self.batch_correction.severity_score, 3),
                "recommended_method": self.batch_correction.recommended_method,
                "alternative_methods": self.batch_correction.alternative_methods,
                "confidence": round(self.batch_correction.confidence, 3),
            }
        else:
            summary["batch_correction"] = {
                "needs_correction": False,
                "severity_score": 0.0,
                "recommended_method": None,
                "note": "No batch_key provided; batch assessment skipped.",
            }

        return summary
