"""
Configuration classes for intelligent preprocessing.

This module provides Pydantic-based configuration for automated
preprocessing parameter selection.
"""

from typing import Any, Dict, List, Literal, Optional

from pydantic import ConfigDict, Field, field_validator

from ...base_config import SclucidBaseConfig


class IntelligentPreprocessConfig(SclucidBaseConfig):
    """
    Configuration for intelligent preprocessing recommendations.

    Controls the statistical methods and heuristics used for
    automatic parameter selection.

    Example:
        >>> config = IntelligentPreprocessConfig(
        ...     variance_explained_threshold=0.85,
        ...     pca_method="elbow",
        ...     batch_effect_threshold=0.25
        ... )
    """

    model_config = ConfigDict(extra="ignore")

    # Variance explanation parameters for n_top_genes
    variance_explained_threshold: float = Field(
        default=0.8,
        ge=0.5,
        le=0.99,
        description="Target cumulative variance explained for HVG selection",
    )
    min_hvg_genes: int = Field(
        default=500, ge=100, description="Minimum HVGs regardless of variance"
    )
    max_hvg_genes: int = Field(
        default=10000, ge=1000, description="Maximum HVGs regardless of variance"
    )
    hvg_search_points: int = Field(
        default=8, ge=5, le=15, description="Number of HVG thresholds to test"
    )

    # PCA dimensionality parameters
    pca_method: Literal["elbow", "cumulative_variance", "knee"] = Field(
        default="elbow", description="Method for determining optimal n_pcs"
    )
    pca_variance_threshold: float = Field(
        default=0.95,
        ge=0.8,
        le=0.999,
        description="Target variance retention for PCA",
    )
    min_pcs: int = Field(default=10, ge=2, description="Minimum number of PCs")
    max_pcs: int = Field(default=100, ge=10, description="Maximum number of PCs")

    # Neighbors optimization parameters (integrating with existing neighbors.py)
    optimize_neighbors: bool = Field(
        default=True, description="Whether to optimize n_neighbors and n_pcs"
    )
    neighbors_search_space: List[int] = Field(
        default_factory=lambda: [5, 10, 15, 20, 30, 50]
    )
    silhouette_sample_size: int = Field(
        default=10000, ge=1000, description="Sample size for silhouette calculation"
    )

    # Resolution optimization parameters
    resolution_search_space: List[float] = Field(
        default_factory=lambda: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
    )
    resolution_stability_n: int = Field(
        default=5, ge=3, le=10, description="Number of runs for stability assessment"
    )

    # Batch effect detection parameters
    batch_effect_threshold: float = Field(
        default=0.2,
        ge=0.0,
        le=1.0,
        description="kBET rejection rate threshold for triggering batch correction",
    )

    # Bootstrap parameters
    n_bootstrap: int = Field(
        default=20, ge=10, le=100, description="Bootstrap iterations for CI"
    )
    confidence_level: float = Field(
        default=0.95, ge=0.8, le=0.99, description="Confidence level for intervals"
    )

    # Data-aware strategy parameters
    sparse_data_threshold: float = Field(
        default=0.9, description="Fraction of zeros to classify data as 'sparse'"
    )
    small_dataset_threshold: int = Field(
        default=1000, description="Cell count threshold for 'small' datasets"
    )
    large_dataset_threshold: int = Field(
        default=50000, description="Cell count threshold for 'large' datasets"
    )

    @field_validator("resolution_search_space")
    @classmethod
    def validate_resolution_search_space(cls, v: List[float]) -> List[float]:
        """Ensure resolution search space is valid."""
        if not v:
            raise ValueError("resolution_search_space cannot be empty")
        if not all(0 < r <= 3.0 for r in v):
            raise ValueError("All resolutions must be in (0, 3.0]")
        return sorted(v)

    @field_validator("max_hvg_genes")
    @classmethod
    def validate_hvg_bounds(cls, v: int, info) -> int:
        """Ensure max_hvg_genes > min_hvg_genes."""
        values = info.data
        if "min_hvg_genes" in values and v < values["min_hvg_genes"]:
            raise ValueError("max_hvg_genes must be >= min_hvg_genes")
        return v
