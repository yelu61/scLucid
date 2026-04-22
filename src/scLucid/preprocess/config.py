"""
Pydantic-based configuration for the scLucid preprocessing module.

Migrates from dataclasses to Pydantic for consistent validation and serialization.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import ConfigDict, Field, field_validator, model_validator

from ..base_config import SclucidBaseConfig, WorkflowConfigBase, apply_config_overrides as _apply_config_overrides

logger = logging.getLogger(__name__)

# Re-export for backward compatibility within preprocess submodules
apply_config_overrides = _apply_config_overrides


class NormalizationConfig(SclucidBaseConfig):
    """Configuration for data normalization."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["standard", "scran", "pearson_residuals", "clr"] = Field(default="standard")
    target_sum: float = Field(default=1e4, gt=0, description="Target sum for normalization")
    exclude_highly_expressed: bool = Field(default=False)
    max_fraction: float = Field(default=0.05, gt=0, lt=1)
    input_layer: str = Field(default="counts", description="Input layer name")
    output_layer: str = Field(default="normalized", description="Output layer name")
    update_X: bool = Field(default=True, description="Update adata.X with normalized data")

    @field_validator("output_layer")
    @classmethod
    def validate_output_layer(cls, v: str) -> str:
        """Prevent reserved names for output layer."""
        reserved = {"X", "raw"}
        if v in reserved:
            raise ValueError(f"output_layer cannot be one of reserved names: {reserved}")
        return v

    @model_validator(mode="after")
    def validate_method_params(self) -> "NormalizationConfig":
        """Validate method-specific constraints."""
        if self.method == "pearson_residuals" and self.target_sum != 1e4:
            logger.warning("Pearson residuals normalization ignores target_sum parameter")
        return self


class HVGConfig(SclucidBaseConfig):
    """Configuration for HVG selection."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["scanpy", "custom", "triku"] = Field(default="scanpy")
    n_top_genes: int = Field(default=2000, ge=100, le=20000, description="Number of HVGs to select")
    flavor: Literal["seurat", "seurat_v3", "cell_ranger"] = Field(default="seurat")
    span: Optional[float] = Field(default=0.3)
    batch_key: Optional[str] = Field(default=None)
    sample_key: str = Field(default="sampleID")
    min_n_samples: int = Field(default=2, ge=1)
    n_highly_expressed_genes: int = Field(default=50, ge=0)
    n_specific_genes: int = Field(default=20, ge=0)
    exclude_gene_types: Optional[List[str]] = Field(
        default_factory=lambda: ["mitochondrial", "ribosomal"]
    )

    @field_validator("span")
    @classmethod
    def validate_span(cls, v: Optional[float]) -> Optional[float]:
        """Validate span parameter range."""
        if v is not None and not (0.01 < v < 1):
            raise ValueError(f"span must be in (0.01, 1), got {v}")
        return v

    @field_validator("n_top_genes")
    @classmethod
    def warn_n_top_genes(cls, v: int) -> int:
        """Warn about unusual n_top_genes values."""
        if v < 500:
            logger.warning(f"n_top_genes={v} is very low. Typical values are 2000-5000.")
        elif v > 10000:
            logger.warning(f"n_top_genes={v} is very high. This may include too much noise.")
        return v


class ScalingConfig(SclucidBaseConfig):
    """Configuration for data scaling and regression."""

    model_config = ConfigDict(extra="ignore")

    vars_to_regress: Optional[List[str]] = Field(
        default_factory=lambda: ["total_counts", "pct_counts_mt"]
    )
    regress_in_scale: bool = Field(default=True)
    vars_to_regress_in_scale: Optional[List[str]] = Field(default=None)
    input_layer_for_regress: str = Field(default="normalized")
    scale_method: Literal["zscore", "robust", "minmax"] = Field(default="zscore")
    max_value: Optional[float] = Field(default=10.0, gt=0)


class IntegrationConfig(SclucidBaseConfig):
    """Configuration for batch correction and data integration."""

    model_config = ConfigDict(extra="ignore")

    method: Optional[Literal["harmony", "scanorama", "scvi", "bbknn", "combat"]] = Field(
        default="harmony"
    )
    batch_key: Optional[Union[str, List[str]]] = Field(default="sampleID")
    use_rep: str = Field(default="X_pca")
    output_key: Optional[str] = Field(default=None)
    harmony_params: Dict[str, Any] = Field(
        default_factory=lambda: {"max_iter_harmony": 20, "theta": 2.0}
    )
    scvi_params: Dict[str, Any] = Field(
        default_factory=lambda: {"n_latent": 30, "max_epochs": 500}
    )
    hvg_key: Optional[str] = Field(default=None, description="For Scanorama")


class NeighborsConfig(SclucidBaseConfig):
    """Configuration for optimizing nearest neighbor and PCA parameters."""

    model_config = ConfigDict(extra="ignore")

    n_neighbors_list: List[int] = Field(default_factory=lambda: [15, 30, 50])
    n_pcs_list: List[int] = Field(default_factory=lambda: [30, 40, 50])
    use_rep: str = Field(default="X_pca")
    clustering_method: Literal["leiden", "louvain"] = Field(default="leiden")
    resolution: float = Field(default=1.0, gt=0)
    subsample: Optional[int] = Field(default=5000, ge=100)
    n_jobs: int = Field(default=-1)

    @field_validator("n_neighbors_list", "n_pcs_list")
    @classmethod
    def validate_positive_list(cls, v: List[int]) -> List[int]:
        """Ensure list contains positive integers."""
        if not all(isinstance(x, int) and x > 0 for x in v):
            raise ValueError("List must contain positive integers")
        return v


class GraphConfig(SclucidBaseConfig):
    """Configuration for the final neighborhood graph and UMAP embedding."""

    model_config = ConfigDict(extra="ignore")

    n_pcs: int = Field(default=50, ge=2, le=100, description="Number of PCs for neighbors/UMAP")
    n_neighbors: int = Field(default=15, ge=3, le=100, description="Number of neighbors")


class PreprocessingWorkflowConfig(WorkflowConfigBase):
    """Master configuration for the entire preprocessing workflow."""

    model_config = ConfigDict(extra="ignore")

    # Layer naming
    counts_layer: str = Field(default="counts")
    normalized_layer: str = Field(default="normalized")
    regressed_layer: str = Field(default="regressed")
    scaled_layer: str = Field(default="scaled")

    # Sub-configs
    normalization: NormalizationConfig = Field(default_factory=NormalizationConfig)
    hvg: HVGConfig = Field(default_factory=HVGConfig)
    scaling: ScalingConfig = Field(default_factory=ScalingConfig)
    integration: IntegrationConfig = Field(default_factory=IntegrationConfig)
    graph: GraphConfig = Field(default_factory=GraphConfig)

    # Workflow control
    run_regression: bool = Field(default=True, description="Run regression step")
    run_scaling: bool = Field(default=True, description="Run scaling step")
    run_pca: bool = Field(default=True, description="Run PCA step")
    run_neighbors: bool = Field(default=True, description="Run neighbors and UMAP")
    run_integration: bool = Field(default=True, description="Run batch correction")
    # Note: save_dir is inherited from SclucidBaseConfig

    @classmethod
    def from_simple_dict(cls, simple_config: Dict[str, Any]) -> "PreprocessingWorkflowConfig":
        """
        Create PreprocessingWorkflowConfig from a simplified flat dictionary.

        This factory method allows users to create complex nested configurations
        from a simple flat dictionary, reducing the need to understand internal
        config structure.

        Args:
            simple_config: Flat dictionary with keys like:
                - normalization_method, normalization_target_sum
                - hvg_n_top_genes, hvg_method
                - run_regression, run_scaling, run_pca
                - save_dir, n_jobs

        Returns:
            PreprocessingWorkflowConfig: Fully configured workflow config

        Example:
            >>> config = PreprocessingWorkflowConfig.from_simple_dict({
            ...     "normalization_method": "scran",
            ...     "hvg_n_top_genes": 3000,
            ...     "run_regression": False,
            ...     "save_dir": "./results"
            ... })
        """
        config_data = dict(simple_config)
        kwargs = {}

        # Extract normalization parameters
        norm_params = {}
        for key in ["method", "target_sum", "exclude_highly_expressed"]:
            config_key = f"normalization_{key}"
            if config_key in config_data:
                norm_params[key] = config_data.pop(config_key)
        if norm_params:
            kwargs["normalization"] = NormalizationConfig(**norm_params)

        # Extract HVG parameters
        hvg_params = {}
        for key in ["method", "n_top_genes", "flavor", "batch_key"]:
            config_key = f"hvg_{key}"
            if config_key in config_data:
                hvg_params[key] = config_data.pop(config_key)
        if hvg_params:
            kwargs["hvg"] = HVGConfig(**hvg_params)

        # Extract scaling parameters
        scaling_params = {}
        for key in ["max_value", "vars_to_regress"]:
            config_key = f"scaling_{key}"
            if config_key in config_data:
                scaling_params[key] = config_data.pop(config_key)
        if scaling_params:
            kwargs["scaling"] = ScalingConfig(**scaling_params)

        # Extract integration parameters
        integration_params = {}
        for key in ["method", "batch_key"]:
            config_key = f"integration_{key}"
            if config_key in config_data:
                integration_params[key] = config_data.pop(config_key)
        if integration_params:
            kwargs["integration"] = IntegrationConfig(**integration_params)

        # Extract graph parameters
        graph_params = {}
        for key in ["n_pcs", "n_neighbors"]:
            config_key = f"graph_{key}"
            if config_key in config_data:
                graph_params[key] = config_data.pop(config_key)
        if graph_params:
            kwargs["graph"] = GraphConfig(**graph_params)

        # Backward compatibility: results_dir -> save_dir
        if "results_dir" in config_data:
            config_data["save_dir"] = config_data.pop("results_dir")

        # Remaining keys go directly to workflow config
        kwargs.update(config_data)

        return cls(**kwargs)

    @classmethod
    def default(cls, **kwargs) -> "PreprocessingWorkflowConfig":
        """
        Default configuration factory for the standard preprocessing path.

        This represents the canonical default pipeline:
        - Normalization (log1p, target_sum=1e4)
        - Regression (total_counts, pct_counts_mt)
        - HVG selection (2000 genes, seurat flavor)
        - Scaling (z-score, max_value=10)
        - PCA (50 components)
        - Batch correction (harmony, if batch_key present)
        - Neighbors + UMAP (15 neighbors, 50 PCs)

        Args:
            **kwargs: Override any default parameter.

        Returns:
            PreprocessingWorkflowConfig: Pre-configured for the standard path.

        Example:
            >>> config = PreprocessingWorkflowConfig.default()
            >>> adata = run_preprocessing(adata, config=config)
        """
        return cls(
            run_regression=True,
            run_scaling=True,
            run_pca=True,
            run_neighbors=True,
            run_integration=True,
            **kwargs
        )

    @classmethod
    def quick(
        cls,
        n_top_genes: int = 2000,
        run_regression: bool = False,
        run_integration: bool = False,
        **kwargs
    ) -> "PreprocessingWorkflowConfig":
        """
        Quick configuration factory for standard analyses.

        Args:
            n_top_genes: Number of highly variable genes to select
            run_regression: Whether to run regression step
            run_integration: Whether to run batch correction
            **kwargs: Additional parameters (species, n_jobs, etc.)

        Returns:
            PreprocessingWorkflowConfig: Pre-configured for standard analysis

        Example:
            >>> config = PreprocessingWorkflowConfig.quick(
            ...     n_top_genes=3000,
            ...     run_regression=True,
            ...     vars_to_regress=["percent_mito"]
            ... )
        """
        return cls(
            hvg=HVGConfig(n_top_genes=n_top_genes),
            run_regression=run_regression,
            run_integration=run_integration,
            **kwargs
        )


# Backward compatibility aliases
WorkflowConfig = PreprocessingWorkflowConfig

__all__ = [
    "NormalizationConfig",
    "HVGConfig",
    "ScalingConfig",
    "IntegrationConfig",
    "NeighborsConfig",
    "GraphConfig",
    "PreprocessingWorkflowConfig",
    "WorkflowConfig",  # Backward compatibility
    "apply_config_overrides",
]
