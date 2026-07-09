"""
Pydantic-based configuration for the scLucid preprocessing module.

Migrates from dataclasses to Pydantic for consistent validation and serialization.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import ConfigDict, Field, field_validator, model_validator

from ..base_config import (
    SclucidBaseConfig,
    WorkflowConfigBase,
)
from ..base_config import (
    apply_config_overrides as _apply_config_overrides,
)

logger = logging.getLogger(__name__)

# Re-export for backward compatibility within preprocess submodules
apply_config_overrides = _apply_config_overrides


class NormalizationConfig(SclucidBaseConfig):
    """Configuration for data normalization."""

    model_config = ConfigDict(extra="ignore")

    method: Literal[
        "standard",
        "scran",
        "pearson_residuals",
        "clr",
        "quality_aware",  # experimental heuristic
        "deconvolution_pool",  # experimental heuristic, not scran
        "quantile_transform",  # non-parametric transform, not regression
        # Backward compatibility alias
        "quantile_regression",
    ] = Field(default="standard")
    target_sum: float = Field(default=1e4, gt=0, description="Target sum for normalization")
    exclude_highly_expressed: bool = Field(default=False)
    max_fraction: float = Field(default=0.05, gt=0, lt=1)
    clr_pseudocount: float = Field(
        default=1.0,
        gt=0,
        description="Pseudocount for CLR: log(x + pseudocount) minus per-cell mean log.",
    )
    input_layer: str = Field(default="counts", description="Input layer name")
    output_layer: str = Field(default="normalized", description="Output layer name")
    update_X: bool = Field(default=False, description="Update adata.X with normalized data")  # noqa: N815
    set_raw: bool = Field(
        default=False,
        description="Store a full-gene copy of the normalized AnnData in adata.raw before downstream filtering.",
    )

    @field_validator("output_layer")
    @classmethod
    def validate_output_layer(cls, v: str) -> str:
        """Prevent reserved names for output layer."""
        reserved = {"X", "raw"}
        if v in reserved:
            raise ValueError(f"output_layer cannot be one of reserved names: {reserved}")
        return v

    @model_validator(mode="after")
    def validate_method_params(self) -> NormalizationConfig:
        """Validate method-specific constraints."""
        if self.method == "pearson_residuals" and self.target_sum != 1e4:
            logger.warning("Pearson residuals normalization ignores target_sum parameter")
        return self


class HVGConfig(SclucidBaseConfig):
    """Configuration for HVG selection."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["scanpy", "custom", "triku", "deviance"] = Field(default="scanpy")
    n_top_genes: int = Field(default=2000, ge=100, le=20000, description="Number of HVGs to select")
    auto_n_top_genes: bool = Field(
        default=False,
        description="Automatically adapt n_top_genes based on dataset size.",
    )
    auto_n_top_genes_method: Literal["linear", "log"] = Field(
        default="log",
        description="Scaling method for auto n_top_genes: 'linear' (conservative) or 'log' (aggressive for large datasets).",
    )
    flavor: Literal["auto", "seurat", "seurat_v3", "cell_ranger"] = Field(
        default="auto",
        description=(
            "HVG flavor. 'auto' uses dependency-light log-normalized HVG selection by "
            "default and reserves seurat_v3 for raw-count inputs when its optional "
            "dependency is available."
        ),
    )
    span: Optional[float] = Field(default=0.3)
    batch_key: Optional[str] = Field(default=None)
    sample_key: str = Field(default="sampleID")
    min_n_samples: int = Field(default=2, ge=1)
    n_highly_expressed_genes: int = Field(default=50, ge=0)
    n_specific_genes: int = Field(default=0, ge=0)
    exclude_gene_types: Optional[List[str]] = Field(
        default_factory=lambda: ["mitochondrial", "ribosomal"]
    )
    protected_gene_presets: List[
        Literal[
            "immune_receptor",
            "cytokine",
            "transcription_factor",
            "pathway",
            "tumor_heterogeneity",
        ]
    ] = Field(
        default_factory=list,
        description="Biology-oriented gene presets to preserve in the final HVG mask.",
    )
    protected_gene_sets: Dict[str, List[str]] = Field(
        default_factory=dict,
        description="Custom protected genes by set name. Matching is case-insensitive.",
    )
    protect_genes: bool = Field(
        default=True,
        description="If true, genes from protected presets/custom sets present in var_names are included in final HVGs.",
    )
    protection_max_extra_genes: Optional[int] = Field(
        default=None,
        ge=1,
        description="Optional cap on protected genes added beyond algorithm-selected HVGs.",
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
        default=None,
        description=(
            "Covariates to regress out. Keep None by default to avoid removing "
            "biological gradients unless the user explicitly opts in."
        ),
    )
    regress_in_scale: bool = Field(default=False)
    vars_to_regress_in_scale: Optional[List[str]] = Field(default=None)
    input_layer_for_regress: str = Field(default="normalized")
    scale_method: Literal["zscore", "robust", "minmax"] = Field(default="zscore")
    max_value: Optional[float] = Field(default=10.0, gt=0)


class IntegrationConfig(SclucidBaseConfig):
    """Configuration for batch correction and data integration."""

    model_config = ConfigDict(extra="ignore")

    method: Optional[
        Literal["harmony", "scanorama", "scvi", "scanvi", "scANVI", "bbknn", "combat"]
    ] = Field(default="harmony")
    batch_key: Optional[Union[str, List[str]]] = Field(default="sampleID")
    use_rep: str = Field(default="X_pca")
    output_key: Optional[str] = Field(default=None)
    harmony_params: Dict[str, Any] = Field(
        default_factory=lambda: {"max_iter_harmony": 50, "theta": 2.0}
    )
    scvi_params: Dict[str, Any] = Field(default_factory=lambda: {"n_latent": 30, "max_epochs": 500})
    scanvi_params: Dict[str, Any] = Field(
        default_factory=lambda: {
            "n_latent": 30,
            "max_epochs": 500,
            "early_stopping_patience": 20,
        }
    )
    scanvi_labels_key: Optional[str] = Field(
        default=None, description="Key in adata.obs for cell type labels (required for scANVI)"
    )
    hvg_key: Optional[str] = Field(default=None, description="For Scanorama")
    method_kwargs: Dict[str, Any] = Field(
        default_factory=dict,
        description="Additional kwargs passed directly to the underlying integration method. "
        "These override any values in harmony_params / scvi_params / scanvi_params.",
    )
    auto_decide: bool = Field(
        default=False,
        description="If True, run integration only when decide_integration(auto) approves it.",
    )
    evaluate: bool = Field(
        default=False,
        description="If True, evaluate integration quality after successful correction.",
    )
    condition_key: Optional[str] = Field(default=None)
    biology_columns: List[str] = Field(default_factory=list)
    label_key: Optional[str] = Field(default=None)
    tumor: bool = Field(default=False)


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


class GeneBiotypeConfig(SclucidBaseConfig):
    """Configuration for optional gene biotype annotation and filtering."""

    model_config = ConfigDict(extra="ignore")

    annotate: bool = Field(
        default=False,
        description="Annotate genes with biotype metadata before normalization.",
    )
    filter: bool = Field(
        default=False,
        description="Filter genes by biotype after annotation and low-detection filtering.",
    )
    species: Literal["human", "mouse", "rat"] = Field(default="human")
    method: Literal["reference", "ensembl", "custom"] = Field(default="reference")
    custom_biotype_path: Optional[str] = Field(
        default=None,
        description="Local CSV/TSV reference with gene_name and biotype columns.",
    )
    keep_biotypes: Optional[List[str]] = Field(
        default=None,
        description="Biotype categories to keep. Overrides use_recommended when provided.",
    )
    use_recommended: bool = Field(
        default=True,
        description="Keep recommended analysis biotypes when keep_biotypes is not provided.",
    )
    filter_stage: Literal["before_normalization", "after_raw"] = Field(
        default="after_raw",
        description=(
            "When to apply biotype filtering. 'after_raw' preserves full-gene normalized "
            "expression in .raw before subsetting analysis features."
        ),
    )
    fuzzy_match: bool = Field(default=True)
    overwrite: bool = Field(default=True)
    allow_download: bool = Field(
        default=False,
        description="Allow downloading a reference if no bundled/cache/custom table is available.",
    )
    prefer_bundled: bool = Field(default=True)
    cache_dir: Optional[str] = Field(default=None)
    fail_on_error: bool = Field(
        default=False,
        description="Raise biotype annotation/filtering errors instead of recording a skipped status.",
    )

    @model_validator(mode="after")
    def validate_gene_biotype_config(self) -> GeneBiotypeConfig:
        """Validate biotype workflow settings."""
        if self.filter and not self.annotate:
            raise ValueError("gene_biotype.filter=True requires gene_biotype.annotate=True")
        if self.method == "custom" and not self.custom_biotype_path:
            raise ValueError("gene_biotype.custom_biotype_path is required when method='custom'")
        return self


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
    gene_biotype: GeneBiotypeConfig = Field(default_factory=GeneBiotypeConfig)

    # PCA auto-selection
    auto_select_n_pcs: bool = Field(
        default=False, description="Automatically select n_pcs using elbow or cumulative method"
    )
    n_pcs_selection_method: Literal["elbow", "cumulative"] = Field(default="elbow")

    # Workflow control
    run_gene_filtering: bool = Field(
        default=True,
        description="Filter genes expressed in too few cells before normalization/HVG selection.",
    )
    min_cells_per_gene: int = Field(
        default=3,
        ge=1,
        description="Minimum number of cells in which a gene must be detected.",
    )
    run_regression: bool = Field(
        default=False,
        description="Run regression step. Disabled by default to preserve biological signal.",
    )
    run_scaling: bool = Field(default=True, description="Run scaling step")
    run_pca: bool = Field(default=True, description="Run PCA step")
    run_neighbors: bool = Field(default=True, description="Run neighbors and UMAP")
    run_integration: bool = Field(
        default=False,
        description="Run batch correction. Disabled by default; enable when batch effects are evident.",
    )
    # Note: save_dir is inherited from SclucidBaseConfig

    @classmethod
    def from_simple_dict(cls, simple_config: Dict[str, Any]) -> PreprocessingWorkflowConfig:
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

        # Extract gene biotype parameters
        gene_biotype_params = {}
        for key in [
            "annotate",
            "filter",
            "species",
            "method",
            "custom_biotype_path",
            "keep_biotypes",
            "use_recommended",
            "filter_stage",
            "allow_download",
            "fail_on_error",
        ]:
            config_key = f"gene_biotype_{key}"
            if config_key in config_data:
                gene_biotype_params[key] = config_data.pop(config_key)
        if gene_biotype_params:
            kwargs["gene_biotype"] = GeneBiotypeConfig(**gene_biotype_params)

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
    def default(cls, **kwargs) -> PreprocessingWorkflowConfig:
        """
        Default configuration factory for the standard preprocessing path.

        This represents the canonical default pipeline:
        - Low-detection gene filtering (min_cells_per_gene=3)
        - Normalization (log1p, target_sum=1e4)
        - Optional regression only when explicitly enabled
        - HVG selection (adaptive 2000-5000 genes, dependency-light auto flavor)
        - Scaling (z-score, max_value=10)
        - PCA (50 components)
        - Optional batch correction (Harmony or other configured method)
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
            run_regression=False,
            run_scaling=True,
            run_pca=True,
            run_neighbors=True,
            run_integration=False,
            **kwargs,
        )

    @classmethod
    def quick(
        cls,
        n_top_genes: int = 2000,
        run_regression: bool = False,
        run_integration: bool = False,
        **kwargs,
    ) -> PreprocessingWorkflowConfig:
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
            **kwargs,
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
    "GeneBiotypeConfig",
    "PreprocessingWorkflowConfig",
    "WorkflowConfig",  # Backward compatibility
    "apply_config_overrides",
]
