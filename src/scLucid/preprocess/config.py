"""
Configuration for the scLucid preprocessing module.

This module defines dataclasses used to configure the preprocessing workflows,
ensuring a reproducible and customizable analysis pipeline.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Dict

log = logging.getLogger(__name__)

__all__ = [
    "NormalizationConfig",
    "HVGConfig",
    "ScalingConfig",
    "IntegrationConfig",
    "GraphConfig",
    "PreprocessingConfig",
]

@dataclass
class NormalizationConfig:
    """
    Configuration for the data normalization step.
    """
    method: Literal["standard", "scran", "pearson_residuals", "clr"] = "standard"
    target_sum: float = 1e4
    exclude_highly_expressed: bool = False
    max_fraction: float = 0.05
    plot_global_distribution: bool = True
    save_dir: Optional[str] = None

@dataclass
class HVGConfig:
    """
    Configuration for Highly Variable Gene (HVG) selection.
    """
    method: Literal["scanpy", "custom", "triku"] = "scanpy"
    layer: Optional[str] = "normalized"
    n_top_genes: int = 2000
    plot: Optional[bool] = None
    report: bool = False
    flavor: Literal["seurat", "seurat_v3", "cell_ranger", "pearson_residuals"] = "seurat_v3"
    batch_key: Optional[str] = None
    sample_key: Optional[str] = "sampleID"
    min_n_samples: int = 2
    n_highly_expressed_genes: int = 50
    n_specific_genes: int = 20
    exclude_gene_types: Optional[List[str]] = None
    plot: Optional[bool] = None
    save_dir: Optional[str] = None
    report: bool = False

@dataclass
class ScalingConfig:
    """
    Configuration for data scaling and regression.
    """
    vars_to_regress: Optional[List[str]] = field(
        default_factory=lambda: ["total_counts", "pct_counts_mt"]
    )
    max_value: Optional[float] = 10.0

@dataclass
class IntegrationConfig:
    """
    Configuration for batch correction and data integration.
    """
    method: Optional[Literal["harmony", "scanorama", "scvi", "bbknn", "combat"]] = "harmony"
    batch_key: Optional[str] = "sampleID"
    harmony_params: Dict = field(default_factory=lambda: {"max_iter_harmony": 20, "theta": 2.0})
    scvi_params: Dict = field(default_factory=lambda: {"n_latent": 30, "max_epochs": 500})

@dataclass
class GraphConfig:
    """
    Configuration for neighborhood graph and UMAP embedding.
    """
    n_pcs: int = 50
    n_neighbors: int = 15

@dataclass
class PreprocessingConfig:
    """
    Master configuration for the entire preprocessing workflow.
    Encapsulates all step-specific configurations.
    """
    counts_layer: str = "counts"
    normalized_layer: str = "log1p_norm"
    scaled_layer: str = "scaled"
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    hvg: HVGConfig = field(default_factory=HVGConfig)
    scaling: ScalingConfig = field(default_factory=ScalingConfig)
    integration: IntegrationConfig = field(default_factory=IntegrationConfig)
    graph: GraphConfig = field(default_factory=GraphConfig)

    def __post_init__(self):
        log.debug(f"Validating {self.__class__.__name__}...")
        if self.hvg.n_top_genes <= 0:
            raise ValueError("hvg.n_top_genes must be a positive integer.")
        if self.graph.n_pcs <= 0:
            raise ValueError("graph.n_pcs must be a positive integer.")
        if self.graph.n_neighbors <= 0:
            raise ValueError("graph.n_neighbors must be a positive integer.")
        if self.integration.batch_key and not self.integration.method:
            log.warning(
                f"A `batch_key` ('{self.integration.batch_key}') is provided, but `integration.method` is None."
            )
        if not self.integration.batch_key and self.integration.method:
            log.warning(
                f"An `integration.method` ('{self.integration.method}') is provided, but `integration.batch_key` is None."
            )
            self.integration.method = None