"""
Configuration for the scLucid preprocessing module.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional

log = logging.getLogger(__name__)

__all__ = [
    "NormalizationConfig",
    "HVGConfig",
    "ScalingConfig",
    "IntegrationConfig",
    "NeighborsConfig",
    "GraphConfig",
    "WorkflowConfig",
]

@dataclass
class NormalizationConfig:
    """Configuration for the data normalization step."""
    method: Literal["standard", "scran", "pearson_residuals", "clr"] = "standard"
    target_sum: float = 1e4
    exclude_highly_expressed: bool = False
    max_fraction: float = 0.05
    input_layer: str = "counts"
    output_layer: str = "normalized"
    update_X: bool = True
    plot: bool = True
    report: bool = True
    save_dir: Optional[str] = None

@dataclass
class HVGConfig:
    """Configuration for Highly Variable Gene (HVG) selection."""
    method: Literal["scanpy", "custom", "triku"] = "scanpy"
    n_top_genes: int = 2000
    flavor: Literal["seurat", "seurat_v3", "cell_ranger"] = "seurat_v3"
    span: Optional[float] = 0.3,
    batch_key: Optional[str] = None
    sample_key: str = "sampleID"
    min_n_samples: int = 2
    n_highly_expressed_genes: int = 50
    n_specific_genes: int = 20
    exclude_gene_types: Optional[List[str]] = field(default_factory=lambda: ["mitochondrial", "ribosomal"])
    plot: bool = True
    report: bool = True
    save_dir: Optional[str] = None

@dataclass
class ScalingConfig:
    """Configuration for data scaling and regression."""
    vars_to_regress: Optional[List[str]] = field(default_factory=lambda: ["total_counts", "pct_counts_mt"])
    scale_method: Literal["zscore", "robust", "minmax"] = "zscore"
    max_value: Optional[float] = 10.0

@dataclass
class IntegrationConfig:
    """Configuration for batch correction and data integration."""
    method: Optional[Literal["harmony", "scanorama", "scvi", "bbknn", "combat"]] = "harmony"
    batch_key: Optional[str] = "sampleID"
    use_rep: str = "X_pca"
    output_key: Optional[str] = None
    plot: bool = True
    save_dir: Optional[str] = None
    harmony_params: Dict = field(default_factory=lambda: {"max_iter_harmony": 20, "theta": 2.0})
    scvi_params: Dict = field(default_factory=lambda: {"n_latent": 30, "max_epochs": 500})
    hvg_key: Optional[str] = None # For Scanorama

@dataclass
class NeighborsConfig:
    """Configuration for optimizing nearest neighbor and PCA parameters."""
    n_neighbors_list: List[int] = field(default_factory=lambda: [15, 30, 50])
    n_pcs_list: List[int] = field(default_factory=lambda: [30, 40, 50])
    use_rep: str = "X_pca"
    clustering_method: Literal["leiden", "louvain"] = "leiden"
    resolution: float = 1.0
    subsample: Optional[int] = 5000
    n_jobs: int = -1
    plot: bool = True
    save_dir: Optional[str] = None

@dataclass
class GraphConfig:
    """Configuration for the final neighborhood graph and UMAP embedding."""
    n_pcs: int = 50
    n_neighbors: int = 15

@dataclass
class WorkflowConfig:
    """Master configuration for the entire preprocessing workflow."""
    counts_layer: str = "counts"
    normalized_layer: str = "normalized"
    regressed_layer: str = "regressed"
    scaled_layer: str = "scaled"
    
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    hvg: HVGConfig = field(default_factory=HVGConfig)
    scaling: ScalingConfig = field(default_factory=ScalingConfig)
    integration: IntegrationConfig = field(default_factory=IntegrationConfig)
    graph: GraphConfig = field(default_factory=GraphConfig)