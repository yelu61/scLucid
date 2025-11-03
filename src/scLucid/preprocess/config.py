"""
Configuration for the scLucid preprocessing module.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Union

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
    """Configuration for data normalization with validation."""
    
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
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        # Validate target_sum
        if self.target_sum <= 0:
            raise ValueError(f"target_sum must be positive, got {self.target_sum}")
        
        if self.target_sum < 1000 or self.target_sum > 1e7:
            log.warning(
                f"target_sum={self.target_sum:.0e} is unusual. "
                "Typical values are 1e4 (10,000) for UMI counts."
            )
        
        # Validate max_fraction
        if not 0 < self.max_fraction < 1:
            raise ValueError(f"max_fraction must be in (0, 1), got {self.max_fraction}")
        
        # Check method-specific constraints
        if self.method == "pearson_residuals":
            if self.target_sum != 1e4:
                log.warning(
                    "Pearson residuals normalization ignores target_sum parameter"
                )
        
        # Validate layer names
        reserved_names = ['X', 'raw']
        if self.output_layer in reserved_names:
            raise ValueError(
                f"output_layer cannot be '{self.output_layer}'. "
                f"Reserved names: {reserved_names}"
            )

@dataclass
class HVGConfig:
    """Configuration for HVG selection with validation."""
    
    method: Literal["scanpy", "custom", "triku"] = "scanpy"
    n_top_genes: int = 2000
    flavor: Literal["seurat", "seurat_v3", "cell_ranger"] = "seurat_v3"
    span: Optional[float] = 0.3
    batch_key: Optional[str] = None
    sample_key: str = "sampleID"
    min_n_samples: int = 2
    n_highly_expressed_genes: int = 50
    n_specific_genes: int = 20
    exclude_gene_types: Optional[List[str]] = field(
        default_factory=lambda: ["mitochondrial", "ribosomal"]
    )
    plot: bool = True
    report: bool = True
    save_dir: Optional[str] = None
    
    def __post_init__(self):
        """Validate HVG configuration."""
        # Validate n_top_genes
        if self.n_top_genes < 100:
            log.warning(
                f"n_top_genes={self.n_top_genes} is very low. "
                "Typical values are 2000-5000."
            )
        
        if self.n_top_genes > 10000:
            log.warning(
                f"n_top_genes={self.n_top_genes} is very high. "
                "This may include too much noise."
            )
        
        # Validate span for seurat flavor
        if self.flavor == "seurat" and self.span is not None:
            if not 0.01 < self.span < 1:
                raise ValueError(f"span must be in (0.01, 1), got {self.span}")
        
        # Validate method-specific parameters
        if self.method == "custom":
            if self.min_n_samples < 1:
                raise ValueError(
                    f"min_n_samples must be >= 1, got {self.min_n_samples}"
                )
            
            if self.min_n_samples > 10:
                log.warning(
                    f"min_n_samples={self.min_n_samples} is quite strict. "
                    "Genes must appear as HVG in many samples."
                )


@dataclass
class ScalingConfig:
    """Configuration for data scaling and regression."""
    vars_to_regress: Optional[List[str]] = field(
        default_factory=lambda: ["total_counts", "pct_counts_mt"]
    )
    
    # Clearer regression settings
    regress_in_scale: bool = True
    vars_to_regress_in_scale: Optional[List[str]] = None  # If None, use vars_to_regress
    
    # Explicit input layer for regression
    input_layer_for_regress: str = "normalized"  # Layer to regress from
    
    scale_method: Literal["zscore", "robust", "minmax"] = "zscore"
    max_value: Optional[float] = 10.0
    

@dataclass
class IntegrationConfig:
    """Configuration for batch correction and data integration."""
    method: Optional[Literal["harmony", "scanorama", "scvi", "bbknn", "combat"]] = "harmony"
    batch_key: Optional[Union[str, List[str]]] = "sampleID"
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