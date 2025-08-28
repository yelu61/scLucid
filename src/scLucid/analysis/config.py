"""
Configuration classes for the analysis module of scLucid.

Defines all dataclasses for clustering, annotation, scoring, differential expression,
enrichment analysis, and advanced analysis steps, ensuring parameter traceability,
validation, and pipeline-level consistency.
"""

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
from ..utils.marker_manager import Manager
__all__ = [
    "ClusteringConfig",
    "ResolutionSearchConfig",
    "AnnotationConfig",
    "ScoringConfig",
    "DifferentialConfig",
    "EnrichmentConfig",
]

# ===================== Clustering Configs =====================

@dataclass
class ClusteringConfig:
    """
    Configuration for a clustering run.
    
    Attributes:
        method: Clustering algorithm ('leiden', 'louvain', 'kmeans', 'hdbscan').
        resolution: Resolution parameter (for leiden/louvain).
        n_clusters: Number of clusters (for kmeans).
        use_rep: Embedding to use (e.g., 'X_pca').
        key_added: Key for cluster assignments in adata.obs.
        random_state: Random seed.
        plot: Whether to plot results.
        save_dir: Directory to save plots.
        extra_params: Additional algorithm-specific parameters.
    """
    method: Literal["leiden", "louvain", "kmeans", "hdbscan"] = "leiden"
    resolution: float = 1.0
    n_clusters: Optional[int] = None
    use_rep: str = "X_pca"
    key_added: Optional[str] = None
    random_state: int = 42
    plot: bool = True
    save_dir: Optional[str] = None
    extra_params: Dict[str, Any] = field(default_factory=dict)

    def validate(self):
        if self.method not in ["leiden", "louvain", "kmeans", "hdbscan"]:
            raise ValueError("Invalid clustering method")
        if self.method in ["leiden", "louvain"] and self.resolution <= 0:
            raise ValueError("Resolution must be positive")
        if self.method == "kmeans" and (self.n_clusters is None or self.n_clusters <= 0):
            raise ValueError("n_clusters must be positive for kmeans")
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "ClusteringConfig":
        return ClusteringConfig(**d)

@dataclass
class ResolutionSearchConfig:
    """
    Configuration for clustering resolution search/optimization.
    
    Attributes:
        resolution_range: (start, end, n_steps) for grid search.
        metric: Metric for optimization ('marker_separation' or 'silhouette').
        marker_config: Marker config for marker_separation (str or Manager).
        use_raw_for_markers: Whether to evaluate markers on adata.raw.
        use_rep: Embedding to use for silhouette.
        plot: Whether to plot optimization curves.
        save_dir: Where to save plots.
    """
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10)
    metric: Literal["marker_separation", "silhouette"] = "silhouette"
    marker_config: Optional[Union[str, Manager]] = None
    use_raw_for_markers: bool = False
    use_rep: str = "X_pca"
    plot: bool = True
    save_dir: Optional[str] = None

    def validate(self):
        if self.metric not in ["marker_separation", "silhouette"]:
            raise ValueError("Unknown optimization metric")
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "ResolutionSearchConfig":
        return ResolutionSearchConfig(**d)

# ===================== Annotation Configs =====================

@dataclass
class AnnotationConfig:
    """
    Configuration for cell type annotation workflow.
    """
    cluster_key: str
    marker_species: str = "human"
    marker_tissue: Optional[str] = None
    run_celltypist: bool = False
    celltypist_model: str = "Immune_All_Low.pkl"
    run_scoring: bool = True
    final_method: Literal["max_score", "enrichment", "combined"] = "combined"
    key_added: str = "cell_type"
    min_confidence: float = 0.1

    def validate(self):
        assert self.final_method in ["max_score", "enrichment", "combined"]

    def to_dict(self): return self.__dict__
    @staticmethod
    def from_dict(d): return AnnotationConfig(**d)

# ===================== Scoring Configs =====================

@dataclass
class ScoringConfig:
    """Gene set/cell type scoring configuration."""
    marker_config: Union[str, Any] = "base_human"
    layer: Optional[str] = "log1p_norm"
    use_raw: bool = False
    min_genes: int = 3
    ctrl_size: int = 50
    score_name_suffix: str = "_score"
    extra_params: Dict[str, Any] = field(default_factory=dict)
    def validate(self): pass
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "ScoringConfig":
        return ScoringConfig(**d)

# ===================== Differential Expression Configs =====================

@dataclass
class DifferentialConfig:
    """Differential expression analysis configuration."""
    method: Literal["wilcoxon", "t-test", "logreg", "cosg"] = "wilcoxon"
    groupby: str = "leiden"
    layer: Optional[str] = None
    use_raw: bool = False
    min_cells: int = 5
    groups: Optional[List[str]] = None
    reference: Optional[str] = "rest"
    fold_change_max: Optional[float] = None
    pval_cutoff: Optional[float] = 0.05
    key_added: Optional[str] = None
    save_dir: Optional[str] = None
    extra_params: Dict[str, Any] = field(default_factory=dict)
    def validate(self):
        if self.method not in ["wilcoxon", "t-test", "logreg", "cosg"]:
            raise ValueError("Unknown DE method")
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "DifferentialConfig":
        return DifferentialConfig(**d)

@dataclass
class FilterMarkersConfig:
    """Marker gene筛选参数设置。"""
    min_log2fc: float = 1.0
    max_padj: float = 0.05
    min_in_group_pct: float = 0.25
    max_out_group_pct: Optional[float] = None
    min_diff_pct: Optional[float] = None
    keep_top_n: Optional[int] = None
    key_added: Optional[str] = None
    def validate(self): pass
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "FilterMarkersConfig":
        return FilterMarkersConfig(**d)

@dataclass
class CompareGroupsConfig:
    """配置：两组直接差异分析参数。"""
    groupby: str = "leiden"
    group1: str = ""
    group2: str = ""
    layer: Optional[str] = None
    use_raw: bool = False
    n_genes: int = 50
    min_log2fc: float = 0.5
    max_padj: float = 0.05
    min_in_group_pct: float = 0.1
    plot: bool = True
    save_path: Optional[str] = None
    key_added: Optional[str] = None
    def validate(self):
        if not self.group1 or not self.group2:
            raise ValueError("group1 and group2 must be specified")
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "CompareGroupsConfig":
        return CompareGroupsConfig(**d)

@dataclass
class CompareConditionsConfig:
    """配置：同一细胞类型不同条件的DE分析。"""
    groupby: str = "cell_type"
    group_name: str = ""
    condition_key: str = "condition"
    condition1: str = ""
    condition2: str = ""
    key_added: Optional[str] = None
    extra_params: Dict[str, Any] = field(default_factory=dict)
    def validate(self):
        if not self.group_name or not self.condition1 or not self.condition2:
            raise ValueError("group_name, condition1, condition2 must be specified")
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "CompareConditionsConfig":
        return CompareConditionsConfig(**d)

# ===================== Enrichment Configs =====================

@dataclass
class EnrichmentConfig:
    """功能富集分析配置。"""
    de_key: str = "rank_genes_groups"
    gene_sets: List[str] = field(default_factory=lambda: ["GO_Biological_Process_2023"])
    organism: str = "Human"
    n_top_genes: int = 100
    min_genes: int = 10
    max_genes: int = 500
    min_enrichment_score: float = 0.0
    max_padj: float = 0.05
    background_genes: Optional[List[str]] = None
    plot: bool = False
    save_path: Optional[str] = None
    key_added: str = "enrichment"
    def validate(self): pass
    def to_dict(self): return asdict(self)
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "EnrichmentConfig":
        return EnrichmentConfig(**d)
