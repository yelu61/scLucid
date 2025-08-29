"""
Configuration classes for the analysis module of scLucid.

Defines all dataclasses for clustering, annotation, scoring, differential expression,
enrichment analysis, and advanced analysis steps, ensuring parameter traceability,
validation, and pipeline-level consistency.
"""

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

from ..utils.marker_manager import Manager

__all__ = [
    "ClusteringConfig",
    "ResolutionSearchConfig",
    "MergeClustersConfig",
    "AnnotationConfig",
    #"ScoringConfig",
    "DifferentialConfig",
    "FilterMarkersConfig",
    "CompareGroupsConfig",
    "CompareConditionsConfig",
    "EnrichmentConfig",
    "ProportionConfig",
    "AnalysisWorkflowConfig",
]


@dataclass
class BaseConfig:
    """Base class for all analysis configuration objects."""

    namespace: str = "sclucid"

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]):
        return cls(**d)


# ===================== Clustering Configs =====================


@dataclass
class ClusteringConfig(BaseConfig):
    """Configuration for a single clustering run."""

    method: Literal["leiden", "louvain"] = "leiden"
    resolution: float = 1.0
    use_rep: str = "X_pca"
    key_added: Optional[str] = None
    random_state: int = 42
    plot: bool = True
    save_dir: Optional[str] = None
    extra_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ResolutionSearchConfig(BaseConfig):
    """Configuration for optimizing clustering resolution."""

    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10)
    metric: Literal["marker_separation", "silhouette"] = "silhouette"
    marker_config: Optional[Union[str, Manager]] = None
    use_raw_for_markers: bool = False
    use_rep: str = "X_pca"
    plot: bool = True
    save_dir: Optional[str] = None


@dataclass
class MergeClustersConfig(BaseConfig):
    """Configuration for merging similar clusters."""

    cluster_key: str = "leiden_clusters"  # This is a required parameter
    similarity_threshold: float = 0.8
    method: Literal["marker_overlap", "expression_correlation"] = "marker_overlap"
    key_added: Optional[str] = None
    de_method_for_markers: str = "wilcoxon"  # For marker_overlap method


# ===================== Annotation Configs =====================


@dataclass
class AnnotationConfig(BaseConfig):
    """Configuration for the cell type annotation workflow."""

    cluster_key: str = "leiden_clusters"
    marker_species: str = "human"
    marker_tissue: Optional[str] = None
    run_celltypist: bool = False
    celltypist_model: str = "Immune_All_Low.pkl"
    run_scoring: bool = True
    final_method: Literal["max_score", "enrichment", "combined"] = "combined"
    key_added: str = "cell_type_auto"
    min_confidence: float = 0.1
    plot: bool = True
    # Parameters for DE and enrichment steps within the annotation workflow
    min_log2fc: float = 0.5
    min_in_group_pct: float = 0.2
    enrichment_gene_sets: List[str] = field(
        default_factory=lambda: ["GO_Biological_Process_2023"]
    )


# ===================== Scoring Configs =====================


#@dataclass
#class ScoringConfig(BaseConfig):
    """Configuration for gene set and cell type scoring."""

#    gene_sets: Dict[str, List[str]]
#    layer: Optional[str] = "normalized"
#    use_raw: bool = False
#    ctrl_size: int = 50
#    score_name_suffix: str = "_score"


# ===================== Differential Expression Configs =====================


@dataclass
class DifferentialConfig(BaseConfig):
    """Configuration for one-vs-rest differential expression analysis."""

    groupby: str = "leiden_clusters"
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon"
    layer: Optional[str] = None
    use_raw: bool = True
    key_added: Optional[str] = "rank_genes_groups"
    groups: Optional[List[str]] = None
    reference: str = "rest"
    pval_cutoff: Optional[float] = None
    fold_change_max: Optional[float] = None


@dataclass
class FilterMarkersConfig(BaseConfig):
    """Configuration for filtering marker gene results."""

    min_log2fc: float = 1.0
    max_padj: float = 0.05
    min_in_group_pct: float = 0.25
    keep_top_n: Optional[int] = 100


#@dataclass
#class CompareGroupsConfig(BaseConfig):
    """Configuration for comparing two specific groups."""

#    groupby: str
#    group1: str
#    group2: str
#    use_raw: bool = True
    # Add other DE parameters as needed


#@dataclass
#class CompareConditionsConfig(BaseConfig):
    """Configuration for comparing the same cell type across conditions."""

#    groupby: str
#    group_name: str
#    condition_key: str
#    condition1: str
#    condition2: str
#    use_raw: bool = True
    # Add other DE parameters as needed


# ===================== Enrichment & Proportion Configs =====================


@dataclass
class EnrichmentConfig(BaseConfig):
    """Configuration for functional enrichment analysis."""

    de_key: str = "rank_genes_groups"
    gene_sets: List[str] = field(
        default_factory=lambda: ["GO_Biological_Process_2023", "KEGG_2021_Human"]
    )
    organism: Optional[str] = "Human"
    n_top_genes: int = 100
    max_padj: float = 0.05
    plot: bool = True
    save_dir: Optional[str] = None


#@dataclass
#class ProportionConfig(BaseConfig):
    """Configuration for cell type proportion analysis."""

#    celltype_col: str
#    condition_col: str
#    group_col: str = "sampleID"
#    plot_types: List[str] = field(default_factory=lambda: ["box", "bar"])
#    test: Literal["t", "wilcoxon", "anova"] = "wilcoxon"
#    save_dir: Optional[str] = None


# ===================== Master Workflow Config =====================


@dataclass
class AnalysisWorkflowConfig(BaseConfig):
    """Master configuration for the entire analysis workflow."""

    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    de: DifferentialConfig = field(default_factory=DifferentialConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    #scoring: Optional[ScoringConfig] = None
    #proportion: Optional[ProportionConfig] = None


# Update __all__ to include all configuration classes
__all__ = [
    name
    for name, obj in locals().items()
    if isinstance(obj, type) and issubclass(obj, BaseConfig)
]
