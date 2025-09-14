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
    #"ProportionConfig",
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

    method: Literal["leiden", "louvain", "kmeans", "hdbscan"] = "leiden"
    resolution: float = 1.0
    n_clusters: Optional[int] = None  # Specifically for methods like kmeans
    use_rep: str = "X_pca"
    key_added: Optional[str] = None
    random_state: int = 42
    plot: bool = True
    save_dir: Optional[str] = None
    extra_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ResolutionSearchConfig(BaseConfig):
    """
    Configuration for optimizing clustering resolution.
    This function now acts as a guide, providing metrics but not setting a final resolution.
    """
    method: Literal["leiden", "louvain"] = "leiden"
    use_rep: str = "X_pca"
    resolution_range: Tuple[float, float, int] = (0.2, 2.0, 10) # Start from a slightly higher resolution
    
    # Metrics to compute for guidance
    compute_silhouette: bool = True
    compute_marker_abundance: bool = True
    compute_stability: bool = True
    
    # Parameters for marker abundance metric
    de_method_for_markers: str = "wilcoxon"
    min_log2fc_for_markers: float = 0.5
    min_pct_for_markers: float = 0.25
    
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


@dataclass
class ScoringConfig(BaseConfig):
    """Configuration for gene set and cell type scoring."""
    gene_sets: Dict[str, List[str]]
    layer: Optional[str] = "normalized"
    use_raw: bool = True
    ctrl_size: int = 50
    score_name_suffix: str = "_score"
    plot_comparison: bool = True
    comparison_groupby: Optional[str] = None


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
    pval_cutoff: Optional[float] = 0.05 #: Post-hoc filter for adjusted p-value. Only genes with pvals_adj <= cutoff will be kept.
    fold_change_max: Optional[float] = None #: Post-hoc clipping for log fold change. Useful for visualization.

@dataclass
class FilterMarkersConfig(BaseConfig):
    """
    Configuration for filtering marker gene results from differential expression analysis.
    
    This config allows for fine-grained control over marker selection based on statistical
    significance, effect size (log fold change), and expression specificity.
    """
    key: str = "rank_genes_groups" #: The key for the source DE results in adata.uns.
    key_added: Optional[str] = None #: The key under which the filtered DataFrame will be stored. If None, defaults to f"{key}_filtered_df".

    # --- Core Filtering Criteria ---
    min_log2fc: float = 1.0 #: Minimum log2 fold change required.
    max_padj: float = 0.05 #: Maximum adjusted p-value (e.g., FDR) allowed.
    min_in_group_pct: float = 0.25 #: Minimum percentage (as a fraction, e.g., 0.25 for 25%) of cells expressing the gene within the target group.

    # --- Specificity Filtering Criteria ---
    max_out_group_pct: Optional[float] = 0.5 #: Maximum percentage of cells expressing the gene outside the target group. Helps find specific markers.
    min_diff_pct: Optional[float] = 0.1 #: Minimum difference between in-group and out-group expression percentages.
    
    # --- Top N Selection ---
    keep_top_n: Optional[int] = 100 #: After all filters, keep only the top N markers for each group, ranked by log2fc.


@dataclass
class CompareGroupsConfig(BaseConfig):
    """Configuration for comparing two specific groups (e.g., two cell types)."""
    groupby: str = "" #: The column in adata.obs that contains the groups to compare.
    group1: str = "" #: The name of the first group.
    group2: str = "" #: The name of the second group (will be used as the reference).
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon"
    layer: Optional[str] = None
    use_raw: bool = False
    key_added: Optional[str] = None #: Key to store the results DataFrame in adata.uns. Defaults to 'compare_{group1}_vs_{group2}'.
    # Filtering parameters for the comparison results
    n_top_genes: int = 50 #: Number of top genes to keep for each group after filtering.
    min_log2fc: float = 0.5
    max_padj: float = 0.05
    min_in_group_pct: float = 0.1
    plot: bool = True #: Whether to generate a volcano plot of the results.
    save_dir: Optional[str] = None #: Directory to save the volcano plot.


@dataclass
class CompareConditionsConfig(BaseConfig):
    """Configuration for comparing the same cell type across two conditions."""
    groupby: str = "" #: The column identifying the cell type/group to analyze.
    group_name: str = "" #: The specific cell type/group to subset for comparison.
    condition_key: str = "" #: The column identifying the conditions (e.g., 'treatment', 'disease_status').
    condition1: str = "" #: The first condition.
    condition2: str = "" #: The second condition (will be used as the reference).
    # Inherits comparison and plotting parameters from CompareGroupsConfig
    comparison_params: CompareGroupsConfig = field(default_factory=CompareGroupsConfig)
    key_added: Optional[str] = None #: Defaults to 'compare_{c1}_vs_{c2}_in_{group}'.


# ===================== Enrichment & Proportion Configs =====================


@dataclass
class EnrichmentConfig(BaseConfig):
    """
    Configuration for functional enrichment analysis.
    Supports both online (Enrichr) and offline (local GMT files) modes.
    """
    de_key: str = "rank_genes_groups_filtered_df" #: Key for the DE results DataFrame in adata.uns.
    mode: Literal["online", "offline"] = "offline" #: Analysis mode. 'offline' is recommended for stability.
    organism: Literal["human", "mouse"] = "human" #: Species for the analysis.
    gmt_version: str = "v2025"
    
    # For 'online' mode, this is a list of Enrichr library names.
    # For 'offline' mode, this is a list of categories (e.g., 'hallmark', 'go_bp')
    # corresponding to local .gmt files in the resources directory.
    gene_sets_online: List[str] = field(default_factory=lambda: ["GO_Biological_Process_2023"])
    gene_sets_offline: List[str] = field(default_factory=lambda: ["hallmark", "go_bp", "reactome"])
    
    custom_gene_sets: Optional[str] = None #: Path to a user-provided .gmt file for offline mode.
    n_top_genes: int = 100 #: Number of top marker genes from each cluster to use for enrichment.
    min_genes_for_enrichment: int = 10 #: Minimum number of genes required to run enrichment for a cluster.
    max_padj: float = 0.05 #: Adjusted p-value cutoff for filtering enrichment results.
    key_added: str = "enrichment" #: Key under which to store the results in adata.uns.
    plot: bool = True #: Whether to generate and save summary plots for each cluster.
    n_plot_terms: int = 15 #: Number of top terms to show in the plot.
    save_dir: Optional[str] = None #: Directory to save plots and results.


@dataclass
class ProportionConfig(BaseConfig):
    """Configuration for cell type proportion analysis."""
    celltype_col: str
    condition_col: str
    group_col: str = "sample"  # e.g., sampleID
    plot_types: List[str] = field(default_factory=lambda: ["box", "bar"])
    test: Literal["t", "wilcoxon", "anova"] = "wilcoxon"
    out_dir: Optional[str] = None

# ===================== Master Workflow Config =====================


@dataclass
class AnalysisWorkflowConfig(BaseConfig):
    """Master configuration for the entire analysis workflow."""
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    de: DifferentialConfig = field(default_factory=DifferentialConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    scoring: Optional[ScoringConfig] = None
    proportion: Optional[ProportionConfig] = None


# Update __all__ to include all configuration classes
__all__ = [
    name
    for name, obj in locals().items()
    if isinstance(obj, type) and issubclass(obj, BaseConfig)
]
