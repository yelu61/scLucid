"""
Configuration classes for the analysis module of scLucid.

Defines all dataclasses for clustering, annotation, scoring, differential expression,
enrichment analysis, and advanced analysis steps, ensuring parameter traceability,
validation, and pipeline-level consistency.
"""

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple

__all__ = [
    "ClusteringConfig",
    "ResolutionSearchConfig",
    "MergeClustersConfig",
    "AnnotationConfig",
    "DifferentialConfig",
    "FilterMarkersConfig",
    "CompareGroupsConfig",
    "CompareConditionsConfig",
    "ConservedMarkersConfig",
    "EnrichmentConfig",
    "ProportionConfig",
    "ScoringConfig",
    "AnalysisWorkflowConfig",
]


@dataclass
class BaseConfig:
    """Base class for all analysis configuration objects."""

    namespace: str = field(init=False, default="sclucid")

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
    resolution_range: Tuple[float, float, int] = (
        0.2,
        2.0,
        10,
    )  # Start from a slightly higher resolution

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


class ScoringConfig(BaseConfig):
    """Configuration for gene set and cell type scoring."""

    gene_sets: Dict[str, List[str]]
    layer: Optional[str] = "normalized"
    use_raw: bool = True
    ctrl_size: int = 50
    score_name_suffix: str = "_score"
    preserve_missing: bool = True
    min_genes_required: int = 1


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
    pval_cutoff: Optional[float] = (
        None  #: Post-hoc filter for adjusted p-value. Only genes with pvals_adj <= cutoff will be kept.
    )
    fold_change_max: Optional[float] = (
        None  #: Post-hoc clipping for log fold change. Useful for visualization.
    )


@dataclass
class FilterMarkersConfig(BaseConfig):
    """
    Configuration for filtering marker gene results from differential expression analysis.

    This config allows for fine-grained control over marker selection based on statistical
    significance, effect size (log fold change), and expression specificity.
    """

    key: str = "rank_genes_groups"  # The key for the source DE results in adata.uns.
    key_added: Optional[str] = (
        None  # The key under which the filtered DataFrame will be stored. If None, defaults to f"{key}_filtered_df".
    )

    # --- Core Filtering Criteria ---
    min_log2fc: float = 1.0  # Minimum log2 fold change required.
    max_padj: float = 0.05  # Maximum adjusted p-value (e.g., FDR) allowed.
    # All percentage thresholds are specified as fractions (e.g., 0.25 == 25%)
    min_in_group_pct: float = (
        0.25  # Minimum fraction of cells expressing the gene within the target group.
    )

    # --- Specificity Filtering Criteria ---
    max_out_group_pct: Optional[float] = (
        0.5  # Maximum fraction outside the target group.
    )
    min_diff_pct: Optional[float] = 0.1  # Minimum in-group minus out-group fraction.

    # --- Top N Selection ---
    keep_top_n: Optional[int] = (
        100  # After all filters, keep only the top N markers per group.
    )

    # --- Behavior Controls ---
    use_abs_log2fc: bool = False  # If True, filter by |log2FC| >= min_log2fc; otherwise require log2FC >= min_log2fc.
    sort_by: Literal["scores", "logfoldchanges", "diff_pct"] = (
        "scores"  # Column used to rank when keeping top N.
    )


@dataclass
class CompareGroupsConfig(BaseConfig):
    """Configuration for comparing two specific groups (e.g., two cell types)."""

    groupby: str = ""  # The column in adata.obs that contains the groups to compare.
    group1: str = ""  # The name of the first group.
    group2: str = ""  # The name of the second group (used as the reference).
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon"
    layer: Optional[str] = None
    use_raw: bool = False
    key_added: Optional[str] = (
        None  # Key to store the results DataFrame in adata.uns. Defaults to 'compare_{group1}_vs_{group2}'.
    )

    # Filtering parameters for the comparison results
    n_top_genes: int = (
        50  # Number of top genes to keep for each direction after filtering.
    )
    min_log2fc: float = 0.5
    max_padj: float = 0.05
    min_in_group_pct: float = 0.1

    plot: bool = True  # Whether to generate a volcano plot of the results.
    save_dir: Optional[str] = None  # Directory to save the volcano plot.

    # Optional: if you want explicit control (currently hardcoded in function as 0–1 fraction logic)
    # use_fraction_scale: bool = True


@dataclass
class CompareConditionsConfig(BaseConfig):
    """Configuration for comparing the same cell type across two conditions."""

    groupby: str = ""  #: The column identifying the cell type/group to analyze.
    group_name: str = ""  #: The specific cell type/group to subset for comparison.
    condition_key: str = ""  #: The column identifying the conditions (e.g., 'treatment', 'disease_status').
    condition1: str = ""  #: The first condition.
    condition2: str = ""  #: The second condition (will be used as the reference).
    # Inherits comparison and plotting parameters from CompareGroupsConfig
    comparison_params: CompareGroupsConfig = field(default_factory=CompareGroupsConfig)
    key_added: Optional[str] = None  #: Defaults to 'compare_{c1}_vs_{c2}_in_{group}'.


@dataclass
class ConservedMarkersConfig(BaseConfig):
    """Configuration for finding conserved markers across conditions."""

    groupby: str
    condition_key: str
    method: str = "wilcoxon"
    min_cells: int = 10
    min_conditions: Optional[int] = None
    min_log2fc: float = 0.5
    max_padj: float = 0.05
    min_in_group_pct: float = 0.25
    layer: Optional[str] = None
    use_raw: bool = False
    key_added: Optional[str] = None

# ===================== Enrichment & Proportion Configs =====================


@dataclass
class EnrichmentConfig(BaseConfig):
    """
    Configuration for functional enrichment analysis.
    Supports both online (Enrichr) and offline (local GMT files) modes.
    """

    de_key: str = "rank_genes_groups_filtered_df"  # Key for the DE results DataFrame in adata.uns.
    method: Literal["ora", "gsea", "both"] = "ora"
    
    mode: Literal["online", "offline"] = (
        "offline"  # Analysis mode. 'offline' is recommended for stability.
    )
    organism: Literal["human", "mouse"] = "human"  # Species for the analysis.
    gmt_version: str = "v2025"
    max_padj: float = 0.05
    key_added: str = "enrichment"
    plot: bool = True
    n_plot_terms: int = 15
    save_dir: Optional[str] = None
    # --- ORA (Enrichr) 特定参数 ---
    n_top_genes_ora: int = 100
    min_genes_for_ora: int = 10
    # --- GSEA (Prerank) 特定参数 ---
    rank_col_gsea: str = "logfoldchanges"  # or "scores"
    gsea_permutations: int = 100
    gsea_min_size: int = 15
    gsea_max_size: int = 500

    # For 'online' mode, this is a list of Enrichr library names.
    # For 'offline' mode, this is a list of categories (e.g., 'hallmark', 'go_bp')
    # corresponding to local .gmt files in the resources directory.
    gene_sets_online: List[str] = field(
        default_factory=lambda: ["GO_Biological_Process_2023"]
    )
    gene_sets_offline: List[str] = field(
        default_factory=lambda: ["hallmark", "go_bp", "reactome"]
    )

    custom_gene_sets: Optional[str] = (
        None  # Path to a user-provided .gmt file for offline mode.
    )
      # --- New: marker ranking preference for enrichment gene list ---
    prefer_score_for_enrichment: bool = (
        True  # If True and 'scores' exists, rank by 'scores' else by 'logfoldchanges'.
    )


@dataclass
class ProportionConfig(BaseConfig):
    """Configuration for cell type proportion analysis."""
    
    celltype_col: str         # .obs column for cell types (e.g., 'cell_type')
    sample_col: str           # .obs column for sample IDs (e.g., 'sampleID')
    condition_col: str        # .obs column for conditions to compare (e.g., 'disease')
    
    test_method: Literal["t-test", "wilcoxon", "anova", "lmm"] = "wilcoxon"
    
    # --- LMM (Linear Mixed Model) specific params ---
    lmm_batch_key: Optional[str] = None # .obs col for batch/patient ID (e.g., 'patient_id')
    lmm_formula: Optional[str] = None   # e.g., "proportion ~ condition" (auto-builds if None)

    plot_types: List[str] = field(default_factory=lambda: ["bar", "box"])
    out_dir: Optional[str] = None
    figsize: Tuple[float, float] = (10, 6)

# ===================== Master Workflow Config =====================


@dataclass
class AnalysisWorkflowConfig(BaseConfig):
    """Master configuration for the entire analysis workflow."""

    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    de: DifferentialConfig = field(default_factory=DifferentialConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    scoring: Optional[ScoringConfig] = None
    proportion: Optional[ProportionConfig] = None


# Update __all__ to include all configuration classes
__all__ = [
    name
    for name, obj in locals().items()
    if isinstance(obj, type) and issubclass(obj, BaseConfig)
]
