"""
Configuration classes for the analysis module of scLucid.

Defines all dataclasses for clustering, annotation, scoring, differential expression,
enrichment analysis, and advanced analysis steps, ensuring parameter traceability,
validation, and pipeline-level consistency.
"""

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
from typing import TYPE_CHECKING
from pathlib import Path
import logging
from importlib import resources

if TYPE_CHECKING:
    from ..analysis.scoring import FunctionalSignatureManager


__all__ = [
    "ClusteringConfig",
    "ResolutionSearchConfig",
    "MergeClustersConfig",
    "DifferentialConfig",
    "FilterMarkersConfig",
    "ComparisonConfig",
    "CompareGroupsConfig",
    "CompareConditionsConfig",
    "ConservedMarkersConfig",
    "EnrichmentConfig",
    "AnnotationConfig",
    "ProportionConfig",
    "ScoringConfig",
    "AnalysisWorkflowConfig",
]
log = logging.getLogger(__name__)

@dataclass
class BaseConfig:
    """Base class for all analysis configuration objects."""

    namespace: str = field(init=False, default="sclucid")
    
    verbose: bool = True
    plot: bool = False
    save_dir: Optional[str] = None

    def __post_init__(self):
        if self.save_dir:
            Path(self.save_dir).mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]):
        return cls(**d)

    def validate(self):
        """Subclasses can override this method for parameter validation"""
        pass


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
class ComparisonConfig(BaseConfig):
    """Base configuration class for differential expression analysis"""
    # Differential expression method
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon"
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg"

    # Data source
    layer: Optional[str] = None
    use_raw: bool = False

    # Result storage
    key_added: str = "rank_genes_groups"

    # Filtering thresholds
    min_log2fc: Optional[float] = 0.5
    max_padj: Optional[float] = 0.05
    min_pct: Optional[float] = 0.1
    use_abs_log2fc: bool = False

    # Advanced parameters
    tie_correct: bool = True
    rankby_abs: bool = False
    pts: bool = True
    n_genes: int = 5000

    def validate(self):
        if self.min_log2fc is not None and self.min_log2fc < 0:
            raise ValueError("min_log2fc must be non-negative")
        if self.max_padj is not None and not (0 < self.max_padj <= 1):
            raise ValueError("max_padj must be in (0, 1]")
        if self.min_pct is not None and not (0 <= self.min_pct <= 1):
            raise ValueError("min_pct must be in [0, 1]")   
        

@dataclass
class CompareGroupsConfig(ComparisonConfig):
    """Configuration for comparing two specific groups (e.g., two cell types)."""
    groupby: str = ""
    group1: str = ""
    group2: str = ""
    n_top_genes: int = 50


@dataclass
class CompareConditionsConfig(ComparisonConfig):
    """Configuration for comparing the same cell type across two conditions."""

    groupby: str = ""  #: The column identifying the cell type/group to analyze.
    group_name: str = ""  #: The specific cell type/group to subset for comparison.
    condition_key: str = ""  #: The column identifying the conditions (e.g., 'treatment', 'disease_status').
    condition1: str = ""  #: The first condition.
    condition2: str = ""  #: The second condition (will be used as the reference).
    key_added: Optional[str] = None  #: Defaults to 'compare_{c1}_vs_{c2}_in_{group}'.


@dataclass
class ConservedMarkersConfig(ComparisonConfig):
    """Configuration for conserved marker gene discovery"""

    groupby: str = ""
    condition_key: str = ""
    min_cells: int = 10
    min_conditions: Optional[int] = None


# ===================== Enrichment Configs =====================


@dataclass
class EnrichmentConfig(BaseConfig):
    """
    Configuration for functional enrichment analysis.
    Supports both online (Enrichr) and offline (local GMT files) modes.
    """
    # Data source
    de_key: str = "rank_genes_groups_filtered_df"  # Key for the DE results DataFrame in adata.uns.
    # Analysis mode
    method: Literal["ora", "gsea", "both"] = "ora"
    mode: Literal["online", "offline"] = (
        "offline"  # Analysis mode. 'offline' is recommended for stability.
    )
    # Organism and gene sets
    organism: Literal["human", "mouse"] = "human"  # Species for the analysis.
    # For 'online' mode, this is a list of Enrichr library names.
    gene_sets_online: List[str] = field(
        default_factory=lambda: ["GO_Biological_Process_2023"]
    )
    # For 'offline' mode, this is a list of categories (e.g., 'hallmark', 'go_bp')
    gene_sets_offline: List[str] = field(
        default_factory=lambda: ["hallmark", "go_bp", "reactome"]
    )
    # corresponding to local .gmt files in the resources directory.
    gmt_version: str = "v2025"
    
    custom_gene_sets: Optional[str] = (
        None  # Path to a user-provided .gmt file for offline mode.
    )
    # General parameters
    max_padj: float = 0.05
    cutoff_pval: float = 0.05
    key_added: str = "enrichment"
    plot: bool = True
    n_plot_terms: int = 15
    save_dir: Optional[str] = None
    # --- ORA (Enrichr) parameters ---
    n_top_genes_ora: int = 100
    min_genes_for_ora: int = 10
    background: Optional[List[str]] = None
    # --- GSEA (Prerank) parameters ---
    rank_col_gsea: str = "logfoldchanges"  # or "scores"
    gsea_permutations: int = 1000
    gsea_min_size: int = 15
    gsea_max_size: int = 500
    prefer_score_for_enrichment: bool = (
        True  # If True and 'scores' exists, rank by 'scores' else by 'logfoldchanges'.
    )
        
    def __post_init__(self):
        super().__post_init__()
        
        if self.mode == "offline" and not self.custom_gene_sets:
            log.info(f"Offline mode: using default GMT files for {self.organism}")
    

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
    

# ===================== Enrichment & Proportion Configs =====================


@dataclass(kw_only=True)
class ProportionConfig(BaseConfig):
    """Configuration for cell type proportion analysis."""
    # Required fields (no default values)
    celltype_col: str  # .obs column for cell types (e.g., 'cell_type')
    sample_col: str  # .obs column for sample IDs (e.g., 'sampleID')
    condition_col: str  # .obs column for conditions to compare (e.g., 'disease')
    
    # Optional fields (with default values)
    pairing_col: Optional[str] = None  
    """Column for paired analysis (e.g., 'patient_id' for before/after studies).
    If specified, enables paired statistical tests (e.g., paired t-test, Wilcoxon signed-rank)."""
    batch_col: Optional[str] = None  
    """Column for batch effect detection (e.g., 'sequencing_batch').
    Used to identify and control for technical variation in proportion analysis."""
    timepoint_col: Optional[str] = None  
    """Column for time-series analysis (e.g., 'day', 'week').
    Enables longitudinal analysis and trajectory visualization."""
    
    auto_configure: bool = True  
    """If True, automatically detect data structure and select appropriate statistical tests."""
    
    test_method: Literal["deseq2", "t-test", "wilcoxon", "anova", "chi-square", "fisher"] = "deseq2"
    correction_scope: str = "per_test"  # or "global" for multiple comparisons correction.

    # --- Plotting Configuration ---
    # Allowed types: 'counts', 'bar', 'bar_composition', 'box', 'alluvial', 'diff'
    plot_types: List[str] = field(
        default_factory=lambda: [
            "counts",
            "bar",
            "bar_composition",
            "box",
            "alluvial",
            "diff",
        ]
    )
    # "group" calls aggregation/comparison plots; "sample" shows individual sample details
    analysis_level: Literal["sample", "group"] = "group"
    # Sample模式下的 counts 图样式:
    # 'stacked' = 堆叠图 (看总数), 'grouped' = 簇状图 (看细胞类型内的样本差异)
    sample_plot_style: Literal["stacked", "grouped"] = "grouped"
    # Volcano plot settings
    volcano_effect_col: str = "cohen_d"  # or "cliff_delta"
    effect_size_threshold: float = 0.5  # For volcano plots

    # --- Palettes separated by type ---
    # Colors for Cell Types (used in Stacked Bars, Alluvial flows)
    ct_palette: Optional[Dict[str, str]] = None
    # Colors for Conditions (used in Boxplots, Diff bars)
    condition_palette: Optional[Dict[str, str]] = None
    # Colors for Samples (used in Sample-Grouped Counts)
    sample_palette: Optional[Dict[str, str]] = None

    # Order of cell types for plotting (e.g., ['T cells', 'B cells', 'NK cells'])
    celltype_order: Optional[List[str]] = None
    condition_order: Optional[List[str]] = None
    timepoint_order: Optional[List[str]] = None  # For longitudinal studies

    out_dir: Optional[str] = None
    figsize: Tuple[float, float] = (8, 6)
    export_data: bool = True
    export_format: str = "csv"  # or "excel", "parquet"


# ===================== Scoring Configs =====================


@dataclass
class ScoringConfig(BaseConfig):
    """Configuration for gene set and cell type scoring."""
    
    # Gene set source: can be dict, path string, or FunctionalSignatureManager
    gene_sets: Optional[Union[Dict[str, List[str]], str, "FunctionalSignatureManager"]] = None
    species: str = "human"
    custom_signatures: Optional[str] = None
    
    # Scoring parameters
    layer: Optional[str] = "log1p_norm"
    use_raw: bool = True
    ctrl_size: int = 50
    score_name_suffix: str = "_score"
    preserve_missing: bool = True
    min_genes_required: int = 1
    
    # Visualization
    plot_heatmap: bool = True
    z_score: bool = True
    heatmap_cmap: str = "RdBu_r"
    
    
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