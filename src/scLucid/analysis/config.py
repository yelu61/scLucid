"""
Pydantic-based configuration classes for the analysis module.

Migrates from dataclasses to Pydantic for consistent validation and serialization.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

from pydantic import ConfigDict, Field, field_validator

from ..base_config import SclucidBaseConfig, WorkflowConfigBase

logger = logging.getLogger(__name__)


class ClusteringConfig(SclucidBaseConfig):
    """Configuration for a single clustering run."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["leiden", "louvain", "kmeans", "hdbscan"] = Field(default="leiden")
    resolution: float = Field(default=1.0, gt=0, description="Resolution for leiden/louvain")
    n_clusters: Optional[int] = Field(default=None, ge=2, description="For kmeans")
    use_rep: str = Field(default="X_pca", description="Embedding to use")
    key_added: Optional[str] = Field(default=None, description="Key for adata.obs")
    random_state: int = Field(default=42)
    plot: bool = Field(default=False)
    save_dir: Optional[str] = Field(default=None)
    extra_params: Dict[str, Any] = Field(default_factory=dict)


class ResolutionSearchConfig(SclucidBaseConfig):
    """Configuration for optimizing clustering resolution."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["leiden", "louvain"] = Field(default="leiden")
    use_rep: str = Field(default="X_pca")
    resolution_range: Tuple[float, float, int] = Field(
        default=(0.2, 2.0, 10), description="(min, max, n_points) for resolution search"
    )

    # Metrics to compute
    compute_silhouette: bool = Field(default=True)
    compute_marker_abundance: bool = Field(default=True)
    compute_stability: bool = Field(default=True)
    plot: bool = Field(default=False)
    save_dir: Optional[str] = Field(default=None)

    # Parameters for marker abundance metric
    de_method_for_markers: str = Field(default="wilcoxon")
    min_log2fc_for_markers: float = Field(default=0.5)
    min_pct_for_markers: float = Field(default=0.25)

    @field_validator("resolution_range")
    @classmethod
    def validate_range(cls, v: Tuple[float, float, int]) -> Tuple[float, float, int]:
        """Validate resolution range."""
        if len(v) != 3:
            raise ValueError("resolution_range must be (min, max, n_points)")
        if v[0] >= v[1]:
            raise ValueError("min must be less than max")
        if v[2] < 2:
            raise ValueError("n_points must be at least 2")
        return v


class MergeClustersConfig(SclucidBaseConfig):
    """Configuration for merging similar clusters."""

    model_config = ConfigDict(extra="ignore")

    cluster_key: str = Field(default="leiden_clusters")
    similarity_threshold: float = Field(default=0.8, ge=0, le=1)
    method: Literal["marker_overlap", "expression_correlation"] = Field(default="marker_overlap")
    key_added: Optional[str] = Field(default=None)
    de_method_for_markers: str = Field(default="wilcoxon")


class DifferentialConfig(SclucidBaseConfig):
    """Configuration for one-vs-rest differential expression analysis."""

    model_config = ConfigDict(extra="ignore")

    groupby: str = Field(default="leiden_clusters")
    method: Literal["wilcoxon", "t-test", "logreg"] = Field(default="wilcoxon")
    layer: Optional[str] = Field(default=None)
    use_raw: bool = Field(default=False)
    key_added: Optional[str] = Field(default="rank_genes_groups")
    groups: Optional[List[str]] = Field(default=None)
    reference: str = Field(default="rest")
    pval_cutoff: Optional[float] = Field(default=None, ge=0, le=1)
    fold_change_max: Optional[float] = Field(default=None)


class FilterMarkersConfig(SclucidBaseConfig):
    """Configuration for filtering marker gene results."""

    model_config = ConfigDict(extra="ignore")

    key: str = Field(default="rank_genes_groups", description="Source DE results key")
    key_added: Optional[str] = Field(default=None)

    # Core filtering criteria
    min_log2fc: float = Field(default=1.0, ge=0)
    max_padj: float = Field(default=0.05, gt=0, le=1)
    min_in_group_pct: float = Field(default=0.25, ge=0, le=1)

    # Specificity filtering
    max_out_group_pct: Optional[float] = Field(default=0.5, ge=0, le=1)
    min_diff_pct: Optional[float] = Field(default=0.1, ge=0, le=1)

    # Top N selection
    keep_top_n: Optional[int] = Field(default=100, ge=1)

    # Behavior controls
    use_abs_log2fc: bool = Field(default=False)
    sort_by: Literal["scores", "logfoldchanges", "diff_pct"] = Field(default="scores")


class ComparisonConfig(SclucidBaseConfig):
    """Base configuration for differential expression comparisons."""

    model_config = ConfigDict(extra="ignore")

    method: Literal["wilcoxon", "t-test", "logreg"] = Field(default="wilcoxon")
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = Field(default="benjamini-hochberg")
    layer: Optional[str] = Field(default=None)
    use_raw: bool = Field(default=False)
    key_added: str = Field(default="rank_genes_groups")

    # Filtering thresholds
    min_log2fc: Optional[float] = Field(default=0.5, ge=0)
    max_padj: Optional[float] = Field(default=0.05, gt=0, le=1)
    min_pct: Optional[float] = Field(default=0.1, ge=0, le=1)
    use_abs_log2fc: bool = Field(default=False)

    # Advanced parameters
    tie_correct: bool = Field(default=True)
    rankby_abs: bool = Field(default=False)
    pts: bool = Field(default=True)
    n_genes: int = Field(default=5000, ge=1)

    # Exploratory-use acknowledgement. Cell-level DE treats cells as independent
    # observations and is not valid for formal condition inference. Set this to
    # True only after reviewing the pseudoreplication warning.
    acknowledge_exploratory: bool = Field(default=False)


class CompareGroupsConfig(ComparisonConfig):
    """Configuration for comparing two specific groups."""

    groupby: str = Field(default="")
    group1: str = Field(default="")
    group2: str = Field(default="")
    n_top_genes: int = Field(default=50, ge=1)


class CompareConditionsConfig(ComparisonConfig):
    """Configuration for comparing same cell type across conditions."""

    groupby: str = Field(default="")
    group_name: str = Field(default="")
    condition_key: str = Field(default="")
    condition1: str = Field(default="")
    condition2: str = Field(default="")
    key_added: Optional[str] = Field(default=None)
    n_top_genes: int = Field(default=50, ge=1)


class PseudobulkDEConfig(SclucidBaseConfig):
    """Configuration for sample-level pseudobulk differential expression."""

    model_config = ConfigDict(extra="ignore")

    sample_col: str = Field(description="Column containing biological sample identifiers")
    condition_key: str = Field(description="Column containing condition labels")
    contrasts: List[Tuple[str, str]] = Field(
        description="List of (condition1, condition2) contrasts; log2FC is condition2 - condition1"
    )
    groupby: Optional[str] = Field(
        default=None,
        description="Optional cell group column; run each contrast within each selected group",
    )
    group_names: Optional[List[str]] = Field(default=None)
    layer: Optional[str] = Field(default=None, description="Layer containing raw counts")
    use_raw: bool = Field(default=False)
    min_cells_per_sample: int = Field(default=5, ge=1)
    min_counts: int = Field(default=10, ge=0)
    min_samples_per_condition: int = Field(default=1, ge=1)
    pseudocount: float = Field(default=0.5, gt=0)
    method: Literal[
        "auto",
        "deseq2",
        "welch_logcpm",
        "linear_model_logcpm",
        "statsmodels_glm",
        "statsmodels_gee",
        "cell_level_fallback",
    ] = Field(default="auto")
    design_covariates: List[str] = Field(
        default_factory=list,
        description=(
            "Sample-level covariates to include in the Python logCPM linear model "
            "path, for example batch or patient identifiers."
        ),
    )
    robust_cov_type: Literal["HC0", "HC1", "HC2", "HC3", "nonrobust"] = Field(
        default="HC3",
        description=(
            "Covariance estimator for linear_model_logcpm. HC3 is the default "
            "heteroscedasticity-robust standard error; use nonrobust to keep "
            "ordinary OLS standard errors."
        ),
    )
    statsmodels_glm_family: Literal["NegativeBinomial", "Gaussian", "Gamma"] = Field(
        default="NegativeBinomial",
        description=(
            "GLM family for the statsmodels_glm backend. NegativeBinomial uses raw "
            "counts; Gaussian uses logCPM with identity link and HC3 sandwich "
            "covariance; Gamma uses counts with log link."
        ),
    )
    hierarchical_correction: bool = Field(
        default=False,
        description=(
            "Apply hierarchical multiple-testing correction: BH within each "
            "(group, contrast) and Benjamini-Yekutieli across all tests."
        ),
    )
    block_col: Optional[str] = Field(
        default=None,
        description=(
            "Optional sample-level blocking column, typically patient or donor. "
            "It is included as a categorical covariate in linear_model_logcpm."
        ),
    )
    experimental_unit_col: Optional[str] = Field(
        default=None,
        description=(
            "Column identifying independent biological units. When omitted, "
            "block_col is used for paired/repeated designs and sample_col otherwise. "
            "Replicate sufficiency is evaluated on this column, never on cells."
        ),
    )
    fallback_to_cell_level: bool = Field(
        default=False,
        description=(
            "Allow exploratory cell-level fallback when sample-level replicates are insufficient. "
            "Disabled by default to avoid pseudoreplication being mistaken for formal inference."
        ),
    )
    single_sample_mode: Literal["descriptive", "skip"] = Field(
        default="descriptive",
        description="Behavior when a contrast has only one biological sample per condition.",
    )
    p_adjust_method: Literal["fdr_bh", "bonferroni"] = Field(default="fdr_bh")
    key_added: str = Field(default="pseudobulk_de")
    n_jobs: int = Field(default=1, ge=1)

    @field_validator("contrasts")
    @classmethod
    def validate_contrasts(cls, value: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
        """Require explicit, directional, non-duplicated two-level contrasts."""
        if not value:
            raise ValueError("contrasts must contain at least one condition comparison")
        normalized = [(str(condition1), str(condition2)) for condition1, condition2 in value]
        if any(condition1 == condition2 for condition1, condition2 in normalized):
            raise ValueError("A pseudobulk contrast must compare two distinct conditions")
        if len(set(normalized)) != len(normalized):
            raise ValueError("Duplicate pseudobulk contrasts are not allowed")
        return normalized


class ConservedMarkersConfig(ComparisonConfig):
    """Configuration for conserved marker gene discovery."""

    groupby: str = Field(default="")
    condition_key: str = Field(default="")
    min_cells: int = Field(default=10, ge=1)
    min_conditions: Optional[int] = Field(default=None)


class EnrichmentConfig(SclucidBaseConfig):
    """Configuration for functional enrichment analysis."""

    model_config = ConfigDict(extra="ignore")

    # Data source
    de_key: str = Field(default="rank_genes_groups_filtered_df")
    method: Literal["ora", "gsea", "both"] = Field(default="ora")
    mode: Literal["online", "offline"] = Field(default="offline")

    # Organism and gene sets
    organism: Literal["human", "mouse"] = Field(default="human")
    gene_sets_online: List[str] = Field(default_factory=lambda: ["GO_Biological_Process_2023"])
    gene_sets_offline: List[str] = Field(default_factory=lambda: ["hallmark", "go_bp", "reactome"])
    gmt_version: str = Field(default="v2025")
    custom_gene_sets: Optional[str] = Field(default=None)

    # General parameters
    max_padj: float = Field(default=0.05, gt=0, le=1)
    cutoff_pval: float = Field(default=0.05, gt=0, le=1)
    key_added: str = Field(default="enrichment")
    n_plot_terms: int = Field(default=15, ge=1)

    # ORA parameters
    n_top_genes_ora: int = Field(default=100, ge=1)
    min_genes_for_ora: int = Field(default=10, ge=1)
    background: Optional[List[str]] = Field(default=None)

    # GSEA parameters
    rank_col_gsea: str = Field(default="logfoldchanges")
    gsea_permutations: int = Field(default=1000, ge=100)
    gsea_min_size: int = Field(default=15, ge=1)
    gsea_max_size: int = Field(default=500, ge=1)
    gsea_seed: int = Field(default=42, description="Random seed for GSEA permutations")
    gsea_processes: int = Field(default=4, ge=1, description="Number of processes for GSEA")
    strict_mode: bool = Field(
        default=False, description="If True, raise on enrichment failures instead of skipping"
    )
    prefer_score_for_enrichment: bool = Field(default=True)


class AnnotationConfig(SclucidBaseConfig):
    """Configuration for the cell type annotation workflow."""

    model_config = ConfigDict(extra="ignore")

    cluster_key: str = Field(default="leiden_clusters")
    marker_species: str = Field(default="human")
    marker_tissue: Optional[str] = Field(default=None)
    marker_states: Optional[List[str]] = Field(default=None)
    lineage_marker_config: Optional[str] = Field(default=None)
    subtype_marker_config: Optional[str] = Field(default=None)
    state_marker_config: Optional[str] = Field(default=None)
    target_lineage: Optional[str] = Field(default=None)
    run_celltypist: bool = Field(default=False)
    celltypist_model: str = Field(default="Immune_All_Low.pkl")
    run_scoring: bool = Field(default=True)
    final_method: Literal[
        "max_score", "enrichment", "combined", "celltypist", "hybrid", "hierarchical"
    ] = Field(default="combined")
    marker_method: Literal["max_score", "enrichment", "combined"] = Field(default="combined")
    key_added: str = Field(default="cell_type_auto")
    lineage_key: str = Field(default="celltype_lineage_auto")
    subtype_key: str = Field(default="celltype_subtype_auto")
    state_key: str = Field(default="celltype_state_auto")
    annotation_basis_key: str = Field(default="celltype_annotation_basis")
    nomenclature_style: Literal["legacy", "modular"] = Field(default="modular")
    min_confidence: float = Field(default=0.1, ge=0, le=1)
    celltypist_confidence_threshold: float = Field(default=0.5, ge=0, le=1)
    state_signature_names: List[str] = Field(default_factory=list)
    state_signature_categories: List[str] = Field(default_factory=list)
    custom_state_signatures: Optional[Dict[str, List[str]]] = Field(default=None)
    custom_state_metadata: Optional[Dict[str, Dict[str, Any]]] = Field(default=None)
    state_score_suffix: str = Field(default="_program")
    min_state_score: float = Field(
        default=0.1,
        ge=0,
        description="Minimum top state signature score required to assign a state.",
    )

    # Reporting and export
    report: bool = Field(
        default=False, description="Whether to export annotation report visualizations."
    )
    save_dir: Optional[str] = Field(default=None)

    # DE and enrichment parameters
    min_log2fc: float = Field(default=0.5, ge=0)
    min_in_group_pct: float = Field(default=0.2, ge=0, le=1)
    enrichment_gene_sets: List[str] = Field(default_factory=lambda: ["GO_Biological_Process_2023"])

    # Compartment mapping for tumor-aware annotation
    compartment_map: Optional[Dict[str, str]] = Field(default=None)


class ProportionConfig(SclucidBaseConfig):
    """Configuration for cell type proportion analysis."""

    model_config = ConfigDict(extra="ignore")

    # Required fields
    celltype_col: str = Field(description="Column for cell types")
    sample_col: str = Field(description="Column for sample IDs")
    condition_col: str = Field(description="Column for conditions")

    # Optional fields
    pairing_col: Optional[str] = Field(default=None)
    experimental_unit_col: Optional[str] = Field(
        default=None,
        description=(
            "Column identifying independent biological units. Defaults to pairing_col "
            "for repeated designs and sample_col otherwise."
        ),
    )
    batch_col: Optional[str] = Field(default=None)
    timepoint_col: Optional[str] = Field(default=None)

    # scCODA reference selection
    reference_cell_type: Optional[str] = Field(default="auto")

    auto_configure: bool = Field(default=True)
    test_method: Literal[
        "deseq2",
        "clr-t-test",
        "clr-wilcoxon",
        "clr-paired-t-test",
        "clr-paired-wilcoxon",
        "clr-ols",
        "ancom-like-clr",
        "t-test",
        "wilcoxon",
        "anova",
        "kruskal",
        "chi-square",
        "paired-t-test",
        "paired-wilcoxon",
    ] = Field(default="clr-t-test")
    composition_transform: Literal["clr", "none"] = Field(default="clr")
    composition_pseudocount: float = Field(default=1e-6, gt=0)
    require_biological_replicates: bool = Field(default=True)
    min_samples_per_condition: int = Field(default=2, ge=1)
    legacy_exploratory: bool = Field(default=False)
    correction_scope: str = Field(default="per_test")

    # Plotting
    plot_types: List[str] = Field(
        default_factory=lambda: ["counts", "bar", "bar_composition", "box", "alluvial", "diff"]
    )
    analysis_level: Literal["sample", "group"] = Field(default="group")
    sample_plot_style: Literal["stacked", "grouped"] = Field(default="grouped")
    volcano_effect_col: str = Field(default="cohen_d")
    effect_size_threshold: float = Field(default=0.5)

    # Palettes and ordering
    ct_palette: Optional[Dict[str, str]] = Field(default=None)
    condition_palette: Optional[Dict[str, str]] = Field(default=None)
    sample_palette: Optional[Dict[str, str]] = Field(default=None)
    celltype_order: Optional[List[str]] = Field(default=None)
    condition_order: Optional[List[str]] = Field(default=None)
    timepoint_order: Optional[List[str]] = Field(default=None)

    out_dir: Optional[str] = Field(default=None)
    figsize: Tuple[float, float] = Field(default=(8, 6))
    export_data: bool = Field(default=True)
    export_format: str = Field(default="csv")


class ScoringConfig(SclucidBaseConfig):
    """Configuration for gene set and cell type scoring."""

    model_config = ConfigDict(extra="ignore")

    gene_sets: Optional[Union[Dict[str, List[str]], str]] = Field(default=None)
    species: str = Field(default="human")
    custom_signatures: Optional[str] = Field(default=None)

    # Scoring parameters
    layer: Optional[str] = Field(default="normalized")
    use_raw: bool = Field(default=False)
    ctrl_size: int = Field(default=50, ge=1)
    score_name_suffix: str = Field(default="_score")
    preserve_missing: bool = Field(default=True)
    min_genes_required: int = Field(default=1, ge=1)

    # Visualization
    plot_heatmap: bool = Field(default=True)
    z_score: bool = Field(default=True)
    heatmap_cmap: str = Field(default="RdBu_r")


class AnalysisWorkflowConfig(WorkflowConfigBase):
    """Master configuration for the entire analysis workflow."""

    model_config = ConfigDict(extra="ignore")

    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    de: DifferentialConfig = Field(default_factory=DifferentialConfig)
    annotation: Optional[AnnotationConfig] = Field(default_factory=AnnotationConfig)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)
    scoring: Optional[ScoringConfig] = Field(default=None)
    proportion: Optional[ProportionConfig] = Field(default=None)

    # Evidence-first analysis controls
    run_clustering_review: bool = Field(
        default=True,
        description=(
            "Run a lightweight clustering-resolution review before final clustering. "
            "Users can disable it for speed with run_clustering_review=False."
        ),
    )
    candidate_resolutions: Optional[List[float]] = Field(default=None)
    use_recommended_resolution: bool = Field(default=False)
    run_annotation_evidence: bool = Field(default=True)
    annotation_methods: Tuple[str, ...] = Field(
        default=("reference", "marker_manager", "data_driven")
    )
    final_annotation_strategy: Literal["consensus", "legacy"] = Field(default="consensus")
    annotation_level: Literal["lineage", "cell_type", "subtype", "state"] = Field(default="lineage")
    llm_annotations: Optional[Union[Dict[str, Any], List[Dict[str, Any]]]] = Field(default=None)
    find_markers: bool = Field(default=True)
    characterize: bool = Field(default=True)
    marker_method: Literal["wilcoxon", "t-test", "logreg"] = Field(default="wilcoxon")
    run_proportion: bool = Field(
        default=False,
        description=(
            "Run cell type proportion analysis as part of the standard workflow. "
            "Requires a ProportionConfig or sample/condition/cell-type columns "
            "matching the default names."
        ),
    )
    proportion_method: Optional[str] = Field(
        default=None,
        description="Optional method override for analyze_celltype_proportion().",
    )

    # Pseudobulk-first workflow routing. When True, run_standard_analysis will
    # run clustering and annotation, then immediately aggregate to sample-level
    # pseudobulk DE per cell type instead of using cell-level condition DE.
    pseudobulk_first: bool = Field(
        default=False,
        description=(
            "After standard clustering and annotation, run sample-level pseudobulk "
            "DE per cell type as the primary differential output."
        ),
    )

    # Deprecated tumor-specific controls kept for backward compatibility.
    # Tumor interpretation should be managed via ``run_tumor_analysis`` or by
    # attaching callables to ``post_analysis_hooks`` instead.
    run_malignancy_interpretation: bool = Field(
        default=False,
        deprecated=True,
        description="Deprecated: use run_tumor_analysis() or post_analysis_hooks.",
    )
    run_cnv_for_malignancy: bool = Field(default=False)
    run_malignancy_score: bool = Field(default=True)
    malignancy_cancer_type: Optional[str] = Field(default=None)
    malignancy_reference_labels: Optional[List[str]] = Field(default=None)
    malignancy_cnv_score_key: Optional[str] = Field(default=None)
    malignancy_key_added: str = Field(default="malignancy_call")
    malignancy_score_key: str = Field(default="malignancy_interpretation_score")
    malignancy_threshold: float = Field(default=0.55, ge=0, le=1)
    malignancy_suspect_threshold: float = Field(default=0.35, ge=0, le=1)

    # Optional post-analysis hooks called after the main analysis workflow completes.
    # Each hook receives (adata, config) and returns adata. Useful for tumor modules
    # that need annotation outputs without creating a hard import dependency.
    post_analysis_hooks: Optional[List[Any]] = Field(
        default=None,
        description="List of callables(adata, config) -> adata invoked after analysis.",
    )

    # Note: save_dir is inherited from SclucidBaseConfig

    @classmethod
    def from_simple_dict(cls, simple_config: Dict[str, Any]) -> AnalysisWorkflowConfig:
        """
        Create AnalysisWorkflowConfig from a simplified flat dictionary.

        Args:
            simple_config: Flat dictionary with keys like:
                - clustering_method, clustering_resolution
                - annotation_method, marker_method
                - save_dir, n_jobs

        Returns:
            AnalysisWorkflowConfig: Fully configured workflow config
        """
        config_data = dict(simple_config)
        kwargs: Dict[str, Any] = {}

        # Extract clustering parameters
        clustering_params = {}
        for key in ["method", "resolution", "n_clusters", "use_rep"]:
            config_key = f"clustering_{key}"
            if config_key in config_data:
                clustering_params[key] = config_data.pop(config_key)
        if clustering_params:
            kwargs["clustering"] = ClusteringConfig(**clustering_params)

        # Extract annotation parameters
        annotation_params = {}
        for key in ["cluster_key", "marker_species", "final_method", "run_scoring"]:
            config_key = f"annotation_{key}"
            if config_key in config_data:
                annotation_params[key] = config_data.pop(config_key)
        if annotation_params:
            kwargs["annotation"] = AnnotationConfig(**annotation_params)

        # Extract DE parameters
        de_params = {}
        for key in ["groupby", "method"]:
            config_key = f"de_{key}"
            if config_key in config_data:
                de_params[key] = config_data.pop(config_key)
        if de_params:
            kwargs["de"] = DifferentialConfig(**de_params)

        # Backward compatibility: results_dir -> save_dir
        if "results_dir" in config_data:
            config_data["save_dir"] = config_data.pop("results_dir")

        # Remaining keys go directly to workflow config
        kwargs.update(config_data)
        return cls(**kwargs)

    @classmethod
    def quick(
        cls,
        clustering_method: str = "leiden",
        resolution: float = 1.0,
        run_annotation: bool = True,
        **kwargs,
    ) -> AnalysisWorkflowConfig:
        """
        Quick configuration factory for standard analyses.

        Args:
            clustering_method: Clustering algorithm
            resolution: Clustering resolution
            run_annotation: Whether to run annotation step
            **kwargs: Additional parameters (save_dir, n_jobs, etc.)

        Returns:
            AnalysisWorkflowConfig: Pre-configured for standard analysis
        """
        return cls(
            clustering=ClusteringConfig(method=clustering_method, resolution=resolution),
            annotation=AnnotationConfig() if run_annotation else None,
            **kwargs,
        )


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
