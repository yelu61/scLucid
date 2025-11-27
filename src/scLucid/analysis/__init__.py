"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, cell type annotation, and cell type proportion statistics.
"""

__version__ = "1.0.0"

# Import and expose key functions from submodules
from .config import (
    ResolutionSearchConfig,
    ClusteringConfig,
    AnnotationConfig,
    ScoringConfig,
    DifferentialConfig,
    FilterMarkersConfig,
    CompareGroupsConfig,
    CompareConditionsConfig,
    EnrichmentConfig,
    ProportionConfig,
    AnalysisWorkflowConfig,
)
from .clustering import (
    find_resolution,
    cluster_cells,
    merge_clusters,
)
from .annotation import (
    score_cell_types,
    annotate_clusters,
    run_celltypist,
    transfer_labels,
    evaluate_annotation,
    summarize_annotation_evidence,
    apply_annotation_mapping,
    remap_labels,
    run_annotation
)
from .de_enrichment import (
    find_markers,
    filter_markers,
    compare_groups,
    compare_conditions,
    get_conserved_markers,
    run_enrichment,
    summarize_markers_and_enrichment,
    characterize_clusters,
    visualize_markers
)
from .scoring import (
    FunctionalSignatureManager,
    score_by_gene_sets,
    calculate_signature_matrix,
    plot_signature_heatmap,
    plot_delta_heatmap,
    batch_plot_delta_heatmap,
    plot_score_violin_with_stats,
    batch_compare_scores,
)
from .proportion import (
    plot_cell_counts,
    plot_proportion_bar,
    plot_diff_stats,
    plot_composition,
    plot_box_summary,
    plot_proportion_shifts,
    plot_individual_boxplots,
    plot_proportion_heatmap,
    plot_proportion_with_ci,
    plot_celltype_correlation,
    run_statistical_test,
    compute_celltype_proportion,
    celltype_proportion_analysis,
)

__all__ = [
    "ResolutionSearchConfig",
    "ClusteringConfig",
    "MergeClustersConfig",
    "AnnotationConfig",
    "ScoringConfig",
    "DifferentialConfig",
    "FilterMarkersConfig",
    "CompareGroupsConfig",
    "CompareConditionsConfig",
    "EnrichmentConfig",
    "ProportionConfig",
    "AnalysisWorkflowConfig",
    # Clustering and dimensionality reduction
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
    # Annotation
    "score_cell_types",
    "annotate_clusters",
    "run_celltypist",
    "transfer_labels",
    "evaluate_annotation",
    "summarize_annotation_evidence",
    "apply_annotation_mapping",
    "remap_labels",
    "run_annotation",
    # Differential expression and enrichment
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    "run_enrichment",
    "summarize_markers_and_enrichment",
    "characterize_clusters",
    "visualize_markers",
    # Scoring
    "FunctionalSignatureManager",
    "score_by_gene_sets",
    "calculate_signature_matrix",
    "plot_signature_heatmap",
    "plot_delta_heatmap",
    "batch_plot_delta_heatmap",
    "plot_score_violin_with_stats",
    "batch_compare_scores",
    # Proportion (cell type composition analysis)
    "celltype_proportion_analysis",
    "compute_celltype_proportion",
    "run_statistical_test",
    "plot_cell_counts",
    "plot_proportion_bar",
    "plot_diff_stats",
    "plot_composition",
    "plot_box_summary",
    "plot_proportion_shifts",
    "plot_individual_boxplots",
    "plot_proportion_heatmap",
    "plot_proportion_with_ci",
    "plot_celltype_correlation",
]