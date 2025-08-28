"""
Analysis module for single-cell RNA-seq data.

This module provides functions for dimensionality reduction, clustering,
differential expression analysis, cell type annotation, and cell type proportion statistics.
"""

__version__ = "1.0.0"

# Import and expose key functions from submodules
from .annotation import (
    score_cell_types,
    annotate_clusters,
    run_celltypist,
    transfer_labels,
    evaluate_annotation,
    summarize_annotation_evidence,
    apply_annotation_mapping,
    run_annotation
)
from .clustering import (
    cluster_cells,
    find_resolution,
    merge_clusters,
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
    score_by_gene_sets,
    compare_scores,
    plot_score_comparison,
    batch_compare_scores,
    batch_plot_score_comparison,
)
from .proportion import (
    compute_celltype_proportion,
    plot_celltype_proportion,
    celltype_proportion_test,
    celltype_proportion_analysis
)

__all__ = [
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
    "score_by_gene_sets",
    "compare_scores",
    "plot_score_comparison",
    "batch_compare_scores",
    "batch_plot_score_comparison",
    # Proportion (cell type composition analysis)
    "compute_celltype_proportion",
    "plot_celltype_proportion",
    "celltype_proportion_test",
    "celltype_proportion_analysis",
]