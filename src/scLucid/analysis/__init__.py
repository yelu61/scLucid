"""Analysis public API for scLucid."""

import warnings
from collections.abc import Iterable
from importlib import import_module

__version__ = "1.0.0"
__all__ = []


def _export(module: str, names: Iterable[str], *, optional: bool = False) -> bool:
    """Import names from a submodule without breaking module import."""
    try:
        loaded = import_module(f"{__name__}.{module}")
    except Exception as exc:
        level = "optional" if optional else "required"
        warnings.warn(
            f"Could not import {level} analysis module '{module}': {exc}",
            ImportWarning,
        )
        return False

    found = False
    for name in names:
        if hasattr(loaded, name):
            globals()[name] = getattr(loaded, name)
            __all__.append(name)
            found = True
    return found


_export(
    "config",
    [
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
    ],
)
_export("clustering", ["find_resolution", "cluster_cells", "merge_clusters"])
_export(
    "annotation",
    [
        "score_cell_types",
        "annotate_clusters",
        "run_celltypist",
        "transfer_labels",
        "evaluate_annotation",
        "summarize_annotation_evidence",
        "run_lineage_state_annotation",
        "filter_marker_table_for_annotation",
        "flag_suspect_clusters",
        "build_annotation_review_table",
        "apply_annotation_mapping",
        "remap_labels",
        "run_annotation",
    ],
)

# Prefer the reorganized DE package, but keep a legacy fallback.
_de_names = [
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    "run_enrichment",
    "export_enrichment_results",
    "batch_celltype_deg_enrichment",
    "summarize_markers_and_enrichment",
    "characterize_clusters",
    "visualize_markers",
    "plot_volcano",
    "plot_multi_cluster_deg",
    "ResultManager",
    "save_results",
    "load_results",
]
if not _export("differential_expression", _de_names, optional=True):
    _export("de_enrichment", _de_names, optional=True)

_export(
    "scoring",
    [
        "FunctionalSignatureManager",
        "score_by_gene_sets",
        "run_module_scoring_workflow",
        "calculate_signature_matrix",
        "plot_signature_heatmap",
        "plot_delta_heatmap",
        "batch_plot_delta_heatmap",
        "plot_score_violin_with_stats",
        "batch_compare_scores",
    ],
)
_export(
    "workflow",
    [
        "run_standard_analysis",
        "run_custom_analysis",
        "compare_clustering_resolutions",
        "AnalysisWorkflowError",
        "PartialAnalysisResult",
    ],
)

# Proportion analysis (reorganized submodule)
_export(
    "proportion",
    [
        "analyze_celltype_proportion",
        "analyze_all_methods",
        "celltype_proportion_analysis",
        "ProportionMethod",
        "recommend_method",
        "compare_methods",
        "ProportionConfig",
        "MethodSelectionConfig",
        "pb_analysis",
        "compute_celltype_proportion",
        "run_statistical_test",
        "export_analysis_data",
        "run_sccoda",
        "plot_cell_counts",
        "plot_proportion_bar",
        "plot_grouped_celltype_counts",
        "plot_grouped_proportion_bar",
        "plot_celltype_alluvial",
        "plot_box_summary",
        "plot_proportion_heatmap",
        "plot_celltype_correlation",
        "plot_effect_size_volcano",
        "plot_proportion_timeseries",
        "plot_batch_effect",
    ],
    optional=True,
)
