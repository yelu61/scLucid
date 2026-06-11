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
        "PseudobulkDEConfig",
        "EnrichmentConfig",
        "ProportionConfig",
        "AnalysisWorkflowConfig",
    ],
)
_export("clustering", ["run_clustering_review", "cluster_cells", "merge_clusters"])
_export(
    "annotation",
    [
        "score_cell_types",
        "annotate_clusters",
        "run_celltypist",
        "transfer_labels",
        "evaluate_annotation",
        "summarize_annotation_evidence",
        "standardize_cluster_marker_table",
        "build_hierarchical_annotation_plan",
        "build_subset_annotation_reconciliation",
        "apply_subset_annotation_reconciliation",
        "run_marker_annotation_evidence",
        "run_program_annotation_evidence",
        "run_subset_annotation_refinement",
        "run_annotation_evidence",
        "build_llm_annotation_bundle",
        "merge_annotation_evidence",
        "build_annotation_consensus",
        "apply_final_annotation",
        "run_lineage_state_annotation",
        "filter_marker_table_for_annotation",
        "flag_suspect_clusters",
        "build_annotation_review_table",
        "apply_annotation_mapping",
        "remap_labels",
        "run_annotation",
    ],
)
_export(
    "trace",
    [
        "ANALYSIS_REQUIRED_REVIEW_SECTIONS",
        "build_posthoc_qc_review_summary",
        "enrich_analysis_review_summary",
        "get_analysis_module_contract",
        "summarize_analysis_review_summary",
        "validate_analysis_module_completeness",
        "validate_analysis_review_summary",
    ],
)
# Backward compatibility: re-export tumor malignancy interpretation with a
# deprecation warning. Users should call this from ``scLucid.tumor`` instead.
try:
    from ..tumor.malignancy import (
        run_malignancy_interpretation as _run_malignancy_interpretation,
    )

    def run_malignancy_interpretation(*args, **kwargs):  # type: ignore[misc]
        warnings.warn(
            "scLucid.analysis.run_malignancy_interpretation is deprecated and will be "
            "removed in a future release. Use scLucid.tumor.run_malignancy_interpretation "
            "or run_tumor_analysis() instead.",
            FutureWarning,
            stacklevel=2,
        )
        return _run_malignancy_interpretation(*args, **kwargs)

    __all__.append("run_malignancy_interpretation")
except Exception as exc:
    warnings.warn(
        f"Could not import tumor malignancy interpretation bridge: {exc}",
        ImportWarning,
    )

# Prefer the reorganized DE package, but keep a legacy fallback.
_de_names = [
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "run_pseudobulk_de",
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
        "plot_composition",
        "plot_diff_stats",
        "plot_individual_boxplots",
        "plot_proportion_shifts",
        "plot_paired_proportion_shifts",
        "plot_proportion_with_ci",
        "plot_celltype_variability",
        "plot_composition_transform_heatmap",
        "plot_composition_pca",
        "transform_composition",
    ],
    optional=True,
)
