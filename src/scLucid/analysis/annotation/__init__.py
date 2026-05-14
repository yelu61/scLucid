"""
Annotation submodule.

The public annotation API is re-exported here while the implementation is split
out of the former monolithic ``annotation.py`` module.
"""

from .core import (
    annotate_clusters,
    apply_annotation_mapping,
    build_annotation_review_table,
    evaluate_annotation,
    filter_marker_table_for_annotation,
    flag_suspect_clusters,
    remap_labels,
    run_annotation,
    run_celltypist,
    run_lineage_state_annotation,
    score_cell_types,
    summarize_annotation_evidence,
    transfer_labels,
)
from .evidence import (
    ANALYSIS_REVIEW_SUMMARY_SCHEMA,
    ANNOTATION_REVIEW_SCHEMA,
    apply_final_annotation,
    build_llm_annotation_bundle,
    merge_annotation_evidence,
    run_marker_annotation_evidence,
    standardize_cluster_marker_table,
)

__all__ = [
    "ANALYSIS_REVIEW_SUMMARY_SCHEMA",
    "ANNOTATION_REVIEW_SCHEMA",
    "score_cell_types",
    "annotate_clusters",
    "run_celltypist",
    "transfer_labels",
    "evaluate_annotation",
    "summarize_annotation_evidence",
    "standardize_cluster_marker_table",
    "run_marker_annotation_evidence",
    "build_llm_annotation_bundle",
    "merge_annotation_evidence",
    "apply_final_annotation",
    "run_lineage_state_annotation",
    "filter_marker_table_for_annotation",
    "flag_suspect_clusters",
    "build_annotation_review_table",
    "apply_annotation_mapping",
    "remap_labels",
    "run_annotation",
]
