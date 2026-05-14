"""Compatibility facade for the annotation subpackage."""

from .cluster import annotate_clusters
from .mapping import apply_annotation_mapping, remap_labels
from .reference import run_celltypist, transfer_labels
from .review import (
    build_annotation_review_table,
    filter_marker_table_for_annotation,
    flag_suspect_clusters,
    summarize_annotation_evidence,
)
from .scoring import score_cell_types
from .workflow import evaluate_annotation, run_annotation, run_lineage_state_annotation

__all__ = [
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
]
