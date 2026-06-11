"""Review-summary helpers for spatial analysis."""

from __future__ import annotations

from typing import Any, Dict

from ...utils import sanitize_for_hdf5


def build_spatial_review_summary(
    diagnostics: Dict[str, Any],
    neighbors: Dict[str, Any],
    svg: Dict[str, Any],
) -> Dict[str, Any]:
    """Build a review-facing summary for a spatial analysis run."""
    summary = {
        "schema_version": "spatial_review_summary_v1",
        "diagnostics": sanitize_for_hdf5(diagnostics),
        "neighbors": sanitize_for_hdf5(neighbors),
        "svg": sanitize_for_hdf5(svg),
        "review_action_items": [],
    }
    if not diagnostics.get("passed"):
        summary["review_action_items"].append(
            {
                "priority": "high",
                "action": "Review spatial data quality warnings before interpreting results.",
                "evidence_key": "spatial.diagnostics.warnings",
            }
        )
    if diagnostics.get("n_duplicate_coords", 0) > 0:
        summary["review_action_items"].append(
            {
                "priority": "medium",
                "action": "Duplicate spatial coordinates detected; verify spot/cell assignment or image alignment.",
                "evidence_key": "spatial.diagnostics.n_duplicate_coords",
            }
        )
    return summary
