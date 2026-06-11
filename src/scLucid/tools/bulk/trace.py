"""Review-summary and audit helpers for bulk analysis."""

from __future__ import annotations

from typing import Any, Dict

from ...utils import sanitize_for_hdf5


def build_bulk_review_summary(
    diagnostics: Dict[str, Any],
    normalization: Dict[str, Any],
    de: Dict[str, Any],
) -> Dict[str, Any]:
    """Build a review-facing summary for a bulk analysis run."""
    summary = {
        "schema_version": "bulk_review_summary_v1",
        "diagnostics": sanitize_for_hdf5(diagnostics),
        "normalization": sanitize_for_hdf5(normalization),
        "de": sanitize_for_hdf5(de),
        "review_action_items": [],
    }
    if not diagnostics.get("passed"):
        summary["review_action_items"].append(
            {
                "priority": "high",
                "action": "Review bulk data quality warnings before interpreting results.",
                "evidence_key": "bulk.diagnostics.warnings",
            }
        )
    if not diagnostics.get("replicate_requirement_met"):
        summary["review_action_items"].append(
            {
                "priority": "high",
                "action": "Results are descriptive only; collect biological replicates for formal inference.",
                "evidence_key": "bulk.diagnostics.replicate_requirement_met",
            }
        )
    return summary
