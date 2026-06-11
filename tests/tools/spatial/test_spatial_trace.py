"""Tests for spatial trace/review summary helpers."""

from scLucid.tools.spatial import build_spatial_review_summary


def test_build_spatial_review_summary_passed():
    diagnostics = {"passed": True, "n_duplicate_coords": 0, "warnings": []}
    neighbors = {"method": "knn", "n_neigh": 6}
    svg = {"n_genes_tested": 100, "n_significant": 10}
    summary = build_spatial_review_summary(diagnostics, neighbors, svg)
    assert summary["schema_version"] == "spatial_review_summary_v1"
    assert not summary["review_action_items"]


def test_build_spatial_review_summary_warnings():
    diagnostics = {"passed": False, "n_duplicate_coords": 5, "warnings": ["Missing coords"]}
    neighbors = {"method": "knn"}
    svg = {"n_genes_tested": 100}
    summary = build_spatial_review_summary(diagnostics, neighbors, svg)
    assert len(summary["review_action_items"]) >= 2
    priorities = {item["priority"] for item in summary["review_action_items"]}
    assert "high" in priorities
    assert "medium" in priorities
