"""Tests for bulk trace/review summary helpers."""

from scLucid.tools.bulk import build_bulk_review_summary


def test_build_bulk_review_summary_passed():
    diagnostics = {"passed": True, "replicate_requirement_met": True, "warnings": []}
    normalization = {"method": "CPM", "size_factors": [1.0, 1.1]}
    de = {"method": "welch", "n_genes_tested": 100}
    summary = build_bulk_review_summary(diagnostics, normalization, de)
    assert summary["schema_version"] == "bulk_review_summary_v1"
    assert summary["diagnostics"]["passed"]
    assert not summary["review_action_items"]


def test_build_bulk_review_summary_warnings():
    diagnostics = {"passed": False, "replicate_requirement_met": False, "warnings": ["Low counts"]}
    normalization = {"method": "CPM"}
    de = {"method": "descriptive"}
    summary = build_bulk_review_summary(diagnostics, normalization, de)
    assert len(summary["review_action_items"]) == 2
    priorities = {item["priority"] for item in summary["review_action_items"]}
    assert "high" in priorities
