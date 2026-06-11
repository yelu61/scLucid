"""Tests for the scLucid resource trust audit layer."""

import pytest

from scLucid.utils import (
    assert_trusted_resources,
    audit_geneset_resources,
    audit_marker_entry_quality,
    audit_marker_resources,
    audit_resource_manifest,
    build_resource_trust_report,
    classify_literature_resource_utility,
    load_marker_curation_literature_index,
    load_reference_index,
)


def test_resource_trust_report_passes_with_review_queue_warning():
    report = build_resource_trust_report()

    assert report["status"] == "warn"
    assert report["n_errors"] == 0
    assert report["n_warnings"] >= 1
    assert report["known_source_ids"] >= 220
    assert "markers" in report["sections"]
    assert "genesets" in report["sections"]


def test_assert_trusted_resources_allows_non_strict_review_queue():
    report = assert_trusted_resources(strict=False)

    assert report["status"] == "warn"
    assert report["n_errors"] == 0


def test_assert_trusted_resources_strict_requires_full_literature_closure():
    with pytest.raises(AssertionError, match="resource trust audit failed"):
        assert_trusted_resources(strict=True)


def test_active_manifest_resources_exist_and_have_workflows():
    manifest_report = audit_resource_manifest()

    assert manifest_report["issues"] == []
    assert {
        "marker_registry_human",
        "marker_registry_mouse",
        "marker_tissue_human",
        "marker_tumor_human",
        "genesets_cancer_signatures",
        "genesets_cancer_hallmarks",
        "marker_aliases",
    }.issubset(set(manifest_report["active_resources"]))


def test_marker_audit_keeps_views_separated():
    marker_report = audit_marker_resources(known_source_ids=set(load_reference_index()))

    assert marker_report["issues"] == []
    views = marker_report["views"]
    assert views["global_annotation"] > 0
    assert views["state_annotation"] > 0
    assert views["program_scoring"] > 0
    assert views["tumor_interpretation"] > 0
    assert views["doublet_detection"] > 0


def test_geneset_audit_resolves_sources_and_categories():
    geneset_report = audit_geneset_resources(known_source_ids=set(load_reference_index()))

    assert geneset_report["issues"] == []
    assert geneset_report["resources"]["genesets_cancer_signatures"]["n_entries"] >= 20
    assert "ITH_Hallmarks" in geneset_report["resources"]["genesets_cancer_signatures"][
        "categories"
    ]
    assert geneset_report["resources"]["genesets_cancer_hallmarks"]["n_entries"] >= 10


def test_literature_index_remains_full_batch_queue():
    rows = load_marker_curation_literature_index()
    queued = [row for row in rows if row.get("curation_status") == "queued"]
    zotero_matched = [row for row in rows if row.get("zotero_item_key")]

    assert len(rows) == 141
    assert len({row["batch_id"] for row in rows}) == 12
    assert len(queued) <= 4
    assert len(zotero_matched) >= 133
    assert all(row.get("resource_utility") for row in rows)
    assert all(row.get("target_resources") for row in rows)
    assert all(row.get("curation_priority") in {"high", "medium", "low"} for row in rows)


def test_literature_utility_classifier_routes_common_source_types():
    naming = classify_literature_resource_utility(
        {
            "title": "Guidelines for T cell nomenclature",
            "source_type": "review",
            "recommended_resource": "naming_reference",
            "curation_status": "reference_registered",
        }
    )
    tumor = classify_literature_resource_utility(
        {
            "title": "Pan-cancer malignant meta-programs",
            "source_type": "pan_cancer_atlas",
            "recommended_resource": "marker_tumor + geneset",
            "curation_status": "reference_registered",
        }
    )

    assert "nomenclature_reference" in naming["resource_utility"]
    assert "marker_aliases" in naming["target_resources"]
    assert "tumor_context" in tumor["resource_utility"]
    assert "geneset_scoring" in tumor["resource_utility"]
    assert "marker_tumor_human" in tumor["target_resources"]


def test_marker_entry_quality_reports_curation_gaps_without_errors():
    quality = audit_marker_entry_quality()

    assert quality["n_gaps"] > 0
    assert len(quality["gaps"]) == quality["n_gaps"]
    assert len(quality["priority_gaps"]) <= len(quality["gaps"])
    assert "thin_marker_set" in quality["gap_counts"]
    assert "missing_effective_negative_markers" in quality["gap_counts"]
    assert "marker_registry_human" in quality["resources"]
    assert "review_status" in quality["resources"]["marker_tissue_human"]
    assert quality["resources"]["marker_tissue_human"]["review_status"]["needs_review"] >= 1


def test_trust_report_includes_quality_and_literature_triage_sections():
    report = build_resource_trust_report()

    assert "marker_quality" in report["sections"]
    curation = report["sections"]["curation_index"]
    assert curation["resource_utility_counts"]["marker_core"] >= 1
    assert curation["priority_counts"]["high"] >= 1
