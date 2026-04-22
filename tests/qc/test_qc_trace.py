"""Tests for unified QC trace contract under adata.uns['sclucid']['qc']."""

import pytest
import scanpy as sc

from scLucid.qc.workflow import run_standard_qc
from scLucid.qc.config import QCWorkflowConfig, MetricsReportingConfig, MarkingConfig, DoubletConfig
from tests.fixtures.data_loader import load_test_data


REQUIRED_QC_KEYS = {
    "context",
    "recommendation",
    "original_config",
    "applied_config",
    "user_overrides",
    "sample_thresholds",
    "filtering_summary",
    "warnings",
    "review_summary",
}


@pytest.fixture
def adata_pbmc():
    return load_test_data("pbmc3k")


def test_qc_trace_has_all_required_keys(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="hierarchical",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config)
    qc_trace = adata_f.uns.get("sclucid", {}).get("qc", {})
    missing = REQUIRED_QC_KEYS - set(qc_trace.keys())
    assert not missing, f"Missing QC trace keys: {missing}"


def test_qc_trace_context_schema(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="hierarchical",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config)
    context = adata_f.uns["sclucid"]["qc"]["context"]["data"]
    assert context["threshold_mode"] == "hierarchical"
    assert context["n_samples"] == 4
    assert context["use_recommendations"] is True


def test_qc_trace_warnings_is_list(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="hierarchical",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config)
    warnings = adata_f.uns["sclucid"]["qc"]["warnings"]["data"]
    assert isinstance(warnings, list)


def test_run_standard_qc_does_not_mutate_input_config_tissue_type(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        tissue_type=None,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes"]},
    )

    _ = run_standard_qc(adata_pbmc, config=config, tissue_type="lung_tumor")

    assert config.tissue_type is None


# ---------------------------------------------------------------------------
# Review-summary contract tests
# ---------------------------------------------------------------------------

REVIEW_SUMMARY_SECTIONS = {
    "recommendation_summary",
    "applied_threshold_summary",
    "user_override_summary",
    "sample_threshold_summary",
    "tumor_aware_summary",
    "filtering_summary",
    "warnings",
}


def test_qc_review_summary_has_expected_sections(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config)
    review = adata_f.uns["sclucid"]["qc"]["review_summary"]["data"]
    missing = REVIEW_SUMMARY_SECTIONS - set(review.keys())
    assert not missing, f"Missing review_summary sections: {missing}"


def test_qc_review_summary_filtering_consistency(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config)
    review = adata_f.uns["sclucid"]["qc"]["review_summary"]["data"]
    stored_filtering = adata_f.uns["sclucid"]["qc"]["filtering_summary"]["data"]
    assert review["filtering_summary"]["initial_cells"] == stored_filtering.get("initial_cells")
    assert review["filtering_summary"]["final_cells"] == stored_filtering.get("final_cells")


def test_qc_review_summary_tumor_aware_flag(adata_pbmc):
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata_pbmc, config=config, tissue_type="lung_tumor")
    review = adata_f.uns["sclucid"]["qc"]["review_summary"]["data"]
    assert review["tumor_aware_summary"]["enabled"] is True
    assert any("tumor" in note.lower() for note in review["tumor_aware_summary"]["notes"])


def test_qc_review_summary_exported_to_disk(adata_pbmc, tmp_path):
    save_dir = str(tmp_path / "qc_output")
    config = QCWorkflowConfig(
        save_dir=save_dir,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    _ = run_standard_qc(adata_pbmc, config=config)
    json_path = tmp_path / "qc_output" / "qc_review_summary.json"
    md_path = tmp_path / "qc_output" / "qc_review_summary.md"
    assert json_path.exists(), "JSON sidecar should be written when save_dir is set"
    assert md_path.exists(), "Markdown sidecar should be written when save_dir is set"
    import json
    loaded = json.loads(json_path.read_text())
    assert "filtering_summary" in loaded
