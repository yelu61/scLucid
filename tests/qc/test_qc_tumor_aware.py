"""Tests for tumor-aware QC behavior: flagging and relaxed MT filtering."""

import pytest

from scLucid.qc.config import DoubletConfig, MarkingConfig, MetricsReportingConfig, QCWorkflowConfig
from scLucid.qc.workflow import run_standard_qc
from tests.fixtures.data_loader import load_test_data


@pytest.fixture
def adata_pbmc():
    return load_test_data("pbmc3k")


def _make_config(tissue_type: str) -> QCWorkflowConfig:
    return QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="pooled",
        tissue_type=tissue_type,
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )


def test_tumor_aware_removes_outlier_mt_from_filtering(adata_pbmc):
    config = _make_config("tumor")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    filtering_summary = adata_f.uns["sclucid"]["qc"]["filtering_summary"]["data"]
    criteria_used = filtering_summary["criteria_used"]
    assert (
        "outlier_mt" not in criteria_used
    ), f"outlier_mt should be removed in tumor-aware mode, got {criteria_used}"


def test_tumor_aware_stores_flags_and_warning(adata_pbmc):
    config = _make_config("lung_tumor")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    warnings = adata_f.uns["sclucid"]["qc"]["warnings"]["data"]
    assert any("outlier_mt excluded" in w for w in warnings)

    tumor_flags = adata_f.uns["sclucid"]["qc"]["tumor_aware_flags"]["data"]
    assert tumor_flags["tumor_aware_enabled"] is True


def test_non_tumor_keeps_outlier_mt_filtering(adata_pbmc):
    config = _make_config("normal")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    filtering_summary = adata_f.uns["sclucid"]["qc"]["filtering_summary"]["data"]
    criteria_used = filtering_summary["criteria_used"]
    assert (
        "outlier_mt" in criteria_used
    ), f"outlier_mt should be retained for non-tumor tissue, got {criteria_used}"
