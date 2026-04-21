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
