"""Tests for multi-sample threshold policy (hierarchical / pooled / independent)."""

import pytest

from scLucid.qc.config import DoubletConfig, MarkingConfig, MetricsReportingConfig, QCWorkflowConfig
from scLucid.qc.workflow import run_standard_qc
from tests.fixtures.data_loader import load_test_data


@pytest.fixture
def adata_pbmc():
    return load_test_data("pbmc3k")


def _make_config(mode: str) -> QCWorkflowConfig:
    return QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode=mode,
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )


def test_hierarchical_produces_per_sample_thresholds(adata_pbmc):
    config = _make_config("hierarchical")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    sample_thresholds = adata_f.uns["sclucid"]["qc"]["sample_thresholds"]["data"]
    assert sample_thresholds, "Expected per-sample thresholds in hierarchical mode"
    for sample_id, thresholds in sample_thresholds.items():
        assert "n_genes_by_counts" in thresholds or "pct_counts_mt" in thresholds
        assert "total_counts" in thresholds
        if "n_genes_by_counts" in thresholds:
            assert thresholds["n_genes_by_counts"]["method"] == "hierarchical"
        assert thresholds["total_counts"]["method"] == "hierarchical"


def test_pooled_has_empty_sample_thresholds(adata_pbmc):
    config = _make_config("pooled")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    sample_thresholds = adata_f.uns["sclucid"]["qc"]["sample_thresholds"]["data"]
    assert sample_thresholds == {}, "Expected empty sample_thresholds in pooled mode"


def test_hierarchical_and_pooled_are_stable_on_homogeneous_data(adata_pbmc):
    config_hier = _make_config("hierarchical")
    config_pooled = _make_config("pooled")
    adata_hier = run_standard_qc(adata_pbmc.copy(), config=config_hier)
    adata_pooled = run_standard_qc(adata_pbmc.copy(), config=config_pooled)
    hier_cells = adata_hier.n_obs
    pooled_cells = adata_pooled.n_obs
    # For a near-homogeneous baseline dataset, both modes should retain a reasonable
    # fraction of cells. The exact count may differ because hierarchical thresholds
    # are sample-specific, but neither mode should produce pathological removal.
    assert (
        hier_cells > adata_pbmc.n_obs * 0.5
    ), f"Hierarchical mode removed too many cells: {hier_cells}/{adata_pbmc.n_obs}"
    assert (
        pooled_cells > adata_pbmc.n_obs * 0.5
    ), f"Pooled mode removed too many cells: {pooled_cells}/{adata_pbmc.n_obs}"


def test_hierarchical_thresholds_are_clipped_to_valid_ranges(adata_pbmc):
    config = _make_config("hierarchical")
    adata_f = run_standard_qc(adata_pbmc, config=config)
    sample_thresholds = adata_f.uns["sclucid"]["qc"]["sample_thresholds"]["data"]
    for sample_id, thresholds in sample_thresholds.items():
        for metric, th in thresholds.items():
            lower = th.get("lower")
            upper = th.get("upper")
            assert lower is not None and upper is not None
            if metric in ("n_genes_by_counts", "total_counts"):
                assert lower >= 0, f"{metric} lower bound for {sample_id} is negative: {lower}"
            elif metric.startswith("pct_counts_") or metric in ("pct_counts_mt", "pct_counts_hb"):
                assert lower >= 0, f"{metric} lower bound for {sample_id} is negative: {lower}"
                assert upper <= 100, f"{metric} upper bound for {sample_id} exceeds 100: {upper}"
