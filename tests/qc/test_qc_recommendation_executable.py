"""Tests that QC recommendations are actually applied to config and affect filtering."""

import pytest

from scLucid.qc.workflow import run_standard_qc
from scLucid.qc.config import QCWorkflowConfig, MetricsReportingConfig, MarkingConfig, DoubletConfig
from tests.fixtures.data_loader import load_test_data


@pytest.fixture
def adata_pbmc():
    return load_test_data("pbmc3k")


def _make_config(use_recommendations: bool) -> QCWorkflowConfig:
    return QCWorkflowConfig(
        save_dir=None,
        use_recommendations=use_recommendations,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )


def test_recommendation_changes_applied_config(adata_pbmc):
    config_with = _make_config(use_recommendations=True)
    config_without = _make_config(use_recommendations=False)
    adata_with = run_standard_qc(adata_pbmc.copy(), config=config_with)
    adata_without = run_standard_qc(adata_pbmc.copy(), config=config_without)

    rec = adata_with.uns["sclucid"]["qc"]["recommendation"]["data"]
    assert rec is not None, "Expected a recommendation when use_recommendations=True"

    # The applied config should reflect the recommendation
    applied = adata_with.uns["sclucid"]["qc"]["applied_config"]["data"]
    rec_min_genes = rec.get("min_genes", {}).get("threshold")
    if rec_min_genes is not None:
        assert applied["marking_config"]["thresholds"]["min_genes"] == rec_min_genes

    # Without recommendations, default min_genes should stay 200
    applied_no_rec = adata_without.uns["sclucid"]["qc"]["applied_config"]["data"]
    assert applied_no_rec["marking_config"]["thresholds"]["min_genes"] == 200


def test_recommendation_affects_cell_counts(adata_pbmc):
    config_with = _make_config(use_recommendations=True)
    config_without = _make_config(use_recommendations=False)
    adata_with = run_standard_qc(adata_pbmc.copy(), config=config_with)
    adata_without = run_standard_qc(adata_pbmc.copy(), config=config_without)
    assert adata_with.n_obs != adata_without.n_obs, (
        "Recommended thresholds should produce different filtering results than defaults"
    )


def test_caller_config_not_mutated(adata_pbmc):
    from scLucid.qc.config import QCThresholds
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="pooled",
        marking_config=MarkingConfig(
            show_plots=False, plot_outliers=False, thresholds=QCThresholds(min_genes=200)
        ),
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    original_min_genes = config.marking_config.thresholds.min_genes
    original_pc_mt = config.marking_config.thresholds.pc_mt
    _ = run_standard_qc(adata_pbmc.copy(), config=config)
    assert config.marking_config.thresholds.min_genes == original_min_genes
    assert config.marking_config.thresholds.pc_mt == original_pc_mt


def test_explicit_user_thresholds_survive_recommendations(adata_pbmc):
    from scLucid.qc.config import QCThresholds
    user_min_genes = 999
    user_pc_mt = 5.0
    config = QCWorkflowConfig(
        save_dir=None,
        use_recommendations=True,
        threshold_mode="pooled",
        marking_config=MarkingConfig(
            show_plots=False,
            plot_outliers=False,
            thresholds=QCThresholds(min_genes=user_min_genes, pc_mt=user_pc_mt),
        ),
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata_pbmc.copy(), config=config)
    applied = adata_f.uns["sclucid"]["qc"]["applied_config"]["data"]
    assert applied["marking_config"]["thresholds"]["min_genes"] == user_min_genes
    assert applied["marking_config"]["thresholds"]["pc_mt"] == user_pc_mt
    # user_overrides should show divergence from recommendation
    overrides = adata_f.uns["sclucid"]["qc"]["user_overrides"]["data"]
    assert "min_genes" in overrides or "max_mt_percent" in overrides


def test_default_thresholds_do_not_appear_as_user_overrides(adata_pbmc):
    config = _make_config(use_recommendations=True)
    adata_f = run_standard_qc(adata_pbmc.copy(), config=config)
    overrides = adata_f.uns["sclucid"]["qc"]["user_overrides"]["data"]
    assert "min_genes" not in overrides
    assert "max_mt_percent" not in overrides
    assert "n_counts" not in overrides
