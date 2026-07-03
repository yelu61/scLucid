"""Public API checks for scLucid.qc."""

import pytest

import scLucid.qc as qc


@pytest.mark.unit
def test_qc_exports_resolve():
    for symbol in qc.__all__:
        assert hasattr(qc, symbol), f"scLucid.qc missing exported symbol: {symbol}"


@pytest.mark.unit
def test_qc_thresholds_success_and_validation_error():
    thresholds = qc.QCThresholds(min_genes=100, max_genes=1000, pc_mt=20.0)
    assert thresholds.min_genes == 100
    assert thresholds.max_genes == 1000

    with pytest.raises(ValueError):
        qc.QCThresholds(min_genes=500, max_genes=100)


@pytest.mark.unit
def test_qc_no_longer_exports_removed_optional_modules():
    removed_symbols = [
        "QCStrategyDecisionTree",
        "recommend_qc_strategy",
        "InteractiveQCExplorer",
        "InteractiveQCPlotter",
        "create_interactive_dashboard",
        "interactive_filter_preview",
    ]
    for symbol in removed_symbols:
        assert not hasattr(qc, symbol), f"scLucid.qc should not export removed symbol: {symbol}"


@pytest.mark.unit
def test_qc_does_not_export_legacy_threshold_entrypoints():
    hidden = [
        "suggest_qc_thresholds",
        "resolve_qc_thresholds",
        "mark_low_quality_cells_adaptive",
        "run_qc_decision_workflow",
    ]
    for symbol in hidden:
        assert not hasattr(qc, symbol), f"legacy threshold API should not be public: {symbol}"

    assert hasattr(qc, "run_advanced_qc")
    assert "run_advanced_qc" not in qc.__all__


@pytest.mark.unit
def test_qc_threshold_chain_entrypoints_are_public():
    for symbol in [
        "recommend_qc_thresholds",
        "decide_qc_thresholds",
        "apply_qc_threshold_decision",
        "run_qc_threshold_decision",
        "filter_cells",
        "build_qc_decisions",
        "evaluate_qc_benchmark",
    ]:
        assert hasattr(qc, symbol)
        assert symbol in qc.__all__


@pytest.mark.unit
def test_filtering_module_only_exposes_final_subsetting_api():
    import scLucid.qc.filtering as filtering

    assert filtering.__all__ == ["filter_cells"]
    assert hasattr(filtering, "filter_cells")
    assert not hasattr(filtering, "mark_low_quality_cell")
    assert not hasattr(filtering, "mark_low_quality_cells_adaptive")
    assert not hasattr(filtering, "AdaptiveThresholdCalculator")


@pytest.mark.unit
def test_qc_policy_entrypoints_are_public():
    for symbol in ["run_qc", "recommend_qc_policy", "apply_qc_policy"]:
        assert hasattr(qc, symbol)
        assert symbol in qc.__all__


@pytest.mark.unit
def test_qc_artifact_contract_and_cleanup_are_public_and_preserve_required_columns():
    import numpy as np
    from anndata import AnnData

    adata = AnnData(X=np.ones((3, 2), dtype=float))
    adata.obs["total_counts"] = [10, 20, 30]
    adata.obs["n_genes_by_counts"] = [2, 2, 2]
    adata.obs["predicted_doublet"] = [False, True, False]
    adata.obs["scrublet_score"] = [0.1, 0.9, 0.2]
    adata.obs["temporary_note"] = ["a", "b", "c"]

    contract = qc.record_qc_artifact_contract(adata)
    assert contract["schema_version"] == qc.QC_ARTIFACT_CONTRACT_SCHEMA_VERSION
    assert "predicted_doublet" in contract["filtering_required_obs_columns"]

    summary = qc.cleanup_qc_intermediates(adata, mode="review")
    assert "scrublet_score" in summary["dropped_obs_columns"]
    assert "predicted_doublet" in adata.obs
    assert "temporary_note" in adata.obs
    assert "scrublet_score" not in adata.obs

    for symbol in [
        "record_threshold_recommendation",
        "record_threshold_decision",
        "record_mark_evidence",
        "record_qc_decision_artifact",
        "record_filter_result",
        "record_benchmark_review",
    ]:
        assert hasattr(qc, symbol)
        assert symbol in qc.__all__


@pytest.mark.unit
def test_generate_qc_report_canonical_source_is_reporting():
    from scLucid.qc import reporting

    assert qc.generate_qc_report is reporting.generate_qc_report
