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
