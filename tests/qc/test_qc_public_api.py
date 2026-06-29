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
def test_qc_compatibility_aliases_remain_importable_but_hidden_from_all():
    hidden = [
        "run_qc_decision_workflow",
        "run_advanced_qc",
    ]
    for symbol in hidden:
        assert hasattr(qc, symbol), f"compatibility symbol missing: {symbol}"
        assert symbol not in qc.__all__


@pytest.mark.unit
def test_qc_policy_entrypoints_are_public():
    for symbol in ["run_qc", "recommend_qc_policy", "apply_qc_policy"]:
        assert hasattr(qc, symbol)
        assert symbol in qc.__all__


@pytest.mark.unit
def test_generate_qc_report_canonical_source_is_reporting():
    from scLucid.qc import reporting
    from scLucid.qc.filtering import generate_qc_report as filtering_generate_qc_report

    assert qc.generate_qc_report is reporting.generate_qc_report
    assert filtering_generate_qc_report is not reporting.generate_qc_report
