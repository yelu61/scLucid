"""Regression test for legacy analysis.bulk imports."""


import pytest


@pytest.mark.parametrize(
    "name",
    [
        "deconvolve_bulk",
        "run_bulk_de",
        "diagnose_bulk_data_quality",
        "normalize_bulk_counts",
        "estimate_size_factors_median_ratio",
        "run_bulk_abundance_test",
        "correlate_abundance_with_clinical",
        "deduplicate_var_names",
        "filter_bulk_genes",
        "build_bulk_review_summary",
        "BulkDiagnosticsConfig",
        "BulkDEConfig",
    ],
)
def test_analysis_bulk_legacy_imports(name):
    """scLucid.analysis.bulk must remain importable after the move to tools."""
    import importlib

    mod = importlib.import_module("scLucid.analysis.bulk")
    assert hasattr(mod, name)


@pytest.mark.parametrize(
    "name",
    [
        "deconvolve_bulk",
        "run_bulk_de",
        "run_bulk_abundance_test",
        "BulkDiagnosticsConfig",
    ],
)
def test_tools_bulk_imports(name):
    """scLucid.tools.bulk must expose the canonical API."""
    import importlib

    mod = importlib.import_module("scLucid.tools.bulk")
    assert hasattr(mod, name)
