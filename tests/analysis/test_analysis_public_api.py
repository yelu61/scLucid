"""Public API checks for scLucid.analysis."""

import pytest

import scLucid.analysis as analysis


@pytest.mark.unit
def test_analysis_exports_resolve():
    for symbol in analysis.__all__:
        assert hasattr(analysis, symbol), f"scLucid.analysis missing exported symbol: {symbol}"


@pytest.mark.unit
def test_clustering_config_success_and_invalid_method():
    cfg = analysis.ClusteringConfig(method="leiden", resolution=1.0)
    assert cfg.method == "leiden"

    with pytest.raises(ValueError):
        analysis.ClusteringConfig(method="not_a_method")
