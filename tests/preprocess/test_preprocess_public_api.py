"""Public API checks for scLucid.preprocess."""

import pytest

import scLucid.preprocess as pp


@pytest.mark.unit
def test_preprocess_exports_resolve():
    for symbol in pp.__all__:
        assert hasattr(pp, symbol), f"scLucid.preprocess missing exported symbol: {symbol}"


@pytest.mark.unit
def test_normalization_config_success_and_reserved_layer_validation():
    if not hasattr(pp, "NormalizationConfig"):
        pytest.skip("NormalizationConfig not available in current environment")

    cfg = pp.NormalizationConfig(target_sum=1e4, output_layer="normalized")
    assert cfg.target_sum == 1e4

    with pytest.raises(ValueError):
        pp.NormalizationConfig(output_layer="X")
