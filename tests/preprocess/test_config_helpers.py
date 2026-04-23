"""Unit tests for preprocess config helper utilities."""

import logging

import pytest

from scLucid.preprocess.config import NormalizationConfig, ScalingConfig, apply_config_overrides


@pytest.mark.unit
class TestApplyConfigOverrides:
    """Tests for apply_config_overrides helper."""

    def test_returns_deep_copy(self):
        cfg = NormalizationConfig(method="standard", target_sum=1e4)
        result = apply_config_overrides(cfg, target_sum=1e5)
        assert result.target_sum == 1e5
        assert cfg.target_sum == 1e4  # Original unchanged

    def test_unknown_kwargs_logged(self, caplog):
        cfg = NormalizationConfig()
        with caplog.at_level(logging.WARNING):
            apply_config_overrides(cfg, unknown_param=42)
        assert "Ignoring unknown parameter: 'unknown_param'" in caplog.text

    def test_ignored_keys_skipped_silently(self, caplog):
        cfg = NormalizationConfig()
        with caplog.at_level(logging.WARNING):
            apply_config_overrides(cfg, ignored_keys={"force"}, force=True)
        assert "force" not in caplog.text

    def test_valid_kwargs_applied(self):
        cfg = ScalingConfig(scale_method="zscore", max_value=10)
        result = apply_config_overrides(cfg, scale_method="robust", max_value=5)
        assert result.scale_method == "robust"
        assert result.max_value == 5
        assert cfg.scale_method == "zscore"  # Original unchanged

    def test_multiple_overrides(self):
        cfg = NormalizationConfig()
        result = apply_config_overrides(cfg, method="clr", target_sum=1, output_layer="custom")
        assert result.method == "clr"
        assert result.target_sum == 1
        assert result.output_layer == "custom"
