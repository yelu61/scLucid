"""Unit tests for scLucid.preprocess.scale."""

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from scLucid.preprocess.config import ScalingConfig
from scLucid.preprocess.scale import scale_data, regress_out, _robust_scale, _minmax_scale


@pytest.mark.unit
class TestScaleDataZScore:
    """Tests for z-score scaling method."""

    def test_zscore_creates_scaled_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers

    def test_zscore_approximate_zero_mean_unit_var(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        scaled = result.layers["scaled"]
        mean = np.mean(scaled, axis=0)
        std = np.std(scaled, axis=0)
        np.testing.assert_allclose(mean, 0, atol=1e-6)
        np.testing.assert_allclose(std, 1, atol=0.5)

    def test_zscore_max_value_clipping(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="zscore", max_value=5, regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        scaled = result.layers["scaled"]
        assert np.max(scaled) <= 5
        assert np.min(scaled) >= -5

    def test_zscore_updates_X(self, minimal_adata):
        adata = minimal_adata.copy()
        original_x = adata.X.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        # X should be updated to scaled values
        with np.testing.assert_raises(AssertionError):
            np.testing.assert_array_equal(result.X, original_x)


@pytest.mark.unit
class TestScaleDataRobust:
    """Tests for robust scaling method."""

    def test_robust_creates_scaled_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="robust", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers

    def test_robust_dense(self):
        x = np.random.default_rng(0).normal(size=(50, 20))
        adata = AnnData(X=x.copy())
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="robust", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers

    def test_robust_sparse(self):
        x = scipy.sparse.random(50, 20, density=0.4, format="csr", random_state=0)
        x.data = np.abs(x.data) * 10
        adata = AnnData(X=x.copy())
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="robust", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers
        assert scipy.sparse.issparse(result.layers["scaled"])

    def test_robust_max_value_clipping(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="robust", max_value=3, regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        scaled = result.layers["scaled"]
        assert np.max(scaled) <= 3
        assert np.min(scaled) >= -3


@pytest.mark.unit
class TestScaleDataMinMax:
    """Tests for minmax scaling method."""

    def test_minmax_creates_scaled_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="minmax", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers

    def test_minmax_range_zero_to_one(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="minmax", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        scaled = result.layers["scaled"]
        assert np.min(scaled) >= 0
        assert np.max(scaled) <= 1

    def test_minmax_sparse(self):
        x = scipy.sparse.random(50, 20, density=0.4, format="csr", random_state=0)
        x.data = np.abs(x.data) * 10
        adata = AnnData(X=x.copy())
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="minmax", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        assert "scaled" in result.layers
        scaled = result.layers["scaled"]
        assert np.min(scaled.data) >= 0
        assert np.max(scaled.data) <= 1


@pytest.mark.unit
class TestRegressOut:
    """Tests for regress_out function."""

    def test_regress_out_creates_regressed_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        # Add required input layer and covariate
        adata.layers["normalized"] = adata.X.copy()
        adata.obs["total_counts"] = adata.X.sum(axis=1)
        result = regress_out(
            adata,
            config=ScalingConfig(vars_to_regress=["total_counts"], plot=False, report=False, verbose=False),
        )
        assert "regressed" in result.layers

    def test_regress_out_skips_when_no_vars(self, minimal_adata):
        adata = minimal_adata.copy()
        result = regress_out(
            adata,
            config=ScalingConfig(vars_to_regress=None, plot=False, report=False, verbose=False),
        )
        # Should not create regressed layer when no vars specified
        assert "regressed" not in result.layers

    def test_regress_out_warns_on_missing_vars(self, minimal_adata, caplog):
        import logging
        adata = minimal_adata.copy()
        with caplog.at_level(logging.WARNING):
            regress_out(
                adata,
                config=ScalingConfig(vars_to_regress=["missing_var"], plot=False, report=False, verbose=False),
            )
        assert "not found" in caplog.text.lower() or "missing" in caplog.text.lower()


@pytest.mark.unit
class TestScaleDataInlineRegression:
    """Tests for inline regression within scale_data."""

    def test_inline_regression_runs_when_vars_present(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["total_counts"] = adata.X.sum(axis=1)
        adata.obs["pct_counts_mt"] = np.random.default_rng(0).random(adata.n_obs)
        # Need normalized layer for regression input
        adata.layers["normalized"] = adata.X.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                regress_in_scale=True,
                vars_to_regress=["total_counts", "pct_counts_mt"],
                input_layer_for_regress="normalized",
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers
        assert "regress_inline" in result.uns["sclucid"]["preprocess"]

    def test_inline_regression_skips_when_vars_missing(self, minimal_adata, caplog):
        import logging
        adata = minimal_adata.copy()
        adata.layers["normalized"] = adata.X.copy()
        with caplog.at_level(logging.INFO):
            result = scale_data(
                adata,
                config=ScalingConfig(
                    regress_in_scale=True,
                    vars_to_regress=["nonexistent"],
                    input_layer_for_regress="normalized",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )
        assert "scaled" in result.layers

    def test_inline_regression_raises_for_missing_input_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["total_counts"] = adata.X.sum(axis=1)
        with pytest.raises(ValueError, match="not found"):
            scale_data(
                adata,
                config=ScalingConfig(
                    regress_in_scale=True,
                    vars_to_regress=["total_counts"],
                    input_layer_for_regress="missing_layer",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )


@pytest.mark.unit
class TestScaleDataMetadata:
    """Tests for metadata storage."""

    def test_scaling_metadata_stored(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False),
        )
        meta = result.uns["sclucid"]["preprocess"]["scaling"]
        assert meta["params"]["scale_method"] == "zscore"
        assert meta["output_layer"] == "scaled"

    def test_config_not_mutated(self, minimal_adata):
        adata = minimal_adata.copy()
        config = ScalingConfig(scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False)
        original_dict = config.to_dict()
        scale_data(adata, config=config)
        assert config.to_dict() == original_dict


@pytest.mark.unit
class TestScaleDataHelpers:
    """Tests for internal helper functions."""

    def test_robust_scale_dense(self):
        x = np.random.default_rng(0).normal(size=(30, 10))
        scaled = _robust_scale(x, max_value=None)
        # Median should be approximately 0 after robust scaling
        medians = np.median(scaled, axis=0)
        np.testing.assert_allclose(medians, 0, atol=1e-6)

    def test_minmax_scale_dense(self):
        x = np.random.default_rng(0).normal(size=(30, 10))
        scaled = _minmax_scale(x)
        assert np.min(scaled) >= 0
        assert np.max(scaled) <= 1

    def test_robust_scale_with_max_value(self):
        x = np.random.default_rng(0).normal(size=(30, 10))
        scaled = _robust_scale(x, max_value=2)
        assert np.max(scaled) <= 2
        assert np.min(scaled) >= -2
