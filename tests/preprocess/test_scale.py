"""Unit tests for scLucid.preprocess.scale."""

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from scLucid.preprocess.config import ScalingConfig
from scLucid.preprocess.scale import (
    _minmax_scale,
    _robust_scale,
    _robust_scale_sparse,
    diagnose_cell_cycle_regression,
    regress_out,
    scale_data,
)


@pytest.mark.unit
class TestScaleDataZScore:
    """Tests for z-score scaling method."""

    def test_zscore_creates_scaled_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="zscore",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers

    def test_zscore_approximate_zero_mean_unit_var(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="zscore",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
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
            config=ScalingConfig(
                scale_method="zscore",
                max_value=5,
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        scaled = result.layers["scaled"]
        assert np.max(scaled) <= 5
        assert np.min(scaled) >= -5

    def test_zscore_updates_X(self, minimal_adata):
        adata = minimal_adata.copy()
        original_x = adata.X.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="zscore",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
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
            config=ScalingConfig(
                scale_method="robust",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers

    def test_robust_dense(self):
        x = np.random.default_rng(0).normal(size=(50, 20))
        adata = AnnData(X=x.copy())
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="robust",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers

    def test_robust_sparse(self):
        x = scipy.sparse.random(50, 20, density=0.4, format="csr", random_state=0)
        x.data = np.abs(x.data) * 10
        adata = AnnData(X=x.copy())
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="robust",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers
        # Values must match the dense reference that includes zeros in median/MAD.
        expected = _robust_scale(x.toarray(), max_value=10.0)
        scaled = result.layers["scaled"].toarray() if scipy.sparse.issparse(result.layers["scaled"]) else result.layers["scaled"]
        np.testing.assert_allclose(scaled, expected)

    def test_robust_max_value_clipping(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="robust",
                max_value=3,
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
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
            config=ScalingConfig(
                scale_method="minmax",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "scaled" in result.layers

    def test_minmax_range_zero_to_one(self, minimal_adata):
        adata = minimal_adata.copy()
        result = scale_data(
            adata,
            config=ScalingConfig(
                scale_method="minmax",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
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
            config=ScalingConfig(
                scale_method="minmax",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
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
            config=ScalingConfig(
                vars_to_regress=["total_counts"], plot=False, report=False, verbose=False
            ),
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
                config=ScalingConfig(
                    vars_to_regress=["missing_var"], plot=False, report=False, verbose=False
                ),
            )
        assert "not found" in caplog.text.lower() or "missing" in caplog.text.lower()

    def test_cell_cycle_regression_diagnostic_flags_condition_confounding(self):
        adata = AnnData(X=np.ones((20, 5)))
        adata.var_names = ["MKI67", "TOP2A", "g1", "g2", "g3"]
        adata.obs["condition"] = ["A"] * 10 + ["B"] * 10
        adata.obs["S_score"] = np.r_[np.repeat(0.1, 10), np.repeat(2.0, 10)]
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = ["G1"] * 10 + ["S"] * 10

        result = diagnose_cell_cycle_regression(
            adata,
            condition_key="condition",
            tumor=True,
        )

        assert result["status"] == "review_required"
        assert result["metrics"]["condition_eta2_cc_score"] > 0.15
        assert result["warnings"]
        assert "sclucid" not in adata.uns or "preprocess" not in adata.uns.get("sclucid", {})

        recorded = diagnose_cell_cycle_regression(
            adata,
            condition_key="condition",
            tumor=True,
            record=True,
        )
        assert recorded["status"] == "review_required"
        assert "cell_cycle_regression_diagnostic" in adata.uns["sclucid"]["preprocess"]

    def test_cell_cycle_regression_diagnostic_identifies_batch_candidate(self):
        adata = AnnData(X=np.ones((20, 5)))
        adata.obs["batch"] = ["b1"] * 10 + ["b2"] * 10
        adata.obs["S_score"] = np.r_[np.repeat(0.1, 10), np.repeat(2.0, 10)]
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = ["G1"] * 10 + ["S"] * 10

        result = diagnose_cell_cycle_regression(adata, batch_key="batch")

        assert result["status"] == "technical_regression_candidate"
        assert result["metrics"]["batch_eta2_cc_score"] > 0.2

    def test_cell_cycle_regression_diagnostic_warns_on_group_imbalance(self):
        adata = AnnData(X=np.ones((30, 5)))
        adata.obs["condition"] = ["A"] * 25 + ["B"] * 5
        adata.obs["S_score"] = np.r_[np.repeat(0.1, 25), np.repeat(1.0, 5)]
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = ["G1"] * 25 + ["S"] * 5

        result = diagnose_cell_cycle_regression(adata, condition_key="condition")

        assert result["metrics"]["condition_group_sizes"] == {"A": 25, "B": 5}
        assert result["metrics"]["condition_group_size_imbalance_ratio"] == 5.0
        assert any("imbalanced" in warning for warning in result["warnings"])


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
            config=ScalingConfig(
                scale_method="zscore",
                regress_in_scale=False,
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        meta = result.uns["sclucid"]["preprocess"]["scaling"]
        assert meta["params"]["scale_method"] == "zscore"
        assert meta["output_layer"] == "scaled"

    def test_config_not_mutated(self, minimal_adata):
        adata = minimal_adata.copy()
        config = ScalingConfig(
            scale_method="zscore", regress_in_scale=False, plot=False, report=False, verbose=False
        )
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

    def test_robust_scale_uses_gaussian_consistent_mad(self):
        x = np.array([[1.0], [2.0], [3.0], [4.0], [100.0]])
        scaled = _robust_scale(x, max_value=None)
        median = np.median(x[:, 0])
        mad = np.median(np.abs(x[:, 0] - median))
        expected = (x[:, 0] - median) / (1.4826 * mad)
        np.testing.assert_allclose(scaled[:, 0], expected)

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

    def test_robust_scale_sparse_matches_dense_reference(self):
        x = scipy.sparse.csr_matrix(
            np.array(
                [
                    [0.0, 2.0, 0.0],
                    [1.0, 4.0, 0.0],
                    [3.0, 0.0, 5.0],
                    [0.0, 8.0, 7.0],
                ]
            )
        )

        scaled = _robust_scale_sparse(x, max_value=2)
        expected = _robust_scale(x.toarray(), max_value=2)

        np.testing.assert_allclose(scaled, expected)

    def test_robust_scale_sparse_zero_median_preserves_sparsity(self):
        # Each column has median 0, so scaled zeros remain zero and output can stay sparse.
        x = scipy.sparse.csr_matrix(
            np.array(
                [
                    [0.0, -2.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 2.0, -5.0],
                    [0.0, 0.0, 5.0],
                ]
            )
        )

        scaled = _robust_scale_sparse(x, max_value=None)
        expected = _robust_scale(x.toarray(), max_value=None)

        assert scipy.sparse.isspmatrix_csr(scaled)
        np.testing.assert_allclose(scaled.toarray(), expected)

    def test_robust_scale_sparse_warns_on_small_matrix(self):
        import warnings

        x = scipy.sparse.csr_matrix(np.array([[0.0, 1.0], [2.0, 3.0]]))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _robust_scale_sparse(x, max_value=None)
        assert any("densifying" in str(warning.message).lower() for warning in w)

    def test_robust_scale_sparse_large_path_matches_dense_reference(self, monkeypatch):
        from scLucid.preprocess import scale as scale_module

        # Force the large-matrix sparse-aware path by lowering the threshold.
        monkeypatch.setattr(
            scale_module, "SPARSE_ROBUST_SCALE_DENSE_THRESHOLD_ELEMENTS", 0
        )
        x = scipy.sparse.csr_matrix(
            np.array(
                [
                    [0.0, 2.0, 0.0],
                    [1.0, 4.0, 0.0],
                    [3.0, 0.0, 5.0],
                    [0.0, 8.0, 7.0],
                ]
            )
        )

        scaled = _robust_scale_sparse(x, max_value=2)
        expected = _robust_scale(x.toarray(), max_value=2)

        np.testing.assert_allclose(scaled, expected)
