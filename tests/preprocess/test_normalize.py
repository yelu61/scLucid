"""Unit tests for scLucid.preprocess.normalize."""

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from scLucid.preprocess.config import NormalizationConfig
from scLucid.preprocess.normalize import normalize_data


@pytest.mark.unit
class TestNormalizeDataStandard:
    """Tests for standard normalization method."""

    def test_standard_creates_normalized_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="standard", plot=False, report=False, verbose=False),
        )
        assert "normalized" in result.layers

    def test_standard_restores_to_target_sum(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="standard", target_sum=1e4, plot=False, report=False, verbose=False
            ),
        )
        normalized = result.layers["normalized"]
        if scipy.sparse.issparse(normalized):
            restored = normalized.copy()
            restored.data = np.expm1(restored.data)
            restored_counts = np.asarray(restored.sum(axis=1)).ravel()
        else:
            restored_counts = np.expm1(normalized).sum(axis=1)
        np.testing.assert_allclose(restored_counts, 1e4, rtol=0.1)

    def test_standard_updates_x_when_configured(self, minimal_adata):
        adata = minimal_adata.copy()
        original_x = adata.X.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="standard", update_X=True, plot=False, report=False, verbose=False
            ),
        )
        assert result.X is not None
        # X should differ from original after normalization
        with np.testing.assert_raises(AssertionError):
            np.testing.assert_array_equal(result.X, original_x)

    def test_standard_does_not_update_x_when_disabled(self, minimal_adata):
        adata = minimal_adata.copy()
        original_x = adata.X.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="standard", update_X=False, plot=False, report=False, verbose=False
            ),
        )
        np.testing.assert_array_equal(result.X, original_x)
        assert "normalized" in result.layers

    def test_custom_output_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="standard",
                output_layer="my_norm",
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "my_norm" in result.layers
        assert "normalized" not in result.layers

    def test_force_overwrites_existing_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.layers["normalized"] = adata.X.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="standard", plot=False, report=False, verbose=False),
            force=True,
        )
        assert "normalized" in result.layers

    def test_no_force_skips_when_layer_exists(self, minimal_adata):
        adata = minimal_adata.copy()
        original_layer = adata.X.copy()
        adata.layers["normalized"] = original_layer.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="standard", plot=False, report=False, verbose=False),
            force=False,
        )
        # Should skip normalization and keep original layer
        np.testing.assert_array_equal(result.layers["normalized"], original_layer)

    def test_stores_metadata_in_uns(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="standard", plot=False, report=False, verbose=False),
        )
        assert "sclucid" in result.uns
        assert "preprocess" in result.uns["sclucid"]
        meta = result.uns["sclucid"]["preprocess"]["normalization"]
        assert meta["params"]["method"] == "standard"
        assert "input_stats" in meta
        assert "output_stats" in meta
        assert meta["log_transformed"] is True

    def test_kwargs_override_config(self, minimal_adata):
        adata = minimal_adata.copy()
        config = NormalizationConfig(
            method="standard", target_sum=1e4, plot=False, report=False, verbose=False
        )
        result = normalize_data(adata, config=config, target_sum=1e5)
        meta = result.uns["sclucid"]["preprocess"]["normalization"]
        assert meta["params"]["target_sum"] == 1e5


@pytest.mark.unit
class TestNormalizeDataCLR:
    """Tests for CLR normalization method."""

    def test_clr_creates_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="clr", plot=False, report=False, verbose=False),
        )
        assert "normalized" in result.layers

    def test_clr_mean_centered_per_cell(self, minimal_adata):
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="clr", plot=False, report=False, verbose=False),
        )
        normalized = result.layers["normalized"]
        if scipy.sparse.issparse(normalized):
            cell_means = np.asarray(normalized.mean(axis=1)).ravel()
        else:
            cell_means = normalized.mean(axis=1)
        # CLR centers log-values per cell, so mean should be near zero
        np.testing.assert_allclose(cell_means, 0, atol=1e-6)


@pytest.mark.unit
class TestNormalizeDataPearson:
    """Tests for Pearson residuals normalization."""

    def test_pearson_creates_layer(self, minimal_adata):
        pytest.importorskip("scanpy.experimental.pp")
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="pearson_residuals", plot=False, report=False, verbose=False
            ),
        )
        assert "normalized" in result.layers

    def test_pearson_metadata_flagged_as_log_transformed(self, minimal_adata):
        pytest.importorskip("scanpy.experimental.pp")
        adata = minimal_adata.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(
                method="pearson_residuals", plot=False, report=False, verbose=False
            ),
        )
        meta = result.uns["sclucid"]["preprocess"]["normalization"]
        assert meta["log_transformed"] is True


@pytest.mark.unit
class TestNormalizeDataValidation:
    """Tests for input validation and error handling."""

    def test_rejects_negative_values(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.layers["counts"] = adata.layers["counts"].copy()
        adata.layers["counts"][0, 0] = -1

        with pytest.raises(ValueError, match="negative values"):
            normalize_data(
                adata,
                config=NormalizationConfig(
                    method="standard", plot=False, report=False, verbose=False
                ),
            )

    def test_rejects_reserved_output_layer_x(self):
        with pytest.raises(ValueError, match="reserved"):
            NormalizationConfig(output_layer="X")

    def test_rejects_reserved_output_layer_raw(self):
        with pytest.raises(ValueError, match="reserved"):
            NormalizationConfig(output_layer="raw")

    def test_rejects_missing_input_layer(self, minimal_adata):
        adata = minimal_adata.copy()
        with pytest.raises(ValueError, match="not found"):
            normalize_data(
                adata,
                config=NormalizationConfig(
                    method="standard",
                    input_layer="nonexistent",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )

    def test_rejects_unknown_method_at_config_level(self):
        with pytest.raises(ValueError, match="standard"):
            NormalizationConfig(method="unknown_method")

    def test_rejects_empty_matrix(self):
        empty = AnnData(X=np.zeros((0, 10)))
        empty.layers["counts"] = empty.X.copy()
        with pytest.raises(ValueError, match="empty"):
            normalize_data(
                empty,
                config=NormalizationConfig(
                    method="standard", plot=False, report=False, verbose=False
                ),
            )

    def test_config_not_mutated(self, minimal_adata):
        adata = minimal_adata.copy()
        config = NormalizationConfig(
            method="standard", target_sum=1e4, plot=False, report=False, verbose=False
        )
        original_dict = config.to_dict()
        normalize_data(adata, config=config, target_sum=1e5)
        assert config.to_dict() == original_dict


@pytest.mark.unit
class TestNormalizeDataSparse:
    """Tests specifically for sparse input handling."""

    def test_sparse_standard_normalization(self):
        x = scipy.sparse.random(50, 100, density=0.3, format="csr")
        x.data = np.abs(x.data) * 100  # Make positive counts
        adata = AnnData(X=x.copy())
        adata.layers["counts"] = x.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="standard", plot=False, report=False, verbose=False),
        )
        assert "normalized" in result.layers
        assert scipy.sparse.issparse(result.layers["normalized"])

    def test_sparse_clr_normalization(self):
        x = scipy.sparse.random(50, 100, density=0.3, format="csr")
        x.data = np.abs(x.data) * 100 + 1  # Avoid zeros for CLR
        adata = AnnData(X=x.copy())
        adata.layers["counts"] = x.copy()
        result = normalize_data(
            adata,
            config=NormalizationConfig(method="clr", plot=False, report=False, verbose=False),
        )
        assert "normalized" in result.layers
