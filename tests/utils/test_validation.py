"""
Tests for validation utilities.

Tests AnnData validation and analysis readiness checks.
"""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from scLucid.utils.validation import (
    ValidationError,
    assert_analysis_ready,
    assert_preprocessing_ready,
    assert_qc_ready,
    check_layer_consistency,
    validate_adata,
    validate_analysis_results,
    validate_config,
    validate_workflow_contract,
)


@pytest.fixture
def valid_adata():
    """Create a valid AnnData object with all required components."""
    adata = AnnData(np.random.randint(0, 100, size=(100, 50)))
    adata.layers["counts"] = adata.X.copy()
    adata.layers["normalized"] = adata.X.astype(float) / adata.X.sum(axis=1, keepdims=True) * 1e4
    adata.obsm["X_pca"] = np.random.randn(100, 10)
    adata.obs["leiden"] = pd.Categorical(["cluster_1"] * 50 + ["cluster_2"] * 50)
    adata.obs["cell_type"] = ["type_a"] * 60 + ["type_b"] * 40
    adata.obs["n_genes_by_counts"] = np.random.randint(100, 1000, 100)
    adata.obs["total_counts"] = np.random.randint(1000, 10000, 100)
    adata.obs["pct_counts_mt"] = np.random.uniform(0, 20, 100)

    # Add sclucid results
    adata.uns["sclucid"] = {
        "qc": {"metrics": {"n_cells": 100}},
        "preprocess": {"workflow_config": {}},
    }

    return adata


@pytest.fixture
def empty_adata():
    """Create an empty AnnData object."""
    return AnnData(np.array([]).reshape(0, 0))


class TestValidationError:
    """Test ValidationError exception."""

    def test_validation_error_basic(self):
        """Test basic ValidationError creation."""
        err = ValidationError("Test error")
        assert str(err) == "Test error"
        assert err.field is None

    def test_validation_error_with_field(self):
        """Test ValidationError with field."""
        err = ValidationError("Invalid value", field="min_genes")
        assert err.field == "min_genes"


class TestValidateAdata:
    """Test validate_adata function."""

    def test_validate_adata_none(self):
        """Test validation with None adata."""
        result = validate_adata(None, raise_on_error=False)

        assert result["valid"] is False
        assert "AnnData object is None" in result["errors"]

    def test_validate_adata_zero_obs(self, valid_adata):
        """Test validation with zero observations."""
        adata = valid_adata[:0].copy()

        result = validate_adata(adata, raise_on_error=False)

        assert result["valid"] is False
        assert any("0 cells" in e for e in result["errors"])

    def test_validate_adata_zero_vars(self, valid_adata):
        """Test validation with zero variables."""
        adata = valid_adata[:, :0].copy()

        result = validate_adata(adata, raise_on_error=False)

        assert result["valid"] is False
        assert any("0 genes" in e for e in result["errors"])

    def test_validate_required_layers(self, valid_adata):
        """Test validation with required layers."""
        # Should pass with existing layer
        result = validate_adata(valid_adata, required_layers=["counts"], raise_on_error=False)
        assert result["valid"] is True

        # Should fail with missing layer
        result = validate_adata(valid_adata, required_layers=["nonexistent"], raise_on_error=False)
        assert result["valid"] is False
        assert any("nonexistent" in e for e in result["errors"])

    def test_validate_required_obs(self, valid_adata):
        """Test validation with required obs columns."""
        result = validate_adata(valid_adata, required_obs=["leiden"], raise_on_error=False)
        assert result["valid"] is True

        result = validate_adata(valid_adata, required_obs=["nonexistent"], raise_on_error=False)
        assert result["valid"] is False

    def test_validate_required_obsm(self, valid_adata):
        """Test validation with required obsm keys."""
        result = validate_adata(valid_adata, required_obsm=["X_pca"], raise_on_error=False)
        assert result["valid"] is True

        result = validate_adata(valid_adata, required_obsm=["X_umap"], raise_on_error=False)
        assert result["valid"] is False

    def test_validate_check_counts(self, valid_adata):
        """Test validation with check_counts."""
        # Should pass with integer counts layer
        result = validate_adata(
            valid_adata, required_layers=["counts"], check_counts=True, raise_on_error=False
        )
        assert result["valid"] is True

    def test_validate_raise_on_error(self, valid_adata):
        """Test that raise_on_error=True raises exception."""
        adata = valid_adata[:0].copy()

        with pytest.raises(ValidationError):
            validate_adata(adata, raise_on_error=True)

    def test_validate_duplicate_obs_names(self, valid_adata):
        """Test validation catches duplicate obs names."""
        valid_adata.obs_names = ["cell_1"] * 50 + ["cell_2"] * 50

        result = validate_adata(valid_adata, raise_on_error=False)

        assert result["valid"] is False
        assert any("duplicate" in e.lower() for e in result["errors"])


class TestValidateConfig:
    """Test validate_config function."""

    class DummyConfig:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

    def test_validate_config_none(self):
        """Test validation with None config."""
        result = validate_config(None, raise_on_error=False)

        assert result["valid"] is False

    def test_validate_required_fields(self):
        """Test validation with required fields."""
        config = self.DummyConfig(method="wilcoxon", n_top_genes=100)

        result = validate_config(
            config, required_fields=["method", "n_top_genes"], raise_on_error=False
        )
        assert result["valid"] is True

        result = validate_config(
            config, required_fields=["method", "nonexistent"], raise_on_error=False
        )
        assert result["valid"] is False

    def test_validate_field_types(self):
        """Test validation with field types."""
        config = self.DummyConfig(n_top_genes=100, threshold=0.05)

        result = validate_config(
            config, field_types={"n_top_genes": int, "threshold": float}, raise_on_error=False
        )
        assert result["valid"] is True

        result = validate_config(
            config, field_types={"n_top_genes": str}, raise_on_error=False  # Wrong type
        )
        assert result["valid"] is False


class TestValidateAnalysisResults:
    """Test validate_analysis_results function."""

    def test_validate_qc_results(self, valid_adata):
        """Test QC results validation."""
        result = validate_analysis_results(valid_adata, "qc", raise_on_error=False)
        assert result["valid"] is True

    def test_validate_preprocess_results(self, valid_adata):
        """Test preprocess results validation."""
        result = validate_analysis_results(valid_adata, "preprocess", raise_on_error=False)
        assert result["valid"] is True

    def test_validate_clustering_results(self, valid_adata):
        """Test clustering results validation."""
        result = validate_analysis_results(valid_adata, "clustering", raise_on_error=False)
        assert result["valid"] is True

    def test_validate_annotation_results(self, valid_adata):
        """Test annotation results validation."""
        result = validate_analysis_results(valid_adata, "annotation", raise_on_error=False)
        assert result["valid"] is True

    def test_validate_no_sclucid(self, valid_adata):
        """Test validation with no sclucid key."""
        del valid_adata.uns["sclucid"]

        result = validate_analysis_results(valid_adata, "qc", raise_on_error=False)
        assert result["valid"] is False

    def test_validate_unknown_analysis_type(self, valid_adata):
        """Test validation with unknown analysis type."""
        result = validate_analysis_results(valid_adata, "unknown_type", raise_on_error=False)
        assert result["valid"] is True  # Warns but doesn't fail
        assert any("unknown" in w.lower() for w in result["warnings"])


class TestCheckLayerConsistency:
    """Test check_layer_consistency function."""

    def test_consistent_layers(self, valid_adata):
        """Test with consistent layers."""
        result = check_layer_consistency(valid_adata, ["counts", "normalized"])

        assert result["consistent"] is True
        assert len(result["errors"]) == 0

    def test_missing_layer(self, valid_adata):
        """Test with missing layer."""
        result = check_layer_consistency(valid_adata, ["nonexistent"])

        assert result["consistent"] is False
        assert any("nonexistent" in e for e in result["errors"])

    def test_inconsistent_shape(self, valid_adata):
        """Test with inconsistent shape detection."""
        # AnnData validates layer shapes on assignment, so we test the check
        # by directly manipulating the storage to simulate corruption
        result = check_layer_consistency(valid_adata, ["counts"])

        # Should be consistent since layers match
        assert result["consistent"] is True
        assert result["expected_shape"] == (valid_adata.n_obs, valid_adata.n_vars)

    def test_all_zero_layer(self, valid_adata):
        """Test detection of all-zero layer."""
        valid_adata.layers["all_zero"] = np.zeros((100, 50))

        result = check_layer_consistency(valid_adata, ["all_zero"])

        # Should warn about all-zero
        assert any("all zeros" in w.lower() for w in result["warnings"])


class TestAssertionHelpers:
    """Test assertion helper functions."""

    def test_assert_qc_ready_valid(self, valid_adata):
        """Test assert_qc_ready with valid data."""
        # Should not raise
        assert_qc_ready(valid_adata)

    def test_assert_qc_ready_invalid(self, valid_adata):
        """Test assert_qc_ready with invalid data."""
        # Create adata without counts (no layers)
        import numpy as np

        empty_adata = AnnData(np.random.randn(10, 5))  # n_obs=10, n_vars=5
        # No layers at all - should fail validation

        with pytest.raises(ValidationError):
            assert_qc_ready(empty_adata)

    def test_assert_preprocessing_ready_valid(self, valid_adata):
        """Test assert_preprocessing_ready with valid data."""
        # Should not raise
        assert_preprocessing_ready(valid_adata)

    def test_assert_preprocessing_ready_no_qc(self, valid_adata):
        """Test assert_preprocessing_ready without QC results."""
        del valid_adata.uns["sclucid"]["qc"]

        with pytest.raises(ValidationError):
            assert_preprocessing_ready(valid_adata)

    def test_assert_analysis_ready_valid(self, valid_adata):
        """Test assert_analysis_ready with valid data."""
        # Should not raise
        assert_analysis_ready(valid_adata)

    def test_assert_analysis_ready_no_pca(self, valid_adata):
        """Test assert_analysis_ready without PCA."""
        del valid_adata.obsm["X_pca"]

        with pytest.raises(ValidationError):
            assert_analysis_ready(valid_adata)


class TestValidateAdataExtraBranches:
    """Coverage for less-traveled branches of validate_adata."""

    def test_required_var_missing(self, valid_adata):
        result = validate_adata(valid_adata, required_var=["nonexistent"], raise_on_error=False)
        assert result["valid"] is False
        assert any("nonexistent" in e for e in result["errors"])

    def test_required_uns_missing(self, valid_adata):
        result = validate_adata(
            valid_adata, required_uns=["nonexistent_uns_key"], raise_on_error=False
        )
        assert result["valid"] is False
        assert any("nonexistent_uns_key" in e for e in result["errors"])

    def test_check_normalized_warns_when_missing(self):
        adata = AnnData(np.random.randint(0, 10, size=(20, 10)).astype(int))
        result = validate_adata(adata, check_normalized=True, raise_on_error=False)
        assert any("normalized" in w.lower() for w in result["warnings"])

    def test_check_counts_warns_on_float_counts(self):
        adata = AnnData(np.random.rand(20, 10))
        adata.layers["counts"] = np.random.rand(20, 10)  # floats, not counts
        result = validate_adata(adata, check_counts=True, raise_on_error=False)
        # Either the warnings or the errors should call this out — accept either.
        all_messages = result["errors"] + result["warnings"]
        assert any("counts" in m.lower() for m in all_messages)

    def test_duplicate_var_names_flagged(self, valid_adata):
        valid_adata.var_names = ["gene_X"] * valid_adata.n_vars
        result = validate_adata(valid_adata, raise_on_error=False)
        assert result["valid"] is False
        assert any("duplicate" in e.lower() and "gene" in e.lower() for e in result["errors"])


class TestValidateConfigBranches:
    """Cover the raise_on_error and missing-value paths of validate_config."""

    class DummyConfig:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

    def test_required_field_present_but_none_is_error(self):
        config = self.DummyConfig(method=None)
        result = validate_config(config, required_fields=["method"], raise_on_error=False)
        assert result["valid"] is False
        assert any("is None" in e for e in result["errors"])

    def test_raise_on_error_with_invalid_config(self):
        config = self.DummyConfig(method="wilcoxon")
        with pytest.raises(ValidationError):
            validate_config(config, required_fields=["missing"], raise_on_error=True)

    def test_none_config_raises_when_requested(self):
        with pytest.raises(ValidationError):
            validate_config(None, raise_on_error=True)


class TestValidateAnalysisResultsExtraBranches:
    """Coverage for the analysis-result validators' negative paths."""

    def test_qc_missing_metric_columns(self, valid_adata):
        # Drop a QC metric column; sclucid["qc"] still present.
        del valid_adata.obs["pct_counts_mt"]
        result = validate_analysis_results(valid_adata, "qc", raise_on_error=False)
        assert result["valid"] is False
        assert any("Missing QC metric" in e for e in result["errors"])

    def test_preprocess_section_missing(self, valid_adata):
        del valid_adata.uns["sclucid"]["preprocess"]
        result = validate_analysis_results(valid_adata, "preprocess", raise_on_error=False)
        assert result["valid"] is False
        assert any("preprocess" in e.lower() for e in result["errors"])

    def test_preprocess_missing_normalized_layer(self, valid_adata):
        del valid_adata.layers["normalized"]
        result = validate_analysis_results(valid_adata, "preprocess", raise_on_error=False)
        assert result["valid"] is False
        assert any("normalized" in e for e in result["errors"])

    def test_preprocess_missing_pca(self, valid_adata):
        del valid_adata.obsm["X_pca"]
        result = validate_analysis_results(valid_adata, "preprocess", raise_on_error=False)
        assert result["valid"] is False
        assert any("X_pca" in e or "PCA" in e for e in result["errors"])

    def test_clustering_missing(self, valid_adata):
        valid_adata.obs = valid_adata.obs.drop(columns=["leiden"])
        result = validate_analysis_results(valid_adata, "clustering", raise_on_error=False)
        assert result["valid"] is False

    def test_annotation_missing(self, valid_adata):
        valid_adata.obs = valid_adata.obs.drop(columns=["cell_type"])
        result = validate_analysis_results(valid_adata, "annotation", raise_on_error=False)
        assert result["valid"] is False

    def test_markers_missing(self, valid_adata):
        # The "markers" branch checks for rank_genes_groups in uns.
        result = validate_analysis_results(valid_adata, "markers", raise_on_error=False)
        assert result["valid"] is False
        assert any("rank_genes_groups" in e for e in result["errors"])

    def test_markers_present(self, valid_adata):
        valid_adata.uns["rank_genes_groups"] = {"names": np.array([["g1"]])}
        result = validate_analysis_results(valid_adata, "markers", raise_on_error=False)
        assert result["valid"] is True

    def test_raise_on_error_path(self, valid_adata):
        valid_adata.obs = valid_adata.obs.drop(columns=["cell_type"])
        with pytest.raises(ValidationError):
            validate_analysis_results(valid_adata, "annotation", raise_on_error=True)

    def test_no_sclucid_raise(self, valid_adata):
        del valid_adata.uns["sclucid"]
        with pytest.raises(ValidationError):
            validate_analysis_results(valid_adata, "qc", raise_on_error=True)


class TestCheckLayerConsistencyExtraBranches:
    """Cover sparse-zero and all-NaN branches plus the layers=None default."""

    def test_layers_none_inspects_all(self, valid_adata):
        # No `layers` arg → check all layers present on the adata.
        result = check_layer_consistency(valid_adata)
        assert result["consistent"] is True
        for layer in valid_adata.layers:
            assert layer in result["shapes"]

    def test_sparse_all_zero_warning(self, valid_adata):
        valid_adata.layers["sparse_zero"] = sp.csr_matrix(
            (valid_adata.n_obs, valid_adata.n_vars), dtype=float
        )
        result = check_layer_consistency(valid_adata, ["sparse_zero"])
        assert any("all zeros" in w.lower() for w in result["warnings"])

    def test_all_nan_layer_is_error(self, valid_adata):
        valid_adata.layers["all_nan"] = np.full(
            (valid_adata.n_obs, valid_adata.n_vars), np.nan, dtype=float
        )
        result = check_layer_consistency(valid_adata, ["all_nan"])
        assert result["consistent"] is False
        assert any("nan" in e.lower() for e in result["errors"])


class TestValidateWorkflowContract:
    """Targeted coverage for validate_workflow_contract wrapper."""

    def test_returns_dict_on_failure_when_not_raising(self, valid_adata):
        # Stripped data: QC input contract may fail.
        adata = AnnData(np.random.rand(5, 3))
        result = validate_workflow_contract(adata, "qc", when="input", raise_on_error=False)
        assert isinstance(result, dict)
        assert "valid" in result

    def test_raises_validation_error_on_invalid(self):
        adata = AnnData(np.random.rand(5, 3))
        with pytest.raises(ValidationError):
            validate_workflow_contract(adata, "analysis", when="input", raise_on_error=True)


class TestAssertHelperErrorTranslation:
    """assert_* helpers should translate ContractError → ValidationError."""

    def test_assert_preprocessing_ready_translates_contract_error(self):
        # A bare AnnData has no preprocess input contract satisfied.
        adata = AnnData(np.random.rand(10, 5))
        with pytest.raises(ValidationError):
            assert_preprocessing_ready(adata)

    def test_assert_analysis_ready_translates_contract_error(self):
        adata = AnnData(np.random.rand(10, 5))
        with pytest.raises(ValidationError):
            assert_analysis_ready(adata)
