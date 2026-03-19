"""
Integration tests for the preprocessing workflow.

Tests the complete preprocessing pipeline from counts to UMAP.
"""

import pytest
import numpy as np
from anndata import AnnData

import sys
sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

import scLucid as scl
from scLucid.preprocess import run_preprocessing
from scLucid.preprocess.config import WorkflowConfig, NormalizationConfig

# Import synthetic data fixtures
from tests.fixtures.synthetic_data import (
    synthetic_generator,
    minimal_adata,
    integration_test_adata,
)


def _workflow_config_for_tests() -> WorkflowConfig:
    """Create a workflow config that avoids interactive plotting/report side effects."""
    config = WorkflowConfig()
    config.integration.method = None
    config.hvg.flavor = "seurat"
    config.hvg.n_top_genes = 100
    for section in [
        config.normalization,
        config.hvg,
        config.scaling,
        config.integration,
        config.graph,
    ]:
        if hasattr(section, "plot"):
            section.plot = False
        if hasattr(section, "report"):
            section.report = False
        if hasattr(section, "verbose"):
            section.verbose = False
    return config


@pytest.mark.integration
class TestPreprocessingWorkflow:
    """Test suite for the preprocessing workflow."""

    def test_workflow_runs_without_errors(self, minimal_adata):
        """Test that the basic workflow runs without errors."""
        config = _workflow_config_for_tests()

        result = run_preprocessing(minimal_adata, config=config)

        assert isinstance(result, AnnData)
        assert result.n_obs == minimal_adata.n_obs
        assert result.n_vars <= minimal_adata.n_vars  # HVG selection reduces genes

    def test_workflow_creates_expected_layers(self, minimal_adata):
        """Test that workflow creates expected layers."""
        config = _workflow_config_for_tests()

        result = run_preprocessing(minimal_adata, config=config)

        # Check required layers exist
        assert "counts" in result.layers
        assert "normalized" in result.layers
        assert "scaled" in result.layers

        # Check .raw is set
        assert result.raw is not None

    def test_workflow_creates_expected_obsm(self, minimal_adata):
        """Test that workflow creates expected obsm entries."""
        config = _workflow_config_for_tests()

        result = run_preprocessing(minimal_adata, config=config)

        # Check PCA and UMAP
        assert "X_pca" in result.obsm
        assert "X_umap" in result.obsm

        # Check dimensions
        assert result.obsm["X_pca"].shape[1] == config.graph.n_pcs
        assert result.obsm["X_umap"].shape[1] == 2

    def test_workflow_creates_expected_uns(self, minimal_adata):
        """Test that workflow stores configuration in uns."""
        config = _workflow_config_for_tests()

        result = run_preprocessing(minimal_adata, config=config)

        # Check sclucid metadata
        assert "sclucid" in result.uns
        assert "preprocess" in result.uns["sclucid"]
        assert "workflow_config" in result.uns["sclucid"]["preprocess"]

    def test_normalization_step(self, minimal_adata):
        """Test that normalization produces correct results."""
        from scLucid.preprocess.normalize import normalize_data
        from scipy import sparse

        config = NormalizationConfig(target_sum=1e4, plot=False, report=False, verbose=False)
        result = normalize_data(minimal_adata.copy(), config=config)

        # Check normalized layer exists
        assert "normalized" in result.layers

        # normalize_data returns log1p-transformed normalized values; invert to validate target_sum.
        normalized = result.layers["normalized"]
        if sparse.issparse(normalized):
            restored = normalized.copy()
            restored.data = np.expm1(restored.data)
            restored_counts = np.asarray(restored.sum(axis=1)).ravel()
        else:
            restored_counts = np.expm1(normalized).sum(axis=1)
        np.testing.assert_allclose(restored_counts, 1e4, rtol=0.1)

    def test_hvg_selection(self, minimal_adata):
        """Test that HVG selection reduces gene count."""
        from scLucid.preprocess.hvg import find_hvgs
        from scLucid.preprocess.config import HVGConfig

        # First normalize
        from scLucid.preprocess.normalize import normalize_data
        adata = normalize_data(
            minimal_adata.copy(),
            config=NormalizationConfig(plot=False, report=False, verbose=False),
        )

        # Then find HVGs
        hvg_config = HVGConfig(n_top_genes=100, plot=False, report=False, verbose=False)
        result = find_hvgs(adata, config=hvg_config, input_layer="normalized")

        # Check HVGs are marked under the configured output key.
        hvg_meta = result.uns["sclucid"]["preprocess"]["hvg"]
        output_key = hvg_meta["output_key"]
        assert output_key in result.var.columns
        n_hvg = int(result.var[output_key].sum())
        assert n_hvg == hvg_meta["n_hvg"]
        assert 0 < n_hvg <= hvg_config.n_top_genes

    def test_scaling_step(self, minimal_adata):
        """Test that scaling produces zero-mean, unit-variance data."""
        from scLucid.preprocess.scale import scale_data
        from scLucid.preprocess.config import ScalingConfig

        config = ScalingConfig(regress_in_scale=False, plot=False, report=False, verbose=False)
        result = scale_data(minimal_adata.copy(), config=config)

        # Check scaled layer
        assert "scaled" in result.layers

        # Check approximate zero mean and unit variance
        # (using a tolerance since we're working with HVGs)
        scaled_data = result.layers["scaled"]
        mean = np.mean(scaled_data, axis=0)
        std = np.std(scaled_data, axis=0)

        np.testing.assert_allclose(mean, 0, atol=1e-6)
        np.testing.assert_allclose(std, 1, atol=0.5)  # Allow some variation

    def test_pca_step(self, minimal_adata):
        """Test that PCA produces correct output dimensions."""
        import scanpy as sc

        # Need scaled data first
        from scLucid.preprocess.scale import scale_data
        scaling_cfg = _workflow_config_for_tests().scaling
        scaling_cfg.regress_in_scale = False
        adata = scale_data(minimal_adata.copy(), config=scaling_cfg)

        n_pcs = 30
        sc.tl.pca(adata, n_comps=n_pcs)

        assert "X_pca" in adata.obsm
        assert adata.obsm["X_pca"].shape == (adata.n_obs, n_pcs)
        assert "PCs" in adata.varm

    def test_neighbors_and_umap(self, minimal_adata):
        """Test that neighbors and UMAP work correctly."""
        import scanpy as sc

        # Prepare data
        config = _workflow_config_for_tests()
        result = run_preprocessing(minimal_adata, config=config)

        # Check neighbors
        assert "neighbors" in result.uns
        assert "connectivities" in result.obsp
        assert "distances" in result.obsp

        # Check UMAP
        assert "X_umap" in result.obsm
        assert result.obsm["X_umap"].shape == (result.n_obs, 2)


@pytest.mark.integration
@pytest.mark.slow
class TestPreprocessingWithBatchEffects:
    """Test preprocessing with batch correction."""

    def test_batch_correction_harmony(self, integration_test_adata):
        """Test Harmony batch correction if available."""
        pytest.importorskip("harmonypy")

        config = _workflow_config_for_tests()
        config.integration.method = "harmony"
        config.integration.batch_key = "batch"

        result = run_preprocessing(integration_test_adata, config=config)

        # Check that batch-corrected representation exists
        assert "X_harmony" in result.obsm

    def test_batch_correction_scanorama(self, integration_test_adata):
        """Test Scanorama batch correction if available."""
        pytest.importorskip("scanorama")

        config = _workflow_config_for_tests()
        config.integration.method = "scanorama"
        config.integration.batch_key = "batch"

        result = run_preprocessing(integration_test_adata, config=config)

        # Check that batch-corrected representation exists
        assert "X_scanorama" in result.obsm

    def test_without_batch_correction(self, integration_test_adata):
        """Test workflow without batch correction."""
        config = _workflow_config_for_tests()

        result = run_preprocessing(integration_test_adata, config=config)

        # Should still complete successfully
        assert "X_umap" in result.obsm
        assert "X_pca" in result.obsm


@pytest.mark.integration
class TestPreprocessingConfigValidation:
    """Test configuration validation."""

    def test_invalid_target_sum(self):
        """Test that negative target_sum raises error."""
        with pytest.raises(ValueError):
            NormalizationConfig(target_sum=-1000)

    def test_invalid_n_top_genes(self):
        """Test that invalid n_top_genes raises validation error."""
        from scLucid.preprocess.config import HVGConfig

        # n_top_genes below minimum should raise validation error
        with pytest.raises(ValueError):
            HVGConfig(n_top_genes=50)

        # Valid n_top_genes should work
        config = HVGConfig(n_top_genes=100)
        assert config.n_top_genes == 100

    def test_layer_naming_conflict(self):
        """Test that reserved layer names are rejected."""
        from scLucid.preprocess.config import NormalizationConfig

        # Should raise error for reserved name
        with pytest.raises(ValueError):
            NormalizationConfig(output_layer="X")


@pytest.mark.integration
class TestPreprocessingIdempotency:
    """Test that running preprocessing multiple times is safe."""

    def test_rerunning_workflow(self, minimal_adata):
        """Test that rerunning workflow doesn't corrupt data."""
        config = _workflow_config_for_tests()

        # Run once
        result1 = run_preprocessing(minimal_adata.copy(), config=config)

        # Run again
        result2 = run_preprocessing(result1, config=config)

        # Should still have valid results
        assert "X_umap" in result2.obsm
        assert "X_pca" in result2.obsm


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
