"""
Integration tests for the preprocessing workflow.

Tests the complete preprocessing pipeline from counts to UMAP.
"""

import sys

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.preprocess import run_preprocessing
from scLucid.preprocess.config import NormalizationConfig, WorkflowConfig
from scLucid.preprocess.gene_biotype import (
    annotate_gene_biotypes,
    apply_gene_biotype_strategy,
    filter_genes_by_biotype,
    list_gene_biotype_resources,
    load_gene_biotypes,
)

# Import synthetic data fixtures


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
        from scipy import sparse

        from scLucid.preprocess.normalize import normalize_data

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
        from scLucid.preprocess.config import HVGConfig
        from scLucid.preprocess.hvg import find_hvgs

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
        from scLucid.preprocess.config import ScalingConfig
        from scLucid.preprocess.scale import scale_data

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

    def test_gene_biotype_annotation_and_filtering_live_under_preprocess(self, minimal_adata):
        """Gene biotype utilities should store metadata under preprocess, not QC."""
        adata = minimal_adata.copy()
        adata.var_names = np.array(
            ["GAPDH", "MALAT1", "RPLP0", "MT-CO1", "IGHG1", "TRAC"]
            + [f"GENE{i}" for i in range(adata.n_vars - 6)]
        )
        custom_df = pd.DataFrame(
            {
                "gene_name": ["GAPDH", "MALAT1", "RPLP0", "MT-CO1", "IGHG1", "TRAC"],
                "biotype": [
                    "protein_coding",
                    "lncRNA",
                    "protein_coding",
                    "protein_coding",
                    "IG_C_gene",
                    "TR_C_gene",
                ],
            }
        )

        adata = annotate_gene_biotypes(adata, biotype_df=custom_df, method="custom")
        assert "preprocess" in adata.uns["sclucid"]
        assert "gene_biotypes" in adata.uns["sclucid"]["preprocess"]
        assert "qc" not in adata.uns["sclucid"] or "gene_biotypes" not in adata.uns["sclucid"].get(
            "qc", {}
        )

        filtered = filter_genes_by_biotype(adata, keep_biotypes=["protein_coding"], copy=True)
        assert filtered is not None
        assert (filtered.var["biotype_category"] == "protein_coding").all()

    def test_apply_gene_biotype_strategy_supports_custom_path(self, minimal_adata, tmp_path):
        adata = minimal_adata[:, :3].copy()
        adata.var_names = np.array(["GAPDH", "MALAT1", "TRAC"])
        custom_path = tmp_path / "custom_biotypes.tsv"
        pd.DataFrame(
            {
                "external_gene_name": ["GAPDH", "MALAT1", "TRAC"],
                "gene_biotype": ["protein_coding", "lncRNA", "TR_C_gene"],
            }
        ).to_csv(custom_path, sep="\t", index=False)

        result = apply_gene_biotype_strategy(
            adata,
            method="custom",
            custom_biotype_path=custom_path,
            keep_biotypes=["protein_coding"],
            copy=True,
        )

        assert result is not None
        assert list(result.var_names) == ["GAPDH"]

    def test_load_gene_biotypes_prefers_cache_without_network(self, tmp_path):
        cache_dir = tmp_path / "gene_annotations"
        cache_dir.mkdir()
        cache_file = cache_dir / "human_reference_latest.csv.gz"
        pd.DataFrame(
            {
                "gene_name": ["GAPDH", "MALAT1"],
                "biotype": ["protein_coding", "lncRNA"],
                "gene_id": ["ENSG1", "ENSG2"],
            }
        ).to_csv(cache_file, index=False, compression="gzip")

        df, meta = load_gene_biotypes(
            species="human",
            allow_download=False,
            cache_dir=cache_dir,
            prefer_bundled=False,
            return_metadata=True,
        )

        assert list(df["gene_name"]) == ["GAPDH", "MALAT1"]
        assert meta["source"].startswith("cache:")

    def test_load_gene_biotypes_uses_bundled_resource(self):
        df, meta = load_gene_biotypes(
            species="human",
            allow_download=False,
            return_metadata=True,
        )

        assert {"gene_name", "biotype"}.issubset(df.columns)
        assert len(df) > 1000
        assert meta["source"].startswith("package:")

    def test_list_gene_biotype_resources_reports_cache(self, tmp_path):
        cache_dir = tmp_path / "gene_annotations"
        cache_dir.mkdir()
        (cache_dir / "mouse_reference_latest.csv").write_text(
            "gene_name,biotype\nGapdh,protein_coding\n"
        )

        resources = list_gene_biotype_resources(species="mouse", cache_dir=cache_dir)

        assert "mouse" in resources
        assert resources["mouse"]["cached"]

    def test_list_gene_biotype_resources_reports_bundled_resources(self):
        resources = list_gene_biotype_resources(species="human")

        assert "human" in resources
        assert resources["human"]["bundled"]

    def test_annotate_gene_biotypes_records_reference_source(self, minimal_adata, tmp_path):
        cache_dir = tmp_path / "gene_annotations"
        cache_dir.mkdir()
        pd.DataFrame(
            {
                "gene_name": list(minimal_adata.var_names[:3]),
                "biotype": ["protein_coding", "lncRNA", "protein_coding"],
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            }
        ).to_csv(cache_dir / "human_ensembl_latest.csv.gz", index=False, compression="gzip")

        adata = minimal_adata[:, :3].copy()
        adata = annotate_gene_biotypes(
            adata,
            species="human",
            method="reference",
            allow_download=False,
            cache_dir=cache_dir,
            prefer_bundled=False,
            overwrite=True,
        )

        meta = adata.uns["sclucid"]["preprocess"]["gene_biotypes"]
        assert meta["reference_source"].startswith("cache:")

    def test_annotate_gene_biotypes_accepts_legacy_ensembl_alias(self, minimal_adata, tmp_path):
        cache_dir = tmp_path / "gene_annotations"
        cache_dir.mkdir()
        pd.DataFrame(
            {
                "gene_name": list(minimal_adata.var_names[:2]),
                "biotype": ["protein_coding", "lncRNA"],
                "gene_id": ["ENSG1", "ENSG2"],
            }
        ).to_csv(cache_dir / "human_reference_latest.csv.gz", index=False, compression="gzip")

        adata = minimal_adata[:, :2].copy()
        adata = annotate_gene_biotypes(
            adata,
            species="human",
            method="ensembl",
            allow_download=False,
            cache_dir=cache_dir,
            prefer_bundled=False,
            overwrite=True,
        )

        meta = adata.uns["sclucid"]["preprocess"]["gene_biotypes"]
        assert meta["reference_source"].startswith("cache:")

    def test_workflow_run_flags_are_honored(self, minimal_adata):
        """Workflow-level run_* flags should disable optional downstream steps."""
        config = _workflow_config_for_tests()
        config.run_regression = False
        config.run_scaling = False
        config.run_pca = False
        config.run_neighbors = False
        config.run_integration = False

        result = run_preprocessing(minimal_adata, config=config)

        preprocess_meta = result.uns["sclucid"]["preprocess"]
        assert preprocess_meta["steps_executed"] == [
            "normalization",
            "set_raw",
            "hvg_selection",
            "subset_hvg",
        ]
        assert config.scaling.vars_to_regress == ["total_counts", "pct_counts_mt"]
        assert "regressed" not in result.layers
        assert "scaled" not in result.layers
        assert "X_pca" not in result.obsm
        assert "X_umap" not in result.obsm
        assert "regress_inline" not in preprocess_meta

    def test_workflow_honors_custom_layer_names(self, minimal_adata):
        """Workflow should use top-level layer naming consistently across steps."""
        config = _workflow_config_for_tests()
        config.run_pca = False
        config.run_neighbors = False
        config.normalized_layer = "lognorm"
        config.regressed_layer = "resid"
        config.scaled_layer = "zscore"

        result = run_preprocessing(minimal_adata, config=config)

        assert "lognorm" in result.layers
        assert "resid" in result.layers
        assert "zscore" in result.layers
        assert result.raw is not None
        assert result.raw.shape[1] >= result.n_vars
        assert result.uns["sclucid"]["preprocess"]["normalization"]["output_layer"] == "lognorm"
        assert result.uns["sclucid"]["preprocess"]["regress"]["output_layer"] == "resid"
        assert result.uns["sclucid"]["preprocess"]["scaling"]["output_layer"] == "zscore"
        assert "regress_inline" not in result.uns["sclucid"]["preprocess"]


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
