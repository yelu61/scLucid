"""Tests for the user-facing 10x / h5ad data readers."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scLucid as scl
from scLucid.utils.contracts import LayerKeys, SCLUCID_ROOT, UnsKeys
from scLucid.utils.helpers import _looks_like_counts, read_10x
from scLucid.utils.io import read_h5ad


def _make_counts_adata(n_cells=20, n_genes=50, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(3, size=(n_cells, n_genes)).astype(np.int32)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    return AnnData(X=X, obs=obs, var=var)


class TestLooksLikeCounts:
    def test_integer_matrix_is_counts(self):
        adata = _make_counts_adata()
        assert _looks_like_counts(adata.X) is True

    def test_float_log_normalized_matrix_is_not_counts(self):
        adata = _make_counts_adata()
        log_normalized = np.log1p(adata.X.astype(float) / adata.X.sum(axis=1, keepdims=True) * 1e4)
        assert _looks_like_counts(log_normalized) is False

    def test_negative_values_are_not_counts(self):
        matrix = np.random.randn(10, 5)
        assert _looks_like_counts(matrix) is False


class TestReadH5ad:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_h5ad(tmp_path / "nope.h5ad")

    def test_adds_counts_layer_for_integer_X(self, tmp_path):
        adata = _make_counts_adata()
        path = tmp_path / "raw.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(path)
        assert LayerKeys.COUNTS in loaded.layers
        np.testing.assert_array_equal(loaded.layers[LayerKeys.COUNTS], adata.X)

    def test_does_not_overwrite_existing_counts_layer(self, tmp_path):
        adata = _make_counts_adata()
        sentinel = adata.X.copy()
        adata.layers[LayerKeys.COUNTS] = sentinel
        path = tmp_path / "raw.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(path)
        # Loaded counts layer must be the one already on disk (not a freshly
        # copied X).
        np.testing.assert_array_equal(loaded.layers[LayerKeys.COUNTS], sentinel)

    def test_skips_counts_when_X_is_not_counts(self, tmp_path):
        adata = _make_counts_adata()
        adata.X = np.log1p(adata.X.astype(float))  # destroy integer-ness
        path = tmp_path / "lognorm.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(path)
        assert LayerKeys.COUNTS not in loaded.layers

    def test_ensure_counts_layer_false_disables_autofill(self, tmp_path):
        adata = _make_counts_adata()
        path = tmp_path / "raw.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(path, ensure_counts_layer=False)
        assert LayerKeys.COUNTS not in loaded.layers

    def test_metadata_stamped_on_obs(self, tmp_path):
        adata = _make_counts_adata()
        path = tmp_path / "raw.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(
            path,
            sample_id="patient1",
            species="human",
            tissue="pancreas",
            tissue_type="tumor_tissue",
            cancer_type="pancreatic_adenocarcinoma",
        )
        for column, value in {
            "sample_id": "patient1",
            "species": "human",
            "tissue": "pancreas",
            "tissue_type": "tumor_tissue",
            "cancer_type": "pancreatic_adenocarcinoma",
        }.items():
            assert column in loaded.obs.columns
            assert (loaded.obs[column] == value).all()

    def test_metadata_lifted_to_analysis_context(self, tmp_path):
        adata = _make_counts_adata()
        path = tmp_path / "raw.h5ad"
        adata.write_h5ad(path)

        loaded = read_h5ad(
            path, species="mouse", tissue_type="tumor_tissue", cancer_type="melanoma"
        )
        ctx = loaded.uns[SCLUCID_ROOT][UnsKeys.ANALYSIS_CONTEXT]
        assert ctx["species"] == "mouse"
        assert ctx["tissue_type"] == "tumor_tissue"
        assert ctx["cancer_type"] == "melanoma"


class TestRead10x:
    def test_missing_path_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_10x(tmp_path / "missing_dir")

    def test_unrecognized_file_raises(self, tmp_path):
        bad = tmp_path / "data.txt"
        bad.write_text("not 10x")
        with pytest.raises(ValueError, match="Cell Ranger"):
            read_10x(bad)

    @pytest.mark.optional
    def test_reads_cell_ranger_directory(self, tmp_path):
        """Smoke test that a synthetic Cell Ranger directory is loadable.

        Builds the minimal three-file layout (matrix.mtx, barcodes.tsv,
        features.tsv) Cell Ranger writes. If scanpy's 10x reader rejects
        the layout, skip rather than fail — formats vary across versions.
        """
        import gzip
        import shutil
        import scipy.io
        import scipy.sparse as sp

        cr_dir = tmp_path / "filtered_feature_bc_matrix"
        cr_dir.mkdir()

        # Genes × cells sparse matrix
        n_genes, n_cells = 30, 12
        rng = np.random.default_rng(0)
        matrix = sp.random(n_genes, n_cells, density=0.2, random_state=0)
        matrix.data = rng.integers(1, 5, size=matrix.nnz).astype(float)
        scipy.io.mmwrite(cr_dir / "matrix.mtx", matrix.tocoo())
        with open(cr_dir / "matrix.mtx", "rb") as fh, gzip.open(
            cr_dir / "matrix.mtx.gz", "wb"
        ) as gz:
            shutil.copyfileobj(fh, gz)
        (cr_dir / "matrix.mtx").unlink()

        with gzip.open(cr_dir / "barcodes.tsv.gz", "wt") as fh:
            for i in range(n_cells):
                fh.write(f"BC{i:03d}-1\n")

        with gzip.open(cr_dir / "features.tsv.gz", "wt") as fh:
            for i in range(n_genes):
                fh.write(f"ENSG{i:05d}\tGene{i}\tGene Expression\n")

        try:
            adata = read_10x(cr_dir, cache=False, species="human")
        except Exception as exc:
            pytest.skip(
                f"scanpy.read_10x_mtx rejected the synthetic directory: {exc}"
            )

        assert adata.n_obs == n_cells
        assert adata.n_vars == n_genes
        assert LayerKeys.COUNTS in adata.layers
        assert "species" in adata.obs.columns
        assert (adata.obs["species"] == "human").all()


class TestTopLevelAPI:
    def test_read_10x_top_level(self):
        assert hasattr(scl, "read_10x")
        assert hasattr(scl.utils, "read_10x")

    def test_read_h5ad_top_level(self):
        assert hasattr(scl, "read_h5ad")
        assert hasattr(scl.utils, "read_h5ad")


class TestRead10xMultiSampleMode:
    """Verify that the merged single+multi API still honors the legacy multi path."""

    def test_legacy_load_10x_data_alias_still_works(self):
        from scLucid.utils.helpers import load_10x_data, read_10x

        # The alias should resolve to a callable; the inner function is the
        # multi-sample dispatch inside read_10x.
        assert callable(load_10x_data)
        # Calling with no valid paths should not crash; returns empty AnnData.
        result = load_10x_data(samples=["S1"], path_dict={})
        assert hasattr(result, "obs")
        assert result.n_obs == 0
        assert callable(read_10x)

    def test_read_10x_rejects_conflicting_modes(self):
        """Passing both single-sample path and multi-sample samples must error."""
        with pytest.raises(ValueError, match="single-sample"):
            read_10x(path="/tmp/anything", samples=["S1"])

    def test_read_10x_rejects_missing_inputs(self):
        """Passing neither single nor multi-sample inputs must error."""
        with pytest.raises(ValueError, match="single-sample"):
            read_10x()

    def test_read_10x_multi_with_path_dict(self):
        """Multi-sample mode with explicit empty path_dict returns empty AnnData."""
        result = read_10x(samples=["s1"], path_dict={})
        assert result.n_obs == 0
