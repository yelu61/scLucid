"""Deep smoke tests for tumor CNV — execute core algorithm paths."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _make_cnv_adata(n_cells=50, n_genes=100):
    """Create AnnData suitable for CNV inference with cell types."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
    adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]
    adata.obs_names = [f"cell_{i:03d}" for i in range(n_cells)]
    adata.layers["counts"] = X.copy()
    # Add normalized layer (required by infer_cnv)
    adata.layers["normalized"] = np.log1p(X)
    # Add cell type labels for reference-based CNV
    adata.obs["cell_type"] = pd.Categorical(
        ["Normal"] * (n_cells // 2) + ["Tumor"] * (n_cells - n_cells // 2)
    )
    return adata


class TestCNVAnalyzerDeep:
    """Test CNVAnalyzer core fit/predict path."""

    def test_analyzer_fit_with_reference(self):
        """CNVAnalyzer.fit() with explicit reference cells."""
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        analyzer = CNVAnalyzer()
        result = analyzer.fit(adata, reference_cells="Normal", reference_key="cell_type")
        assert result is analyzer
        assert analyzer.cnv_matrix_ is not None
        assert analyzer.cnv_matrix_.shape == (adata.n_obs, adata.n_vars)
        assert analyzer.tumor_scores_ is not None
        assert len(analyzer.tumor_scores_) == adata.n_obs

    def test_analyzer_fit_no_reference(self):
        """CNVAnalyzer.fit() without reference (uses mean)."""
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        analyzer = CNVAnalyzer()
        analyzer.fit(adata)
        assert analyzer.cnv_matrix_ is not None
        assert analyzer.tumor_scores_ is not None

    def test_predict_tumor_cells(self):
        """predict_tumor_cells after fit."""
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        analyzer = CNVAnalyzer()
        analyzer.fit(adata, reference_cells="Normal")
        predictions = analyzer.predict_tumor_cells(threshold=0.5)
        assert len(predictions) == adata.n_obs
        assert predictions.dtype == bool

    def test_predict_without_fit_raises(self):
        """predict_tumor_cells without fit raises ValueError."""
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        analyzer = CNVAnalyzer()
        with pytest.raises(ValueError, match="Must call fit"):
            analyzer.predict_tumor_cells()

    def test_fit_with_gene_order(self):
        """CNVAnalyzer.fit() with chromosome gene order."""
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        adata = _make_cnv_adata(n_cells=30, n_genes=100)
        # Create fake gene order with chromosomes
        gene_order = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 50 + ["chr2"] * 50,
                "start": list(range(100)),
                "end": list(range(100, 200)),
            },
            index=adata.var_names,
        )
        analyzer = CNVAnalyzer(gene_order=gene_order, window_size=2)
        analyzer.fit(adata, reference_cells="Normal")
        assert analyzer.cnv_matrix_ is not None


class TestInferCNVDeep:
    """Test infer_cnv end-to-end."""

    def test_infer_cnv_basic(self):
        """infer_cnv with default parameters."""
        from scLucid.tumor.cnv.infercnv import infer_cnv

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        result = infer_cnv(adata, reference_cells="Normal")
        assert "X_cnv" in result.obsm
        assert "cnv_score" in result.obs
        assert result.obsm["X_cnv"].shape == (adata.n_obs, adata.n_vars)

    def test_infer_cnv_custom_key(self):
        """infer_cnv with custom key_added."""
        from scLucid.tumor.cnv.infercnv import infer_cnv

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        result = infer_cnv(adata, reference_cells="Normal", key_added="my_cnv")
        assert "X_my_cnv" in result.obsm
        assert "my_cnv_score" in result.obs

    def test_infer_cnv_copy(self):
        """infer_cnv with copy=True does not modify original."""
        from scLucid.tumor.cnv.infercnv import infer_cnv

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        original_keys = set(adata.obsm.keys())
        result = infer_cnv(adata, reference_cells="Normal", copy=True)
        assert set(adata.obsm.keys()) == original_keys
        assert "X_cnv" in result.obsm

    def test_infer_cnv_without_normalized_layer(self):
        """infer_cnv falls back to .X when normalized layer absent."""
        from scLucid.tumor.cnv.infercnv import infer_cnv

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        del adata.layers["normalized"]
        result = infer_cnv(adata, reference_cells="Normal")
        assert "X_cnv" in result.obsm
        assert "cnv_score" in result.obs


class TestFindTumorCellsDeep:
    """Test find_tumor_cells with pre-computed CNV."""

    def test_find_by_cnv_score(self):
        """find_tumor_cells with cnv_score method."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, find_tumor_cells

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = find_tumor_cells(adata, method="cnv_score", threshold=0.5)
        assert len(result) == adata.n_obs
        assert result.dtype == bool

    def test_find_by_clustering(self):
        """find_tumor_cells with clustering method."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, find_tumor_cells

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = find_tumor_cells(adata, method="clustering")
        assert len(result) == adata.n_obs
        assert result.dtype == bool


class TestIdentifyClonesDeep:
    """Test identify_clones with pre-computed CNV."""

    def test_identify_clones_hierarchical(self):
        """identify_clones with hierarchical clustering."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, identify_clones

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = identify_clones(adata, n_clusters=3, method="hierarchical")
        assert len(result) == adata.n_obs
        assert result.nunique() == 3

    def test_identify_clones_kmeans(self):
        """identify_clones with kmeans clustering."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, identify_clones

        adata = _make_cnv_adata(n_cells=40, n_genes=80)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = identify_clones(adata, n_clusters=2, method="kmeans")
        assert len(result) == adata.n_obs
        assert result.nunique() == 2


class TestCalculateCNVScoreDeep:
    """Test calculate_cnv_score with all methods."""

    def test_cnv_score_mean_absolute(self):
        """calculate_cnv_score with mean_absolute method."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, calculate_cnv_score

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = calculate_cnv_score(adata, method="mean_absolute")
        assert len(result) == adata.n_obs
        assert (result >= 0).all()

    def test_cnv_score_variance(self):
        """calculate_cnv_score with variance method."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, calculate_cnv_score

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = calculate_cnv_score(adata, method="variance")
        assert len(result) == adata.n_obs
        assert (result >= 0).all()

    def test_cnv_score_gini(self):
        """calculate_cnv_score with gini method."""
        from scLucid.tumor.cnv.infercnv import infer_cnv, calculate_cnv_score

        adata = _make_cnv_adata(n_cells=30, n_genes=60)
        adata = infer_cnv(adata, reference_cells="Normal")
        result = calculate_cnv_score(adata, method="gini")
        assert len(result) == adata.n_obs
        assert (result >= 0).all()
        assert (result <= 1).all()


class TestGiniCoefficient:
    """Unit tests for _gini_coefficient helper."""

    def test_gini_perfect_equality(self):
        from scLucid.tumor.cnv.infercnv import _gini_coefficient

        assert _gini_coefficient(np.ones(10)) == 0.0

    def test_gini_perfect_inequality(self):
        from scLucid.tumor.cnv.infercnv import _gini_coefficient

        result = _gini_coefficient(np.array([0, 0, 0, 0, 100]))
        assert 0 <= result <= 1
        assert result > 0.5

    def test_gini_single_element(self):
        from scLucid.tumor.cnv.infercnv import _gini_coefficient

        result = _gini_coefficient(np.array([5.0]))
        assert result == 0.0

    def test_gini_zeros(self):
        from scLucid.tumor.cnv.infercnv import _gini_coefficient

        result = _gini_coefficient(np.zeros(10))
        assert result == 0.0
