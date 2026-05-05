"""Deep smoke tests for tumor malignancy — execute core algorithm paths."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _make_adata_with_genes(genes, n_cells=50):
    """Create a minimal AnnData with specific gene names."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, len(genes))).astype(np.float32)
    adata = AnnData(X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i:03d}" for i in range(n_cells)]
    adata.layers["counts"] = X.copy()
    return adata


class TestMalignancyScorerDeep:
    """Test MalignancyScorer core path with real gene names."""

    def test_scorer_fit_with_default_genes(self):
        """MalignancyScorer.fit() with default gene sets present."""
        from scLucid.tumor.malignancy.scoring import MalignancyScorer

        genes = (
            ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1"]
            + ["MYC", "KRAS", "EGFR", "BRAF", "PIK3CA"]
            + ["TP53", "PTEN", "RB1", "CDKN2A", "APC"]
            + ["other_gene_1", "other_gene_2"]
        )
        adata = _make_adata_with_genes(genes)
        scorer = MalignancyScorer()
        result = scorer.fit(adata)
        assert result is scorer
        assert scorer.scores_ is not None
        assert len(scorer.scores_) == adata.n_obs
        assert scorer.scores_.min() >= 0
        assert scorer.scores_.max() <= 1

    def test_scorer_fit_no_matching_genes(self):
        """Graceful fallback when no default genes present."""
        from scLucid.tumor.malignancy.scoring import MalignancyScorer

        adata = _make_adata_with_genes(["gene_a", "gene_b", "gene_c"])
        scorer = MalignancyScorer()
        scorer.fit(adata)
        assert scorer.scores_ is not None
        # When no genes match, proliferation and oncogene scores are 0,
        # but tumor suppressor contributes (1 - 0) * 0.25 = 0.25
        assert scorer.scores_.notna().all()
        np.testing.assert_array_almost_equal(
            scorer.scores_.values, np.full(adata.n_obs, 0.25)
        )

    def test_score_malignancy_integration(self):
        """score_malignancy end-to-end with real gene names."""
        from scLucid.tumor.malignancy.scoring import score_malignancy

        genes = ["MKI67", "PCNA", "MYC", "KRAS", "TP53", "PTEN", "gene_x"]
        adata = _make_adata_with_genes(genes)
        result = score_malignancy(adata, key_added="malignancy")
        assert "malignancy" in result.obs.columns
        assert result.obs["malignancy"].notna().all()

    def test_proliferation_index_with_match(self):
        """calculate_proliferation_index with matching genes."""
        from scLucid.tumor.malignancy.scoring import calculate_proliferation_index

        genes = ["MKI67", "PCNA", "TOP2A", "AURKA", "CDK1", "gene_y"]
        adata = _make_adata_with_genes(genes)
        result = calculate_proliferation_index(adata, gene_set="classic")
        assert len(result) == adata.n_obs
        assert result.notna().all()

    def test_proliferation_index_no_match_raises(self):
        """calculate_proliferation_index raises when no genes match."""
        from scLucid.tumor.malignancy.scoring import calculate_proliferation_index

        adata = _make_adata_with_genes(["gene_a", "gene_b"])
        with pytest.raises(ValueError, match="No proliferation genes"):
            calculate_proliferation_index(adata)

    def test_metastatic_potential_with_emt_genes(self):
        """estimate_metastatic_potential with EMT genes present."""
        from scLucid.tumor.malignancy.scoring import estimate_metastatic_potential

        genes = ["VIM", "CDH2", "FN1", "SNAI1", "MMP2", "MMP9", "gene_z"]
        adata = _make_adata_with_genes(genes)
        result = estimate_metastatic_potential(adata)
        assert len(result) == adata.n_obs
        assert result.notna().all()

    def test_metastatic_potential_no_match(self):
        """estimate_metastatic_potential returns zeros when no genes match."""
        from scLucid.tumor.malignancy.scoring import estimate_metastatic_potential

        adata = _make_adata_with_genes(["gene_a", "gene_b"])
        result = estimate_metastatic_potential(adata)
        assert len(result) == adata.n_obs
        np.testing.assert_array_almost_equal(result.values, np.zeros(adata.n_obs))


class TestMalignancyClassificationDeep:
    def test_classify_malignant_cells(self):
        from scLucid.tumor.malignancy.classification import classify_malignant_cells

        genes = ["MKI67", "PCNA", "MYC", "TP53"]
        adata = _make_adata_with_genes(genes)
        result = classify_malignant_cells(adata)
        assert result is not None

    def test_score_malignancy_potential(self):
        from scLucid.tumor.malignancy.classification import score_malignancy_potential

        genes = ["MKI67", "PCNA", "MYC", "TP53"]
        adata = _make_adata_with_genes(genes)
        result = score_malignancy_potential(adata)
        assert result is not None


class TestStemnessDeep:
    def test_calculate_stemness_score(self):
        from scLucid.tumor.malignancy.stemness import calculate_stemness_score

        genes = ["PROM1", "NANOG", "SOX2", "ALDH1A1"]
        adata = _make_adata_with_genes(genes)
        result = calculate_stemness_score(adata)
        assert result is not None
