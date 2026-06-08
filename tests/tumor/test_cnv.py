"""Smoke tests for CNV analysis."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


class TestCNV:
    def test_import(self):
        from scLucid.tumor.cnv.infercnv import (
            CNVAnalyzer,
            calculate_cnv_score,
            find_tumor_cells,
            identify_clones,
            infer_cnv,
            _gini_coefficient,
        )
        assert callable(infer_cnv)
        assert callable(calculate_cnv_score)
        assert callable(find_tumor_cells)
        assert callable(identify_clones)

    def test_gini_coefficient(self):
        from scLucid.tumor.cnv.infercnv import _gini_coefficient

        # Perfect equality
        assert _gini_coefficient(np.ones(10)) == 0.0
        # Perfect inequality
        result = _gini_coefficient(np.array([0, 0, 0, 0, 100]))
        assert 0 <= result <= 1
        assert result > 0.5

    def test_cnv_analyzer_init(self):
        from scLucid.tumor.cnv.infercnv import CNVAnalyzer

        analyzer = CNVAnalyzer()
        assert analyzer is not None

    def test_calculate_cnv_score_smoke(self, qc_test_adata):
        from scLucid.tumor.cnv.infercnv import calculate_cnv_score

        try:
            result = calculate_cnv_score(qc_test_adata)
            assert result is not None
        except (ValueError, KeyError):
            # Expected when reference cells not specified or genes missing
            pass

    def test_clone_analysis_import(self):
        from scLucid.tumor.cnv.clone_analysis import CloneAnalyzer

        assert CloneAnalyzer

    def test_cnv_signature_import(self):
        from scLucid.tumor.cnv.cnv_signature import CNVSigExtractor

        assert CNVSigExtractor

    def test_infer_cnv_records_input_quality(self):
        from scLucid.tumor.cnv.infercnv import infer_cnv

        rng = np.random.default_rng(42)
        n_cells, n_genes = 100, 80
        X = rng.integers(0, 20, size=(n_cells, n_genes)).astype(np.float32)
        adata = AnnData(X)
        adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]
        adata.obs_names = [f"cell_{i:03d}" for i in range(n_cells)]
        adata.layers["counts"] = X.copy()
        adata.obs["cell_type"] = pd.Categorical(["Normal"] * 50 + ["Tumor"] * 50)

        result = infer_cnv(adata, reference_cells="Normal")
        assert "cnv_summary" in result.uns
        quality = result.uns["cnv_summary"]["input_quality"]
        assert quality["n_cells_input"] == n_cells
        assert quality["reference_cells_used"] == 50
        assert quality["has_genomic_coordinates"] is False

    def test_infer_cnv_records_genomic_coordinates(self):
        from scLucid.tumor.cnv.infercnv import infer_cnv

        rng = np.random.default_rng(42)
        n_cells, n_genes = 100, 80
        X = rng.integers(0, 20, size=(n_cells, n_genes)).astype(np.float32)
        adata = AnnData(X)
        adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]
        adata.obs_names = [f"cell_{i:03d}" for i in range(n_cells)]
        adata.layers["counts"] = X.copy()
        adata.obs["cell_type"] = pd.Categorical(["Normal"] * 50 + ["Tumor"] * 50)
        adata.var["chromosome"] = ["1"] * 40 + ["2"] * 40
        adata.var["start"] = list(range(n_genes))
        adata.var["end"] = list(range(1, n_genes + 1))

        result = infer_cnv(adata, reference_cells="Normal")
        quality = result.uns["cnv_summary"]["input_quality"]
        assert quality["has_genomic_coordinates"] is True
        assert "1" in quality["chromosomes_used"]
        assert "2" in quality["chromosomes_used"]
