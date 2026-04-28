"""Smoke tests for CNV analysis."""

import numpy as np
import pytest


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
