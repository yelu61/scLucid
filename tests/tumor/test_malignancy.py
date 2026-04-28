"""Smoke tests for tumor malignancy scoring."""

import numpy as np
import pytest


class TestMalignancyScoring:
    def test_import(self):
        from scLucid.tumor.malignancy.scoring import (
            MalignancyScorer,
            calculate_proliferation_index,
            estimate_metastatic_potential,
            score_malignancy,
        )
        assert callable(score_malignancy)
        assert callable(calculate_proliferation_index)
        assert callable(estimate_metastatic_potential)

    def test_score_malignancy_smoke(self, qc_test_adata):
        from scLucid.tumor.malignancy.scoring import score_malignancy

        result = score_malignancy(qc_test_adata)
        assert result is not None

    def test_malignancy_scorer_init(self):
        from scLucid.tumor.malignancy.scoring import MalignancyScorer

        scorer = MalignancyScorer()
        assert scorer is not None

    def test_proliferation_index_graceful_fallback(self, qc_test_adata):
        from scLucid.tumor.malignancy.scoring import calculate_proliferation_index

        try:
            result = calculate_proliferation_index(qc_test_adata)
            assert result is not None
        except ValueError:
            # Expected when no proliferation genes found in synthetic data
            pass

    def test_metastatic_potential(self, qc_test_adata):
        from scLucid.tumor.malignancy.scoring import estimate_metastatic_potential

        result = estimate_metastatic_potential(qc_test_adata)
        assert result is not None

    def test_classification_import(self):
        from scLucid.tumor.malignancy.classification import (
            MalignancyClassifier,
            classify_malignant_cells,
            score_malignancy_potential,
        )
        assert callable(classify_malignant_cells)
        assert callable(score_malignancy_potential)

    def test_stemness_import(self):
        from scLucid.tumor.malignancy.stemness import (
            calculate_stemness_score,
        )
        assert callable(calculate_stemness_score)
