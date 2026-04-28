"""Smoke tests for therapy response prediction."""

import pytest


class TestTherapy:
    def test_import(self):
        from scLucid.tumor.therapy.prediction import (
            ResponsePredictor,
            evaluate_biomarker,
            predict_therapy_response,
            stratify_patients,
        )
        assert callable(predict_therapy_response)
        assert callable(stratify_patients)
        assert callable(evaluate_biomarker)

    def test_response_predictor_init(self):
        from scLucid.tumor.therapy.prediction import ResponsePredictor

        predictor = ResponsePredictor()
        assert predictor is not None

    def test_predict_therapy_response_smoke(self, qc_test_adata):
        from scLucid.tumor.therapy.prediction import predict_therapy_response

        try:
            result = predict_therapy_response(qc_test_adata)
            assert result is not None
        except (ValueError, KeyError):
            # Expected when cancer type or biomarkers not specified
            pass

    def test_resistance_import(self):
        from scLucid.tumor.therapy.resistance import (
            ResistanceAnalyzer,
        )
        assert ResistanceAnalyzer

    def test_target_import(self):
        from scLucid.tumor.therapy.target import TargetDiscovery

        assert TargetDiscovery
