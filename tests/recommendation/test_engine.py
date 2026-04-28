"""Tests for the unified RecommendationEngine."""

from __future__ import annotations

import importlib
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.qc.intelligent_qc import (
    QCRecommendation,
    StrategyType,
    ThresholdRecommendation,
)
from scLucid.recommendation.config import RecommendationConfig
from scLucid.recommendation.engine import (
    RecommendationEngine,
    recommend_analysis_parameters,
)
from scLucid.recommendation.schema import ParameterRecommendation, RecommendationSection


def _make_qc_recommendation() -> QCRecommendation:
    return QCRecommendation(
        min_genes=ThresholdRecommendation(
            threshold=200, ci_lower=180, ci_upper=220,
            method="gmm", confidence=0.85,
            evidence={"bic": 450},
        ),
        max_mt_percent=ThresholdRecommendation(
            threshold=20.0, ci_lower=15.0, ci_upper=25.0,
            method="lognormal_fit", confidence=0.9,
            evidence={"shape": 1.2},
        ),
        doublet_threshold=ThresholdRecommendation(
            threshold=0.3, ci_lower=0.2, ci_upper=0.4,
            method="distribution_fit", confidence=0.75,
            evidence={"n_cells": 500},
        ),
        n_counts=ThresholdRecommendation(
            threshold=500.0, ci_lower=400.0, ci_upper=600.0,
            method="lognormal", confidence=0.8,
            evidence={"mean": 520},
        ),
        overall_strategy=StrategyType.STANDARD,
        overall_confidence=0.83,
        data_quality_score=85.0,
        concerns=[],
        tumor_specific_considerations=[],
    )


def _make_preprocess_mock_strategy():
    """Create a mock PreprocessingStrategy with required attributes."""
    from scLucid.preprocess.intelligent.data_classes import (
        BatchCorrectionRecommendation,
        DataProfile,
        HVGRecommendation,
        NeighborsRecommendation,
        PCARecommendation,
        PreprocessingStrategy,
        ResolutionRecommendation,
    )
    profile = DataProfile(
        n_cells=500, n_genes=2000, sparsity=0.85,
        median_counts_per_cell=2500, median_genes_per_cell=800,
        is_sparse=False, is_small_dataset=True,
        is_medium_dataset=False, is_large_dataset=False,
        has_batch_info=False, data_quality_score=75.0,
        strategy_type="standard",
    )
    return PreprocessingStrategy(
        data_profile=profile,
        hvg=HVGRecommendation(
            n_top_genes=1500, variance_explained=0.65,
            ci_lower=1200, ci_upper=1800,
            method="elbow", confidence=0.8,
        ),
        pca=PCARecommendation(
            n_pcs=30, variance_explained=0.85,
            ci_lower=20, ci_upper=40,
            method="elbow", confidence=0.85,
        ),
        neighbors=NeighborsRecommendation(
            n_neighbors=15, n_pcs=30, silhouette_score=0.45,
            ci_lower_neighbors=10, ci_upper_neighbors=20,
            ci_lower_pcs=20, ci_upper_pcs=40,
            method="grid_search", confidence=0.8,
        ),
        resolution=ResolutionRecommendation(
            resolution=1.2, n_clusters=8, stability_score=0.75,
            ci_lower=0.8, ci_upper=1.5,
            method="stability_analysis", confidence=0.78,
        ),
        batch_correction=None,
        overall_confidence=0.82,
        concerns=[],
        recommendations=[],
    )


class TestAdaptQC:
    def test_adapt_qc(self):
        engine = RecommendationEngine()
        qc = _make_qc_recommendation()
        section = engine._adapt_qc(qc)

        assert isinstance(section, RecommendationSection)
        assert section.name == "qc"
        assert section.confidence == qc.overall_confidence
        assert len(section.parameters) == 4

        param_names = {p.name for p in section.parameters}
        assert param_names == {"min_genes", "max_mt_percent", "doublet_threshold", "n_counts"}

    def test_adapt_qc_parameter_values(self):
        engine = RecommendationEngine()
        qc = _make_qc_recommendation()
        section = engine._adapt_qc(qc)

        min_genes = section.get_parameter("min_genes")
        assert min_genes.value == 200
        assert min_genes.ci_lower == 180
        assert min_genes.ci_upper == 220

    def test_adapt_qc_metadata(self):
        engine = RecommendationEngine()
        qc = _make_qc_recommendation()
        section = engine._adapt_qc(qc)
        assert section.metadata["data_quality_score"] == 85.0
        assert section.metadata["strategy"] == "standard"


class TestAdaptPreprocess:
    @mock.patch.object(RecommendationEngine, "_recommend_preprocess")
    def test_adapt_preprocess_no_batch(self, mock_rec):
        mock_strategy = _make_preprocess_mock_strategy()
        mock_rec.return_value = mock_strategy

        engine = RecommendationEngine()
        section = engine._adapt_preprocess(mock_strategy)

        assert section.name == "preprocess"
        assert section.confidence == mock_strategy.overall_confidence
        assert section.metadata["strategy_type"] == "standard"

        param_names = {p.name for p in section.parameters}
        assert "n_top_genes" in param_names
        assert "n_pcs" in param_names
        assert "n_neighbors" in param_names
        assert "graph_n_pcs" in param_names
        # No batch_correction param when batch_correction is None
        assert "batch_correction_method" not in param_names

    def test_adapt_preprocess_with_batch(self):
        from scLucid.preprocess.intelligent.data_classes import (
            BatchCorrectionRecommendation,
            DataProfile,
            HVGRecommendation,
            NeighborsRecommendation,
            PCARecommendation,
            PreprocessingStrategy,
            ResolutionRecommendation,
        )

        batch = BatchCorrectionRecommendation(
            needs_correction=True,
            severity_score=0.7,
            recommended_method="harmony",
            alternative_methods=["bbknn"],
            method_scores={"harmony": 0.8},
            confidence=0.85,
        )
        profile = DataProfile(
            n_cells=1000, n_genes=2000, sparsity=0.85,
            median_counts_per_cell=2500, median_genes_per_cell=800,
            is_sparse=False, is_small_dataset=False,
            is_medium_dataset=True, is_large_dataset=False,
            has_batch_info=True, n_batches=2,
            data_quality_score=75.0, strategy_type="standard",
        )
        strategy = PreprocessingStrategy(
            data_profile=profile,
            hvg=HVGRecommendation(
                n_top_genes=1500, variance_explained=0.65,
                ci_lower=1200, ci_upper=1800,
                method="elbow", confidence=0.8,
            ),
            pca=PCARecommendation(
                n_pcs=30, variance_explained=0.85,
                ci_lower=20, ci_upper=40,
                method="elbow", confidence=0.85,
            ),
            neighbors=NeighborsRecommendation(
                n_neighbors=15, n_pcs=30, silhouette_score=0.45,
                ci_lower_neighbors=10, ci_upper_neighbors=20,
                ci_lower_pcs=20, ci_upper_pcs=40,
                method="grid_search", confidence=0.8,
            ),
            resolution=ResolutionRecommendation(
                resolution=1.2, n_clusters=8, stability_score=0.75,
                ci_lower=0.8, ci_upper=1.5,
                method="stability_analysis", confidence=0.78,
            ),
            batch_correction=batch,
            overall_confidence=0.82,
        )

        engine = RecommendationEngine()
        section = engine._adapt_preprocess(strategy)

        param_names = {p.name for p in section.parameters}
        assert "batch_correction_method" in param_names

        batch_param = section.get_parameter("batch_correction_method")
        assert batch_param.value == "harmony"
        assert section.metadata["batch_effect_severity"] == 0.7


class TestRecommend:
    """Integration-style tests exercising the full recommend() flow with mocks."""

    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_recommend_qc_only(self, mock_qc, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc"])
        )
        result = engine.recommend(minimal_adata, tissue_type="normal")

        assert "qc" in result.sections
        assert result.sections["qc"].confidence == 0.83
        assert result.overall_confidence > 0
        mock_qc.assert_called_once()

    @mock.patch("scLucid.recommendation.engine.core.IntelligentPreprocessRecommender")
    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_recommend_qc_and_preprocess(self, mock_qc, mock_preprocess, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        mock_instance = mock.MagicMock()
        mock_instance.recommend.return_value = _make_preprocess_mock_strategy()
        mock_preprocess.return_value = mock_instance

        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc", "preprocess"])
        )
        result = engine.recommend(minimal_adata, tissue_type="normal")

        assert "qc" in result.sections
        assert "preprocess" in result.sections

    @mock.patch("scLucid.recommendation.engine.core.IntelligentPreprocessRecommender")
    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_recommend_all_modules(self, mock_qc, mock_preprocess, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        mock_instance = mock.MagicMock()
        mock_instance.recommend.return_value = _make_preprocess_mock_strategy()
        mock_preprocess.return_value = mock_instance

        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc", "preprocess", "clustering", "annotation", "tumor"])
        )
        result = engine.recommend(minimal_adata, tissue_type="normal")

        assert "qc" in result.sections
        assert "preprocess" in result.sections
        assert "clustering" in result.sections
        assert "annotation" in result.sections
        assert "tumor" in result.sections

    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_overall_confidence_average(self, mock_qc, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc"])
        )
        result = engine.recommend(minimal_adata, tissue_type="normal")
        # Single section: overall = section confidence
        assert result.overall_confidence == pytest.approx(0.83, abs=0.01)

    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_concerns_deduplication(self, mock_qc, minimal_adata):
        qc = _make_qc_recommendation()
        # Add same concern to both
        qc.concerns = ["shared concern"]
        mock_qc.return_value = qc

        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc"])
        )
        result = engine.recommend(minimal_adata, tissue_type="normal")
        # Check concerns are deduplicated (since we only have qc, just verify it's present)
        assert "shared concern" in result.concerns

    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_recommend_with_tissue_type_tumor(self, mock_qc, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        engine = RecommendationEngine(
            config=RecommendationConfig(modules=["qc"])
        )
        result = engine.recommend(minimal_adata, tissue_type="tumor")
        mock_qc.assert_called_once()
        assert mock_qc.call_args[1]["tissue_type"] == "tumor"


class TestBuildContext:
    def test_build_context(self, minimal_adata):
        engine = RecommendationEngine()
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="normal_tissue")
        result = engine._build_context(minimal_adata, analysis_context=ctx, batch_key=None)

        assert result["n_cells"] == minimal_adata.n_obs
        assert result["n_genes"] == minimal_adata.n_vars


class TestAdaptClusteringSearch:
    def test_empty_dataframe(self):
        engine = RecommendationEngine()
        df = pd.DataFrame(columns=["resolution", "n_clusters", "silhouette", "stability"])
        section = engine._adapt_clustering_search(df, None)

        assert section.name == "clustering"
        assert section.confidence == 0.0
        assert len(section.concerns) == 1

    def test_with_data(self):
        engine = RecommendationEngine()
        eval_data = {
            "resolution": [0.5, 1.0, 1.5],
            "n_clusters": [5, 10, 15],
            "silhouette": [0.3, 0.5, 0.4],
            "stability": [0.6, 0.8, 0.7],
        }
        df = pd.DataFrame(eval_data)
        section = engine._adapt_clustering_search(df, 1.0)

        assert section.name == "clustering"
        assert section.confidence > 0
        resolution_param = section.get_parameter("resolution")
        assert resolution_param is not None
        assert resolution_param.value == 1.0

    def test_with_marker_abundance(self):
        engine = RecommendationEngine()
        eval_data = {
            "resolution": [0.5, 1.0, 1.5],
            "n_clusters": [5, 10, 15],
            "silhouette": [0.3, 0.5, 0.4],
            "marker_abundance": [5, 12, 8],
        }
        df = pd.DataFrame(eval_data)
        section = engine._adapt_clustering_search(df, 1.0)
        assert section.confidence > 0


class TestConvenienceFunction:
    @mock.patch("scLucid.recommendation.engine.core.recommend_intelligent_qc")
    def test_recommend_analysis_parameters(self, mock_qc, minimal_adata):
        mock_qc.return_value = _make_qc_recommendation()
        result = recommend_analysis_parameters(
            minimal_adata,
            tissue_type="normal",
            config=RecommendationConfig(modules=["qc"]),
        )
        assert "qc" in result.sections


class TestPrepareRepresentation:
    def test_prepare_representation(self, minimal_adata):
        engine = RecommendationEngine()
        result = engine._prepare_representation(minimal_adata.copy(), use_rep="X_pca")

        assert result is not None
        assert "X_pca" in result.obsm


class TestCelltypistDetection:
    def test_celltypist_available(self):
        engine = RecommendationEngine()
        # This should not crash; result depends on environment
        available = engine._celltypist_available()
        assert isinstance(available, bool)

    @mock.patch.object(RecommendationEngine, "_celltypist_available", return_value=True)
    def test_assess_marker_support(self, mock_ct, minimal_adata):
        engine = RecommendationEngine()
        result = engine._assess_marker_support(
            minimal_adata,
            marker_species="human",
            marker_tissue="blood",
        )
        assert "total_types" in result
        assert "eligible_types" in result
        assert "eligible_ratio" in result
        assert 0 <= result["eligible_ratio"] <= 1

    @mock.patch.object(RecommendationEngine, "_celltypist_available", return_value=True)
    def test_evaluate_existing_celltypist_no_columns(self, mock_ct, minimal_adata):
        engine = RecommendationEngine()
        result = engine._evaluate_existing_celltypist_evidence(
            minimal_adata, cluster_key="nonexistent"
        )
        assert result["mean_confidence"] == 0.0
        assert result["label_source"] is None


class TestAnnotationStrategySelection:
    def test_select_annotation_strategy_no_celltypist(self):
        engine = RecommendationEngine()
        result = engine._select_annotation_strategy(
            marker_evidence={
                "eligible_ratio": 0.5,
                "cluster_marker_signal": 0.4,
            },
            celltypist_evidence={"mean_confidence": 0.0, "cluster_purity": 0.0},
            celltypist_available=False,
            existing_celltypist=False,
            expected_n_clusters=10,
        )
        assert result["run_celltypist"] is False
        assert result["final_method"] in ("max_score", "combined", "enrichment")

    def test_select_annotation_strategy_with_celltypist_available(self):
        engine = RecommendationEngine()
        result = engine._select_annotation_strategy(
            marker_evidence={
                "eligible_ratio": 0.5,
                "cluster_marker_signal": 0.35,
            },
            celltypist_evidence={"mean_confidence": 0.0, "cluster_purity": 0.0},
            celltypist_available=True,
            existing_celltypist=False,
            expected_n_clusters=10,
        )
        assert result["run_celltypist"] is True
        assert result["final_method"] in ("hybrid", "celltypist")

    def test_select_annotation_strategy_existing_celltypist_good(self):
        engine = RecommendationEngine()
        result = engine._select_annotation_strategy(
            marker_evidence={
                "eligible_ratio": 0.6,
                "cluster_marker_signal": 0.5,
                "cluster_best_labels": {},
            },
            celltypist_evidence={
                "mean_confidence": 0.85,
                "cluster_purity": 0.8,
                "cluster_labels": {},
            },
            celltypist_available=True,
            existing_celltypist=True,
            expected_n_clusters=8,
        )
        assert result["run_celltypist"] is True

    def test_select_annotation_strategy_existing_celltypist_poor(self):
        engine = RecommendationEngine()
        result = engine._select_annotation_strategy(
            marker_evidence={
                "eligible_ratio": 0.7,
                "cluster_marker_signal": 0.5,
                "cluster_best_labels": {"0": "T cell"},
            },
            celltypist_evidence={
                "mean_confidence": 0.5,
                "cluster_purity": 0.5,
                "cluster_labels": {"0": "B cell"},
            },
            celltypist_available=True,
            existing_celltypist=True,
            expected_n_clusters=8,
        )
        if result["run_celltypist"]:
            # With poor CellTypist but good marker signal, should favor marker-based
            assert result["final_method"] in ("max_score", "combined", "hybrid", "marker_method")
