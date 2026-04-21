"""Targeted semantic tests for preprocess robustness and traceability."""

import numpy as np
import pytest

from scLucid.preprocess.config import (
    HVGConfig,
    IntegrationConfig,
    NormalizationConfig,
    WorkflowConfig,
)
from scLucid.preprocess.hvg import find_hvgs
from scLucid.preprocess.integrate import batch_correction
from scLucid.preprocess.intelligent.data_classes import (
    BatchCorrectionRecommendation,
    DataProfile,
    HVGRecommendation,
    NeighborsRecommendation,
    PCARecommendation,
    PreprocessingStrategy,
    ResolutionRecommendation,
)
from scLucid.preprocess.intelligent.recommender import run_intelligent_preprocessing
from scLucid.preprocess.normalize import normalize_data
from scLucid.preprocess.workflow import run_preprocessing


@pytest.mark.unit
def test_normalization_rejects_negative_values(minimal_adata):
    adata = minimal_adata.copy()
    adata.layers["counts"] = adata.layers["counts"].copy()
    adata.layers["counts"][0, 0] = -1

    with pytest.raises(ValueError, match="negative values"):
        normalize_data(
            adata,
            config=NormalizationConfig(
                method="standard",
                plot=False,
                report=False,
                verbose=False,
            ),
        )


@pytest.mark.unit
def test_hvg_raises_for_missing_input_layer(minimal_adata):
    with pytest.raises(KeyError, match="missing_layer"):
        find_hvgs(
            minimal_adata.copy(),
            config=HVGConfig(n_top_genes=100, plot=False, report=False, verbose=False),
            input_layer="missing_layer",
        )


@pytest.mark.integration
def test_workflow_hvg_selection_uses_normalized_layer(minimal_adata):
    config = WorkflowConfig()
    config.integration.method = None
    config.hvg.n_top_genes = 100
    config.hvg.plot = False
    config.hvg.report = False
    config.normalization.plot = False
    config.normalization.report = False
    config.scaling.plot = False
    config.scaling.report = False
    config.graph.plot = False
    config.graph.report = False
    config.run_pca = False
    config.run_neighbors = False

    result = run_preprocessing(minimal_adata, config=config)

    assert result.uns["sclucid"]["preprocess"]["hvg"]["input_layer"] == config.normalized_layer


@pytest.mark.unit
def test_batch_correction_preserves_method_metadata(monkeypatch, minimal_adata):
    import scLucid.preprocess.integrate as integrate_module

    adata = minimal_adata.copy()
    adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
    adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

    def fake_harmony(adata, covariate_keys, basis, embedding_key, **kwargs):
        adata.obsm[embedding_key] = adata.obsm[basis].copy()
        adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
            "integration", {}
        )["harmony"] = {
            "covariate_keys": covariate_keys,
            "output_dims": adata.obsm[embedding_key].shape[1],
        }
        return adata

    monkeypatch.setattr(integrate_module, "_integrate_harmony", fake_harmony)

    result = batch_correction(
        adata,
        config=IntegrationConfig(
            method="harmony",
            batch_key="batch",
            use_rep="X_pca",
            plot=False,
            report=False,
            verbose=False,
        ),
    )

    integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
    assert "harmony" in integration_meta
    assert "workflow" in integration_meta
    assert integration_meta["workflow"]["output_key"] == "X_harmony"
    assert integration_meta["workflow"]["method"] == "harmony"


@pytest.mark.unit
def test_strategy_to_config_preserves_pydantic_base_and_batch_key():
    strategy = PreprocessingStrategy(
        data_profile=DataProfile(
            n_cells=100,
            n_genes=2000,
            sparsity=0.9,
            median_counts_per_cell=1500,
            median_genes_per_cell=700,
            is_sparse=True,
            is_small_dataset=True,
            is_medium_dataset=False,
            is_large_dataset=False,
        ),
        hvg=HVGRecommendation(
            n_top_genes=1500,
            variance_explained=0.8,
            ci_lower=1200,
            ci_upper=1800,
            method="test",
            confidence=0.9,
        ),
        pca=PCARecommendation(
            n_pcs=25,
            variance_explained=0.75,
            ci_lower=20,
            ci_upper=30,
            method="test",
            confidence=0.8,
        ),
        neighbors=NeighborsRecommendation(
            n_neighbors=20,
            n_pcs=25,
            silhouette_score=0.4,
            ci_lower_neighbors=15,
            ci_upper_neighbors=25,
            ci_lower_pcs=20,
            ci_upper_pcs=30,
            method="test",
            confidence=0.8,
        ),
        resolution=ResolutionRecommendation(
            resolution=0.8,
            n_clusters=8,
            stability_score=0.7,
            ci_lower=0.6,
            ci_upper=1.0,
            method="test",
            confidence=0.8,
        ),
        batch_correction=BatchCorrectionRecommendation(
            needs_correction=True,
            severity_score=0.6,
            recommended_method="scanorama",
            confidence=0.9,
            evidence={"batch_key": "batch_col"},
        ),
    )

    base_config = WorkflowConfig()
    base_config.hvg.flavor = "seurat_v3"
    base_config.integration.batch_key = "original_batch"

    applied = strategy.to_config(base_config=base_config)

    assert base_config.hvg.n_top_genes == 2000
    assert base_config.integration.batch_key == "original_batch"
    assert applied.hvg.n_top_genes == 1500
    assert applied.hvg.flavor == "seurat_v3"
    assert applied.graph.n_neighbors == 20
    assert applied.integration.method == "scanorama"
    assert applied.integration.batch_key == "batch_col"


@pytest.mark.unit
def test_run_intelligent_preprocessing_stores_trace(monkeypatch, minimal_adata):
    import scLucid.preprocess.intelligent.recommender as recommender_module

    strategy = PreprocessingStrategy(
        data_profile=DataProfile(
            n_cells=100,
            n_genes=2000,
            sparsity=0.9,
            median_counts_per_cell=1500,
            median_genes_per_cell=700,
            is_sparse=True,
            is_small_dataset=True,
            is_medium_dataset=False,
            is_large_dataset=False,
        ),
        hvg=HVGRecommendation(1000, 0.7, 800, 1200, "test", 0.8),
        pca=PCARecommendation(20, 0.7, 15, 25, "test", 0.8),
        neighbors=NeighborsRecommendation(15, 20, 0.3, 10, 20, 15, 25, "test", 0.8),
        resolution=ResolutionRecommendation(1.0, 5, 0.6, 0.8, 1.2, "test", 0.8),
    )

    def fake_recommend(*args, **kwargs):
        return strategy

    def fake_run(adata, config=None, results_dir=None):
        result = adata.copy()
        result.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["workflow_config"] = (
            config.to_dict() if hasattr(config, "to_dict") else {}
        )
        return result

    monkeypatch.setattr(recommender_module, "recommend_intelligent_preprocessing", fake_recommend)
    monkeypatch.setattr("scLucid.preprocess.workflow.run_preprocessing", fake_run)

    result, returned_strategy = run_intelligent_preprocessing(
        minimal_adata.copy(),
        batch_key="sampleID",
        apply_recommendations=True,
    )

    assert returned_strategy is strategy
    intelligent_meta = result.uns["sclucid"]["preprocess"]["intelligent_recommendation"]
    assert intelligent_meta["batch_key"] == "sampleID"
    assert intelligent_meta["strategy"]["hvg"]["n_top_genes"] == 1000
    assert intelligent_meta["applied_config"]["hvg"]["n_top_genes"] == 1000
