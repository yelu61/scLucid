"""Tests for the batch-correction diagnostic benchmark runner."""

from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

from validation.preprocess_analysis.run_batch_correction_diagnostic_benchmark import (
    _run_method_comparison,
)


@pytest.fixture
def tiny_batch_adata():
    """Small synthetic data with batch labels and a PCA embedding."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 50
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(float)
    adata = AnnData(X)
    adata.obs["batch"] = ["b1"] * 30 + ["b2"] * 30
    adata.obs["cell_type"] = ["a"] * 15 + ["b"] * 15 + ["a"] * 15 + ["b"] * 15
    adata.layers["counts"] = X.copy()
    # Add minimal QC metrics needed by the benchmark path.
    adata.obs["total_counts"] = X.sum(axis=1)
    adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    adata.obs["pct_counts_mt"] = rng.uniform(0, 10, size=n_cells)
    # Add a PCA embedding so diagnose_integration_risk has before/after reps.
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10))
    return adata


class TestBatchCorrectionDiagnosticBenchmark:
    def test_run_method_comparison_records_production_risk(self, tmp_path, tiny_batch_adata):
        result = _run_method_comparison(
            tiny_batch_adata,
            dataset="demo",
            batch_key="batch",
            label_key="cell_type",
            method="no_correction",
            max_epochs=2,
        )
        assert result["method_status"] == "ok"
        assert result["production_risk_level"] in {"low", "moderate", "high", "unknown"}
        assert "production_risk_score" in result
        assert "production_warnings" in result

    def test_run_method_comparison_rejects_unknown_method(self, tmp_path, tiny_batch_adata):
        result = _run_method_comparison(
            tiny_batch_adata,
            dataset="demo",
            batch_key="batch",
            label_key="cell_type",
            method="unknown_method",
            max_epochs=2,
        )
        assert result["method_status"] == "failed"
        assert result["production_risk_level"] == "unknown"

    def test_harmony_uses_deeper_params(self, tmp_path, tiny_batch_adata, monkeypatch):
        import scLucid.preprocess.integrate as integrate_module

        captured = {}

        def fake_harmony(adata, covariate_keys, basis, embedding_key, **kwargs):
            captured["kwargs"] = kwargs
            adata.obsm[embedding_key] = adata.obsm[basis].copy()
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_harmony", fake_harmony)

        _run_method_comparison(
            tiny_batch_adata,
            dataset="demo",
            batch_key="batch",
            label_key="cell_type",
            method="harmony",
            max_epochs=2,
        )
        assert captured["kwargs"]["max_iter_harmony"] == 50
        assert captured["kwargs"]["theta"] == 2.0
        assert captured["kwargs"]["lambda_val"] == 1.0

    def test_scvi_model_save_path_is_recorded(self, tmp_path, tiny_batch_adata, monkeypatch):
        import scLucid.preprocess.integrate as integrate_module

        def fake_scvi(adata, batch_key, embedding_key, **kwargs):
            adata.obsm[embedding_key] = np.random.default_rng(0).normal(
                size=(adata.n_obs, kwargs.get("n_latent", 15))
            )
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
                "integration", {}
            )["scvi"] = {
                "model_saved": kwargs.get("save_model", False),
                "model_path": kwargs.get("model_path", ""),
            }
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_scvi", fake_scvi)

        model_save_dir = tmp_path / "models"
        result = _run_method_comparison(
            tiny_batch_adata,
            dataset="demo",
            batch_key="batch",
            label_key="cell_type",
            method="scvi",
            max_epochs=2,
            model_save_dir=model_save_dir,
        )
        assert result["method_status"] == "ok"
        assert result["scvi_model_saved"] is True
        assert result["scvi_model_path"]
        assert Path(result["scvi_model_path"]).parent == model_save_dir
