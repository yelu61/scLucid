"""Unit tests for scLucid.preprocess.integrate."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.preprocess.config import IntegrationConfig
from scLucid.preprocess.integrate import batch_correction, evaluate_integration


@pytest.mark.unit
class TestBatchCorrection:
    """Tests for batch_correction function."""

    def test_no_method_returns_adata_unchanged(self, minimal_adata):
        adata = minimal_adata.copy()
        result = batch_correction(
            adata,
            config=IntegrationConfig(method=None, plot=False, report=False, verbose=False),
        )
        assert result is adata or isinstance(result, AnnData)

    def test_harmony_integration_mock(self, monkeypatch, minimal_adata):
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

        assert "X_harmony" in result.obsm
        integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
        assert "harmony" in integration_meta
        assert integration_meta["workflow"]["method"] == "harmony"
        assert integration_meta["workflow"]["output_key"] == "X_harmony"

    def test_scanorama_integration_mock(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        def fake_scanorama(adata, batch_key, embedding_key, **kwargs):
            adata.obsm[embedding_key] = adata.obsm["X_pca"].copy()
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
                "integration", {}
            )["scanorama"] = {"batch_key": batch_key}
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_scanorama", fake_scanorama)

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="scanorama",
                batch_key="batch",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
        )

        assert "X_scanorama" in result.obsm
        integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
        assert "scanorama" in integration_meta

    def test_combat_integration(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.layers["normalized"] = adata.X.copy()

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="combat",
                batch_key="batch",
                plot=False,
                report=False,
                verbose=False,
            ),
        )

        assert (
            "X_combat" in result.obsm
            or "combat" in result.layers
            or "combat"
            in result.uns.get("sclucid", {}).get("preprocess", {}).get("integration", {})
        )

    def test_no_method_returns_early_without_error(self, minimal_adata):
        adata = minimal_adata.copy()
        result = batch_correction(
            adata,
            config=IntegrationConfig(method=None, plot=False, report=False, verbose=False),
        )
        # When method is None, function returns early without storing metadata
        assert isinstance(result, AnnData)
        assert result.n_obs == minimal_adata.n_obs

    def test_config_not_mutated(self, minimal_adata):
        adata = minimal_adata.copy()
        config = IntegrationConfig(method=None, plot=False, report=False, verbose=False)
        original_dict = config.to_dict()
        batch_correction(adata, config=config)
        assert config.to_dict() == original_dict


@pytest.mark.unit
class TestEvaluateIntegration:
    """Tests for evaluate_integration function."""

    def test_evaluate_integration_basic(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        result = evaluate_integration(
            adata,
            batch_key="batch",
            use_rep="X_pca",
            methods=["silhouette"],
        )

        assert isinstance(result, dict)
        assert "method" in result
        assert "n_batches" in result
        assert result["n_batches"] == 2

    def test_evaluate_integration_missing_batch_key(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        with pytest.raises(ValueError, match="batch"):
            evaluate_integration(adata, batch_key="missing", use_rep="X_pca")

    def test_evaluate_integration_explicit_missing_rep_raises(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)

        # When an explicit use_rep does not exist, the function raises KeyError
        with pytest.raises(KeyError):
            evaluate_integration(adata, batch_key="batch", use_rep="X_missing")
