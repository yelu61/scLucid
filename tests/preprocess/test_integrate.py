"""Unit tests for scLucid.preprocess.integrate."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.preprocess.config import IntegrationConfig
from scLucid.preprocess.integrate import (
    batch_correction,
    decide_integration,
    diagnose_integration_risk,
    evaluate_integration,
)


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

    def test_diagnose_integration_risk_flags_batch_condition_confounding(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["b1"] * (adata.n_obs // 2) + ["b2"] * (
            adata.n_obs - adata.n_obs // 2
        )
        adata.obs["condition"] = ["ctrl"] * (adata.n_obs // 2) + ["treat"] * (
            adata.n_obs - adata.n_obs // 2
        )

        result = diagnose_integration_risk(
            adata,
            batch_key="batch",
            condition_key="condition",
            tumor=True,
        )

        assert result["risk_level"] in {"moderate", "high"}
        assert result["metrics"]["batch_condition_cramers_v"] >= 0.8
        assert result["warnings"]
        assert "integration_risk" in adata.uns["sclucid"]["preprocess"]["integration"]

    def test_diagnose_integration_risk_detects_biology_confounding(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["b1", "b2"] * (adata.n_obs // 2)
        adata.obs["group"] = ["g1", "g2"] * (adata.n_obs // 2)

        result = diagnose_integration_risk(
            adata,
            batch_key="batch",
            biology_columns=["group"],
        )

        assert result["risk_level"] in {"moderate", "high"}
        assert any("one-to-one" in w for w in result["warnings"])

    def test_decide_integration_auto_skips_when_confounded(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["b1", "b2"] * (adata.n_obs // 2)
        adata.obs["group"] = ["g1", "g2"] * (adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        run, warnings, risk = decide_integration(
            adata,
            batch_key="batch",
            run_integration="auto",
            biology_columns=["group"],
        )

        assert run is False
        assert any("one-to-one" in w for w in warnings)

    def test_decide_integration_forced_true(self, minimal_adata):
        adata = minimal_adata.copy()
        run, warnings, risk = decide_integration(
            adata,
            batch_key="batch",
            run_integration=True,
        )
        assert run is True
        assert warnings == []
        assert risk is None

    def test_decide_integration_forced_false(self, minimal_adata):
        adata = minimal_adata.copy()
        run, warnings, risk = decide_integration(
            adata,
            batch_key="batch",
            run_integration=False,
        )
        assert run is False
        assert "disabled by user" in warnings[0]
        assert risk is None

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

    def test_auto_decide_success_preserves_decision_trace(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        def fake_harmony(adata, covariate_keys, basis, embedding_key, **kwargs):
            adata.obsm[embedding_key] = adata.obsm[basis].copy()
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_harmony", fake_harmony)

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="harmony",
                batch_key="batch",
                use_rep="X_pca",
                auto_decide=True,
                plot=False,
                report=False,
                verbose=False,
            ),
        )

        workflow = result.uns["sclucid"]["preprocess"]["integration"]["workflow"]
        assert workflow["auto_decide"] is True
        assert workflow["decision"]["run"] is True
        assert workflow["decision"]["risk"]["risk_level"] == "low"
        assert workflow["risk"] == workflow["decision"]["risk"]

    def test_harmony_low_level_default_max_iter_matches_config(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))
        seen = {}

        class FakeHarmony:
            Z_corr = adata.obsm["X_pca"].T
            objective_history = [2.0, 1.0]

        def fake_run_harmony(**kwargs):
            seen["max_iter_harmony"] = kwargs["max_iter_harmony"]
            return FakeHarmony()

        monkeypatch.setattr(integrate_module.hm, "run_harmony", fake_run_harmony)

        integrate_module._integrate_harmony(
            adata,
            covariate_keys="batch",
            basis="X_pca",
            check_convergence=False,
        )

        assert seen["max_iter_harmony"] == 50

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

    def test_scvi_integration_mock(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        def fake_scvi(adata, batch_key, embedding_key, **kwargs):
            adata.obsm[embedding_key] = adata.obsm["X_pca"].copy()
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
                "integration", {}
            )["scvi"] = {"batch_key": batch_key}
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_scvi", fake_scvi)

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="scvi",
                batch_key="batch",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert "X_scvi" in result.obsm
        integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
        assert integration_meta["workflow"]["method"] == "scvi"

    def test_scanvi_method_is_case_normalized(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obs["cell_label"] = ["T"] * adata.n_obs
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        def fake_scanvi(adata, batch_key, labels_key, embedding_key, **kwargs):
            adata.obsm[embedding_key] = adata.obsm["X_pca"].copy()
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
                "integration", {}
            )["scanvi"] = {"batch_key": batch_key, "labels_key": labels_key}
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_scanvi", fake_scanvi)

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="scANVI",
                batch_key="batch",
                scanvi_labels_key="cell_label",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
        )

        assert "X_scanvi" in result.obsm
        integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
        assert integration_meta["workflow"]["method"] == "scanvi"

    def test_bbknn_integration_mock(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        def fake_bbknn(adata, batch_key, use_rep, **kwargs):
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
                "integration", {}
            )["bbknn"] = {"batch_key": batch_key}
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_bbknn", fake_bbknn)

        result = batch_correction(
            adata,
            config=IntegrationConfig(
                method="bbknn",
                batch_key="batch",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        integration_meta = result.uns["sclucid"]["preprocess"]["integration"]
        assert integration_meta["workflow"]["method"] == "bbknn"

    def test_unknown_method_raises_validation_error(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        # Pydantic's Literal type catches unknown methods at config validation time
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="method"):
            IntegrationConfig(
                method="nonexistent_method",
                batch_key="batch",
                plot=False,
                report=False,
                verbose=False,
            )

    def test_missing_batch_key_raises_valueerror(self, minimal_adata):
        adata = minimal_adata.copy()
        with pytest.raises(ValueError, match="not found in adata.obs"):
            batch_correction(
                adata,
                config=IntegrationConfig(
                    method="harmony",
                    batch_key="missing_batch_column",
                    use_rep="X_pca",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )

    def test_force_rerun_when_output_exists(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))
        adata.obsm["X_harmony"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        call_count = [0]

        def fake_harmony(adata, covariate_keys, basis, embedding_key, **kwargs):
            call_count[0] += 1
            adata.obsm[embedding_key] = adata.obsm[basis].copy()
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_harmony", fake_harmony)

        # Without force: should return early (harmony already exists)
        batch_correction(
            adata,
            config=IntegrationConfig(
                method="harmony",
                batch_key="batch",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
            force=False,
        )
        assert call_count[0] == 0  # never called

        # With force: should re-run
        batch_correction(
            adata,
            config=IntegrationConfig(
                method="harmony",
                batch_key="batch",
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
            force=True,
        )
        assert call_count[0] == 1  # called once

    def test_batch_key_list_for_harmony(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch1"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obs["batch2"] = "c"
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 5))

        received_keys = []

        def fake_harmony(adata, covariate_keys, basis, embedding_key, **kwargs):
            received_keys.append(list(covariate_keys) if isinstance(covariate_keys, list) else covariate_keys)
            adata.obsm[embedding_key] = adata.obsm[basis].copy()
            return adata

        monkeypatch.setattr(integrate_module, "_integrate_harmony", fake_harmony)

        batch_correction(
            adata,
            config=IntegrationConfig(
                method="harmony",
                batch_key=["batch1", "batch2"],
                use_rep="X_pca",
                plot=False,
                report=False,
                verbose=False,
            ),
        )
        assert received_keys[0] == ["batch1", "batch2"]


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
        assert result["interpretation"]["status"] in {"acceptable", "review"}
        assert "evaluation" in adata.uns["sclucid"]["preprocess"]["integration"]

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
