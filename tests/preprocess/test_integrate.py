"""Unit tests for scLucid.preprocess.integrate."""

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from scLucid.preprocess.config import IntegrationConfig
from scLucid.preprocess.integrate import (
    _compute_kbet_reference,
    _compute_kbet_score,
    _compute_local_batch_chi2_acceptance,
    batch_correction,
    decide_integration,
    diagnose_integration_risk,
    evaluate_integration,
)


@pytest.mark.unit
class TestBatchCorrection:
    """Tests for batch_correction function."""

    def test_kbet_metadata_marks_chi_square_approximation(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(40, 5))
        batch_labels = pytest.importorskip("pandas").Series(["a", "b"] * 20)

        result = _compute_local_batch_chi2_acceptance(
            x, batch_labels, n_neighbors=10, alpha=0.05, n_sample_cells=40
        )

        assert result["n_neighbors"] == 10
        assert result["alpha"] == 0.05
        assert result["model_type"] == "chi_square_local_batch_approximation"
        assert result["claim_level"] == "batch_mixing_diagnostic_heuristic"
        assert "approximated" in result["review_note"]
        assert "kbet_score" not in result

    def test_kbet_deprecated_alias_warns(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(40, 5))
        batch_labels = pytest.importorskip("pandas").Series(["a", "b"] * 20)

        with pytest.warns(FutureWarning, match="deprecated"):
            result = _compute_kbet_score(
                x, batch_labels, n_neighbors=10, alpha=0.05, n_sample_cells=40
            )

        assert "acceptance_rate" in result
        assert "rejection_rate" in result

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
        assert "sclucid" not in adata.uns or "preprocess" not in adata.uns.get("sclucid", {})

        recorded = diagnose_integration_risk(
            adata,
            batch_key="batch",
            condition_key="condition",
            tumor=True,
            record=True,
        )
        assert recorded["risk_level"] in {"moderate", "high"}
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

    def test_decide_integration_consumes_precomputed_risk(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["b1", "b2"] * (adata.n_obs // 2)

        def fail_if_called(*args, **kwargs):
            raise AssertionError("diagnose_integration_risk should not be called")

        monkeypatch.setattr(integrate_module, "diagnose_integration_risk", fail_if_called)

        risk = {"risk_level": "low", "warnings": [], "metrics": {}}
        run, warnings, returned_risk = decide_integration(
            adata,
            batch_key="batch",
            run_integration="auto",
            risk=risk,
        )

        assert run is True
        assert warnings == []
        assert returned_risk is risk

    def test_decide_integration_precomputed_risk_can_skip(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["b1", "b2"] * (adata.n_obs // 2)

        risk = {
            "risk_level": "moderate",
            "warnings": ["batch and biology are partially confounded"],
            "metrics": {},
        }
        run, warnings, returned_risk = decide_integration(
            adata,
            batch_key="batch",
            run_integration="auto",
            risk=risk,
        )

        assert run is False
        assert returned_risk is risk
        assert "partially confounded" in warnings[0]
        assert any("low_risk_only" in w for w in warnings)

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

    def test_combat_warns_on_sparse_densification(self, minimal_adata, caplog):
        import logging

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.X = scipy.sparse.csr_matrix(adata.X.copy())

        with caplog.at_level(logging.WARNING):
            batch_correction(
                adata,
                config=IntegrationConfig(
                    method="combat",
                    batch_key="batch",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )

        messages = [record.message.lower() for record in caplog.records]
        assert any("densifying" in m and ("harmony" in m or "scvi" in m) for m in messages)

    def test_combat_raises_memory_error_on_large_sparse_without_force_dense(self):
        import scipy.sparse

        # Create a sparse matrix above the 50M threshold.
        large_sparse = scipy.sparse.random(10000, 6000, density=0.01, format="csr", random_state=0)
        adata = AnnData(X=large_sparse)
        adata.obs["batch"] = ["a"] * 5000 + ["b"] * 5000

        with pytest.raises(MemoryError):
            batch_correction(
                adata,
                config=IntegrationConfig(
                    method="combat",
                    batch_key="batch",
                    plot=False,
                    report=False,
                    verbose=False,
                ),
            )

    def test_combat_allows_large_sparse_with_force_dense(self, monkeypatch):
        import scipy.sparse

        import scLucid.preprocess.integrate as integrate_module

        large_sparse = scipy.sparse.random(10000, 6000, density=0.01, format="csr", random_state=0)
        adata = AnnData(X=large_sparse)
        adata.obs["batch"] = ["a"] * 5000 + ["b"] * 5000

        combat_called = [False]

        def fake_combat(temp_adata, *args, **kwargs):
            combat_called[0] = True
            temp_adata.X = (
                temp_adata.X.toarray() if hasattr(temp_adata.X, "toarray") else temp_adata.X
            )

        monkeypatch.setattr(integrate_module.sc.pp, "combat", fake_combat)

        # Without force_dense, a large sparse matrix raises MemoryError.
        with pytest.raises(MemoryError):
            integrate_module._integrate_combat(adata, batch_key="batch", force_dense=False)

        # With force_dense=True, the call proceeds past the guard.
        result = integrate_module._integrate_combat(
            adata, batch_key="batch", force_dense=True
        )

        assert combat_called[0] is True
        assert "combat" in result.uns["sclucid"]["preprocess"]["integration"]

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

    def test_evaluate_integration_exposes_local_batch_chi2_semantics(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        result = evaluate_integration(
            adata,
            batch_key="batch",
            use_rep="X_pca",
            n_neighbors=10,
            methods=["kbet"],
            plot=False,
        )

        assert result["local_batch_chi2_model_type"] == "chi_square_local_batch_approximation"
        assert result["local_batch_chi2_claim_level"] == "batch_mixing_diagnostic_heuristic"
        assert "approximated" in result["local_batch_chi2_review_note"]
        stored = adata.uns["sclucid"]["preprocess"]["integration"]["evaluation"]
        assert stored["local_batch_chi2_model_type"] == "chi_square_local_batch_approximation"

    def test_evaluate_integration_backward_compatible_kbet_aliases(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        with pytest.warns(FutureWarning, match="deprecated"):
            result = evaluate_integration(
                adata,
                batch_key="batch",
                use_rep="X_pca",
                n_neighbors=10,
                methods=["kbet"],
                plot=False,
            )

        assert result["kbet_acceptance"] == result["local_batch_chi2_acceptance"]
        assert result["kbet_rejection_rate"] == result["local_batch_chi2_rejection_rate"]

    def test_evaluate_integration_reference_kbet_fallback(self, monkeypatch, minimal_adata):
        import scLucid.preprocess.integrate as integrate_module

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        def fake_kbet_reference(*args, **kwargs):
            return {
                "acceptance_rate": 0.85,
                "rejection_rate": 0.15,
                "model_type": "scib_kbet_reference",
                "claim_level": "batch_mixing_diagnostic_reference",
                "review_note": "Reference kBET used.",
                "n_neighbors": kwargs.get("n_neighbors", 25),
                "reference_backend": "scib",
            }

        monkeypatch.setattr(integrate_module, "_compute_kbet_reference", fake_kbet_reference)

        result = evaluate_integration(
            adata,
            batch_key="batch",
            use_rep="X_pca",
            n_neighbors=10,
            methods=["kbet"],
            plot=False,
            use_reference_kbet=True,
        )

        assert result["local_batch_chi2_acceptance"] == 0.85
        assert result["local_batch_chi2_rejection_rate"] == 0.15
        assert result["local_batch_chi2_reference_backend"] == "scib"

    def test_compute_kbet_reference_fallback_when_scib_missing(self, monkeypatch, minimal_adata):
        import builtins

        adata = minimal_adata.copy()
        adata.obs["batch"] = ["a"] * (adata.n_obs // 2) + ["b"] * (adata.n_obs - adata.n_obs // 2)
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        real_import = builtins.__import__

        def fail_scib_import(name, *args, **kwargs):
            if name == "scib":
                raise ImportError("No module named 'scib'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fail_scib_import)

        result = _compute_kbet_reference(
            adata,
            batch_key="batch",
            use_rep="X_pca",
            n_neighbors=10,
            n_sample_cells=adata.n_obs,
        )

        assert result["reference_backend"] == "chi_square_approximation"
        assert "acceptance_rate" in result
        assert "rejection_rate" in result

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
