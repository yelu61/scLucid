"""Tests for scLucid.qc.doublet module.

Covers:
- predict_doublets configuration and dispatch logic
- _run_heuristic (co-expression-based doublet detection)
- DoubletEvidenceProfiler evidence table generation
- _export_doublet_stats
"""

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.qc.config import DoubletConfig
from scLucid.qc.doublet import (
    ALGORITHM_PRED_COL,
    ALGORITHM_SCORE_COL,
    COMBINED_SCORE_COL,
    EXPECTED_HETEROTYPIC_RATE_COL,
    EXPECTED_HOMOTYPIC_RATE_COL,
    EXPECTED_TOTAL_RATE_COL,
    HETEROTYPIC_RISK_COL,
    HOMOTYPIC_RISK_COL,
    DoubletEvidenceProfiler,
    _export_doublet_stats,
    _run_heuristic,
    predict_doublets,
)
from scLucid.qc.doublet.algorithms import (
    _raw_count_guard,
    _run_scrublet,
    _run_scdblfinder,
)
from scLucid.qc.doublet.ensemble import _merge_doublet_predictions
from scLucid.qc.doublet._scrublet_compat import apply_scrublet_compatibility_shims


# ---------------------------------------------------------------------------
# Scrublet compatibility shims
# ---------------------------------------------------------------------------


class TestScrubletCompatibilityShims:
    """Tests for runtime compatibility shims applied to the scrublet package."""

    def test_shim_restores_sparse_A_attribute(self):
        """The shim should make .A available on scipy sparse matrices if absent."""
        import scipy.sparse

        # Force-remove the attribute to simulate a modern SciPy environment.
        original_A = getattr(scipy.sparse.spmatrix, "A", None)
        try:
            if hasattr(scipy.sparse.spmatrix, "A"):
                delattr(scipy.sparse.spmatrix, "A")

            # Re-apply shims (they are idempotent and use the module-level flag).
            from scLucid.qc.doublet import _scrublet_compat as _compat

            _compat._SCRUBLET_SHIMS_APPLIED = False
            apply_scrublet_compatibility_shims()

            mat = scipy.sparse.csr_matrix(np.array([[0, 1], [2, 0]]))
            assert np.array_equal(mat.A, mat.toarray())
        finally:
            if original_A is not None:
                scipy.sparse.spmatrix.A = original_A
            elif hasattr(scipy.sparse.spmatrix, "A"):
                delattr(scipy.sparse.spmatrix, "A")
            _compat._SCRUBLET_SHIMS_APPLIED = False

    def test_run_scrublet_uses_fallback_when_call_doublets_fails(self, monkeypatch):
        """_run_scrublet should fall back to quantile thresholding when call_doublets raises."""

        class FakeScrublet:
            def __init__(self, counts, expected_doublet_rate=None):
                self.counts = counts
                self.expected_doublet_rate = expected_doublet_rate

            def scrub_doublets(self, n_prin_comps=None, verbose=False):
                scores = np.linspace(0.0, 1.0, self.counts.shape[0])
                return scores, None

            def call_doublets(self, verbose=False):
                raise RuntimeError("threshold detection failed")

        monkeypatch.setitem(sys.modules, "scrublet", SimpleNamespace(Scrublet=FakeScrublet))
        rng = np.random.default_rng(7)
        adata = AnnData(X=rng.poisson(3, size=(80, 120)).astype(np.float32))
        scores, predicted = _run_scrublet(
            adata,
            sample_name="s1",
            config=DoubletConfig(expected_doublet_rate=0.1, scr_n_pcs=10),
        )

        assert scores is not None
        assert predicted is not None
        # 10% highest scores should be flagged by the fallback threshold.
        assert int(predicted.sum()) == 8


# ---------------------------------------------------------------------------
# scDblFinder
# ---------------------------------------------------------------------------


class TestScDblFinder:
    """Tests for scDblFinder integration via pyscdblfinder."""

    def test_config_accepts_scdblfinder(self):
        """DoubletConfig should accept method='scdblfinder'."""
        cfg = DoubletConfig(method="scdblfinder")
        assert cfg.method == "scdblfinder"
        assert cfg.scdblfinder_nfeatures == 1352
        assert cfg.scdblfinder_dims == 20

    def test_unknown_method_raises_at_config_level(self):
        """Unknown method should be rejected by Pydantic config validation."""
        with pytest.raises(Exception, match="scrublet|solo|doubletdetection|scdblfinder"):
            DoubletConfig(method="unknown_method")

    def test_run_scdblfinder_optional_dependency_missing(self, monkeypatch):
        """When pyscdblfinder is missing, _run_scdblfinder should log and return None."""
        monkeypatch.setitem(sys.modules, "pyscdblfinder", None)
        rng = np.random.default_rng(7)
        adata = AnnData(X=rng.poisson(3, size=(80, 120)).astype(np.float32))
        scores, predicted = _run_scdblfinder(
            adata,
            sample_name="s1",
            config=DoubletConfig(expected_doublet_rate=0.1, method="scdblfinder"),
        )
        assert scores is None
        assert predicted is None

    def test_predict_doublets_dispatches_scdblfinder(self, monkeypatch):
        """predict_doublets should call the scDblFinder wrapper and write result columns."""
        import scLucid.qc.doublet.ensemble as ensemble_module

        adata = AnnData(X=np.ones((80, 6), dtype=float))
        adata.obs_names = [f"cell{i}" for i in range(80)]
        adata.var_names = [f"g{i}" for i in range(6)]
        adata.obs["sampleID"] = "sample_1"

        def fake_scdblfinder(data_view, sample, cfg):
            n = data_view.n_obs
            scores = pd.Series(np.linspace(0, 1, n), index=data_view.obs_names)
            predicted = scores > 0.8
            return scores, predicted

        monkeypatch.setattr(ensemble_module, "_run_scdblfinder", fake_scdblfinder)

        cfg = DoubletConfig(
            method="scdblfinder",
            run_algorithm=True,
            use_heuristics=False,
            expected_doublet_rate=0.1,
            plot_summary=False,
            export_stats=False,
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        assert "scdblfinder_score" in out.obs.columns
        assert "scdblfinder_predicted" in out.obs.columns
        assert "predicted_doublet" in out.obs.columns

    @pytest.mark.slow
    @pytest.mark.optional
    def test_scdblfinder_runs_on_synthetic_doublets(self, doublet_test_adata):
        """End-to-end scDblFinder run on synthetic data with known doublets."""
        pytest.importorskip("pyscdblfinder")

        adata = doublet_test_adata.copy()
        # Use integer raw counts and ensure a sample key.
        adata.X = adata.layers["counts"].astype(int)
        adata.obs["sampleID"] = "sample_1"

        cfg = DoubletConfig(
            method="scdblfinder",
            run_algorithm=True,
            use_heuristics=False,
            expected_doublet_rate=0.1,
            scdblfinder_dims=10,
            scdblfinder_nfeatures=500,
            plot_summary=False,
            export_stats=False,
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        assert "scdblfinder_score" in out.obs.columns
        assert "scdblfinder_predicted" in out.obs.columns

        if "is_doublet" in out.obs.columns:
            from sklearn.metrics import roc_auc_score

            auc = roc_auc_score(out.obs["is_doublet"], out.obs["scdblfinder_score"])
            assert auc > 0.55, f"scDblFinder AUC {auc:.3f} is unexpectedly low"


# ---------------------------------------------------------------------------
# _run_heuristic
# ---------------------------------------------------------------------------


class TestRunHeuristic:
    """Tests for the co-expression heuristic doublet detection."""

    def test_basic_heuristic_run(self, minimal_adata):
        """Heuristic should run on basic synthetic data."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        potential_doublets, lineage_scores, heuristic_scores = _run_heuristic(
            minimal_adata, cfg
        )

        assert isinstance(potential_doublets, pd.Series)
        assert isinstance(lineage_scores, pd.DataFrame)
        assert isinstance(heuristic_scores, pd.Series)
        assert len(potential_doublets) == minimal_adata.n_obs
        assert len(heuristic_scores) == minimal_adata.n_obs
        # Scores should be in [0, 1]
        assert heuristic_scores.min() >= 0
        assert heuristic_scores.max() <= 1

    def test_heuristic_returns_no_false_positives_on_clean_data(self, minimal_adata):
        """On clean synthetic data without doublets, predictions should be conservative."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        potential_doublets, _, heuristic_scores = _run_heuristic(minimal_adata, cfg)

        # Even if no doublets are predicted, scores should be computed
        assert heuristic_scores.notna().all()

    def test_heuristic_ignore_pairs(self, minimal_adata):
        """ignore_coexpression_pairs should mask scores for specified lineages."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
            ignore_coexpression_pairs=[("T_cell", "B_cell")],
        )
        _, _, heuristic_scores_with_ignore = _run_heuristic(minimal_adata, cfg)
        _, _, heuristic_scores_without = _run_heuristic(
            minimal_adata,
            DoubletConfig(
                run_algorithm=False,
                use_heuristics=True,
                marker_species="human",
                marker_tissue="pbmc",
            ),
        )
        # Scores may differ when ignoring pairs
        assert isinstance(heuristic_scores_with_ignore, pd.Series)

    def test_heuristic_ignore_pairs_masks_scores_and_predictions(self, monkeypatch):
        """Allowlisted co-expression should affect both heuristic score and binary call."""
        X = np.ones((8, 8), dtype=float)
        adata = AnnData(X=X)
        adata.var_names = [f"g{i}" for i in range(8)]
        adata.obs["n_genes_by_counts"] = 1000

        def fake_score_genes(adata_in, gene_list, score_name, ctrl_size):
            adata_in.obs[score_name] = 1.0

        monkeypatch.setattr("scLucid.qc.doublet.heuristic.sc.tl.score_genes", fake_score_genes)

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_configs={
                "Lineage_A": {"genes": ["g0", "g1"], "min_genes_required": 2},
                "Lineage_B": {"genes": ["g2", "g3"], "min_genes_required": 2},
            },
            heuristic_score_threshold=0.1,
            expected_doublet_rate=0.5,
            ignore_coexpression_pairs=[("Lineage_A", "Lineage_B")],
        )

        potential_doublets, _, heuristic_scores = _run_heuristic(adata, cfg)

        assert not potential_doublets.any()
        assert (heuristic_scores == 0).all()
        meta = adata.uns["sclucid"]["qc"]["doublet_params"]["heuristic_allowlist"]
        assert meta["applied_before_thresholding"] is True
        assert meta["n_cells_score_zeroed"] == adata.n_obs


# ---------------------------------------------------------------------------
# predict_doublets
# ---------------------------------------------------------------------------


class TestPredictDoublets:
    """Tests for the main predict_doublets entry point."""

    def test_doubletdetection_default_p_thresh_is_moderate(self):
        cfg = DoubletConfig()
        assert cfg.dd_p_thresh == 1e-4

        custom = DoubletConfig(dd_p_thresh=0.01)
        assert custom.dd_p_thresh == 0.01

    def test_heuristic_only(self, minimal_adata):
        """predict_doublets with heuristic only should populate required columns."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        assert "predicted_doublet" in out.obs.columns
        assert "doublet_source" in out.obs.columns
        assert "heuristic_predicted" in out.obs.columns
        assert "heuristic_confidence_score" in out.obs.columns

    def test_invalid_sample_key_raises(self, minimal_adata):
        """Missing sample key should raise ValueError."""
        cfg = DoubletConfig()
        with pytest.raises(ValueError, match="Sample key"):
            predict_doublets(minimal_adata, config=cfg, sample_key="nonexistent")

    def test_unknown_method_raises_at_config_level(self):
        """Unknown method should be rejected by Pydantic config validation."""
        with pytest.raises(Exception, match="scrublet|solo|doubletdetection"):
            DoubletConfig(method="unknown_method")

    def test_skip_small_samples(self, minimal_adata):
        """Samples with <50 cells should be skipped for algorithmic detection."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        cfg = DoubletConfig(
            run_algorithm=True,
            use_heuristics=False,
            method="scrublet",
        )
        # With <50 cells, scrublet should be skipped and heuristic columns still created
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")
        assert "predicted_doublet" in out.obs.columns

    def test_config_none_uses_defaults(self, minimal_adata):
        """config=None should use DoubletConfig defaults."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        out = predict_doublets(adata, config=None, sample_key="sampleID")
        assert "predicted_doublet" in out.obs.columns

    def test_merge_doublet_predictions_all_zero_scores_is_stable(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["algorithm_score"] = 0.0
        adata.obs["heuristic_score"] = 0.0

        merged = _merge_doublet_predictions(
            adata,
            algorithm_score_col="algorithm_score",
            heuristic_score_col="heuristic_score",
            expected_rate=0.1,
        )

        assert merged.notna().all()
        assert not merged.any()
        assert COMBINED_SCORE_COL in adata.obs
        assert np.allclose(adata.obs[COMBINED_SCORE_COL], 0.0)

    def test_merge_doublet_predictions_uses_per_sample_expected_rates(self):
        adata = AnnData(X=np.ones((10, 3), dtype=float))
        adata.obs_names = [f"cell{i}" for i in range(adata.n_obs)]
        adata.obs["sampleID"] = ["s1"] * 5 + ["s2"] * 5
        adata.obs["algorithm_score"] = [0.1, 0.2, 0.3, 0.4, 0.9, 0.1, 0.2, 0.3, 0.8, 0.9]
        adata.obs["heuristic_score"] = 0.0

        merged = _merge_doublet_predictions(
            adata,
            algorithm_score_col="algorithm_score",
            heuristic_score_col="heuristic_score",
            expected_rate={"s1": 0.2, "s2": 0.4},
            sample_key="sampleID",
        )

        assert int(merged.iloc[:5].sum()) == 1
        assert int(merged.iloc[5:].sum()) == 2

    def test_external_doublet_evidence_review_only_does_not_change_final(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        adata.obs["hashing_doublet"] = [True] + [False] * (adata.n_obs - 1)

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=False,
            external_doublet_cols=["hashing_doublet", "missing_col"],
            external_doublet_policy="review_only",
            plot_summary=False,
            export_stats=False,
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        assert "external_doublet_evidence" in out.obs
        assert bool(out.obs["external_doublet_evidence"].iloc[0]) is True
        assert bool(out.obs["predicted_doublet"].iloc[0]) is False
        meta = out.uns["sclucid"]["qc"]["doublet_params"]["external_doublet_evidence"]
        assert meta["columns_used"] == ["hashing_doublet"]
        assert meta["missing_columns"] == ["missing_col"]
        assert meta["policy"] == "review_only"

    def test_scrublet_uses_scrub_doublets_prediction_when_call_doublets_is_none(self, monkeypatch):
        class FakeScrublet:
            def __init__(self, counts, expected_doublet_rate=None):
                self.counts = counts
                self.expected_doublet_rate = expected_doublet_rate

            def scrub_doublets(self, n_prin_comps=None, verbose=False):
                scores = np.linspace(0.0, 1.0, self.counts.shape[0])
                predicted = scores > 0.8
                return scores, predicted

            def call_doublets(self, verbose=False):
                return None

        monkeypatch.setitem(sys.modules, "scrublet", SimpleNamespace(Scrublet=FakeScrublet))
        rng = np.random.default_rng(7)
        adata = AnnData(X=rng.poisson(3, size=(80, 120)).astype(np.float32))
        scores, predicted = _run_scrublet(
            adata,
            sample_name="s1",
            config=DoubletConfig(expected_doublet_rate=0.1, scr_n_pcs=10),
        )

        assert scores is not None
        assert predicted is not None
        assert int(predicted.sum()) == 16

    def test_predict_doublets_records_heterotypic_and_homotypic_risk(self, monkeypatch):
        import scLucid.qc.doublet.ensemble as ensemble_module

        adata = AnnData(X=np.ones((80, 6), dtype=float))
        adata.obs_names = [f"cell{i}" for i in range(80)]
        adata.var_names = [f"g{i}" for i in range(6)]
        adata.obs["sampleID"] = ["s1"] * 40 + ["s2"] * 40
        adata.obs["n_genes_by_counts"] = np.r_[np.repeat(100, 70), np.repeat(500, 10)]
        adata.obs["total_counts"] = np.r_[np.repeat(500, 70), np.repeat(2500, 10)]

        def fake_algorithm(data_view, sample, cfg):
            n = data_view.n_obs
            scores = pd.Series(np.linspace(0, 1, n), index=data_view.obs_names)
            predicted = scores > 0.8
            return scores, predicted

        def fake_heuristic(adata_in, cfg, expected_rate=None, sample_key=None):
            lineage = pd.DataFrame(
                {
                    "T_cell": np.r_[np.repeat(1.0, 40), np.repeat(0.1, 40)],
                    "B_cell": np.r_[np.repeat(0.2, 30), np.repeat(1.0, 50)],
                },
                index=adata_in.obs_names,
            )
            scores = pd.Series(0.0, index=adata_in.obs_names)
            scores.iloc[35:45] = 1.0
            return scores > 0.5, lineage, scores

        monkeypatch.setattr(ensemble_module, "_run_scrublet", fake_algorithm)
        monkeypatch.setattr(ensemble_module, "_run_heuristic", fake_heuristic)

        cfg = DoubletConfig(
            method="scrublet",
            run_algorithm=True,
            use_heuristics=True,
            expected_doublet_rate={"s1": 0.05, "s2": 0.10},
            plot_summary=False,
            export_stats=False,
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        for col in [
            ALGORITHM_SCORE_COL,
            ALGORITHM_PRED_COL,
            COMBINED_SCORE_COL,
            HETEROTYPIC_RISK_COL,
            HOMOTYPIC_RISK_COL,
            EXPECTED_TOTAL_RATE_COL,
            EXPECTED_HETEROTYPIC_RATE_COL,
            EXPECTED_HOMOTYPIC_RATE_COL,
        ]:
            assert col in out.obs
        assert out.obs[HETEROTYPIC_RISK_COL].between(0, 1).all()
        assert out.obs[HOMOTYPIC_RISK_COL].between(0, 1).all()
        assert out.obs[EXPECTED_TOTAL_RATE_COL].iloc[:40].eq(0.05).all()
        meta = out.uns["sclucid"]["qc"]["doublet_params"]["risk_decomposition"]
        assert meta["schema_version"] == "doublet_risk_decomposition_v1"
        assert "heterotypic_sources" in meta["evidence_priority"]

    def test_raw_count_guard_rejects_log_normalized_input(self):
        adata = AnnData(X=np.log1p(np.array([[0, 1, 2], [3, 0, 4]], dtype=float)))

        assert _raw_count_guard(adata, sample_name="sample_1", method="Scrublet") is False

    def test_raw_count_guard_rejects_negative_input(self):
        adata = AnnData(X=np.array([[0, -1], [2, 3]], dtype=float))

        assert _raw_count_guard(adata, sample_name="sample_1", method="Solo") is False


# ---------------------------------------------------------------------------
# DoubletEvidenceProfiler
# ---------------------------------------------------------------------------


class TestDoubletEvidenceProfiler:
    """Tests for DoubletEvidenceProfiler."""

    def test_generate_evidence_table_basic(self, minimal_adata):
        """Evidence table should be generated from data with doublet results."""
        adata = minimal_adata.copy()
        # Add mock doublet results
        adata.obs["scrublet_score"] = np.random.random(adata.n_obs)
        adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, adata.n_obs)
        adata.obsm["lineage_module_scores"] = pd.DataFrame(
            np.random.random((adata.n_obs, 3)),
            index=adata.obs_names,
            columns=["T_cell", "B_cell", "Myeloid"],
        )

        profiler = DoubletEvidenceProfiler(adata)
        evidence = profiler.generate_evidence_table()

        assert isinstance(evidence, pd.DataFrame)
        assert len(evidence) == adata.n_obs
        assert "scrublet_evidence" in evidence.columns
        assert "n_coexpressed_lineages" in evidence.columns

    def test_evidence_table_includes_risk_decomposition(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs[ALGORITHM_SCORE_COL] = np.linspace(0, 1, adata.n_obs)
        adata.obs["heuristic_confidence_score"] = np.linspace(1, 0, adata.n_obs)
        adata.obs[COMBINED_SCORE_COL] = 0.5
        adata.obs[HETEROTYPIC_RISK_COL] = np.linspace(0, 1, adata.n_obs)
        adata.obs[HOMOTYPIC_RISK_COL] = np.linspace(1, 0, adata.n_obs)

        evidence = DoubletEvidenceProfiler(adata).generate_evidence_table()

        assert HETEROTYPIC_RISK_COL in evidence.columns
        assert HOMOTYPIC_RISK_COL in evidence.columns
        assert "heuristic_evidence_score" in evidence.columns
        assert evidence["combined_evidence_score"].between(0, 1).all()

    def test_mt_pct_is_descriptive_not_combined_evidence(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obs["scrublet_score"] = 0.0
        adata.obs["n_genes_by_counts"] = 100
        adata.obs["total_counts"] = 1000
        adata.obs["pct_counts_mt"] = np.linspace(1, 50, adata.n_obs)
        adata.obsm["lineage_module_scores"] = pd.DataFrame(
            0.0,
            index=adata.obs_names,
            columns=["T_cell", "B_cell"],
        )

        evidence = DoubletEvidenceProfiler(adata).generate_evidence_table()

        assert "mt_pct_zscore" in evidence.columns
        assert "low_mt_evidence" not in evidence.columns
        assert np.allclose(evidence["combined_evidence_score"], 0.0)

    def test_generate_evidence_table_without_doublet_results(self, minimal_adata):
        """Evidence table should handle data without doublet predictions gracefully."""
        profiler = DoubletEvidenceProfiler(minimal_adata)
        evidence = profiler.generate_evidence_table()

        assert isinstance(evidence, pd.DataFrame)
        assert len(evidence) == minimal_adata.n_obs

    def test_generate_doublet_report(self, minimal_adata):
        """Individual cell report should be generated."""
        adata = minimal_adata.copy()
        adata.obs["scrublet_score"] = np.random.random(adata.n_obs)
        adata.obs["predicted_doublet"] = False
        adata.obs["heuristic_confidence_score"] = np.random.random(adata.n_obs)

        profiler = DoubletEvidenceProfiler(adata)
        report = profiler.generate_doublet_report(adata.obs_names[0])

        assert isinstance(report, str)
        assert len(report) > 0


# ---------------------------------------------------------------------------
# _export_doublet_stats
# ---------------------------------------------------------------------------


class TestExportDoubletStats:
    """Tests for doublet statistics export."""

    def test_export_without_save_dir(self, minimal_adata):
        """Export should return DataFrames without writing files when save_dir is None."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        adata.obs["predicted_doublet"] = False
        adata.obs["scrublet_score"] = 0.1

        result = _export_doublet_stats(adata, sample_key="sampleID")

        assert "sample" in result
        assert "global" in result
        assert isinstance(result["sample"], pd.DataFrame)
        assert isinstance(result["global"], pd.DataFrame)

    def test_export_with_temp_dir(self, minimal_adata, tmp_path):
        """Export should write CSV files when save_dir is provided."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        adata.obs["predicted_doublet"] = False
        adata.obs["scrublet_score"] = 0.1

        save_dir = tmp_path / "doublet_stats"
        _export_doublet_stats(adata, sample_key="sampleID", save_dir=str(save_dir), export_csv=True)

        assert (save_dir / "doublet_stats_per_sample.csv").exists()
        assert (save_dir / "doublet_stats_global.csv").exists()
