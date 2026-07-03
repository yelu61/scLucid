"""Tests for scLucid.qc.policy.adaptive_threshold module."""

import numpy as np
import pytest

from scLucid.qc.policy.adaptive_threshold import (
    AdaptiveThresholdLearner,
    MultiMetricAdaptiveLearner,
    THRESHOLD_RESULT_SCHEMA_VERSION,
    compute_kdistance_eps,
    fit_count_mixture_threshold_model,
    infer_qc_metric_type,
    recommended_threshold_methods,
    _zinb_loglik,
    _zinb_loglik_logparams,
)


class TestCountMixtureThresholdModel:
    """Tests for count-aware threshold fitting of n_genes-like metrics."""

    def test_nb_data_fit_and_threshold(self):
        """NB-generated data should be fit by NB or ZINB with a sensible threshold."""
        rng = np.random.default_rng(42)
        mu_true, alpha_true = 1200.0, 0.3
        r = 1.0 / alpha_true
        p = 1.0 / (1.0 + alpha_true * mu_true)
        y = rng.negative_binomial(r, p, size=5000)
        result = fit_count_mixture_threshold_model(
            y, direction="lower", percentile=10, model="auto", random_state=42
        )
        assert result["is_success"]
        assert result["model"] in {"nb", "zinb", "poisson"}
        # 10th percentile threshold should be below the mean for this right-skewed
        # distribution, and it should be a non-negative finite count.
        assert 0 < result["threshold"] < float(np.mean(y))
        assert "aic" in result
        assert result["aic"] is not None

    def test_zinb_threshold_above_nb_threshold(self):
        """ZINB lower-tail threshold should account for structural zeros."""
        rng = np.random.default_rng(7)
        mu_true, alpha_true = 1500.0, 0.2
        r = 1.0 / alpha_true
        p = 1.0 / (1.0 + alpha_true * mu_true)
        nb = rng.negative_binomial(r, p, size=5000)
        # Inject 20% structural zeros
        mask = rng.random(size=5000) < 0.2
        y = nb.copy()
        y[mask] = 0

        nb_result = fit_count_mixture_threshold_model(
            y, direction="lower", percentile=10, model="nb", random_state=7
        )
        zinb_result = fit_count_mixture_threshold_model(
            y, direction="lower", percentile=10, model="zinb", random_state=7
        )
        assert zinb_result["is_success"]
        # ZINB must recognise the extra zeros and not under-shoot.
        assert zinb_result["threshold"] >= nb_result["threshold"]

    def test_poisson_selected_for_equidispersed_data(self):
        """Under/equidispersed count data should prefer Poisson."""
        rng = np.random.default_rng(13)
        y = rng.poisson(lam=800, size=2000)
        result = fit_count_mixture_threshold_model(
            y, direction="lower", percentile=10, model="auto", random_state=13
        )
        assert result["is_success"]
        assert result["model"] == "poisson"

    def test_fallback_for_empty_input(self):
        """Empty input falls back to GMM when fallback=True."""
        result = fit_count_mixture_threshold_model(
            np.array([]), direction="lower", percentile=10, fallback=True
        )
        assert result["is_success"]
        assert result["fallback_used"] is True
        assert np.isnan(result["threshold"])

    def test_zinb_logparams_treats_first_parameter_as_logit(self):
        """Fallback ZINB optimizer should interpret the first parameter as logit(pi)."""
        y = np.array([0, 0, 1, 2, 4, 8], dtype=float)
        pi, mu, alpha = 0.2, 3.0, 0.4
        logit_pi = np.log(pi / (1.0 - pi))

        direct = _zinb_loglik(np.array([pi, mu, alpha]), y)
        transformed = _zinb_loglik_logparams(np.array([logit_pi, np.log(mu), np.log(alpha)]), y)

        assert transformed == pytest.approx(direct)


class TestDBSCANLargeSample:
    """Tests for DBSCAN threshold learning on large datasets."""

    def test_dbscan_detects_outliers_at_large_n(self):
        """DBSCAN should still find outliers when N is large."""
        rng = np.random.default_rng(99)
        normal = rng.normal(loc=1000.0, scale=50.0, size=100_000)
        outliers = rng.uniform(low=1500, high=1600, size=100)
        values = np.concatenate([normal, outliers])
        learner = AdaptiveThresholdLearner(method="dbscan", random_state=99)
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)
        # The detected outlier boundary should separate most normal cells from
        # the outlier cloud.
        assert threshold > np.percentile(normal, 95)
        assert threshold < np.max(outliers)

    def test_compute_kdistance_eps_subsampling(self):
        """k-distance eps estimation should subsample large inputs."""
        rng = np.random.default_rng(101)
        values = rng.normal(loc=0.0, scale=1.0, size=50_000)
        eps_full = compute_kdistance_eps(values, k=5, max_samples=50_000)
        eps_sub = compute_kdistance_eps(values, k=5, max_samples=5_000, random_state=101)
        assert np.isfinite(eps_full)
        assert np.isfinite(eps_sub)
        # Subsampled estimate should be within a reasonable factor of the full one.
        assert 0.05 < eps_sub / eps_full < 20.0


class TestAdaptiveThresholdLearner:
    """Tests for single-metric adaptive threshold learning."""

    def test_percentile_method_upper(self):
        """Percentile method should return sensible upper threshold."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)
        assert threshold > values.min()

    def test_percentile_method_lower(self):
        """Percentile method should return sensible lower threshold."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="lower"
        )["threshold"]
        assert not np.isnan(threshold)
        assert threshold < values.max()

    def test_mad_method(self):
        """MAD method should detect outliers."""
        learner = AdaptiveThresholdLearner(method="mad")
        # Data with clear outliers
        values = np.concatenate([
            np.random.normal(100, 5, 95),
            np.array([200, 210, 220]),  # outliers
        ])
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)
        # Threshold should catch the outliers
        assert threshold < 200

    def test_gmm_method(self):
        """GMM method should learn threshold from bimodal data."""
        learner = AdaptiveThresholdLearner(method="gmm")
        # Bimodal distribution
        values = np.concatenate([
            np.random.normal(50, 5, 50),
            np.random.normal(150, 10, 50),
        ])
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)
        # Threshold should separate the two modes
        assert 50 < threshold < 150

    def test_kde_method(self):
        """KDE method should learn threshold."""
        learner = AdaptiveThresholdLearner(method="kde")
        values = np.random.normal(100, 20, 100)
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)

    def test_dbscan_method(self):
        """DBSCAN method should detect outliers."""
        learner = AdaptiveThresholdLearner(method="dbscan")
        values = np.concatenate([
            np.random.normal(100, 5, 90),
            np.array([200, 205, 210]),  # clear outliers
        ])
        threshold = learner.learn_threshold_result(
            values, "test_metric", direction="upper"
        )["threshold"]
        assert not np.isnan(threshold)

    def test_empty_values_returns_nan(self):
        """Empty input should return NaN."""
        learner = AdaptiveThresholdLearner(method="percentile")
        threshold = learner.learn_threshold_result(np.array([]), "test_metric")["threshold"]
        assert np.isnan(threshold)

    def test_nan_values_handled(self):
        """NaN values should be ignored."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, np.nan, 4, 5, np.nan])
        threshold = learner.learn_threshold_result(values, "test_metric")["threshold"]
        assert not np.isnan(threshold)

    def test_unknown_method_raises(self):
        """Unknown method should raise ValueError."""
        learner = AdaptiveThresholdLearner(method="unknown")
        with pytest.raises(ValueError, match="Unknown method"):
            learner.learn_threshold_result(np.array([1, 2, 3]), "test_metric")

    def test_learned_thresholds_stored(self):
        """Learned thresholds should be stored in _learned_thresholds."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5])
        learner.learn_threshold_result(values, "metric_a", direction="upper")
        learner.learn_threshold_result(values, "metric_b", direction="lower")

        assert "metric_a" in learner._learned_thresholds
        assert "metric_b" in learner._learned_thresholds

    def test_numeric_compatibility_wrapper_warns(self):
        """Numeric-only wrapper should remain visibly deprecated."""
        learner = AdaptiveThresholdLearner(method="percentile")
        with pytest.warns(FutureWarning, match="learn_threshold_result"):
            threshold = learner.learn_threshold(
                np.array([1, 2, 3, 4, 5]), "test_metric", direction="upper"
            )
        assert isinstance(threshold, float)

    def test_auto_threshold_result_has_audit_schema_for_count_metric(self):
        """Structured threshold results should expose provenance and removal estimates."""
        rng = np.random.default_rng(123)
        values = rng.negative_binomial(n=8, p=0.01, size=500).astype(float)
        learner = AdaptiveThresholdLearner(method="auto", random_state=123)

        result = learner.learn_threshold_result(
            values, "n_genes_by_counts", direction="lower"
        )

        assert result["schema_version"] == THRESHOLD_RESULT_SCHEMA_VERSION
        assert result["metric_name"] == "n_genes_by_counts"
        assert result["metric_type"] == "count"
        assert result["direction"] == "lower"
        assert result["method"] == "count_mixture"
        assert result["lower"] == pytest.approx(result["threshold"])
        assert result["upper"] is None
        assert result["n_cells"] == len(values)
        assert 0.0 <= result["removed_fraction_estimate"] <= 1.0
        assert 0.0 <= result["confidence"] <= 1.0
        assert "count_mixture" in result["recommended_methods"]
        assert result["score_semantics"] == (
            "Threshold recommendation, not an automatic removal decision."
        )

    def test_learn_all_threshold_results_preserves_numeric_api(self, minimal_adata):
        """Batch schema API should coexist with legacy float threshold output."""
        from tests.fixtures.synthetic_data import SyntheticDataGenerator

        gen = SyntheticDataGenerator()
        adata = gen.generate_adata(n_cells=120, n_genes=200)
        adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, adata.n_obs)
        adata.obs["pct_counts_mt"] = np.random.uniform(0, 20, adata.n_obs)

        metrics = {"n_genes_by_counts": "lower", "pct_counts_mt": "upper"}
        learner = AdaptiveThresholdLearner(method="percentile")
        with pytest.warns(FutureWarning, match="learn_all_threshold_results"):
            numeric = learner.learn_all_thresholds(adata, metrics=metrics)
        structured = learner.learn_all_threshold_results(adata, metrics=metrics)

        assert set(numeric) == set(metrics)
        assert set(structured) == set(metrics)
        assert all(isinstance(value, float) for value in numeric.values())
        assert all(
            item["schema_version"] == THRESHOLD_RESULT_SCHEMA_VERSION
            for item in structured.values()
        )

    def test_review_only_metric_result_is_not_hard_filter_semantics(self):
        """Review evidence metrics should be explicitly labelled as non-default filters."""
        values = np.linspace(0, 1, 200)
        learner = AdaptiveThresholdLearner(method="auto")

        result = learner.learn_threshold_result(values, "ambient_score", direction="upper")

        assert result["metric_type"] == "review_evidence"
        assert result["method"] == "percentile"
        assert "review evidence" in result["review_note"]
        assert "not an automatic removal decision" in result["score_semantics"]

    def test_metric_type_and_recommended_methods_are_explicit(self):
        """Metric family inference should drive default threshold method selection."""
        assert infer_qc_metric_type("n_genes_by_counts") == "count"
        assert infer_qc_metric_type("log1p_total_counts") == "log_count"
        assert infer_qc_metric_type("pct_counts_mt") == "percentage"
        assert infer_qc_metric_type("ambient_score") == "review_evidence"

        assert recommended_threshold_methods("n_genes_by_counts", n_cells=500)[0] == (
            "count_mixture"
        )
        assert recommended_threshold_methods("pct_counts_mt", n_cells=500)[0] == "mad"
        assert recommended_threshold_methods("ambient_score", n_cells=500)[0] == (
            "percentile"
        )


class TestMultiMetricAdaptiveLearner:
    """Tests for multi-metric adaptive threshold learning."""

    def test_fit_isolation_forest(self, minimal_adata):
        """Should fit isolation forest on QC metrics."""
        from tests.fixtures.synthetic_data import SyntheticDataGenerator

        gen = SyntheticDataGenerator()
        adata = gen.generate_adata(n_cells=100, n_genes=200)
        adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, adata.n_obs)
        adata.obs["pct_counts_mt"] = np.random.uniform(0, 20, adata.n_obs)

        learner = MultiMetricAdaptiveLearner(method="isolation_forest")
        learner.fit(adata, metrics=["n_genes_by_counts", "pct_counts_mt"])
        assert learner._model is not None

    def test_predict_after_fit(self, minimal_adata):
        """Should predict outlier labels after fitting."""
        from tests.fixtures.synthetic_data import SyntheticDataGenerator

        gen = SyntheticDataGenerator()
        adata = gen.generate_adata(n_cells=100, n_genes=200)
        adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, adata.n_obs)
        adata.obs["pct_counts_mt"] = np.random.uniform(0, 20, adata.n_obs)

        learner = MultiMetricAdaptiveLearner(method="isolation_forest")
        learner.fit(adata, metrics=["n_genes_by_counts", "pct_counts_mt"])
        predictions = learner.predict(adata, metrics=["n_genes_by_counts", "pct_counts_mt"])

        assert len(predictions) == adata.n_obs
        assert set(np.unique(predictions)).issubset({False, True})

    def test_unknown_method_raises(self):
        """Unknown method should raise ValueError."""
        learner = MultiMetricAdaptiveLearner(method="unknown")
        with pytest.raises(ValueError, match="Unknown method"):
            import numpy as np
            from anndata import AnnData
            adata = AnnData(np.random.random((10, 5)))
            adata.obs["m1"] = np.random.random(10)
            learner.fit(adata, metrics=["m1"])
