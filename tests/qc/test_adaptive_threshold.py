"""Tests for scLucid.qc.adaptive_threshold module."""

import numpy as np
import pytest

from scLucid.qc.adaptive_threshold import AdaptiveThresholdLearner, MultiMetricAdaptiveLearner


class TestAdaptiveThresholdLearner:
    """Tests for single-metric adaptive threshold learning."""

    def test_percentile_method_upper(self):
        """Percentile method should return sensible upper threshold."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        threshold = learner.learn_threshold(values, "test_metric", direction="upper")
        assert not np.isnan(threshold)
        assert threshold > values.min()

    def test_percentile_method_lower(self):
        """Percentile method should return sensible lower threshold."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        threshold = learner.learn_threshold(values, "test_metric", direction="lower")
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
        threshold = learner.learn_threshold(values, "test_metric", direction="upper")
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
        threshold = learner.learn_threshold(values, "test_metric", direction="upper")
        assert not np.isnan(threshold)
        # Threshold should separate the two modes
        assert 50 < threshold < 150

    def test_kde_method(self):
        """KDE method should learn threshold."""
        learner = AdaptiveThresholdLearner(method="kde")
        values = np.random.normal(100, 20, 100)
        threshold = learner.learn_threshold(values, "test_metric", direction="upper")
        assert not np.isnan(threshold)

    def test_dbscan_method(self):
        """DBSCAN method should detect outliers."""
        learner = AdaptiveThresholdLearner(method="dbscan")
        values = np.concatenate([
            np.random.normal(100, 5, 90),
            np.array([200, 205, 210]),  # clear outliers
        ])
        threshold = learner.learn_threshold(values, "test_metric", direction="upper")
        assert not np.isnan(threshold)

    def test_empty_values_returns_nan(self):
        """Empty input should return NaN."""
        learner = AdaptiveThresholdLearner(method="percentile")
        threshold = learner.learn_threshold(np.array([]), "test_metric")
        assert np.isnan(threshold)

    def test_nan_values_handled(self):
        """NaN values should be ignored."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, np.nan, 4, 5, np.nan])
        threshold = learner.learn_threshold(values, "test_metric")
        assert not np.isnan(threshold)

    def test_unknown_method_raises(self):
        """Unknown method should raise ValueError."""
        learner = AdaptiveThresholdLearner(method="unknown")
        with pytest.raises(ValueError, match="Unknown method"):
            learner.learn_threshold(np.array([1, 2, 3]), "test_metric")

    def test_learned_thresholds_stored(self):
        """Learned thresholds should be stored in _learned_thresholds."""
        learner = AdaptiveThresholdLearner(method="percentile")
        values = np.array([1, 2, 3, 4, 5])
        learner.learn_threshold(values, "metric_a", direction="upper")
        learner.learn_threshold(values, "metric_b", direction="lower")

        assert "metric_a" in learner._learned_thresholds
        assert "metric_b" in learner._learned_thresholds


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
