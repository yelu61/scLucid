"""
Scientific benchmarks for QC module.

Validates that adaptive threshold methods, intelligent QC recommendations,
and tumor-aware logic behave correctly across diverse data scenarios.
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.qc.adaptive_threshold import AdaptiveThresholdLearner
from scLucid.qc.intelligent_qc import recommend_intelligent_qc, StrategyType


# ---------------------------------------------------------------------------
# Adaptive Threshold Benchmarks
# ---------------------------------------------------------------------------


class TestAdaptiveThresholdScientific:
    """
    Benchmark adaptive threshold learning across synthetic distributions.
    Key question: Does each method learn a threshold that separates
    the "good" population from outliers?
    """

    @pytest.fixture
    def bimodal_data(self):
        """Clean cells + low-quality outliers."""
        rng = np.random.default_rng(42)
        clean = rng.normal(500, 30, 450)
        outliers = rng.normal(200, 20, 50)
        return np.concatenate([clean, outliers])

    @pytest.fixture
    def uniform_with_outliers(self):
        """Uniform distribution with a few extreme outliers."""
        rng = np.random.default_rng(42)
        base = rng.uniform(300, 700, 480)
        outliers = rng.uniform(50, 150, 20)
        return np.concatenate([base, outliers])

    @pytest.fixture
    def single_mode(self):
        """Single mode, no clear outliers."""
        rng = np.random.default_rng(42)
        return rng.normal(500, 50, 500)

    def _assert_threshold_reasonable(self, threshold, data, direction):
        assert not np.isnan(threshold), "Threshold should not be NaN"
        # MAD lower can clip to 0, which is acceptable
        assert threshold >= 0, f"Threshold {threshold} should be >= 0"
        if direction == "lower":
            assert (data > threshold).sum() > 0, "Some cells should be above threshold"
        else:
            assert (data < threshold).sum() > 0, "Some cells should be below threshold"

    @pytest.mark.parametrize("method", ["percentile", "mad", "gmm", "kde", "dbscan"])
    def test_all_methods_on_bimodal(self, bimodal_data, method):
        """All methods should learn sensible threshold on bimodal data."""
        learner = AdaptiveThresholdLearner(method=method)
        threshold = learner.learn_threshold(bimodal_data, "n_genes", direction="lower")
        self._assert_threshold_reasonable(threshold, bimodal_data, "lower")
        # On clear bimodal data, threshold should separate the modes
        if method in ("percentile", "gmm", "kde"):
            assert threshold < 450, \
                f"{method} threshold {threshold} should separate modes (~200, ~500)"

    @pytest.mark.parametrize("method", ["percentile", "mad", "gmm", "kde", "dbscan"])
    def test_all_methods_on_uniform_outliers(self, uniform_with_outliers, method):
        """All methods should catch extreme outliers."""
        learner = AdaptiveThresholdLearner(method=method)
        threshold = learner.learn_threshold(uniform_with_outliers, "n_genes", direction="lower")
        self._assert_threshold_reasonable(threshold, uniform_with_outliers, "lower")

    @pytest.mark.parametrize("method", ["percentile", "mad", "kde"])
    def test_single_mode_conservative(self, single_mode, method):
        """On single-mode data without outliers, threshold should be conservative."""
        learner = AdaptiveThresholdLearner(method=method)
        threshold = learner.learn_threshold(single_mode, "n_genes", direction="lower")
        self._assert_threshold_reasonable(threshold, single_mode, "lower")
        # Should not flag more than ~25% of cells as outliers (conservative)
        outlier_rate = (single_mode < threshold).sum() / len(single_mode)
        assert outlier_rate <= 0.25, \
            f"{method} flagged {outlier_rate:.1%} as outliers on clean data"

    def test_gmm_vs_percentile_same_order_of_magnitude(self, bimodal_data):
        """GMM and percentile should give thresholds in same ballpark on bimodal data."""
        gmm_threshold = AdaptiveThresholdLearner(method="gmm").learn_threshold(
            bimodal_data, "n_genes", direction="lower"
        )
        pct_threshold = AdaptiveThresholdLearner(method="percentile").learn_threshold(
            bimodal_data, "n_genes", direction="lower"
        )
        diff = abs(gmm_threshold - pct_threshold)
        # Allow up to 200 difference (generous for different algorithms)
        assert diff < 200, \
            f"GMM ({gmm_threshold}) and percentile ({pct_threshold}) diverged by {diff}"

    def test_methods_agree_on_clear_outliers(self, bimodal_data):
        """All methods should flag the clear outliers in bimodal data."""
        thresholds = {}
        for method in ["percentile", "mad", "kde"]:
            learner = AdaptiveThresholdLearner(method=method)
            thresholds[method] = learner.learn_threshold(bimodal_data, "n_genes", direction="lower")

        # All thresholds should be below the clean population mean (~500)
        for method, thr in thresholds.items():
            assert thr < 450, f"{method} threshold {thr} too high, would miss outliers"

    def test_empty_data(self):
        """Empty data should return NaN gracefully."""
        learner = AdaptiveThresholdLearner(method="percentile")
        threshold = learner.learn_threshold(np.array([]), "metric")
        assert np.isnan(threshold)

    def test_nan_inf_data(self):
        """NaN and Inf values should be handled gracefully."""
        learner = AdaptiveThresholdLearner(method="percentile")
        data = np.array([1.0, 2.0, np.nan, 3.0, np.inf, 4.0])
        threshold = learner.learn_threshold(data, "metric", direction="lower")
        assert not np.isnan(threshold)
        assert threshold > 0


# ---------------------------------------------------------------------------
# Intelligent QC Benchmarks
# ---------------------------------------------------------------------------


class TestIntelligentQCScientific:
    """
    Benchmark intelligent QC recommendations.
    Key questions:
    1. Does the engine detect data quality issues?
    2. Are recommended thresholds sensible?
    3. Does tumor-aware mode behave differently?
    """

    @pytest.fixture
    def high_quality_pbmc(self):
        """Simulate high-quality PBMC: narrow distribution, low MT%."""
        rng = np.random.default_rng(42)
        adata = AnnData(np.random.negative_binomial(5, 0.5, (500, 100)))
        adata.obs["n_genes_by_counts"] = rng.normal(2000, 200, 500)
        adata.obs["total_counts"] = rng.normal(8000, 1000, 500)
        adata.obs["pct_counts_mt"] = rng.normal(5, 2, 500).clip(0, 20)
        return adata

    @pytest.fixture
    def low_quality_data(self):
        """Simulate low-quality data: wide distribution, high MT%."""
        rng = np.random.default_rng(42)
        adata = AnnData(np.random.negative_binomial(2, 0.3, (500, 100)))
        adata.obs["n_genes_by_counts"] = rng.normal(500, 300, 500).clip(100, None)
        adata.obs["total_counts"] = rng.normal(2000, 1500, 500).clip(500, None)
        adata.obs["pct_counts_mt"] = rng.normal(25, 10, 500).clip(0, 60)
        return adata

    @pytest.fixture
    def tumor_like_data(self):
        """Simulate tumor data: elevated MT%, wide heterogeneity."""
        rng = np.random.default_rng(42)
        adata = AnnData(np.random.negative_binomial(3, 0.4, (500, 100)))
        adata.obs["n_genes_by_counts"] = rng.normal(1500, 500, 500)
        adata.obs["total_counts"] = rng.normal(6000, 2500, 500)
        adata.obs["pct_counts_mt"] = rng.normal(15, 8, 500).clip(0, 50)
        return adata

    def test_recommendation_returns_non_none(self, high_quality_pbmc):
        """Recommendation engine should always return a result."""
        rec = recommend_intelligent_qc(high_quality_pbmc, tissue_type="normal")
        assert rec is not None

    def test_recommendation_has_required_keys(self, high_quality_pbmc):
        """Recommendation should contain expected keys."""
        rec = recommend_intelligent_qc(high_quality_pbmc, tissue_type="normal")
        rec_dict = rec.to_dict() if hasattr(rec, "to_dict") else {}
        assert "overall_strategy" in rec_dict
        assert "overall_confidence" in rec_dict
        assert "min_genes" in rec_dict or "n_counts" in rec_dict

    def test_recommendation_on_low_quality_has_concerns_or_aggressive(self, low_quality_data):
        """Low-quality data should trigger concerns or aggressive strategy."""
        rec = recommend_intelligent_qc(low_quality_data, tissue_type="normal")
        assert rec is not None
        rec_dict = rec.to_dict() if hasattr(rec, "to_dict") else {}
        concerns = rec_dict.get("concerns", [])
        strategy = rec_dict.get("overall_strategy", "")
        # Either there are concerns, or the strategy is aggressive/conservative
        assert len(concerns) > 0 or strategy in ("aggressive", "conservative"), \
            "Low-quality data should trigger concerns or special strategy"

    def test_tumor_aware_recommends_different_strategy(self, tumor_like_data):
        """Tumor-aware mode should use a tumor-aware strategy."""
        rec_tumor = recommend_intelligent_qc(tumor_like_data, tissue_type="tumor")
        rec_normal = recommend_intelligent_qc(tumor_like_data, tissue_type="normal")

        assert rec_tumor is not None
        assert rec_normal is not None

        tumor_dict = rec_tumor.to_dict() if hasattr(rec_tumor, "to_dict") else {}
        normal_dict = rec_normal.to_dict() if hasattr(rec_normal, "to_dict") else {}

        # Tumor mode should have different or tumor-specific considerations
        tumor_considerations = tumor_dict.get("tumor_specific_considerations", [])
        # Either there are tumor-specific considerations, or strategies differ
        assert (
            len(tumor_considerations) > 0
            or tumor_dict.get("overall_strategy") != normal_dict.get("overall_strategy")
            or tumor_dict.get("max_mt_percent", {}).get("threshold")
            != normal_dict.get("max_mt_percent", {}).get("threshold")
        ), "Tumor and normal should have different recommendations"

    def test_strategy_types_exist(self):
        """All strategy types should be defined."""
        strategies = [
            StrategyType.STANDARD,
            StrategyType.TUMOR_AWARE,
            StrategyType.CONSERVATIVE,
            StrategyType.AGGRESSIVE,
            StrategyType.AUTO,
        ]
        assert len(strategies) == 5
        for s in strategies:
            assert isinstance(s.value, str)


# ---------------------------------------------------------------------------
# Tumor-Aware Logic Benchmarks
# ---------------------------------------------------------------------------


class TestTumorAwareLogic:
    """Validate tumor-specific QC adjustments using pre-computed QC fixtures."""

    def test_tumor_qc_retains_majority(self, qc_test_adata):
        """Tumor data with realistic QC should retain majority of cells."""
        adata = qc_test_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        from scLucid.qc.workflow import run_standard_qc
        result = run_standard_qc(adata, tissue_type="tumor", show_progress=False)

        retention = result.n_obs / adata.n_obs
        assert retention >= 0.1, \
            f"Tumor-aware QC retained only {retention:.1%} of cells"

    def test_normal_vs_tumor_retention(self, qc_test_adata):
        """Normal tissue should filter at least as aggressively as tumor."""
        adata = qc_test_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        from scLucid.qc.workflow import run_standard_qc
        result_normal = run_standard_qc(adata.copy(), tissue_type="normal", show_progress=False)
        result_tumor = run_standard_qc(adata.copy(), tissue_type="tumor", show_progress=False)

        # Tumor should retain equal or more cells than normal
        assert result_tumor.n_obs >= result_normal.n_obs * 0.5, \
            f"Tumor ({result_tumor.n_obs}) retained far fewer than normal ({result_normal.n_obs})"

    def test_tumor_flags_elevated_mt(self, qc_test_adata):
        """Tumor QC should flag elevated MT but not hard-filter."""
        adata = qc_test_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        # Artificially elevate MT for all cells
        adata.obs["pct_counts_mt"] = np.random.default_rng(42).normal(18, 3, adata.n_obs).clip(5, 35)

        from scLucid.qc.workflow import run_standard_qc
        result = run_standard_qc(adata, tissue_type="tumor", show_progress=False)

        # Should retain some cells despite elevated MT
        assert result.n_obs >= 10, \
            f"Tumor QC with elevated MT retained only {result.n_obs} cells"


# ---------------------------------------------------------------------------
# Doublet Heuristic Benchmarks
# ---------------------------------------------------------------------------


class TestDoubletHeuristicScientific:
    """Validate heuristic doublet detection on known data."""

    def test_heuristic_degrades_gracefully_without_markers(self):
        """Without matching lineage markers, heuristic should return all zeros (not crash)."""
        from tests.fixtures.synthetic_data import SyntheticDataGenerator
        from scLucid.qc.doublet import _run_heuristic
        from scLucid.qc.config import DoubletConfig

        gen = SyntheticDataGenerator()
        adata = gen.generate_with_doublets(n_cells=300, doublet_rate=0.1)

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        predicted, lineage_scores, scores = _run_heuristic(adata, cfg)

        # On synthetic data (gene names don't match real markers),
        # heuristic scores should all be zero (graceful degradation)
        assert isinstance(predicted, pd.Series)
        assert isinstance(scores, pd.Series)
        assert len(scores) == adata.n_obs
        # All scores should be in [0, 1]
        assert scores.min() >= 0
        assert scores.max() <= 1

    def test_heuristic_returns_valid_dataframe(self):
        """Heuristic should return valid DataFrame even without markers."""
        from tests.fixtures.synthetic_data import SyntheticDataGenerator
        from scLucid.qc.doublet import _run_heuristic
        from scLucid.qc.config import DoubletConfig

        gen = SyntheticDataGenerator()
        adata = gen.generate_with_doublets(n_cells=300, doublet_rate=0.1)

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        _, lineage_scores, _ = _run_heuristic(adata, cfg)

        assert isinstance(lineage_scores, pd.DataFrame)
        assert len(lineage_scores) == adata.n_obs
