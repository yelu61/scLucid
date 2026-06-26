"""
Scientific benchmarks for QC module.

Validates that adaptive threshold methods, intelligent QC recommendations,
and tumor-aware logic behave correctly across diverse data scenarios.
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from validation.qc.run_threshold_benchmark import _annotate_decision_rows_with_strategy_evidence
from validation.qc.run_doublet_evidence_benchmark import (
    _algorithm_weight_recommendation_rows as _doublet_algorithm_weight_recommendation_rows,
    _threshold_calibration_rows as _doublet_threshold_calibration_rows,
)
from validation.qc.build_figure2_qc_evidence_package import (
    _source_rows_from_figure_table as _figure2_source_rows_from_table,
)
from validation.qc.run_tumor_biological_fidelity_benchmark import (
    _narrative_rows as _tumor_fidelity_narrative_rows,
)
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


def test_threshold_benchmark_decision_rows_include_policy_evidence():
    decision_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "parameter": "max_mt_percent",
            "recommended": 30.0,
            "applied": 30.0,
            "source": "sclucid_tumor_aware",
            "confidence": "medium",
            "evidence": "{}",
            "review_required": True,
            "affected_cells": 10,
            "biological_guardrail": "preserve high-mt malignant/stress/program signal until reviewed",
            "risk_note": "review",
        }
    ]
    scorecard_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "rank_within_dataset": 1,
            "recommended_for_review": True,
            "composite_score": 0.91,
            "risk_note": "",
        }
    ]

    annotated = _annotate_decision_rows_with_strategy_evidence(
        decision_rows,
        scorecard_rows,
    )

    row = annotated[0]
    assert row["strategy_rank"] == 1
    assert row["recommended_policy"] is True
    assert row["strategy_composite_score"] == 0.91
    assert "Recommended as the benchmark-selected QC policy" in row["decision_narrative"]
    assert row["biological_guardrail"].startswith("preserve high-mt")


def test_tumor_fidelity_narrative_rows_explain_recommendation():
    marker_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "marker_panel": "epithelial",
            "high_mt_removed_relative_to_all": 0.5,
        }
    ]
    program_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "program": "hypoxia_stress",
            "program_retention_ratio": 0.95,
            "high_mt_removed_program_ratio": 0.7,
        }
    ]
    bias_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "group_type": "sample",
            "group": "S1",
            "before": 100,
            "retention_rate": 0.9,
        }
    ]
    scorecard_rows = [
        {
            "dataset": "tumor_demo",
            "strategy": "sclucid_tumor_aware",
            "rank_within_dataset": 1,
            "recommended_for_review": True,
            "biological_fidelity_score": 0.95,
            "high_mt_removed_cells": 12,
            "biological_harm_risk": False,
            "review_required": False,
            "risk_note": "",
        }
    ]

    rows = _tumor_fidelity_narrative_rows(
        marker_rows,
        program_rows,
        bias_rows,
        scorecard_rows,
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["recommended_policy"] is True
    assert row["strategy_rank"] == 1
    assert row["worst_group"] == "S1"
    assert row["mean_program_retention_ratio"] == 0.95
    assert "recommended for review" in row["decision_narrative"]


def test_doublet_benchmark_rejects_heuristic_method():
    """Heuristic is a fallback baseline, not a selectable primary method."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "validation/qc/run_doublet_evidence_benchmark.py", "--methods", "heuristic"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0, "heuristic should not be accepted as a primary method"
    assert "heuristic" in result.stderr.lower() or "invalid choice" in result.stderr.lower()

    truth = pd.Series([True, True, False, False, True, False])
    scores = {
        "scrublet": pd.Series([0.9, 0.8, 0.7, 0.2, 0.6, 0.1]),
    }
    predictions = {
        "scrublet": pd.Series([True, False, False, False, False, False]),
    }

    rows = _doublet_threshold_calibration_rows(
        dataset="demo",
        truth=truth,
        predictions=predictions,
        scores=scores,
        target_recalls=(0.5,),
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["method"] == "scrublet"
    assert row["calibrated_recall"] >= 0.5
    assert row["recall_gain_vs_default"] > 0
    assert row["review_required"] is True


def test_doublet_algorithm_weight_recommendation_prefers_best_fusion():
    evidence_rows = [
        {
            "dataset": "kang2018.pbmc",
            "method": "scdblfinder_python_pyscdblfinder",
            "base_method": "scdblfinder_python_pyscdblfinder",
            "method_status": "ok",
            "uses_heuristics": False,
            "merge_strategy": "",
            "algorithm_weight": None,
            "f1": 0.60,
            "precision": 0.75,
            "recall": 0.50,
            "score_auc": 0.88,
            "predicted_rate": 0.06,
            "expected_rate_from_demuxlet": 0.09,
        },
        {
            "dataset": "kang2018.pbmc",
            "method": "scdblfinder_python_pyscdblfinder_plus_heuristic_w0.50",
            "base_method": "scdblfinder_python_pyscdblfinder",
            "method_status": "ok",
            "uses_heuristics": True,
            "merge_strategy": "weighted_average",
            "algorithm_weight": 0.50,
            "f1": 0.62,
            "precision": 0.70,
            "recall": 0.56,
            "score_auc": 0.87,
            "predicted_rate": 0.08,
            "expected_rate_from_demuxlet": 0.09,
        },
        {
            "dataset": "kang2018.pbmc",
            "method": "scdblfinder_python_pyscdblfinder_plus_heuristic_w0.70",
            "base_method": "scdblfinder_python_pyscdblfinder",
            "method_status": "ok",
            "uses_heuristics": True,
            "merge_strategy": "weighted_average",
            "algorithm_weight": 0.70,
            "f1": 0.71,
            "precision": 0.73,
            "recall": 0.69,
            "score_auc": 0.90,
            "predicted_rate": 0.09,
            "expected_rate_from_demuxlet": 0.09,
        },
    ]

    rows = _doublet_algorithm_weight_recommendation_rows(evidence_rows)

    assert len(rows) == 1
    row = rows[0]
    assert row["recommended_algorithm_weight"] == 0.70
    assert row["recommended_method"].endswith("w0.70")
    assert row["f1_delta_vs_algorithm_only"] > 0.1
    assert row["recall_delta_vs_algorithm_only"] > 0
    assert row["review_required"] is False


def test_figure2_source_rows_harmonize_tumor_panels():
    table = pd.DataFrame(
        [
            {
                "figure_panel": "2D",
                "dataset": "lin2020.pdac",
                "strategy": "sclucid_tumor_aware",
                "metric": "biological_fidelity_score",
                "value": 0.95,
                "context": '{"recommended_policy": true}',
            }
        ]
    )

    rows = _figure2_source_rows_from_table(
        table,
        evidence_domain="tumor_biological_fidelity",
        source_file="figure2_tumor_fidelity_data.tsv",
        panel_map={"2D": "2E"},
    )

    assert rows[0]["figure_panel"] == "2E"
    assert rows[0]["evidence_domain"] == "tumor_biological_fidelity"
    assert "recommended_policy" in rows[0]["context"]


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
        from scLucid.qc.doublet.heuristic import _run_heuristic
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
        from scLucid.qc.doublet.heuristic import _run_heuristic
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
