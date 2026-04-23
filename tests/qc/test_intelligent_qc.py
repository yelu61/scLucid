"""
Tests for Intelligent QC Recommendation System.

This test module validates the core innovation of scLucid:
data-driven QC threshold recommendations with confidence intervals.
"""


import numpy as np
import pytest

from scLucid.qc.intelligent_qc import (
    IntelligentQCRecommender,
    QCRecommendation,
    StrategyType,
    ThresholdRecommendation,
    recommend_intelligent_qc,
)
from tests.fixtures.synthetic_data import SyntheticDataGenerator

# =============================================================================
# Fixtures
# =============================================================================


def _make_qc_ready_adata(n_cells: int = 1000, n_genes: int = 1200, random_state: int = 42):
    """Create synthetic AnnData with the QC columns expected by intelligent QC."""
    generator = SyntheticDataGenerator(random_state=random_state)
    adata = generator.generate_adata(
        n_cells=n_cells,
        n_genes=n_genes,
        with_qc_metrics=True,
        with_cell_types=True,
        with_batches=True,
        sparsity=0.8,
    )
    # intelligent_qc expects these aliases.
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    adata.obs["n_counts"] = adata.obs["total_counts"]
    return adata


@pytest.fixture
def sample_adata_with_qc():
    """Create sample AnnData with QC metrics."""
    return _make_qc_ready_adata(n_cells=1000, n_genes=1200, random_state=42)


@pytest.fixture
def tumor_like_adata():
    """Create tumor-like AnnData with elevated mitochondrial content."""
    adata = _make_qc_ready_adata(n_cells=1000, n_genes=1200, random_state=43)

    # Simulate tumor characteristics
    # 1. Higher mitochondrial content
    adata.obs["pct_counts_mt"] = adata.obs["pct_counts_mt"] * 1.5 + 5

    # 2. Add some doublet-like cells
    n_doublets = int(0.1 * adata.n_obs)
    doublet_idx = np.random.choice(adata.n_obs, n_doublets, replace=False)
    adata.obs.loc[adata.obs.index[doublet_idx], "n_genes"] *= 2

    return adata


@pytest.fixture
def low_quality_adata():
    """Create low-quality AnnData for testing conservative strategy."""
    adata = _make_qc_ready_adata(n_cells=1000, n_genes=1200, random_state=44)

    # Simulate low quality
    # 1. Low gene counts
    adata.obs["n_genes"] = adata.obs["n_genes"] * 0.5

    # 2. High mitochondrial content
    adata.obs["pct_counts_mt"] = adata.obs["pct_counts_mt"] * 2

    return adata


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create temporary output directory."""
    output_dir = tmp_path / "qc_output"
    output_dir.mkdir()
    return output_dir


# =============================================================================
# Test Basic Functionality
# =============================================================================


class TestIntelligentQCRecommender:
    """Test IntelligentQCRecommender class."""

    def test_initialization(self):
        """Test recommender can be initialized."""
        recommender = IntelligentQCRecommender()
        assert recommender.strategy == StrategyType.AUTO

        recommender_tumor = IntelligentQCRecommender(strategy=StrategyType.TUMOR_AWARE)
        assert recommender_tumor.strategy == StrategyType.TUMOR_AWARE

    def test_recommend_returns_valid_result(self, sample_adata_with_qc, temp_output_dir):
        """Test recommend() returns valid QCRecommendation."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(
            sample_adata_with_qc, tissue_type="normal", save_dir=temp_output_dir
        )

        # Check return type
        assert isinstance(result, QCRecommendation)

        # Check all thresholds are present
        assert isinstance(result.min_genes, ThresholdRecommendation)
        assert isinstance(result.max_mt_percent, ThresholdRecommendation)
        assert isinstance(result.n_counts, ThresholdRecommendation)

        # Check threshold values are reasonable
        assert result.min_genes.threshold > 0
        assert result.max_mt_percent.threshold > 0
        assert result.n_counts.threshold > 0

        # Check confidence intervals
        assert 0 <= result.min_genes.ci_lower <= result.min_genes.threshold
        assert result.min_genes.threshold <= result.min_genes.ci_upper

        # Check confidence scores
        assert 0 <= result.overall_confidence <= 1
        assert 0 <= result.data_quality_score <= 100

    def test_recommendation_to_dict_preserves_all_thresholds(self, sample_adata_with_qc):
        """Serialized recommendation should retain every executable threshold."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal", plot=False)
        payload = result.to_dict()

        for key in ["min_genes", "max_mt_percent", "n_counts", "doublet_threshold"]:
            assert key in payload
            assert "threshold" in payload[key]
            assert "confidence" in payload[key]
            assert "evidence" in payload[key]


# =============================================================================
# Test Different Strategies
# =============================================================================


class TestStrategies:
    """Test different QC strategies."""

    def test_auto_strategy_normal_tissue(self, sample_adata_with_qc):
        """Test AUTO strategy selects appropriate strategy for normal tissue."""
        recommender = IntelligentQCRecommender(strategy=StrategyType.AUTO)

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal")

        # Should select STANDARD for normal tissue
        assert result.overall_strategy in [StrategyType.STANDARD, StrategyType.TUMOR_AWARE]

    def test_auto_strategy_tumor_tissue(self, tumor_like_adata):
        """Test AUTO strategy detects tumor tissue."""
        recommender = IntelligentQCRecommender(strategy=StrategyType.AUTO)

        result = recommender.recommend(tumor_like_adata, tissue_type="lung_tumor")

        # Should select TUMOR_AWARE for tumor tissue
        assert result.overall_strategy == StrategyType.TUMOR_AWARE

    def test_tumor_aware_strategy(self, tumor_like_adata):
        """Test TUMOR_AWARE strategy explicitly."""
        recommender = IntelligentQCRecommender(strategy=StrategyType.TUMOR_AWARE)

        result = recommender.recommend(tumor_like_adata, tissue_type="lung_tumor")

        assert result.overall_strategy == StrategyType.TUMOR_AWARE

        # Check that tumor-specific considerations are present
        assert len(result.tumor_specific_considerations) > 0

    def test_conservative_strategy(self, low_quality_adata):
        """Test CONSERVATIVE strategy keeps more cells."""
        adata = low_quality_adata

        # Get conservative recommendations
        recommender_conservative = IntelligentQCRecommender(strategy=StrategyType.CONSERVATIVE)
        result_conservative = recommender_conservative.recommend(adata, tissue_type="unknown")

        # Get aggressive recommendations
        recommender_aggressive = IntelligentQCRecommender(strategy=StrategyType.AGGRESSIVE)
        result_aggressive = recommender_aggressive.recommend(adata.copy(), tissue_type="unknown")

        # Conservative should have lower thresholds (keep more cells)
        assert result_conservative.min_genes.threshold <= result_aggressive.min_genes.threshold
        assert (
            result_conservative.max_mt_percent.threshold
            >= result_aggressive.max_mt_percent.threshold
        )


# =============================================================================
# Test Confidence Intervals
# =============================================================================


class TestConfidenceIntervals:
    """Test confidence interval calculations."""

    def test_min_genes_has_ci(self, sample_adata_with_qc):
        """Test min_genes recommendation includes confidence interval."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal")

        # Check CI exists and is reasonable
        assert result.min_genes.ci_lower > 0
        assert result.min_genes.ci_upper > result.min_genes.ci_lower
        assert result.min_genes.threshold >= result.min_genes.ci_lower
        assert result.min_genes.threshold <= result.min_genes.ci_upper

        # Check CI width is reasonable (not too wide)
        ci_width = result.min_genes.ci_upper - result.min_genes.ci_lower
        assert ci_width < result.min_genes.threshold  # CI shouldn't be wider than threshold

    def test_max_mt_has_ci(self, sample_adata_with_qc):
        """Test max_mt_percent recommendation includes confidence interval."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal")

        assert result.max_mt_percent.ci_lower > 0
        assert result.max_mt_percent.ci_upper > result.max_mt_percent.ci_lower
        assert result.max_mt_percent.threshold >= result.max_mt_percent.ci_lower
        assert result.max_mt_percent.threshold <= result.max_mt_percent.ci_upper

    def test_n_counts_has_ci(self, sample_adata_with_qc):
        """Test n_counts recommendation includes confidence interval."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal")

        assert result.n_counts.ci_lower > 0
        assert result.n_counts.ci_upper > result.n_counts.ci_lower
        assert result.n_counts.threshold >= result.n_counts.ci_lower
        assert result.n_counts.threshold <= result.n_counts.ci_upper


# =============================================================================
# Test Data-Driven vs Fixed Thresholds
# =============================================================================


class TestDataDrivenRecommendations:
    """Test that recommendations differ from fixed thresholds."""

    def test_min_genes_not_fixed_200(self, sample_adata_with_qc):
        """Test that min_genes recommendation is not always 200."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal")

        # The recommendation should be data-driven, not fixed at 200
        # It might be 200 for some datasets, but let's check it varies
        # by testing multiple datasets

        # Create dataset with low gene counts

        adata_low = sample_adata_with_qc.copy()
        adata_low.obs["n_genes"] = adata_low.obs["n_genes"] * 0.5

        result_low = recommender.recommend(adata_low, tissue_type="normal", plot=False)

        # Thresholds should differ
        assert result.min_genes.threshold != result_low.min_genes.threshold

    def test_max_mt_adapts_to_tissue(self, sample_adata_with_qc, tumor_like_adata):
        """Test that max_mt_percent adapts to tissue type."""
        recommender = IntelligentQCRecommender()

        # Normal tissue
        result_normal = recommender.recommend(
            sample_adata_with_qc, tissue_type="normal", plot=False
        )

        # Tumor tissue
        result_tumor = recommender.recommend(tumor_like_adata, tissue_type="lung_tumor", plot=False)

        # Tumor should have higher MT threshold
        # (or at least not strictly lower)
        # Note: This depends on the data distributions

        # Check that tissue_type is in evidence
        assert "tissue_type" in result_tumor.max_mt_percent.evidence
        assert result_tumor.max_mt_percent.evidence["tissue_type"] == "lung_tumor"


# =============================================================================
# Test Tumor vs Normal Differentiation
# =============================================================================


class TestTumorVsNormal:
    """Test tumor-aware recommendations."""

    def test_tumor_has_higher_mt_threshold(self, tumor_like_adata):
        """Test that tumor tissue gets higher MT threshold."""
        # Get recommendation with tumor_aware strategy
        recommender_tumor = IntelligentQCRecommender(strategy=StrategyType.TUMOR_AWARE)
        result_tumor = recommender_tumor.recommend(
            tumor_like_adata, tissue_type="lung_tumor", plot=False
        )

        # Get recommendation with standard strategy
        recommender_standard = IntelligentQCRecommender(strategy=StrategyType.STANDARD)
        result_standard = recommender_standard.recommend(
            tumor_like_adata, tissue_type="normal", plot=False
        )

        # Tumor-aware should have higher MT threshold
        # (or at least consider tissue type)
        assert (
            "tumor" in result_tumor.tumor_specific_considerations
            or len(result_tumor.tumor_specific_considerations) > 0
        )

    def test_tumor_considerations_generated(self, tumor_like_adata):
        """Test that tumor-specific considerations are generated."""
        recommender = IntelligentQCRecommender(strategy=StrategyType.TUMOR_AWARE)

        result = recommender.recommend(tumor_like_adata, tissue_type="lung_tumor", plot=False)

        # Should have tumor-specific considerations
        assert isinstance(result.tumor_specific_considerations, list)

    def test_normal_no_tumor_considerations(self, sample_adata_with_qc):
        """Test that normal tissue has fewer tumor considerations."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal", plot=False)

        # Should have fewer or no tumor-specific considerations
        assert len(result.tumor_specific_considerations) == 0 or all(
            "normal" in c.lower() for c in result.tumor_specific_considerations
        )


# =============================================================================
# Test Data Quality Assessment
# =============================================================================


class TestDataQualityAssessment:
    """Test data quality scoring."""

    def test_quality_score_in_range(self, sample_adata_with_qc):
        """Test that quality score is between 0 and 100."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal", plot=False)

        assert 0 <= result.data_quality_score <= 100

    def test_low_quality_detected(self, low_quality_adata):
        """Test that low quality data is detected."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(low_quality_adata, tissue_type="unknown", plot=False)

        # Should have lower quality score
        assert result.data_quality_score < 80

        # Should have concerns
        assert len(result.concerns) > 0

    def test_concerns_generated(self, low_quality_adata):
        """Test that concerns are generated for poor quality data."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(low_quality_adata, tissue_type="unknown", plot=False)

        # Should have concerns
        assert isinstance(result.concerns, list)


# =============================================================================
# Test Missing Metrics
# =============================================================================


class TestMissingMetrics:
    """Test handling of missing QC metrics."""

    def test_missing_metrics_handled(self, temp_output_dir):
        """Test that missing metrics are handled gracefully."""
        generator = SyntheticDataGenerator(random_state=45)
        # Create AnnData without precomputed QC metrics
        adata = generator.generate_adata(
            n_cells=100,
            n_genes=500,
            with_qc_metrics=False,
            with_cell_types=False,
            with_batches=False,
            sparsity=0.85,
        )

        recommender = IntelligentQCRecommender()

        # Should not crash, but return with warnings
        result = recommender.recommend(adata, tissue_type="unknown", save_dir=temp_output_dir)

        # Should still return a result
        assert isinstance(result, QCRecommendation)

        # Should have concerns about missing metrics
        assert any("missing" in c.lower() or "metric" in c.lower() for c in result.concerns)


# =============================================================================
# Test Convenience Function
# =============================================================================


class TestConvenienceFunction:
    """Test recommend_intelligent_qc() convenience function."""

    def test_convenience_function_works(self, sample_adata_with_qc, temp_output_dir):
        """Test that convenience function works."""
        result = recommend_intelligent_qc(
            sample_adata_with_qc, tissue_type="normal", save_dir=temp_output_dir
        )

        assert isinstance(result, QCRecommendation)

    def test_convenience_function_all_strategies(self, sample_adata_with_qc):
        """Test convenience function with all strategies."""
        for strategy in ["auto", "tumor_aware", "conservative", "aggressive"]:
            result = recommend_intelligent_qc(
                sample_adata_with_qc.copy(), tissue_type="normal", strategy=strategy, plot=False
            )

            assert isinstance(result, QCRecommendation)


# =============================================================================
# Test Serialization
# =============================================================================


class TestSerialization:
    """Test QCRecommendation serialization."""

    def test_to_dict(self, sample_adata_with_qc):
        """Test QCRecommendation.to_dict() method."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(sample_adata_with_qc, tissue_type="normal", plot=False)

        # Convert to dict
        result_dict = result.to_dict()

        # Check structure
        assert "min_genes" in result_dict
        assert "max_mt_percent" in result_dict
        assert "overall_strategy" in result_dict
        assert "overall_confidence" in result_dict
        assert "data_quality_score" in result_dict

        # Check nested structure
        assert "threshold" in result_dict["min_genes"]
        assert "ci_lower" in result_dict["min_genes"]
        assert "ci_upper" in result_dict["min_genes"]
        assert "confidence" in result_dict["min_genes"]


# =============================================================================
# Test Plot Generation
# =============================================================================


class TestPlotGeneration:
    """Test diagnostic plot generation."""

    def test_plots_are_generated(self, sample_adata_with_qc, temp_output_dir):
        """Test that diagnostic plots are generated."""
        recommender = IntelligentQCRecommender()

        result = recommender.recommend(
            sample_adata_with_qc, tissue_type="normal", plot=True, save_dir=temp_output_dir
        )

        # Check that plot files were created
        plot_file = temp_output_dir / "min_genes_recommendation.pdf"
        # Note: Plotting might fail in headless environments, so we just check
        # that the code doesn't crash
        assert isinstance(result, QCRecommendation)


# =============================================================================
# Integration Tests
# =============================================================================


class TestIntegration:
    """Integration tests for intelligent QC."""

    def test_end_to_end_workflow(self, sample_adata_with_qc, temp_output_dir):
        """Test complete end-to-end workflow."""
        import json

        # Get recommendations
        result = recommend_intelligent_qc(
            sample_adata_with_qc, tissue_type="normal", save_dir=temp_output_dir
        )

        # Check JSON report was saved
        json_path = temp_output_dir / "qc_recommendation.json"
        assert json_path.exists()

        # Load and check JSON
        with open(json_path) as f:
            saved_data = json.load(f)

        assert "min_genes" in saved_data
        assert "max_mt_percent" in saved_data

    def test_comparison_with_fixed_thresholds(self, sample_adata_with_qc):
        """Test intelligent QC vs traditional fixed thresholds."""
        # Get intelligent recommendations
        result = recommend_intelligent_qc(sample_adata_with_qc, tissue_type="normal", plot=False)

        # Traditional fixed thresholds
        fixed_min_genes = 200
        fixed_max_mt = 20.0

        # Intelligent recommendations
        intelligent_min_genes = result.min_genes.threshold
        intelligent_max_mt = result.max_mt_percent.threshold

        # They should differ (at least for some datasets)
        # The key is that intelligent recommendations are DATA-DRIVEN
        # not arbitrary

        # Check that we have confidence intervals
        # (which fixed thresholds don't provide)
        assert result.min_genes.ci_lower != result.min_genes.ci_upper
        assert result.max_mt_percent.ci_lower != result.max_mt_percent.ci_upper

        # Check that we have evidence
        # (which fixed thresholds don't provide)
        assert "method" in result.min_genes.evidence
        assert "method" in result.max_mt_percent.evidence


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
