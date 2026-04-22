"""Tests for analysis configuration classes."""

import pytest
from scLucid.analysis.config import (
    AnalysisWorkflowConfig,
    ClusteringConfig,
    AnnotationConfig,
    DifferentialConfig,
    EnrichmentConfig,
)


class TestAnalysisWorkflowConfig:
    """Test AnalysisWorkflowConfig factory methods."""

    def test_from_simple_dict_clustering_params(self):
        """Extract clustering_* prefixed keys into nested config."""
        config = AnalysisWorkflowConfig.from_simple_dict({
            "clustering_method": "louvain",
            "clustering_resolution": 0.8,
            "clustering_n_clusters": 5,
        })
        assert config.clustering.method == "louvain"
        assert config.clustering.resolution == 0.8
        assert config.clustering.n_clusters == 5

    def test_from_simple_dict_annotation_params(self):
        """Extract annotation_* prefixed keys into nested config."""
        config = AnalysisWorkflowConfig.from_simple_dict({
            "annotation_cluster_key": "clusters",
            "annotation_final_method": "celltypist",
            "annotation_run_scoring": False,
        })
        assert config.annotation.cluster_key == "clusters"
        assert config.annotation.final_method == "celltypist"
        assert config.annotation.run_scoring is False

    def test_from_simple_dict_de_params(self):
        """Extract de_* prefixed keys into nested config."""
        config = AnalysisWorkflowConfig.from_simple_dict({
            "de_groupby": "clusters",
            "de_method": "t-test",
        })
        assert config.de.groupby == "clusters"
        assert config.de.method == "t-test"

    def test_from_simple_dict_results_dir_backward_compat(self):
        """results_dir is mapped to save_dir for backward compatibility."""
        config = AnalysisWorkflowConfig.from_simple_dict({
            "results_dir": "./old_results",
        })
        assert config.save_dir == "./old_results"

    def test_from_simple_dict_unknown_keys_warned(self):
        """Unknown keys should be warned but not break construction."""
        config = AnalysisWorkflowConfig.from_simple_dict({
            "clustering_method": "leiden",
            "unknown_param": 123,
        })
        assert config.clustering.method == "leiden"

    def test_quick_defaults(self):
        """quick() factory uses sensible defaults."""
        config = AnalysisWorkflowConfig.quick()
        assert config.clustering.method == "leiden"
        assert config.clustering.resolution == 1.0
        assert config.annotation is not None

    def test_quick_custom_clustering(self):
        """quick() accepts clustering overrides."""
        config = AnalysisWorkflowConfig.quick(
            clustering_method="kmeans",
            resolution=2.0,
        )
        assert config.clustering.method == "kmeans"
        assert config.clustering.resolution == 2.0

    def test_quick_no_annotation(self):
        """quick() can skip annotation step."""
        config = AnalysisWorkflowConfig.quick(run_annotation=False)
        assert config.annotation is None

    def test_quick_extra_kwargs(self):
        """quick() passes extra kwargs to workflow config."""
        config = AnalysisWorkflowConfig.quick(save_dir="./out", n_jobs=4)
        assert config.save_dir == "./out"
        assert config.n_jobs == 4

    def test_to_dict_roundtrip(self):
        """Serialization preserves nested structure."""
        config = AnalysisWorkflowConfig.quick()
        d = config.to_dict()
        assert "clustering" in d
        assert "annotation" in d
        assert d["clustering"]["method"] == "leiden"


class TestDifferentialConfig:
    """Test DifferentialConfig validation."""

    def test_default_groupby(self):
        config = DifferentialConfig()
        assert config.groupby == "leiden_clusters"

    def test_invalid_method_rejected(self):
        with pytest.raises(ValueError):
            DifferentialConfig(method="invalid")

    def test_pval_cutoff_range(self):
        with pytest.raises(ValueError):
            DifferentialConfig(pval_cutoff=1.5)


class TestFilterMarkersConfig:
    """Test FilterMarkersConfig validation."""

    def test_defaults(self):
        from scLucid.analysis.config import FilterMarkersConfig
        config = FilterMarkersConfig()
        assert config.min_log2fc == 1.0
        assert config.max_padj == 0.05
        assert config.keep_top_n == 100

    def test_invalid_sort_by(self):
        from scLucid.analysis.config import FilterMarkersConfig
        with pytest.raises(ValueError):
            FilterMarkersConfig(sort_by="invalid")
