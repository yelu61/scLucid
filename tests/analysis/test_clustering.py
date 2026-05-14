"""
Tests for the analysis clustering module.

Tests clustering, practical resolution review, and cluster merging.
"""

import sys

import pytest
import scanpy as sc

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.analysis.clustering import (
    cluster_cells,
    merge_clusters,
    run_clustering_review,
)
from scLucid.analysis.config import ClusteringConfig, MergeClustersConfig


@pytest.fixture
def preprocessed_adata(minimal_adata):
    """Provide preprocessed data for clustering tests."""
    adata = minimal_adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("clustering", {})
    adata.uns["sclucid"]["analysis"]["clustering"].setdefault("resolution_search", {})
    return adata


@pytest.mark.integration
class TestClustering:
    """Test suite for clustering functionality."""

    def test_leiden_clustering(self, preprocessed_adata):
        """Test basic Leiden clustering."""
        config = ClusteringConfig(method="leiden", resolution=1.0, plot=False)
        result = cluster_cells(preprocessed_adata, config=config)

        # Check that clusters were assigned
        assert "leiden_clusters" in result.obs.columns
        assert result.obs["leiden_clusters"].nunique() > 1

    def test_louvain_clustering(self, preprocessed_adata):
        """Test Louvain clustering."""
        pytest.importorskip("louvain")
        config = ClusteringConfig(method="louvain", resolution=0.8, plot=False)
        result = cluster_cells(preprocessed_adata, config=config)

        assert "louvain_clusters" in result.obs.columns
        assert result.obs["louvain_clusters"].nunique() > 1

    def test_kmeans_clustering(self, preprocessed_adata):
        """Test K-means clustering."""
        config = ClusteringConfig(method="kmeans", n_clusters=4, plot=False)
        result = cluster_cells(preprocessed_adata, config=config)

        assert "kmeans_clusters" in result.obs.columns
        assert result.obs["kmeans_clusters"].nunique() == 4

    def test_hdbscan_clustering(self, preprocessed_adata):
        """Test HDBSCAN clustering if available."""
        pytest.importorskip("hdbscan")

        config = ClusteringConfig(method="hdbscan", plot=False)
        result = cluster_cells(preprocessed_adata, config=config)

        assert "hdbscan_clusters" in result.obs.columns
        # HDBSCAN may assign -1 to noise points
        unique_labels = result.obs["hdbscan_clusters"].unique()
        assert len(unique_labels[unique_labels != -1]) > 0

    def test_clustering_with_custom_key(self, preprocessed_adata):
        """Test clustering with custom key_added."""
        config = ClusteringConfig(method="leiden", key_added="my_clusters", plot=False)
        result = cluster_cells(preprocessed_adata, config=config)

        assert "my_clusters" in result.obs.columns

    def test_clustering_different_resolutions(self, preprocessed_adata):
        """Test that different resolutions give different numbers of clusters."""
        config_low = ClusteringConfig(method="leiden", resolution=0.2, plot=False)
        config_high = ClusteringConfig(method="leiden", resolution=2.0, plot=False)

        result_low = cluster_cells(preprocessed_adata.copy(), config=config_low)
        result_high = cluster_cells(preprocessed_adata.copy(), config=config_high)

        n_low = result_low.obs["leiden_clusters"].nunique()
        n_high = result_high.obs["leiden_clusters"].nunique()

        # Higher resolution should give more clusters
        assert n_high >= n_low

    def test_run_clustering_review_stores_practical_evidence(self, preprocessed_adata):
        """Practical resolution review should summarize candidates and store artifacts."""
        review = run_clustering_review(
            preprocessed_adata,
            resolutions=[0.5, 0.8],
            n_top_markers=5,
            min_cluster_cells=2,
        )

        assert {"resolution", "cluster_key", "n_clusters", "interpretability_score"}.issubset(
            review.columns
        )
        clustering_ns = preprocessed_adata.uns["sclucid"]["analysis"]["clustering"]
        assert "clustering_review_summary" in clustering_ns
        assert clustering_ns["clustering_review_summary"]["recommended_resolution"] in [0.5, 0.8]


@pytest.mark.integration
class TestClusterMerging:
    """Test cluster merging functionality."""

    def test_merge_clusters_marker_overlap(self, preprocessed_adata):
        """Test merging based on marker overlap."""
        # First cluster
        cluster_config = ClusteringConfig(method="leiden", resolution=1.5, plot=False)
        adata = cluster_cells(preprocessed_adata, config=cluster_config)

        # Then merge
        merge_config = MergeClustersConfig(
            cluster_key="leiden_clusters",
            method="marker_overlap",
            similarity_threshold=0.5,
        )

        result = merge_clusters(adata, config=merge_config)

        # Check that merged clusters column was added
        assert "leiden_clusters_merged" in result.obs.columns

    def test_merge_clusters_expression_correlation(self, preprocessed_adata):
        """Test merging based on expression correlation."""
        # First cluster
        cluster_config = ClusteringConfig(method="leiden", resolution=1.5, plot=False)
        adata = cluster_cells(preprocessed_adata, config=cluster_config)

        # Then merge
        merge_config = MergeClustersConfig(
            cluster_key="leiden_clusters",
            method="expression_correlation",
            similarity_threshold=0.8,
        )

        result = merge_clusters(adata, config=merge_config)

        assert "leiden_clusters_merged" in result.obs.columns


@pytest.mark.integration
class TestClusteringConfigValidation:
    """Test configuration validation."""

    def test_invalid_resolution(self):
        """Test that negative resolution raises error."""
        with pytest.raises(ValueError):
            ClusteringConfig(resolution=-0.5)

    def test_invalid_n_clusters(self):
        """Test that invalid n_clusters raises error."""
        with pytest.raises(ValueError):
            ClusteringConfig(method="kmeans", n_clusters=1)  # Need at least 2

    def test_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError):
            ClusteringConfig(method="invalid_method")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
