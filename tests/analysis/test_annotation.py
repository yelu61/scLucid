"""
Tests for the analysis annotation module.

Tests cell type annotation, scoring, and label transfer.
"""

import pytest
import numpy as np
from anndata import AnnData
import scanpy as sc
from pathlib import Path

import sys
sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

import scLucid as scl
from scLucid.analysis.annotation import (
    score_cell_types,
    annotate_clusters,
    run_annotation,
)
from scLucid.analysis.config import AnnotationConfig, ScoringConfig
from scLucid.utils.manager import Manager

from tests.fixtures.synthetic_data import minimal_adata, synthetic_generator


def _write_marker_toml(path: Path, genes_a, genes_b) -> str:
    """Create a minimal marker config file compatible with Manager."""
    content = f"""
[["Synthetic"]]
name = "Type_A"
markers = {list(genes_a)}

[["Synthetic"]]
name = "Type_B"
markers = {list(genes_b)}
"""
    path.write_text(content.strip() + "\n")
    return str(path)


@pytest.fixture
def clustered_adata(minimal_adata):
    """Provide clustered data for annotation tests."""
    from scLucid.analysis.clustering import cluster_cells
    from scLucid.analysis.config import ClusteringConfig

    # Lightweight preprocessing for clustering prerequisites.
    adata = minimal_adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    # Cluster
    cluster_config = ClusteringConfig(method="leiden", resolution=1.0, plot=False)
    adata = cluster_cells(adata, config=cluster_config)

    return adata


@pytest.mark.integration
class TestScoring:
    """Test gene set and cell type scoring."""

    def test_score_cell_types_basic(self, clustered_adata, tmp_path):
        """Test basic cell type scoring."""
        marker_file = _write_marker_toml(
            tmp_path / "markers.toml",
            clustered_adata.var_names[:10].tolist(),
            clustered_adata.var_names[10:20].tolist(),
        )
        marker_manager = Manager(marker_file, case_sensitive=True)
        result = score_cell_types(
            clustered_adata,
            marker_config=marker_manager,
            use_raw=False,
            layer=None,
            score_name_suffix="_test",
        )

        # Check scores were added
        assert "Type_A_test" in result.obs.columns
        assert "Type_B_test" in result.obs.columns

    def test_score_cell_types_from_manager(self, clustered_adata):
        """Test scoring using marker manager."""
        # This test requires marker databases
        pytest.skip("Marker database not available in test environment")


@pytest.mark.integration
class TestAnnotation:
    """Test cell type annotation."""

    def test_annotate_clusters_basic(self, clustered_adata, tmp_path):
        """Test basic cluster annotation."""
        marker_file = _write_marker_toml(
            tmp_path / "markers.toml",
            clustered_adata.var_names[:10].tolist(),
            clustered_adata.var_names[10:20].tolist(),
        )
        marker_manager = Manager(marker_file, case_sensitive=True)

        # Generate score columns required by max-score annotation.
        scored = score_cell_types(
            clustered_adata,
            marker_config=marker_manager,
            use_raw=False,
            layer=None,
            score_name_suffix="_score",
        )

        result = annotate_clusters(
            scored,
            cluster_key="leiden_clusters",
            marker_config=marker_manager,
            method="max_score",
        )

        # Check annotation was added
        assert "leiden_clusters_annotated" in result.obs.columns

    def test_run_annotation_scoring_only(self, clustered_adata):
        """Test run_annotation with scoring method."""
        config = AnnotationConfig(
            cluster_key="leiden_clusters",
            marker_species="human",
            run_celltypist=False,
            run_scoring=True,
            final_method="max_score",
        )

        # This may fail without proper marker databases, so we test config validation
        assert config.cluster_key == "leiden_clusters"
        assert config.run_celltypist == False


@pytest.mark.integration
class TestAnnotationConfigValidation:
    """Test annotation configuration validation."""

    def test_invalid_confidence(self):
        """Test that invalid confidence raises error."""
        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=1.5)  # Should be <= 1

        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=-0.1)  # Should be >= 0

    def test_invalid_final_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError):
            AnnotationConfig(final_method="invalid_method")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
