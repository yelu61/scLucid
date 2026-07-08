"""Tests for the graph stability benchmark runner."""

import numpy as np
import pytest
from anndata import AnnData

from validation.preprocess.run_graph_stability_benchmark import (
    _cluster_leiden,
    _neighbor_overlap,
    _pca_embedding,
)


class TestGraphStabilityBenchmark:
    def test_cluster_leiden_does_not_require_n_clusters(self):
        embedding = np.random.default_rng(0).normal(size=(50, 10))
        labels = _cluster_leiden(embedding, seed=0)
        assert labels.shape == (50,)
        assert len(np.unique(labels)) >= 1

    def test_cluster_leiden_seed_stability(self):
        embedding = np.random.default_rng(0).normal(size=(50, 10))
        labels_a = _cluster_leiden(embedding, seed=42)
        labels_b = _cluster_leiden(embedding, seed=42)
        np.testing.assert_array_equal(labels_a, labels_b)

    def test_neighbor_overlap_between_identical_embeddings(self):
        embedding = np.random.default_rng(0).normal(size=(30, 5))
        overlap = _neighbor_overlap(embedding, embedding, n_neighbors=5)
        assert overlap == pytest.approx(1.0, abs=1e-6)

    def test_pca_embedding_truncates_to_available_dims(self):
        X = np.random.default_rng(0).poisson(1.0, size=(20, 8)).astype(float)
        embedding, ratio = _pca_embedding(X, n_pcs=50, seed=0)
        assert embedding.shape[1] <= 7  # min(n_obs - 1, n_vars - 1)
        assert ratio.shape[0] == embedding.shape[1]
