"""Tests for recommendation schema data structures."""

import pytest

from scLucid.recommendation.schema import (
    ParameterRecommendation,
    RecommendationSection,
    WorkflowRecommendations,
)


class TestParameterRecommendation:
    def test_create_basic(self):
        p = ParameterRecommendation(
            name="min_genes",
            value=200,
            method="gmm",
            confidence=0.85,
            rationale="GMM separation threshold.",
        )
        assert p.name == "min_genes"
        assert p.value == 200
        assert p.ci_lower is None
        assert p.ci_upper is None
        assert p.evidence == {}
        assert p.alternatives == []

    def test_create_with_ci(self):
        p = ParameterRecommendation(
            name="max_mt",
            value=20.5,
            method="lognormal_fit",
            confidence=0.9,
            rationale="Log-normal fit cutoff.",
            ci_lower=15.0,
            ci_upper=25.0,
            evidence={"n_outliers": 12},
            alternatives=[18.0, 22.0],
        )
        assert p.ci_lower == 15.0
        assert p.ci_upper == 25.0
        assert p.evidence["n_outliers"] == 12
        assert len(p.alternatives) == 2

    def test_to_dict_roundtrip(self):
        p = ParameterRecommendation(
            name="n_top_genes",
            value=2000,
            method="elbow",
            confidence=0.78,
            rationale="Elbow point of variance explained.",
            ci_lower=1500,
            ci_upper=2500,
            evidence={"variance_explained": 0.87},
            alternatives=[1800, 2200],
        )
        d = p.to_dict()
        assert d["name"] == "n_top_genes"
        assert d["value"] == 2000
        assert d["ci_lower"] == 1500
        assert d["ci_upper"] == 2500
        assert d["evidence"]["variance_explained"] == 0.87


class TestRecommendationSection:
    def test_create_empty(self):
        section = RecommendationSection(name="qc", summary="Empty", confidence=0.0)
        assert section.name == "qc"
        assert section.parameters == []
        assert section.concerns == []
        assert section.notes == []
        assert section.metadata == {}

    def test_get_parameter_exists(self):
        p1 = ParameterRecommendation(
            name="resolution", value=1.2, method="grid", confidence=0.8, rationale="test"
        )
        p2 = ParameterRecommendation(
            name="n_pcs", value=30, method="elbow", confidence=0.9, rationale="test"
        )
        section = RecommendationSection(
            name="clustering",
            summary="test",
            confidence=0.8,
            parameters=[p1, p2],
        )
        found = section.get_parameter("resolution")
        assert found is not None
        assert found.value == 1.2

    def test_get_parameter_missing(self):
        section = RecommendationSection(name="qc", summary="Empty", confidence=0.0)
        assert section.get_parameter("nonexistent") is None

    def test_to_dict(self):
        p = ParameterRecommendation(
            name="min_genes", value=200, method="gmm", confidence=0.85, rationale="test"
        )
        section = RecommendationSection(
            name="qc",
            summary="QC summary",
            confidence=0.85,
            parameters=[p],
            concerns=["Small dataset"],
            notes=["note 1"],
            metadata={"strategy": "standard"},
        )
        d = section.to_dict()
        assert d["name"] == "qc"
        assert len(d["parameters"]) == 1
        assert d["parameters"][0]["name"] == "min_genes"
        assert d["concerns"] == ["Small dataset"]
        assert d["metadata"]["strategy"] == "standard"


class TestWorkflowRecommendations:
    def _make_qc_section(self):
        p = ParameterRecommendation(
            name="min_genes", value=200, method="gmm", confidence=0.85, rationale="test"
        )
        return RecommendationSection(
            name="qc", summary="QC", confidence=0.85, parameters=[p]
        )

    def _make_clustering_section(self, resolution=1.2, n_clusters=8):
        p_res = ParameterRecommendation(
            name="resolution", value=resolution, method="grid", confidence=0.8, rationale="test"
        )
        p_n = ParameterRecommendation(
            name="n_clusters", value=n_clusters, method="derived", confidence=0.8, rationale="test"
        )
        return RecommendationSection(
            name="clustering",
            summary="Clustering",
            confidence=0.8,
            parameters=[p_res, p_n],
        )

    def test_create_and_get_section(self):
        qc = self._make_qc_section()
        wr = WorkflowRecommendations(sections={"qc": qc}, overall_confidence=0.85)
        assert wr.get_section("qc") is qc
        assert wr.get_section("missing") is None

    def test_overall_confidence(self):
        wr = WorkflowRecommendations(sections={}, overall_confidence=0.5)
        assert wr.overall_confidence == 0.5

    def test_to_qc_thresholds(self):
        p_genes = ParameterRecommendation(
            name="min_genes", value=200, method="gmm", confidence=0.85, rationale="test"
        )
        p_mt = ParameterRecommendation(
            name="max_mt_percent", value=20.0, method="lognormal", confidence=0.9, rationale="test"
        )
        section = RecommendationSection(
            name="qc",
            summary="QC",
            confidence=0.85,
            parameters=[p_genes, p_mt],
        )
        wr = WorkflowRecommendations(
            sections={"qc": section}, overall_confidence=0.85
        )
        thresholds = wr.to_qc_thresholds()
        assert thresholds == {"min_genes": 200, "max_mt_percent": 20.0}

    def test_to_qc_thresholds_no_section(self):
        wr = WorkflowRecommendations(sections={}, overall_confidence=0.5)
        assert wr.to_qc_thresholds() == {}

    def test_to_clustering_config_no_section(self):
        wr = WorkflowRecommendations(sections={}, overall_confidence=0.5)
        assert wr.to_clustering_config() is None

    def test_to_clustering_config(self):
        section = self._make_clustering_section(resolution=1.5, n_clusters=10)
        wr = WorkflowRecommendations(
            sections={"clustering": section}, overall_confidence=0.8
        )
        config = wr.to_clustering_config(method="leiden", use_rep="X_pca")
        assert config is not None
        assert config.resolution == 1.5
        assert config.method == "leiden"
        assert config.use_rep == "X_pca"

    def test_to_clustering_config_kmeans(self):
        section = self._make_clustering_section(resolution=1.0, n_clusters=12)
        wr = WorkflowRecommendations(
            sections={"clustering": section}, overall_confidence=0.8
        )
        config = wr.to_clustering_config(method="kmeans", use_rep="X_pca")
        assert config is not None
        assert config.n_clusters == 12

    def test_to_annotation_config_no_section(self):
        wr = WorkflowRecommendations(sections={}, overall_confidence=0.5)
        assert wr.to_annotation_config() is None

    def test_to_tumor_config_no_section(self):
        wr = WorkflowRecommendations(sections={}, overall_confidence=0.5)
        assert wr.to_tumor_config() is None

    def test_to_dict(self):
        qc = self._make_qc_section()
        wr = WorkflowRecommendations(
            sections={"qc": qc},
            overall_confidence=0.85,
            context={"n_cells": 100},
            concerns=["small dataset"],
        )
        d = wr.to_dict()
        assert d["overall_confidence"] == 0.85
        assert "qc" in d["sections"]
        assert d["context"]["n_cells"] == 100
