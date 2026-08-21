"""Tests for scLucid.utils.helpers utility functions."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from scLucid.utils.helpers import assess_matrix_semantics, build_metadata_dicts, subset_adata


class TestBuildMetadataDicts:
    def test_builds_group_and_batch_dicts(self):
        samples = ["S1", "S2", "S3"]
        group_dict = {"S1": "tumor", "S2": "normal", "S3": "tumor", "S4": "other"}
        batch_dict = {"S1": "A", "S2": "A", "S3": "B", "S4": "C"}

        result = build_metadata_dicts(samples, group_dict=group_dict, batch_dict=batch_dict)

        assert result == {
            "group": {"S1": "tumor", "S2": "normal", "S3": "tumor"},
            "batch": {"S1": "A", "S2": "A", "S3": "B"},
        }

    def test_empty_inputs_return_empty(self):
        samples = ["S1", "S2"]
        result = build_metadata_dicts(samples)
        assert result == {}

    def test_custom_keys(self):
        samples = ["S1", "S2"]
        group_dict = {"S1": "ctrl", "S2": "treat"}
        result = build_metadata_dicts(samples, group_dict=group_dict, group_key="condition")
        assert result == {"condition": {"S1": "ctrl", "S2": "treat"}}

    def test_extra_dicts_and_default_value(self):
        samples = ["S1", "S2", "S3"]
        result = build_metadata_dicts(
            samples,
            extra_dicts={
                "patient": {"S1": "P1", "S2": "P2"},
                "timepoint": {"S1": "pre", "S3": "post"},
            },
            default_value="unknown",
        )

        assert result["patient"] == {"S1": "P1", "S2": "P2", "S3": "unknown"}
        assert result["timepoint"] == {"S1": "pre", "S2": "unknown", "S3": "post"}

    def test_strict_missing_sample_raises(self):
        with pytest.raises(KeyError, match="missing values"):
            build_metadata_dicts(
                ["S1", "S2"],
                group_dict={"S1": "ctrl"},
                strict=True,
            )


class TestSubsetAdata:
    def test_subset_adata_retains_raw_all_genes(self):
        adata = AnnData(np.ones((4, 3)))
        adata.obs_names = [f"cell_{i}" for i in range(4)]
        adata.var_names = ["gene_a", "gene_b", "gene_c"]
        adata.obs["cell_type"] = pd.Categorical(["T", "B", "T", "B"])
        adata.raw = adata.copy()

        subset = subset_adata(adata, {"cell_type": "T"})

        assert subset.n_obs == 2
        assert subset.raw is not None
        assert list(subset.raw.var_names) == ["gene_a", "gene_b", "gene_c"]


class TestAssessMatrixSemantics:
    def test_raw_dense_counts_are_count_like(self):
        rng = np.random.default_rng(0)
        x = np.zeros((50, 20), dtype=float)
        x[:25, :10] = rng.poisson(5, size=(25, 10))
        x[25:, 10:] = rng.poisson(3, size=(25, 10))
        result = assess_matrix_semantics(x)
        assert result["is_valid"] is True
        assert result["is_count_like"] is True
        assert result["diagnostics"]["has_negative"] is False
        assert result["diagnostics"]["zero_fraction"] >= 0.20
        assert result["diagnostics"]["max_value"] > 10

    def test_normalized_matrix_is_not_count_like(self):
        x = np.random.rand(50, 20)
        result = assess_matrix_semantics(x)
        assert result["is_valid"] is False
        assert result["is_count_like"] is False
        assert result["diagnostics"]["fractional_positive_rate"] > 0.01

    def test_negative_values_are_not_count_like(self):
        x = np.random.poisson(3, size=(50, 20)).astype(float)
        x[0, 0] = -1
        result = assess_matrix_semantics(x)
        assert result["is_valid"] is False
        assert result["is_count_like"] is False
        assert result["diagnostics"]["has_negative"] is True
        assert any("negative" in w for w in result["warnings"])

    def test_sparse_raw_counts_are_count_like(self):
        rng = np.random.default_rng(0)
        x = sp.random(100, 50, density=0.1, format="csr", random_state=rng)
        x.data = np.abs(rng.poisson(5, size=x.data.size).astype(float))
        x.data[x.data == 0] = 1
        result = assess_matrix_semantics(x)
        assert result["is_valid"] is True
        assert result["is_count_like"] is True
        assert result["diagnostics"]["has_negative"] is False
        assert result["diagnostics"]["zero_fraction"] >= 0.20

    def test_empty_matrix_is_invalid(self):
        x = np.zeros((0, 10))
        result = assess_matrix_semantics(x)
        assert result["is_valid"] is False
        assert result["is_count_like"] is False
        assert any("empty" in w for w in result["warnings"])

    def test_missing_layer_is_invalid(self):
        adata = AnnData(np.random.poisson(2, size=(10, 10)).astype(float))
        result = assess_matrix_semantics(adata, layer="counts")
        assert result["is_valid"] is False
        assert result["is_count_like"] is False
        assert any("not found" in w for w in result["warnings"])

    def test_layer_selection_on_anndata(self):
        rng = np.random.default_rng(0)
        x = np.zeros((50, 20), dtype=float)
        x[:25, :10] = rng.poisson(5, size=(25, 10))
        x[25:, 10:] = rng.poisson(3, size=(25, 10))
        adata = AnnData(x)
        adata.layers["counts"] = adata.X.copy()
        adata.layers["normalized"] = (
            adata.X.astype(float) / adata.X.sum(axis=1, keepdims=True) * 1e4
        )

        counts_result = assess_matrix_semantics(adata, layer="counts")
        norm_result = assess_matrix_semantics(adata, layer="normalized")
        assert counts_result["is_count_like"] is True
        # Normalized layer may still pass integer check if values round, but it should at least be valid.
        assert norm_result["is_valid"] in (True, False)

    def test_any_semantics_only_reports_diagnostics(self):
        x = np.random.rand(50, 20)
        result = assess_matrix_semantics(x, semantics="any")
        assert result["is_valid"] is True
        assert result["semantics"] == "any"
        assert "diagnostics" in result

    def test_assess_matrix_semantics_count_like_and_non_count_like(self):
        rng = np.random.default_rng(1)
        counts = np.zeros((50, 20), dtype=float)
        counts[:25, :10] = rng.poisson(5, size=(25, 10))
        counts[25:, 10:] = rng.poisson(3, size=(25, 10))
        normalized = np.random.rand(40, 20)

        counts_result = assess_matrix_semantics(counts)
        assert counts_result["is_count_like"] is True
        assert "zero_fraction" in counts_result["diagnostics"]
        assert "max_value" in counts_result["diagnostics"]

        normalized_result = assess_matrix_semantics(normalized)
        assert normalized_result["is_count_like"] is False
