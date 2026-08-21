from __future__ import annotations

import warnings

import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy import sparse

import scLucid as scl
from validation.preprocess.run_controlled_preprocess_contract_benchmark import (
    _controlled_umi_adata,
)


def _context() -> scl.ProjectContext:
    return scl.ProjectContext(
        dataset_type="cell_line",
        assay="scrna",
        input_provenance="filtered_counts",
        sample_key="sample",
        batch_key="batch",
        condition_key="condition",
        cell_type_key="cell_type",
        is_multi_sample=True,
    )


def test_sparse_public_policy_records_bounded_pca_densification_without_warning():
    adata = _controlled_umi_adata(seed=29)
    card = scl.recommend_preprocess_policy(adata, _context(), consumer="exploration")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        evidence = scl.apply_preprocess_policy(adata, card.policy)

    implicit_warnings = [
        str(item.message)
        for item in caught
        if "zero-centering a sparse array/matrix densifies it" in str(item.message)
    ]
    assert implicit_warnings == []

    result = evidence.adata
    assert sparse.issparse(result.layers["counts"])
    assert sparse.issparse(result.layers["normalized_full"])
    assert sparse.issparse(result.X)
    assert result.raw is not None and sparse.issparse(result.raw.X)

    contract = evidence.result["discovery_temporary_contract"]
    assert contract["densification_occurred"] is True
    assert contract["scope"] == "temporary_discovery_feature_matrix"
    assert contract["temporary_shape"] == [
        result.n_obs,
        int(result.var["discovery_feature"].sum()),
    ]
    assert contract["estimated_peak_bytes"] >= contract["dense_matrix_bytes"] > 0
    assert contract["persistent"] is False
    assert contract["consumer"] == "PCA_and_neighbor_graph_only"
    assert contract["expression_inference_eligible"] is False


def test_explicit_bounded_dense_pca_preserves_legacy_geometry():
    adata = _controlled_umi_adata(seed=29)
    card = scl.recommend_preprocess_policy(adata, _context(), consumer="exploration")
    result = scl.apply_preprocess_policy(adata, card.policy).adata

    mask = result.var["discovery_feature"].to_numpy(bool)
    legacy = AnnData(
        X=result.layers["normalized_full"][:, mask].toarray(),
        obs=result.obs.copy(),
        var=result.var.loc[mask].copy(),
    )
    sc.pp.scale(legacy, max_value=10.0, zero_center=True)
    n_comps = min(50, legacy.n_obs - 1, legacy.n_vars - 1)
    sc.tl.pca(legacy, n_comps=n_comps, svd_solver="arpack")

    observed_gram = np.asarray(result.obsm["X_pca"]) @ np.asarray(result.obsm["X_pca"]).T
    legacy_gram = np.asarray(legacy.obsm["X_pca"]) @ np.asarray(legacy.obsm["X_pca"]).T
    assert np.allclose(observed_gram, legacy_gram, atol=1e-5, rtol=1e-5)
