from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

import scLucid as scl


def _count_adata(*, n_cells: int = 80, n_genes: int = 160) -> AnnData:
    rng = np.random.default_rng(11)
    adata = AnnData(rng.poisson(1.2, size=(n_cells, n_genes)).astype(np.int32))
    adata.obs_names = [f"cell_{idx}" for idx in range(n_cells)]
    adata.var_names = [f"G{idx}" for idx in range(n_genes)]
    return adata


def test_qc_review_is_read_only_and_returns_public_policy_types():
    adata = _count_adata(n_cells=120, n_genes=240)
    adata.obs["sample"] = ["S1"] * 60 + ["S2"] * 60
    before_obs = adata.obs.copy(deep=True)
    before_uns = dict(adata.uns)

    card = scl.recommend_qc_policy(
        adata,
        scl.ProjectContext(
            sample_key="sample",
            is_multi_sample=True,
            input_provenance="filtered_counts",
        ),
    )

    assert isinstance(card, scl.DecisionCard)
    assert isinstance(card.policy, scl.QCPolicy)
    pd.testing.assert_frame_equal(adata.obs, before_obs)
    assert adata.uns == before_uns
    assert card.policy.evidence_heads["ambient"]["status"] == "NOT_EVALUABLE"
    assert "data_quality_score" not in adata.obs


def test_qc_blocks_jointly_catastrophic_sample_and_refuses_application():
    adata = _count_adata(n_cells=80, n_genes=240)
    adata.obs["sample"] = ["good"] * 40 + ["lin_like_failure"] * 40
    adata.obs["n_genes_by_counts"] = [1200.0] * 40 + [45.0] * 40
    adata.obs["total_counts"] = [4000.0] * 40 + [80.0] * 40
    adata.obs["pct_counts_mt"] = [5.0] * 40 + [72.0] * 40
    adata.obs["pct_counts_in_top_20_genes"] = [35.0] * 40 + [97.0] * 40
    card = scl.recommend_qc_policy(
        adata,
        context=scl.ProjectContext(sample_key="sample", is_multi_sample=True),
    )

    assert card.status == "BLOCKED"
    assert card.affected["blocked_samples"] == ["lin_like_failure"]
    with pytest.raises(RuntimeError, match="Blocked QC policies"):
        scl.apply_qc_policy(adata, card.policy)


def test_qc_fails_closed_when_multisample_key_is_not_resolvable():
    adata = _count_adata()
    card = scl.recommend_qc_policy(
        adata,
        context=scl.ProjectContext(is_multi_sample=True),
    )
    assert card.status == "BLOCKED"
    assert any("sample_key" in item for item in card.policy.blockers)


def test_relative_sample_outlier_requires_review_but_is_not_catastrophic():
    adata = _count_adata(n_cells=100, n_genes=240)
    adata.obs["sample"] = np.repeat([f"S{i}" for i in range(10)], 10)
    adata.obs["n_genes_by_counts"] = 1200.0
    adata.obs["total_counts"] = 4000.0
    adata.obs["pct_counts_mt"] = 5.0
    adata.obs["pct_counts_in_top_20_genes"] = 35.0
    adata.obs.loc[adata.obs["sample"] == "S9", "pct_counts_mt"] = 10.0
    adata.obs.loc[adata.obs["sample"] == "S9", "pct_counts_in_top_20_genes"] = 45.0

    card = scl.recommend_qc_policy(
        adata,
        context=scl.ProjectContext(sample_key="sample", is_multi_sample=True),
    )

    assert card.status == "REVIEW"
    assert card.affected["blocked_samples"] == []
    assert card.affected["review_samples"] == ["S9"]


def test_qc_application_requires_unchanged_input_fingerprint():
    adata = _count_adata(n_cells=120, n_genes=240)
    card = scl.recommend_qc_policy(adata, context=scl.ProjectContext())
    adata.X[0, 0] += 1
    with pytest.raises(ValueError, match="fingerprint"):
        scl.apply_qc_policy(adata, card.policy)


def test_preprocess_review_blocks_transformed_input_without_counts():
    adata = _count_adata()
    adata.X = np.log1p(np.asarray(adata.X, dtype=float))
    card = scl.recommend_preprocess_policy(adata, scl.ProjectContext())
    assert card.status == "BLOCKED"
    assert "count matrix" in card.policy.blockers[0]


def test_preprocess_integration_is_blocked_by_condition_confounding():
    adata = _count_adata()
    adata.obs["batch"] = ["B1"] * 40 + ["B2"] * 40
    adata.obs["condition"] = ["control"] * 40 + ["treated"] * 40
    context = scl.ProjectContext(batch_key="batch", condition_key="condition")

    card = scl.recommend_preprocess_policy(adata, context, consumer="integration")

    assert card.status == "BLOCKED"
    assert card.details["integration_review"]["status"] == "BLOCKED"
    assert card.policy.run_integration is False


def test_preprocess_application_preserves_four_space_contract(tmp_path):
    adata = _count_adata()
    adata.obs["sample"] = ["S1"] * 40 + ["S2"] * 40
    context = scl.ProjectContext(sample_key="sample", is_multi_sample=True)
    card = scl.recommend_preprocess_policy(adata, context)

    evidence = scl.apply_preprocess_policy(adata, card.policy)
    result = evidence.adata

    assert isinstance(evidence, scl.RunEvidence)
    assert "counts" in result.layers
    assert "normalized_full" in result.layers
    assert result.raw is not None and result.raw.n_vars == adata.n_vars
    assert "discovery_feature" in result.var
    assert "X_pca" in result.obsm
    assert "X_umap" in result.obsm
    assert result.n_vars == adata.n_vars
    contract = result.uns["sclucid"]["preprocess"]["representation_contract"]
    assert contract["formal_count_model_source"] == "layers[counts]"
    assert contract["marker_program_source"] == "layers[normalized_full]"
    assert contract["integrated_rep"] == "not_selected"
    output = tmp_path / "policy_result.h5ad"
    result.write_h5ad(output)
    roundtrip = read_h5ad(output)
    assert "normalized_full" in roundtrip.layers
    assert "X_pca" in roundtrip.obsm


def test_legacy_qc_recommendation_shape_remains_available():
    adata = _count_adata(n_cells=40, n_genes=120)
    legacy = scl.recommend_qc_policy(adata)
    assert isinstance(legacy, dict)
    assert legacy["schema_version"] == "qc_policy_bundle_v1"
