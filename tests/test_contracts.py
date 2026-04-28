"""Tests for canonical scLucid workflow contracts."""

import numpy as np
from anndata import AnnData


def _adata(n_obs=20, n_vars=30):
    counts = np.random.default_rng(0).poisson(1, size=(n_obs, n_vars)).astype(int)
    adata = AnnData(X=counts)
    adata.layers["counts"] = counts.copy()
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    return adata


def test_contract_constants_define_core_keys():
    from scLucid.utils.contracts import (
        LayerKeys,
        Modules,
        ObsKeys,
        ObsmKeys,
        STAGE_ORDER,
        UnsKeys,
    )

    assert LayerKeys.COUNTS == "counts"
    assert LayerKeys.NORMALIZED == "normalized"
    assert ObsmKeys.PCA == "X_pca"
    assert ObsKeys.QC_N_GENES == "n_genes_by_counts"
    assert Modules.PREPROCESS == "preprocess"
    assert UnsKeys.REVIEW_SUMMARY == "review_summary"
    assert STAGE_ORDER == ("qc", "preprocess", "analysis")


def test_contract_spec_is_serializable_and_documents_stages():
    from scLucid.utils.contracts import get_contract_spec, get_stage_contract

    spec = get_contract_spec()
    preprocess = get_stage_contract("preprocess")

    assert spec["schema_version"] == "1.0"
    assert spec["storage_root"] == "sclucid"
    assert "review_summary" in spec
    assert spec["canonical_keys"]["layers"]["counts"] == "counts"
    assert spec["stages"]["qc"]["name"] == "qc"
    assert preprocess["input_layers"] == ["counts"]
    assert preprocess["output_obsm"] == ["X_pca"]


def test_stage_contract_reports_missing_preprocess_inputs():
    from scLucid.utils.contracts import validate_stage_contract

    adata = _adata()
    result = validate_stage_contract(adata, "preprocess", when="input")

    assert result.valid is False
    assert any('adata.uns["sclucid"]["qc"]' in error for error in result.errors)


def test_stage_contract_accepts_preprocess_output():
    from scLucid.utils.contracts import normalize_review_summary, validate_stage_contract

    adata = _adata()
    adata.layers["normalized"] = adata.X.copy()
    adata.obsm["X_pca"] = np.random.default_rng(1).normal(size=(adata.n_obs, 5))
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["review_summary"] = (
        normalize_review_summary(
            {"steps_executed": ["normalization", "pca"]},
            module="preprocess",
            workflow_name="standard",
            adata=adata,
            steps_executed=["normalization", "pca"],
        )
    )

    result = validate_stage_contract(adata, "preprocess", when="output")

    assert result.valid is True


def test_review_summary_schema_is_backwards_compatible():
    from scLucid.utils.contracts import normalize_review_summary, validate_review_summary_schema

    adata = _adata()
    summary = normalize_review_summary(
        {"steps_executed": ["qc_metrics"], "recommendation_summary": {"available": False}},
        module="qc",
        workflow_name="standard",
        adata=adata,
        steps_executed=["qc_metrics"],
        config={"sample_key": "sampleID"},
    )

    assert summary["recommendation_summary"]["available"] is False
    assert summary["schema_version"] == "1.0"
    assert summary["module"] == "qc"
    assert summary["data_shape"] == {"n_cells": adata.n_obs, "n_genes": adata.n_vars}
    assert validate_review_summary_schema(summary, module="qc").valid is True


def test_review_summary_schema_rejects_invalid_core_types():
    from scLucid.utils.contracts import validate_review_summary_schema

    result = validate_review_summary_schema(
        {
            "schema_version": "1.0",
            "module": "qc",
            "workflow_name": "standard",
            "steps_executed": "qc_metrics",
            "data_shape": {"n_cells": "20", "n_genes": 30},
        },
        module="qc",
    )

    assert result.valid is False
    assert any("steps_executed" in error for error in result.errors)
    assert any("n_cells" in error for error in result.errors)


def test_validate_all_stage_contracts_returns_per_stage_results():
    from scLucid.utils.contracts import validate_all_stage_contracts

    adata = _adata()
    adata.obs["n_genes_by_counts"] = 10
    adata.obs["total_counts"] = 100
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["review_summary"] = {
        "schema_version": "1.0",
        "module": "qc",
        "workflow_name": "standard",
        "steps_executed": ["qc_metrics"],
        "data_shape": {"n_cells": adata.n_obs, "n_genes": adata.n_vars},
    }

    results = validate_all_stage_contracts(adata, stages=["qc", "preprocess"], when="output")

    assert results["qc"].valid is True
    assert results["preprocess"].valid is False


def test_run_pipeline_records_stage_contract(monkeypatch):
    import scLucid as scl

    adata = _adata()

    def fake_qc(input_adata, **kwargs):
        input_adata.obs["n_genes_by_counts"] = 10
        input_adata.obs["total_counts"] = 100
        input_adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["review_summary"] = {
            "schema_version": "1.0",
            "module": "qc",
            "workflow_name": "standard",
            "steps_executed": ["qc_metrics"],
            "data_shape": {"n_cells": input_adata.n_obs, "n_genes": input_adata.n_vars},
        }
        return input_adata

    monkeypatch.setattr(scl, "run_standard_qc", fake_qc)

    out = scl.run_pipeline(adata, stages=["qc"], show_progress=False)

    contract = out.uns["sclucid"]["qc"]["contract"]["output_validation"]
    assert contract["valid"] is True
    assert out.uns["sclucid"]["pipeline_context"]["config_lineage"]["precedence"][0] == (
        "explicit stage config"
    )
