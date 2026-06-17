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
        API_LAYER_ORDER,
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
    assert UnsKeys.CONFIG_LINEAGE == "config_lineage"
    assert UnsKeys.ARTIFACTS == "artifacts"
    assert STAGE_ORDER == ("qc", "preprocess", "analysis")
    assert API_LAYER_ORDER == ("workflow", "simple_api", "advanced")


def test_contract_spec_is_serializable_and_documents_stages():
    from scLucid.utils.contracts import (
        get_api_layer_spec,
        get_contract_spec,
        get_minimal_workflow_contract,
        get_stage_contract,
    )

    spec = get_contract_spec()
    preprocess = get_stage_contract("preprocess")
    workflow_layer = get_api_layer_spec("workflow")
    workflow_contract = get_minimal_workflow_contract()

    assert spec["schema_version"] == "1.0"
    assert spec["storage_root"] == "sclucid"
    assert spec["api_layers"]["layer_order"] == ["workflow", "simple_api", "advanced"]
    assert spec["minimal_workflow"]["stage_order"] == ["qc", "preprocess", "analysis"]
    assert "review_summary" in spec
    assert spec["canonical_keys"]["uns"]["namespace_metadata"] == "_metadata"
    assert spec["canonical_keys"]["uns"]["config_lineage"] == "config_lineage"
    assert spec["canonical_keys"]["layers"]["counts"] == "counts"
    assert spec["canonical_keys"]["modalities"]["spatial"] == "spatial"
    assert spec["canonical_keys"]["layer_semantics"]["raw_counts"] == "raw_counts"
    assert spec["stages"]["qc"]["name"] == "qc"
    assert preprocess["input_layers"] == ["counts"]
    assert preprocess["output_obsm"] == ["X_pca"]
    assert workflow_layer["primary_entrypoints"][0] == "scLucid.run_pipeline"
    assert workflow_contract["required_stage_namespace_keys"] == [
        "workflow_config",
        "steps_executed",
        "review_summary",
    ]


def test_normalize_review_summary_preserves_review_lists():
    from scLucid.utils.contracts import normalize_review_summary

    review = normalize_review_summary(
        {
            "step_evidence_summary": {"steps": [{"step": "normalization"}]},
            "evidence_bundle": {"evidence_chain": [{"name": "step_evidence_summary"}]},
            "review_action_items": [{"evidence_key": "x"}],
        },
        module="preprocess",
        workflow_name="run_preprocessing",
        steps_executed=["normalization"],
    )

    assert isinstance(review["step_evidence_summary"]["steps"], list)
    assert review["step_evidence_summary"]["steps"][0]["step"] == "normalization"
    assert isinstance(review["evidence_bundle"]["evidence_chain"], list)
    assert isinstance(review["review_action_items"], list)
    assert isinstance(review["data"]["step_evidence_summary"]["steps"], list)


def test_frozen_api_layer_entrypoints_resolve():
    import importlib

    from scLucid.utils.contracts import get_api_layer_spec

    def resolve_dotted(name: str):
        parts = name.split(".")
        module = importlib.import_module(parts[0])
        current = module
        for part in parts[1:]:
            current = getattr(current, part)
        return current

    spec = get_api_layer_spec()
    for layer in ("workflow", "simple_api"):
        for entrypoint in spec["layers"][layer]["primary_entrypoints"]:
            assert callable(resolve_dotted(entrypoint)), entrypoint


def test_frozen_api_layer_examples_exist():
    from pathlib import Path

    from scLucid.utils.contracts import get_api_layer_spec

    repo_root = Path(__file__).resolve().parents[1]
    spec = get_api_layer_spec()
    for layer in spec["layers"].values():
        for artifact in layer["example_artifacts"]:
            assert (repo_root / artifact).exists(), artifact


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
    assert summary["data"]["recommendation_summary"]["available"] is False
    assert summary["data"]["schema_version"] == "1.0"
    assert validate_review_summary_schema(summary, module="qc").valid is True


def test_review_summary_list_of_dict_evidence_is_h5ad_safe(tmp_path):
    from scLucid.utils.contracts import normalize_review_summary

    adata = _adata()
    summary = normalize_review_summary(
        {
            "empty_droplet_summary": {
                "top_background_genes": [
                    {"gene": "mt-Co3", "background_fraction": 0.1},
                    {"gene": "mt-Atp6", "background_fraction": 0.2},
                ]
            }
        },
        module="qc",
        workflow_name="standard",
        adata=adata,
        steps_executed=["diagnose_empty_droplets"],
    )
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["review_summary"] = summary

    path = tmp_path / "review_summary_safe.h5ad"
    adata.write_h5ad(path)

    assert path.exists()


def test_namespace_helpers_add_metadata_and_stage_records():
    from scLucid.utils.contracts import (
        UnsKeys,
        build_config_lineage,
        ensure_sclucid_namespace,
        module_namespace,
        record_artifact,
        record_config_lineage,
    )

    adata = _adata()
    root = ensure_sclucid_namespace(adata)
    qc_ns = module_namespace(adata, "qc")
    lineage = build_config_lineage(
        global_config={"default_species": "human"},
        inherited={"dataset_type": "pbmc_or_blood"},
        stage_config={"sample_key": "sampleID"},
        effective_config={"sample_key": "sampleID", "species": "human"},
    )

    qc_ns[UnsKeys.REVIEW_SUMMARY] = {
        "schema_version": "1.0",
        "module": "qc",
        "workflow_name": "standard",
        "steps_executed": [],
        "data_shape": {"n_cells": adata.n_obs, "n_genes": adata.n_vars},
    }
    record_config_lineage(adata, "qc", lineage)
    record_artifact(adata, "qc", "summary_json", "/tmp/qc_review_summary.json")

    assert root[UnsKeys.NAMESPACE_METADATA]["schema_version"] == "1.0"
    assert qc_ns[UnsKeys.NAMESPACE_METADATA]["module"] == "qc"
    assert qc_ns[UnsKeys.CONFIG_LINEAGE]["effective"]["species"] == "human"
    assert qc_ns[UnsKeys.REVIEW_SUMMARY][UnsKeys.CONFIG_LINEAGE]["stage"]["sample_key"] == (
        "sampleID"
    )
    assert qc_ns[UnsKeys.ARTIFACTS]["summary_json"]["path"].endswith(
        "qc_review_summary.json"
    )


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


def test_register_and_validate_spatial_input_semantics():
    from scLucid.utils.contracts import (
        LayerSemanticKeys,
        ModalityKeys,
        register_anndata_semantics,
        validate_modality_contract,
    )

    adata = _adata()
    adata.obsm["spatial"] = np.random.default_rng(1).normal(size=(adata.n_obs, 2))
    register_anndata_semantics(
        adata,
        modality=ModalityKeys.SPATIAL,
        assay="visium",
        layer_semantics={"counts": LayerSemanticKeys.RAW_COUNTS},
    )

    result = validate_modality_contract(adata, require_spatial=True, require_raw_counts=True)

    assert result.valid is True
    assert result.modality == "spatial"
    assert adata.uns["sclucid"]["input_semantics"]["modality"] == "spatial"
    assert adata.uns["sclucid"]["input_semantics_validation"]["valid"] is True


def test_validate_multimodal_protein_and_atac_contract():
    from scLucid.utils.contracts import (
        ModalityKeys,
        register_anndata_semantics,
        validate_modality_contract,
    )

    adata = _adata()
    adata.obsm["protein_counts"] = np.random.default_rng(2).poisson(1, size=(adata.n_obs, 3))
    adata.layers["atac_counts"] = np.random.default_rng(3).poisson(1, size=adata.shape)
    adata.var["feature_types"] = "Gene Expression"
    register_anndata_semantics(
        adata,
        modality=ModalityKeys.MULTIOME,
        assay="multiome",
        protein_obsm_key="protein_counts",
        atac_layer="atac_counts",
    )

    result = validate_modality_contract(
        adata,
        require_protein=True,
        require_atac=True,
    )

    assert result.valid is True
    assert result.detected["protein_obsm_key"] == "protein_counts"
    assert result.detected["atac_layer"] == "atac_counts"


def test_validate_modality_contract_reports_missing_spatial_coordinates():
    from scLucid.utils.contracts import ModalityKeys, register_anndata_semantics, validate_modality_contract

    adata = _adata()
    register_anndata_semantics(adata, modality=ModalityKeys.SPATIAL, assay="visium")

    result = validate_modality_contract(adata, require_spatial=True)

    assert result.valid is False
    assert any("Spatial modality" in error for error in result.errors)


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
