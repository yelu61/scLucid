"""Tests for the lightweight QC/preprocess validation scaffold."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.utils.contracts import Modules, SCHEMA_VERSION, SCLUCID_ROOT, UnsKeys
from scLucid.utils.validation_scaffold import (
    COMPARATIVE_READINESS_LABEL,
    build_qc_preprocess_validation,
    validation_table_to_dataframe,
    write_validation_outputs,
)


def _review(module: str, n_cells: int, n_genes: int) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "module": module,
        "workflow_name": f"test_{module}",
        "steps_executed": ["step_1"],
        "data_shape": {"n_cells": n_cells, "n_genes": n_genes},
        "generated_at": "2026-05-14T00:00:00Z",
        "warnings": [],
        "config": {},
        "contract": {"valid": True},
        "config_lineage": {},
        "artifacts": {},
    }


def _validated_adata() -> AnnData:
    rng = np.random.default_rng(0)
    adata = AnnData(rng.poisson(2, size=(30, 50)).astype(float))
    adata.layers["counts"] = adata.X.copy()
    adata.layers["normalized"] = np.log1p(adata.X)
    adata.raw = adata.copy()
    adata.var["highly_variable"] = [True] * 10 + [False] * 40
    adata.obs["low_quality"] = [True] * 3 + [False] * 27
    adata.obs["predicted_doublet"] = [False] * 28 + [True] * 2
    adata.obsm["X_pca"] = rng.normal(size=(30, 5))
    adata.obsm["X_umap"] = rng.normal(size=(30, 2))
    adata.uns["neighbors"] = {"params": {"n_neighbors": 10}}
    adata.uns[SCLUCID_ROOT] = {
        Modules.QC: {
            UnsKeys.REVIEW_SUMMARY: _review(Modules.QC, adata.n_obs, adata.n_vars),
            "warnings": {"data": ["example warning"]},
        },
        Modules.PREPROCESS: {
            UnsKeys.REVIEW_SUMMARY: _review(
                Modules.PREPROCESS, adata.n_obs, adata.n_vars
            ),
            "hvg_stability": {"available": True},
        },
    }
    return adata


def test_build_qc_preprocess_validation_ready():
    adata = _validated_adata()
    validation = build_qc_preprocess_validation(
        adata,
        run_manifest={
            "workflow": "unit_golden_path",
            "input_shape": {"n_cells": 40, "n_genes": 50},
            "retention_fraction": 0.75,
        },
        dataset_role="unit_baseline",
    )

    assert validation["schema_version"] == "0.1"
    assert validation["scope"] == "qc_preprocess_lightweight"
    assert validation["readiness_status"] == COMPARATIVE_READINESS_LABEL
    assert validation["ready_for_comparative_validation"] is True
    assert validation["qc_metrics"]["retention_fraction"] == 0.75
    assert validation["qc_metrics"]["warning_count"] == 1
    assert validation["preprocess_metrics"]["hvg_count"] == 10
    assert validation["preprocess_metrics"]["hvg_key"] == "highly_variable"
    assert validation["preprocess_metrics"]["canonical_hvg_available"] is True
    assert validation["preprocess_metrics"]["layer_contract"]["counts_present"] is True
    assert validation["preprocess_metrics"]["representation_contract"]["umap_present"] is True
    assert "does not claim scientific superiority" in validation["claim_boundary"]


def test_validation_table_and_outputs(tmp_path):
    validation = build_qc_preprocess_validation(
        _validated_adata(),
        run_manifest={"input_shape": {"n_cells": 30, "n_genes": 50}},
    )
    frame = validation_table_to_dataframe(validation)
    assert {"metric", "value", "status", "interpretation"}.issubset(frame.columns)
    assert "hvg_count" in set(frame["metric"])

    artifacts = write_validation_outputs(validation, tmp_path)
    json_path = tmp_path / "qc_preprocess_validation.json"
    table_path = tmp_path / "qc_preprocess_validation_table.csv"
    assert artifacts == {"json": str(json_path), "table_csv": str(table_path)}
    assert json.loads(json_path.read_text())["scope"] == "qc_preprocess_lightweight"
    table = pd.read_csv(table_path)
    assert "qc_review_summary_valid" in set(table["metric"])


def test_missing_contracts_are_blocking_failures():
    adata = AnnData(np.ones((5, 6)))
    validation = build_qc_preprocess_validation(adata)

    assert validation["ready_for_comparative_validation"] is False
    failed = {row["metric"] for row in validation["blocking_failures"]}
    assert "counts_layer_present" in failed
    assert "qc_review_summary_valid" in failed
    assert "preprocess_review_summary_valid" in failed


def test_hvg_selected_fallback_is_comparison_ready():
    adata = _validated_adata()
    adata.var["highly_variable_selected"] = adata.var["highly_variable"]
    del adata.var["highly_variable"]

    validation = build_qc_preprocess_validation(
        adata,
        run_manifest={"input_shape": {"n_cells": 30, "n_genes": 50}},
    )

    assert validation["preprocess_metrics"]["hvg_count"] == 10
    assert validation["preprocess_metrics"]["hvg_key"] == "highly_variable_selected"
    assert validation["preprocess_metrics"]["canonical_hvg_available"] is False
    assert validation["ready_for_comparative_validation"] is True
