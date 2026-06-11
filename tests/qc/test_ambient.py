"""Tests for Python-native ambient RNA diagnostics."""

import numpy as np
from anndata import AnnData

from scLucid.qc import (
    diagnose_ambient_rna,
    diagnose_empty_droplets,
    record_ambient_correction_status,
    register_external_ambient_result,
)


def test_diagnose_ambient_rna_returns_risk_summary():
    X = np.ones((100, 30), dtype=float)
    X[:10, :] = 0
    X[:10, 0] = 50
    X[:10, 1] = 30
    adata = AnnData(X=X)
    adata.var_names = [f"Gene{i}" for i in range(30)]

    summary = diagnose_ambient_rna(adata, low_count_quantile=0.1, top_n_genes=5)

    assert summary["available"] is True
    assert summary["diagnostic_only"] is True
    assert summary["method"] == "python_heuristic_low_count_enrichment"
    assert "no expression correction" in summary["method_note"]
    assert summary["calibration_status"] == "heuristic_unvalidated"
    assert summary["risk_score_weights"] == {
        "top_gene_dominance": 0.4,
        "low_count_enrichment": 0.4,
        "enriched_gene_breadth": 0.2,
    }
    assert "not calibrated" in summary["risk_score_note"]
    assert summary["risk_level"] in {"low", "moderate", "high"}
    assert "top_genes" in summary
    assert summary["correction_status"]["corrected"] is False


def test_record_ambient_correction_status_is_used_by_diagnostic():
    adata = AnnData(X=np.ones((50, 10), dtype=float))
    record_ambient_correction_status(
        adata,
        corrected=True,
        backend="cellbender",
        output_layer="cellbender_corrected",
        details={"source": "external_cli"},
    )

    summary = diagnose_ambient_rna(adata, top_n_genes=3)

    assert summary["correction_status"]["corrected"] is True
    assert summary["correction_status"]["backend"] == "cellbender"


def test_register_external_ambient_result_copies_matching_matrix_to_layer():
    adata = AnnData(X=np.ones((5, 4), dtype=float))
    corrected = AnnData(X=np.full((5, 4), 2.0, dtype=float))
    adata.obs_names = corrected.obs_names = [f"cell{i}" for i in range(5)]
    adata.var_names = corrected.var_names = [f"gene{i}" for i in range(4)]

    status = register_external_ambient_result(
        adata,
        backend="cellbender",
        corrected_adata=corrected,
        output_layer="cellbender_corrected",
    )

    assert status["corrected"] is True
    assert status["backend"] == "cellbender"
    assert "cellbender_corrected" in adata.layers
    assert np.asarray(adata.layers["cellbender_corrected"]).mean() == 2.0


def test_diagnose_empty_droplets_records_background_profile():
    X = np.ones((120, 20), dtype=float)
    X[:20, :] = 0
    X[:20, 0] = 10
    X[:20, 1] = 8
    adata = AnnData(X=X)
    adata.var_names = [f"Gene{i}" for i in range(20)]

    summary = diagnose_empty_droplets(adata, min_barcodes=50, top_n_genes=5)

    assert summary["available"] is True
    assert summary["diagnostic_only"] is True
    assert summary["method"] == "python_barcode_rank_background_profile"
    assert "not an EmptyDrops replacement" in summary["method_note"]
    assert "top_background_genes" in summary
    assert "empty_droplet_summary" in adata.uns["sclucid"]["qc"]
