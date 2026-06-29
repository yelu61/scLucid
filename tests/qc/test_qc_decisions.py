"""Tests for evidence-based QC decision schema."""

import numpy as np
from anndata import AnnData

from scLucid.qc.decisions import (
    QC_DECISION_SCHEMA_VERSION,
    build_qc_decisions,
    score_qc_gene_panels,
)


def _decision_adata() -> AnnData:
    genes = ["HBA1", "HBB", "PPBP", "PF4", "FOS", "JUN", "BAX", "CASP3", "GAPDH"]
    X = np.ones((5, len(genes)), dtype=float)
    X[1, [0, 1]] = 30.0
    X[2, [2, 3]] = 25.0
    X[3, [4, 5]] = 35.0
    X[4, [6, 7]] = 40.0
    adata = AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.obs["pct_counts_mt"] = [5.0, 8.0, 10.0, 30.0, 45.0]
    adata.obs["outlier_count"] = [False, False, False, False, True]
    adata.obs["outlier_min_genes"] = [False, False, False, False, True]
    adata.obs["outlier_mt"] = [False, False, False, True, True]
    adata.obs["outlier_qc_metrics"] = [False, False, False, False, True]
    adata.obs["predicted_doublet"] = [False, False, True, False, False]
    return adata


def test_score_qc_gene_panels_writes_contamination_and_stress_scores():
    adata = _decision_adata()

    summary = score_qc_gene_panels(adata)

    assert "hemoglobin_score" in adata.obs
    assert "platelet_score" in adata.obs
    assert "stress_score" in adata.obs
    assert summary["schema_version"] == "qc_gene_panel_scores_v1"
    assert adata.obs.loc["cell_1", "hemoglobin_score"] > adata.obs.loc["cell_0", "hemoglobin_score"]
    assert adata.obs.loc["cell_3", "stress_score"] > adata.obs.loc["cell_0", "stress_score"]


def test_build_qc_decisions_uses_joint_evidence_and_preserves_tumor_high_mt_for_review():
    adata = _decision_adata()

    summary = build_qc_decisions(adata, tissue_type="pdac_tumor", policy="conservative")

    assert summary["schema_version"] == QC_DECISION_SCHEMA_VERSION
    assert {"qc_decision", "qc_reason", "qc_confidence", "qc_review_required"}.issubset(
        adata.obs.columns
    )
    assert "qc_remove" in adata.obs
    assert adata.obs.loc["cell_3", "qc_decision"] == "review"
    assert "high_mt" in adata.obs.loc["cell_3", "qc_reason"]
    assert adata.obs.loc["cell_4", "qc_decision"] == "remove"
    assert bool(adata.obs.loc["cell_4", "qc_remove"])
    assert summary["decision_counts"]["remove"] == 1


def test_build_qc_decisions_can_run_without_gene_scoring():
    adata = _decision_adata()

    summary = build_qc_decisions(
        adata,
        tissue_type="blood",
        policy="screening",
        score_panels=False,
    )

    assert summary["n_cells"] == adata.n_obs
    assert "qc_decision_summary" in adata.uns["sclucid"]["qc"]
