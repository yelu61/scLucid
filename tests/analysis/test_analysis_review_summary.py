"""Tests for benchmark-grade analysis review summaries."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from scLucid.analysis import (
    ANALYSIS_REQUIRED_REVIEW_SECTIONS,
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    build_annotation_consensus,
    build_posthoc_qc_review_summary,
    get_analysis_module_contract,
    run_annotation_evidence,
    run_malignancy_interpretation,
    run_standard_analysis,
    summarize_analysis_review_summary,
    validate_analysis_module_completeness,
    validate_analysis_review_summary,
)
from scLucid.analysis.trace import (
    build_analysis_readiness_assessment,
    build_analysis_review_action_items,
)


def _rows(value):
    return list(value.values()) if isinstance(value, dict) else value


def _make_preprocessed_adata(n_obs=120, n_vars=80):
    import anndata

    rng = np.random.default_rng(7)
    counts = rng.poisson(3, size=(n_obs, n_vars)).astype(np.float32)
    counts[: n_obs // 2, :6] += 8
    counts[n_obs // 2 :, 6:12] += 8
    adata = anndata.AnnData(X=counts)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    adata.obs["sampleID"] = np.where(np.arange(n_obs) % 2 == 0, "S1", "S2")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=40, flavor="seurat")
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=20)
    sc.pp.neighbors(adata)
    adata.raw = adata
    return adata


def _write_marker_toml(path: Path, genes_a, genes_b) -> str:
    content = f"""
[["Synthetic"]]
name = "Type_A"
markers = {list(genes_a)}

[["Synthetic"]]
name = "Type_B"
markers = {list(genes_b)}
"""
    path.write_text(content.strip() + "\n")
    return str(path)


def test_analysis_annotation_evidence_and_consensus_wrappers(tmp_path):
    adata = _make_preprocessed_adata()
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_clusters", random_state=42)
    marker_file = _write_marker_toml(
        tmp_path / "markers.toml",
        ["gene_0", "gene_1", "gene_2"],
        ["gene_6", "gene_7", "gene_8"],
    )
    clusters = adata.obs["leiden_clusters"].astype(str).drop_duplicates().tolist()
    markers_df = pd.DataFrame(
        {
            "group": [clusters[0], clusters[0], clusters[-1], clusters[-1]],
            "names": ["gene_0", "gene_1", "gene_6", "gene_7"],
            "scores": [8.0, 7.0, 8.0, 7.0],
            "logfoldchanges": [2.0, 1.8, 2.0, 1.8],
            "pvals_adj": [0.001, 0.001, 0.001, 0.001],
        }
    )

    review = run_annotation_evidence(
        adata,
        "leiden_clusters",
        markers_df=markers_df,
        marker_config=marker_file,
        llm_annotations={
            clusters[0]: {"llm_label": "Type_A", "llm_confidence": 0.8},
            clusters[-1]: {"llm_label": "Type_B", "llm_confidence": 0.8},
        },
    )
    assert {"final_label", "annotation_confidence", "needs_review"}.issubset(review.columns)

    consensus = build_annotation_consensus(
        adata,
        "leiden_clusters",
        review,
        key_added="cell_type_auto",
        lineage_key="celltype_lineage_auto",
    )
    assert consensus is review
    assert "cell_type_auto" in adata.obs
    assert "celltype_lineage_auto" in adata.obs
    assert "annotation_consensus_table" in adata.uns["sclucid"]["analysis"]["annotation"]


def test_run_standard_analysis_creates_analysis_maturity_review_summary(tmp_path):
    adata = _make_preprocessed_adata()
    marker_file = _write_marker_toml(
        tmp_path / "markers_workflow.toml",
        ["gene_0", "gene_1", "gene_2"],
        ["gene_6", "gene_7", "gene_8"],
    )
    config = AnalysisWorkflowConfig(
        clustering=ClusteringConfig(resolution=0.5, key_added="leiden_clusters"),
        annotation=AnnotationConfig(
            cluster_key="leiden_clusters",
            key_added="cell_type_auto",
            lineage_key="celltype_lineage_auto",
            lineage_marker_config=marker_file,
            run_scoring=False,
            final_method="celltypist",
        ),
        run_clustering_review=True,
        candidate_resolutions=[0.5],
        use_recommended_resolution=True,
        run_annotation_evidence=True,
        characterize=False,
    )

    result = run_standard_analysis(
        adata,
        config=config,
        steps=[
            "clustering_review",
            "clustering",
            "markers",
            "annotation_evidence",
            "annotation_consensus",
        ],
        show_progress=False,
    )
    review = result.uns["sclucid"]["analysis"]["review_summary"]

    assert validate_analysis_review_summary(review) == []
    assert ANALYSIS_REQUIRED_REVIEW_SECTIONS.issubset(review)
    assert review["module_maturity"]["module"] == "analysis"
    assert review["analysis_readiness"]["status"] in {"ready", "review_required"}
    policy = review["analysis_inference_policy"]
    assert policy["clustering_review"]["recommended"] is True
    assert policy["clustering_review"]["executed"] is True
    assert policy["condition_de"]["recommended_primary_method"] == "sample_level_pseudobulk"
    assert policy["condition_de"]["cell_level_compare_policy"] == "exploratory_only"
    assert policy["cell_level_compare"]["inference_level"] == "exploratory_cell_level"
    output_contract = review["analysis_output_contract"]
    assert output_contract["schema_version"] == "analysis_output_contract_v1"
    assert output_contract["cluster_key"] == "leiden_clusters"
    assert output_contract["annotation_key"] == "cell_type_auto"
    assert "cell_type" in output_contract["canonical_annotation_aliases"]
    assert {row["stage"] for row in _rows(output_contract["stage_contracts"])} >= {
        "preprocess_handoff",
        "clustering",
        "marker_discovery",
        "annotation",
        "condition_de",
        "posthoc_qc",
    }
    decision_summary = review["analysis_decision_summary"]
    assert decision_summary["schema_version"] == "analysis_decision_summary_v1"
    assert decision_summary["primary_cluster_key"] == "leiden_clusters"
    assert decision_summary["primary_annotation_key"] == "cell_type_auto"
    decision_rows = {row["step"]: row for row in _rows(decision_summary["decisions"])}
    assert {"clustering_review", "annotation_consensus", "condition_de"}.issubset(
        decision_rows
    )
    assert decision_rows["condition_de"]["decision"] in {
        "prefer_pseudobulk",
        "use_pseudobulk_results",
    }
    reviewer_table = review["analysis_reviewer_table"]
    reviewer_table_rows = _rows(reviewer_table)
    reviewer_rows = {row["item"]: row for row in reviewer_table_rows}
    assert {"clustering_review", "annotation_consensus", "cell_level_compare"}.issubset(
        reviewer_rows
    )
    required_reviewer_columns = {
        "recommended_value",
        "applied_value",
        "source",
        "confidence",
        "affected_output",
        "analysis_decision",
        "inference_level",
        "biological_risk_note",
        "review_required",
    }
    for row in reviewer_table_rows:
        assert required_reviewer_columns.issubset(row)
    assert review["clustering_evidence_summary"]["n_clusters"] > 0
    assert review["clustering_evidence_summary"]["recommendation_claim_level"] == (
        "heuristic_review_recommendation"
    )
    assert review["annotation_evidence_summary"]["review_table_rows"] > 0
    assert review["annotation_evidence_summary"]["claim_level"] == (
        "evidence_consensus_not_formal_truth"
    )
    assert review["annotation_evidence_summary"]["cluster_evidence_rows"] > 0
    cluster_evidence = _rows(
        review["annotation_evidence_summary"]["cluster_evidence_table"]
    )
    required_cluster_evidence_cols = {
        "cluster",
        "predicted_label",
        "annotation_confidence",
        "claim_level",
        "evidence_status",
        "positive_marker_support",
        "contradictory_labels",
        "reference_model_label",
        "requires_manual_review",
        "manual_review_recommendation",
    }
    for row in cluster_evidence:
        assert required_cluster_evidence_cols.issubset(row)
    assert review["annotation_consensus_summary"]["final_obs_present"] is True
    claim_summary = review["analysis_claim_level_summary"]
    assert claim_summary["schema_version"] == "analysis_claim_level_summary_v1"
    outputs = {row["output"]: row for row in _rows(claim_summary["outputs"])}
    assert outputs["cluster_marker_discovery"]["claim_level"] == (
        "exploratory_marker_screen"
    )
    assert outputs["condition_de"]["claim_level"] == "not_formal_until_pseudobulk"
    assert outputs["cell_level_compare"]["claim_level"] == (
        "exploratory_hypothesis_generation"
    )
    assert outputs["celltype_proportion"]["not_allowed_claim"]
    assert "cell_type_auto" in result.obs

    validation = validate_analysis_module_completeness(result)
    assert validation["valid"] is True
    compact = summarize_analysis_review_summary(review)
    assert compact["module"] == "analysis"
    assert compact["n_clusters"] == review["clustering_evidence_summary"]["n_clusters"]
    assert compact["condition_de_primary_method"] == "sample_level_pseudobulk"
    assert compact["cell_level_compare_policy"] == "exploratory_only"
    assert compact["global_claim_boundary"] == (
        "heuristic_and_exploratory_until_evidence_review"
    )
    assert compact["recommended_resolution_claim_level"] == (
        "heuristic_review_recommendation"
    )
    assert compact["annotation_cluster_evidence_rows"] > 0
    assert compact["primary_annotation_key"] == "cell_type_auto"
    assert compact["analysis_decision_counts"]


def test_analysis_module_contract_is_public():
    contract = get_analysis_module_contract()
    assert contract["module"] == "analysis"
    assert "scLucid.analysis.run_standard_analysis" in contract["stable_entrypoints"]
    assert "analysis_inference_policy" in contract["required_review_sections"]
    assert "analysis_claim_level_summary" in contract["required_review_sections"]
    assert contract["inference_policy_key"] == "analysis_inference_policy"
    assert contract["claim_level_key"] == "analysis_claim_level_summary"
    assert "analysis_output_contract" in contract["required_review_sections"]
    assert "analysis_decision_summary" in contract["required_review_sections"]
    assert "analysis_reviewer_table" in contract["required_review_sections"]
    assert contract["output_contract_key"] == "analysis_output_contract"
    assert contract["decision_summary_key"] == "analysis_decision_summary"
    assert contract["reviewer_table_key"] == "analysis_reviewer_table"
    assert "clustering_evidence_summary" in contract["required_review_sections"]
    assert "annotation_consensus_summary" in contract["required_review_sections"]
    assert "posthoc_qc_review_summary" in contract["required_review_sections"]
    assert "malignancy_interpretation_summary" in contract["required_review_sections"]


def test_posthoc_qc_review_summary_flags_doublet_heavy_clusters():
    adata = _make_preprocessed_adata(n_obs=30, n_vars=40)
    adata.obs["leiden_clusters"] = ["0"] * 15 + ["1"] * 15
    adata.obs["predicted_doublet"] = [True] * 12 + [False] * 18
    adata.obs["pct_counts_mt"] = [4.0] * 15 + [25.0] * 15

    summary = build_posthoc_qc_review_summary(adata, cluster_key="leiden_clusters")

    assert summary["review_required"] is True
    assert summary["n_doublet_heavy_clusters"] == 1
    assert summary["doublet_heavy_clusters"] == ["0"]
    assert summary["n_high_mitochondrial_clusters"] == 1
    assert summary["high_mitochondrial_clusters"] == ["1"]


def test_analysis_readiness_flags_low_confidence_annotation_clusters():
    adata = _make_preprocessed_adata(n_obs=30, n_vars=40)
    adata.obs["leiden_clusters"] = ["0"] * 15 + ["1"] * 15
    annotation_summary = {
        "needs_review_clusters": 0,
        "low_confidence_clusters": 1,
    }
    consensus_summary = {
        "final_obs_present": True,
        "needs_review_cells": 0,
    }

    readiness = build_analysis_readiness_assessment(
        adata=adata,
        successful_steps=["markers", "annotation_consensus"],
        cluster_key="leiden_clusters",
        preprocess_context={"pca_present": True},
        clustering_summary={},
        annotation_summary=annotation_summary,
        consensus_summary=consensus_summary,
        posthoc_qc_summary={},
        malignancy_summary={},
    )
    actions = build_analysis_review_action_items(
        readiness=readiness,
        clustering_summary={},
        annotation_summary=annotation_summary,
        consensus_summary=consensus_summary,
        posthoc_qc_summary={},
        malignancy_summary={},
    )

    assert readiness["status"] == "review_required"
    assert "annotation_low_confidence_clusters_present" in readiness["review_reasons"]
    assert any("low-confidence annotation clusters" in item["action"] for item in actions)


def test_malignancy_interpretation_bridge_adds_reviewable_outputs():
    import anndata

    genes = [
        "EPCAM",
        "KRT8",
        "KRT18",
        "MUC1",
        "MKI67",
        "TOP2A",
        "PTPRC",
        "CD3D",
        "COL1A1",
        "VWF",
    ]
    X = np.ones((12, len(genes)), dtype=np.float32)
    X[:6, :6] += 6
    X[6:, 6:] += 6
    adata = anndata.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.obs["leiden_clusters"] = ["0"] * 6 + ["1"] * 6
    adata.obs["cell_type_auto"] = ["Epithelial tumor identity"] * 6 + ["T cells"] * 6

    table = run_malignancy_interpretation(
        adata,
        annotation_key="cell_type_auto",
        cluster_key="leiden_clusters",
        run_cnv=False,
        run_malignancy_score=True,
    )

    assert "malignancy_call" in adata.obs
    assert "malignancy_interpretation_score" in adata.obs
    assert table.shape[0] == 2
    summary = adata.uns["sclucid"]["analysis"]["malignancy"][
        "malignancy_interpretation_summary"
    ]
    assert summary["available"] is True
    assert summary["n_malignant"] > 0
    assert 0 < summary["tumor_purity_estimate"] <= 1
    assert summary["low_tumor_purity_warning"] is False


def test_malignancy_interpretation_preserves_unit_interval_external_scores():
    import anndata

    adata = anndata.AnnData(X=np.ones((6, 2), dtype=np.float32))
    adata.var_names = ["GeneA", "GeneB"]
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.obs["cell_type_auto"] = ["Epithelial cells"] * adata.n_obs
    adata.obs["external_malignancy_score"] = [1.0] * adata.n_obs

    run_malignancy_interpretation(
        adata,
        annotation_key="cell_type_auto",
        run_cnv=False,
        run_malignancy_score=False,
        cnv_score_key=None,
        malignancy_score_key="external_malignancy_score",
    )

    assert (adata.obs["malignancy_interpretation_score"] > 0.5).all()
    assert set(adata.obs["malignancy_call"].astype(str)) == {"malignant"}


def test_malignancy_interpretation_validates_threshold_order():
    import anndata

    adata = anndata.AnnData(X=np.ones((2, 1), dtype=np.float32))
    adata.var_names = ["GeneA"]
    adata.obs["cell_type_auto"] = ["T cells", "T cells"]

    with pytest.raises(ValueError, match="suspect_threshold"):
        run_malignancy_interpretation(
            adata,
            annotation_key="cell_type_auto",
            threshold=0.3,
            suspect_threshold=0.5,
        )
