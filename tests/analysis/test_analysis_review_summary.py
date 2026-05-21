"""Tests for benchmark-grade analysis review summaries."""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from scLucid.analysis import (
    ANALYSIS_REQUIRED_REVIEW_SECTIONS,
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    build_annotation_consensus,
    get_analysis_module_contract,
    run_annotation_evidence,
    run_malignancy_interpretation,
    run_standard_analysis,
    summarize_analysis_review_summary,
    validate_analysis_module_completeness,
    validate_analysis_review_summary,
)


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
    assert review["clustering_evidence_summary"]["n_clusters"] > 0
    assert review["annotation_evidence_summary"]["review_table_rows"] > 0
    assert review["annotation_consensus_summary"]["final_obs_present"] is True
    assert "cell_type_auto" in result.obs

    validation = validate_analysis_module_completeness(result)
    assert validation["valid"] is True
    compact = summarize_analysis_review_summary(review)
    assert compact["module"] == "analysis"
    assert compact["n_clusters"] == review["clustering_evidence_summary"]["n_clusters"]


def test_analysis_module_contract_is_public():
    contract = get_analysis_module_contract()
    assert contract["module"] == "analysis"
    assert "scLucid.analysis.run_standard_analysis" in contract["stable_entrypoints"]
    assert "clustering_evidence_summary" in contract["required_review_sections"]
    assert "annotation_consensus_summary" in contract["required_review_sections"]
    assert "malignancy_interpretation_summary" in contract["required_review_sections"]


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
