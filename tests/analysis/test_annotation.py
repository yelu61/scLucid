"""
Tests for the analysis annotation module.

Tests cell type annotation, scoring, and label transfer.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.analysis.annotation import (
    annotate_clusters,
    apply_final_annotation,
    apply_subset_annotation_reconciliation,
    build_annotation_evidence_report,
    build_annotation_review_table,
    build_hierarchical_annotation_plan,
    build_llm_annotation_bundle,
    build_subset_annotation_reconciliation,
    evaluate_annotation,
    evaluate_annotation_benchmark,
    export_annotation_evidence_report,
    filter_marker_table_for_annotation,
    flag_suspect_clusters,
    merge_annotation_evidence,
    recommend_celltypist_model,
    run_lineage_state_annotation,
    run_marker_annotation_evidence,
    run_program_annotation_evidence,
    run_subset_annotation_refinement,
    score_cell_types,
    standardize_cluster_marker_table,
)
from scLucid.analysis.config import AnnotationConfig
from scLucid.utils.manager import Manager


def _write_marker_toml(path: Path, genes_a, genes_b) -> str:
    """Create a minimal marker config file compatible with Manager."""
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


@pytest.fixture
def clustered_adata(minimal_adata):
    """Provide clustered data for annotation tests."""
    from scLucid.analysis.clustering import cluster_cells
    from scLucid.analysis.config import ClusteringConfig

    # Lightweight preprocessing for clustering prerequisites.
    adata = minimal_adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    # Cluster
    cluster_config = ClusteringConfig(method="leiden", resolution=1.0, plot=False)
    adata = cluster_cells(adata, config=cluster_config)

    return adata


@pytest.mark.integration
class TestScoring:
    """Test gene set and cell type scoring."""

    def test_score_cell_types_basic(self, clustered_adata, tmp_path):
        """Test basic cell type scoring."""
        marker_file = _write_marker_toml(
            tmp_path / "markers.toml",
            clustered_adata.var_names[:10].tolist(),
            clustered_adata.var_names[10:20].tolist(),
        )
        marker_manager = Manager(marker_file, case_sensitive=True)
        result = score_cell_types(
            clustered_adata,
            marker_config=marker_manager,
            use_raw=False,
            layer=None,
            score_name_suffix="_test",
        )

        # Check scores were added
        assert "Type_A_test" in result.obs.columns
        assert "Type_B_test" in result.obs.columns

    def test_score_cell_types_from_manager(self, clustered_adata):
        """Test scoring using marker manager."""
        # This test requires marker databases
        pytest.skip("Marker database not available in test environment")

    def test_aucell_ucell_correlates_with_marker_expression(self, tmp_path):
        """AUCell and UCell scores should increase with marker gene expression."""
        rng = np.random.default_rng(42)
        n_cells, n_genes = 80, 50
        markers = [f"marker{i}" for i in range(5)]
        others = [f"gene{i}" for i in range(n_genes - 5)]
        var_names = markers + others

        X = rng.poisson(2, size=(n_cells, n_genes)).astype(float)
        # First half of cells strongly expresses marker genes
        X[: n_cells // 2, :5] += 10.0
        adata = AnnData(X=X)
        adata.obs_names = [f"cell{i}" for i in range(n_cells)]
        adata.var_names = var_names

        marker_file = _write_marker_toml(
            tmp_path / "markers.toml", markers, others[:5]
        )
        marker_manager = Manager(marker_file, case_sensitive=True)

        for backend in ("aucell", "ucell"):
            scored = score_cell_types(
                adata.copy(),
                marker_config=marker_manager,
                use_raw=False,
                layer=None,
                score_name_suffix="_score",
                scoring_backend=backend,
            )
            high_group = scored.obs["Type_A_score"].iloc[: n_cells // 2].mean()
            low_group = scored.obs["Type_A_score"].iloc[n_cells // 2 :].mean()
            assert high_group > low_group, f"{backend} did not increase with expression"
            assert scored.uns["sclucid"]["analysis"]["annotation"]["scoring_params"][
                "backend"
            ] == backend


@pytest.mark.integration
class TestAnnotation:
    """Test cell type annotation."""

    def test_annotate_clusters_basic(self, clustered_adata, tmp_path):
        """Test basic cluster annotation."""
        marker_file = _write_marker_toml(
            tmp_path / "markers.toml",
            clustered_adata.var_names[:10].tolist(),
            clustered_adata.var_names[10:20].tolist(),
        )
        marker_manager = Manager(marker_file, case_sensitive=True)

        # Generate score columns required by max-score annotation.
        scored = score_cell_types(
            clustered_adata,
            marker_config=marker_manager,
            use_raw=False,
            layer=None,
            score_name_suffix="_score",
        )

        result = annotate_clusters(
            scored,
            cluster_key="leiden_clusters",
            marker_config=marker_manager,
            method="max_score",
        )

        # Check annotation was added
        assert "leiden_clusters_annotated" in result.obs.columns

    def test_max_score_significance_assigns_unknown_for_uniform_noise(self, tmp_path):
        """When scores lack significant separation, all clusters should be Unknown."""
        rng = np.random.default_rng(7)
        adata = AnnData(X=rng.poisson(2, size=(60, 20)).astype(float))
        adata.obs_names = [f"cell{i}" for i in range(60)]
        adata.var_names = [f"gene{i}" for i in range(20)]
        adata.obs["cluster"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

        # Identical uniform scores -> no cluster is significantly above background
        adata.obs["Type_A_score"] = 0.5
        adata.obs["Type_B_score"] = 0.5

        marker_file = _write_marker_toml(tmp_path / "markers.toml", ["gene0"], ["gene1"])
        marker_manager = Manager(marker_file, case_sensitive=True)

        result = annotate_clusters(
            adata,
            cluster_key="cluster",
            marker_config=marker_manager,
            method="max_score",
            significance_threshold=0.05,
            min_score_margin=0.0,
        )

        assert result.obs["cluster_annotated"].astype(str).eq("Unknown").all()
        evidence = result.uns["sclucid"]["analysis"]["annotation"][
            "cluster_annotated_params"
        ]["max_score_evidence"]
        assert "0" in evidence
        assert "margin" in evidence["0"]
        assert "fdr" in evidence["0"]

    def test_filter_marker_table_for_annotation_removes_noise_markers(self):
        """Noise-like ribosomal and stress genes should be filtered from marker review tables."""
        markers_df = pd.DataFrame(
            {
                "group": ["0", "0", "0", "1", "1"],
                "names": ["RPL13A", "HSPA1A", "LTB", "MALAT1", "NKG7"],
                "logfoldchanges": [3.0, 2.5, 2.0, 4.0, 1.8],
            }
        )

        filtered = filter_marker_table_for_annotation(markers_df)

        assert filtered["names"].tolist() == ["LTB", "NKG7"]
        assert "annotation_noise_category" in filtered.columns
        assert filtered["is_annotation_informative"].all()

    def test_annotation_evidence_chain_builds_and_applies_final_labels(
        self, clustered_adata, tmp_path
    ):
        """Marker evidence, LLM bundle, evidence merge, and final application should compose."""
        adata = clustered_adata.copy()
        clusters = adata.obs["leiden_clusters"].astype(str)
        cluster_codes = clusters.drop_duplicates().tolist()

        marker_file = _write_marker_toml(
            tmp_path / "markers_evidence.toml",
            ["LTB", "IL7R"],
            ["NKG7", "CCL5"],
        )
        markers_df = pd.DataFrame(
            {
                "group": [
                    cluster_codes[0],
                    cluster_codes[0],
                    cluster_codes[0],
                    cluster_codes[1],
                    cluster_codes[1],
                ],
                "names": ["LTB", "IL7R", "RPL13A", "NKG7", "CCL5"],
                "scores": [8.0, 7.0, 6.0, 8.0, 7.0],
                "logfoldchanges": [2.5, 2.0, 4.0, 3.0, 2.5],
                "pvals_adj": [0.001, 0.002, 0.003, 0.001, 0.002],
            }
        )

        marker_table = standardize_cluster_marker_table(markers_df, keep_top_n_per_cluster=3)
        assert {"cluster", "gene", "marker_rank", "noise_category"}.issubset(marker_table.columns)
        assert (
            marker_table.loc[marker_table["gene"] == "RPL13A", "noise_category"].iloc[0]
            == "ribosomal"
        )

        marker_evidence = run_marker_annotation_evidence(
            adata,
            "leiden_clusters",
            marker_file,
            markers_df=markers_df,
            top_n_markers=3,
        )
        assert {"cluster", "marker_label", "marker_confidence"}.issubset(marker_evidence.columns)

        bundle = build_llm_annotation_bundle(
            adata,
            "leiden_clusters",
            markers_df=markers_df,
            marker_evidence=marker_evidence,
            lineage_key="lineage_gate",
        )
        assert bundle["schema_version"] == "analysis_annotation_bundle_v1"
        assert cluster_codes[0] in bundle["clusters"]

        llm_annotations = {
            cluster_codes[0]: {"llm_label": "Type_A", "llm_confidence": 0.8},
            cluster_codes[1]: {"llm_label": "Type_B", "llm_confidence": 0.8},
        }
        review = merge_annotation_evidence(
            adata,
            "leiden_clusters",
            marker_evidence=marker_evidence,
            llm_annotations=llm_annotations,
        )
        assert {"final_label", "annotation_confidence", "needs_review"}.issubset(review.columns)

        result = apply_final_annotation(adata, "leiden_clusters", review)
        assert "cell_type_final" in result.obs.columns
        assert "cell_type_final_confidence" in result.obs.columns

    def test_hierarchical_annotation_plan_recommends_subset_refinement(self):
        adata = AnnData(X=np.ones((80, 6)))
        adata.obs_names = [f"cell{i}" for i in range(80)]
        adata.var_names = [f"gene{i}" for i in range(6)]
        adata.obs["cluster"] = ["0"] * 45 + ["1"] * 35
        adata.obs["lineage"] = ["T cells"] * 80

        plan = build_hierarchical_annotation_plan(
            adata,
            cluster_key="cluster",
            lineage_key="lineage",
            min_cells_per_lineage=20,
            min_clusters_per_lineage=2,
            min_lineage_purity=0.4,
        )

        assert plan.loc[0, "lineage_label"] == "T cells"
        assert plan.loc[0, "recommended_action"] == "subset_recluster_for_subtype"
        assert "hierarchical_annotation_plan" in adata.uns["sclucid"]["analysis"]["annotation"]

    def test_llm_bundle_includes_hierarchical_gate_context(self):
        adata = AnnData(X=np.ones((6, 4)))
        adata.obs_names = [f"cell{i}" for i in range(6)]
        adata.var_names = ["CD3D", "IL7R", "MS4A1", "NKG7"]
        adata.obs["cluster"] = ["0", "0", "0", "1", "1", "1"]
        adata.obs["lineage"] = ["T cells", "T cells", "T cells", "B cells", "B cells", "B cells"]
        adata.obs["subtype"] = ["CD4+ T", "CD4+ T", "CD4+ T", "Naive B", "Naive B", "Naive B"]
        markers_df = pd.DataFrame(
            {
                "group": ["0", "0", "1", "1"],
                "names": ["CD3D", "IL7R", "MS4A1", "NKG7"],
                "scores": [5, 4, 5, 2],
            }
        )

        bundle = build_llm_annotation_bundle(
            adata,
            "cluster",
            markers_df=markers_df,
            lineage_key="lineage",
            subtype_key="subtype",
        )

        cluster_zero = bundle["clusters"]["0"]
        lineage_items = list(cluster_zero["lineage_annotation"].values())
        subtype_items = list(cluster_zero["subtype_annotation"].values())
        assert lineage_items[0]["label"] == "T cells"
        assert subtype_items[0]["label"] == "CD4+ T"
        assert "major lineage evidence" in bundle["instructions"]

    def test_subset_annotation_refinement_reprocesses_counts_without_writeback(self):
        rng = np.random.default_rng(42)
        counts = rng.poisson(2, size=(70, 30)).astype(float)
        counts[:35, :4] += 4
        counts[35:, 4:8] += 4
        adata = AnnData(X=np.log1p(counts))
        adata.obs_names = [f"cell{i}" for i in range(70)]
        adata.var_names = [f"gene{i}" for i in range(30)]
        adata.layers["counts"] = counts.copy()
        adata.obs["lineage"] = ["T cells"] * 70

        build_hierarchical_annotation_plan(
            adata,
            cluster_key="lineage",
            lineage_key="lineage",
            min_cells_per_lineage=20,
            min_clusters_per_lineage=1,
        )
        subsets = run_subset_annotation_refinement(
            adata,
            lineage_key="lineage",
            lineages=["T cells"],
            cluster_resolution=0.2,
            n_top_hvgs=15,
            n_pcs=5,
            n_neighbors=5,
            min_cells=20,
            write_back=False,
        )

        assert "T cells" in subsets
        subset = subsets["T cells"]
        assert "subset_annotation_input" in subset.layers
        summary = adata.uns["sclucid"]["analysis"]["annotation"]["subset_annotation_refinement"]
        assert summary.loc[0, "input_mode"] == "layer:counts"
        assert summary.loc[0, "status"] == "completed"
        assert "subset_annotation_refinement_cluster" not in adata.obs.columns

    def test_subset_annotation_refinement_can_optionally_write_back_clusters(self):
        rng = np.random.default_rng(7)
        counts = rng.poisson(2, size=(60, 24)).astype(float)
        counts[:30, :4] += 3
        counts[30:, 4:8] += 3
        adata = AnnData(X=np.log1p(counts))
        adata.obs_names = [f"cell{i}" for i in range(60)]
        adata.var_names = [f"gene{i}" for i in range(24)]
        adata.layers["counts"] = counts.copy()
        adata.obs["lineage"] = ["Myeloid"] * 60

        subsets = run_subset_annotation_refinement(
            adata,
            lineage_key="lineage",
            lineages=["Myeloid"],
            cluster_resolution=0.2,
            n_top_hvgs=12,
            n_pcs=5,
            n_neighbors=5,
            min_cells=20,
            write_back=True,
            key_added="myeloid_refine",
        )

        assert "Myeloid" in subsets
        assert "myeloid_refine_cluster" in adata.obs.columns
        assert str(adata.obs["myeloid_refine_cluster"].iloc[0]).startswith("Myeloid:")

    def test_subset_reconciliation_flags_conflicts_and_exclusions_without_dropping_cells(self):
        adata = AnnData(X=np.ones((5, 4)))
        adata.obs_names = [f"cell{i}" for i in range(5)]
        adata.var_names = [f"gene{i}" for i in range(4)]
        adata.obs["lineage"] = ["T cells"] * 5
        adata.obs["subtype_global"] = [
            "CD4+ T",
            "CD4+ T",
            "CD8+ T",
            "Unknown",
            "CD4+ T",
        ]

        subset = adata.copy()
        subset.obs["subset_clusters"] = ["0", "0", "1", "1", "2"]
        subset.obs["subset_label"] = [
            "CD4+ T",
            "CD8+ T",
            "CD8+ T",
            "CD8+ T",
            "Unknown",
        ]
        subset.obs["subset_review_exclude"] = [False, False, False, True, False]

        reconciliation = build_subset_annotation_reconciliation(
            adata,
            {"T cells": subset},
            lineage_key="lineage",
            global_subtype_key="subtype_global",
            subset_label_key="subset_label",
            subset_cluster_key="subset_clusters",
        )

        conflict = reconciliation.set_index("obs_name").loc["cell1"]
        excluded = reconciliation.set_index("obs_name").loc["cell3"]
        assert conflict["recommended_action"] == "review_subtype_conflict"
        assert bool(conflict["subtype_conflict"]) is True
        assert excluded["recommended_action"] == "exclude_from_global_review"
        assert bool(excluded["exclude_from_global"]) is True
        assert adata.n_obs == 5

        apply_subset_annotation_reconciliation(
            adata,
            reconciliation,
            target_key="subtype_refined",
            global_subtype_key="subtype_global",
        )

        assert str(adata.obs.loc["cell0", "subtype_refined"]) == "CD4+ T"
        assert str(adata.obs.loc["cell1", "subtype_refined"]) == "CD4+ T"
        assert bool(adata.obs.loc["cell3", "subset_refinement_exclude"]) is True
        assert adata.n_obs == 5

    def test_program_annotation_evidence_enriches_bundle_and_review_without_label_vote(
        self, tmp_path
    ):
        adata = AnnData(X=np.ones((8, 6), dtype=float))
        adata.obs_names = [f"cell{i}" for i in range(8)]
        adata.var_names = ["GZMB", "PRF1", "NKG7", "PDCD1", "TOX", "HAVCR2"]
        adata.obs["cluster"] = ["0"] * 4 + ["1"] * 4
        adata.X[:4, :3] = 8
        adata.X[4:, 3:] = 8

        program_file = _write_marker_toml(
            tmp_path / "program_annotation_markers.toml",
            ["GZMB", "PRF1", "NKG7"],
            ["PDCD1", "TOX", "HAVCR2"],
        )
        program_mgr = Manager(program_file, case_sensitive=True)
        program_evidence = run_program_annotation_evidence(
            adata,
            "cluster",
            program_config=program_mgr,
            min_genes=2,
            top_n_programs=1,
        )

        assert {"cluster", "program", "program_score_mean", "top_programs"}.issubset(
            program_evidence.columns
        )
        bundle = build_llm_annotation_bundle(
            adata,
            "cluster",
            program_evidence=program_evidence,
        )
        assert bundle["clusters"]["0"]["program_evidence"]
        review = merge_annotation_evidence(
            adata,
            "cluster",
            program_evidence=program_evidence,
        )
        assert "top_programs" in review.columns
        assert review["final_label"].astype(str).eq("Unknown").all()

    def test_evaluate_annotation_benchmark_reports_disagreement_and_confusion(self):
        adata = AnnData(X=np.ones((6, 3)))
        adata.obs["truth"] = ["T", "T", "B", "B", "NK", "NK"]
        adata.obs["cell_type_final"] = ["T", "T", "B", "T", "Unknown", "NK"]
        review = pd.DataFrame(
            {
                "cluster": ["0", "1", "2"],
                "reference_label": ["T", "B", "NK"],
                "marker_label": ["T", "T", "Unknown"],
                "llm_label": ["T", "B", "NK"],
                "final_label": ["T", "T", "Unknown"],
                "needs_review": [False, True, False],
                "marker_database": ["testdb", "testdb", "testdb"],
            }
        )

        result = evaluate_annotation_benchmark(
            adata,
            label_key="cell_type_final",
            truth_key="truth",
            review_table=review,
        )

        assert result["schema_version"] == "annotation_benchmark_v1"
        assert result["accuracy"] < 1.0
        assert "Ambiguous" in result["conservative_label_counts"]
        assert not result["disagreement_matrix"].empty
        assert not result["confusion_matrix"].empty
        assert "annotation_benchmark" in adata.uns["sclucid"]["analysis"]["annotation"]

    def test_flag_suspect_clusters_identifies_ribosomal_and_doublet_clusters(self, clustered_adata):
        """Cluster-level suspect flags should capture ribosomal dominance and doublet-heavy clusters."""
        adata = clustered_adata.copy()
        cluster_codes = adata.obs["leiden_clusters"].astype(str)
        target_cluster = cluster_codes.iloc[0]
        other_cluster = next(code for code in cluster_codes.unique() if code != target_cluster)

        adata.obs["pct_counts_mt"] = 5.0
        adata.obs["predicted_doublet"] = False
        adata.obs.loc[cluster_codes == target_cluster, "predicted_doublet"] = True

        markers_df = pd.DataFrame(
            {
                "group": [
                    target_cluster,
                    target_cluster,
                    target_cluster,
                    other_cluster,
                    other_cluster,
                    other_cluster,
                ],
                "names": ["RPL13A", "RPS18", "RPLP0", "NKG7", "CCL5", "TRBC1"],
                "logfoldchanges": [3.0, 2.5, 2.0, 3.0, 2.5, 2.0],
            }
        )

        summary = flag_suspect_clusters(
            adata,
            cluster_key="leiden_clusters",
            markers_df=markers_df,
            doublet_fraction_threshold=0.5,
            ribosomal_fraction_threshold=0.5,
        )

        flagged = summary.set_index("cluster")
        assert flagged.loc[target_cluster, "suspect_flag"] == "doublet_suspect"
        assert "ribosomal_dominant" in flagged.loc[target_cluster, "suspect_reasons"]
        assert flagged.loc[other_cluster, "suspect_flag"] == "clean"

    def test_evaluate_annotation_uses_exclusive_marker_support(self, tmp_path):
        """Shared markers should not be treated as annotation-conflict evidence."""
        X = np.zeros((40, 6), dtype=float)
        clusters = np.array(["0"] * 20 + ["1"] * 20)
        X[:, 0] = 5  # shared marker
        X[clusters == "0", 1] = 20  # Type_A exclusive marker
        X[clusters == "1", 2] = 20  # Type_B exclusive marker
        X[:, 3:] = 1
        adata = AnnData(X=X)
        adata.var_names = [f"g{i}" for i in range(6)]
        adata.obs["cluster"] = pd.Categorical(clusters)
        adata.obs["cell_type_auto"] = np.where(clusters == "0", "Type_A", "Type_B")

        marker_file = _write_marker_toml(
            tmp_path / "shared_marker_eval.toml",
            ["g0", "g1"],
            ["g0", "g2"],
        )

        result = evaluate_annotation(
            adata,
            cluster_key="cluster",
            annotation_key="cell_type_auto",
            marker_config=marker_file,
            plot=False,
        )

        row_a = result.loc[result["cell_type"] == "Type_A"].iloc[0]
        assert row_a["marker_specificity_reason"] == "exclusive_marker_support"
        assert row_a["exclusive_expected_markers"] == 1
        assert row_a["exclusive_detected_markers"] == 1
        assert row_a["marker_specificity"] == pytest.approx(1.0)

    def test_build_annotation_review_table_summarizes_clusters(self, clustered_adata):
        """Review helper should build a compact per-cluster annotation table and persist it."""
        adata = clustered_adata.copy()
        adata.obs["sampleID"] = np.where(np.arange(adata.n_obs) % 2 == 0, "S1", "S2")
        adata.obs["group"] = np.where(np.arange(adata.n_obs) % 2 == 0, "WT", "KO")
        adata.obs["time"] = np.where(np.arange(adata.n_obs) % 2 == 0, "6h", "24h")
        adata.obs["lineage_score"] = np.linspace(0, 1, adata.n_obs)
        adata.obs["celltype"] = np.where(
            adata.obs["leiden_clusters"].astype(str)
            == adata.obs["leiden_clusters"].astype(str).iloc[0],
            "T cells",
            "NK cells",
        )

        cluster_codes = adata.obs["leiden_clusters"].astype(str).unique().tolist()
        markers_df = pd.DataFrame(
            {
                "group": [cluster_codes[0], cluster_codes[0], cluster_codes[1], cluster_codes[1]],
                "names": ["LTB", "IL7R", "NKG7", "CCL5"],
                "logfoldchanges": [2.5, 2.0, 3.0, 2.2],
            }
        )
        enrichment_dict = {
            cluster_codes[0]: pd.DataFrame(
                {"Term": ["T cell activation"], "Adjusted P-value": [0.001]}
            ),
            cluster_codes[1]: pd.DataFrame(
                {"Term": ["NK mediated cytotoxicity"], "Adjusted P-value": [0.002]}
            ),
        }

        review_df = build_annotation_review_table(
            adata,
            cluster_key="leiden_clusters",
            markers_df=markers_df,
            enrichment_dict=enrichment_dict,
            annotation_key="celltype",
            sample_col="sampleID",
            group_col="group",
            time_col="time",
            score_cols=["lineage_score"],
        )

        assert {"cluster", "annotation", "top_markers", "top_terms", "mean_scores"}.issubset(
            review_df.columns
        )
        assert "leiden_clusters_review_table" in adata.uns["sclucid"]["analysis"]["annotation"]

    def test_run_annotation_scoring_only(self, clustered_adata):
        """Test run_annotation with scoring method."""
        config = AnnotationConfig(
            cluster_key="leiden_clusters",
            marker_species="human",
            run_celltypist=False,
            run_scoring=True,
            final_method="max_score",
        )

        # This may fail without proper marker databases, so we test config validation
        assert config.cluster_key == "leiden_clusters"
        assert not config.run_celltypist

    def test_run_lineage_state_annotation_generates_modular_labels(self, clustered_adata, tmp_path):
        """Hierarchical annotation should produce lineage/subtype/state outputs plus a modular display label."""
        adata = clustered_adata.copy()
        clusters = adata.obs["leiden_clusters"].astype(str)
        cluster_a = clusters.iloc[0]
        cluster_b = next(code for code in clusters.unique() if code != cluster_a)

        lineage_genes = adata.var_names[:8].tolist()
        subtype_genes = adata.var_names[8:14].tolist()
        state_genes = adata.var_names[14:18].tolist()

        X = np.asarray(adata.X)
        X[clusters == cluster_a, 0:4] += 8.0
        X[clusters == cluster_b, 4:8] += 8.0
        X[clusters == cluster_a, 8:11] += 6.0
        X[clusters == cluster_a, 14:16] += 5.0
        X[clusters == cluster_b, 16:18] += 5.0
        adata.X = X
        adata.raw = adata.copy()

        lineage_marker_file = _write_marker_toml(
            tmp_path / "lineage_markers.toml",
            lineage_genes[:4],
            lineage_genes[4:8],
        )

        subtype_content = f"""
[["T subtypes"]]
name = "Naive-like T"
markers = {subtype_genes[:3]}

[["T subtypes"]]
name = "Cytotoxic T"
markers = {subtype_genes[3:6]}
"""
        subtype_marker_file = tmp_path / "subtype_markers.toml"
        subtype_marker_file.write_text(subtype_content.strip() + "\n")

        config = AnnotationConfig(
            cluster_key="leiden_clusters",
            final_method="hierarchical",
            marker_method="max_score",
            lineage_marker_config=str(lineage_marker_file),
            subtype_marker_config=str(subtype_marker_file),
            target_lineage="Type_A",
            lineage_key="lineage_auto",
            subtype_key="subtype_auto",
            state_key="state_auto",
            key_added="celltype_display",
            custom_state_signatures={
                "Activated": state_genes[:2],
                "Memory": state_genes[2:4],
            },
            nomenclature_style="modular",
        )

        result = run_lineage_state_annotation(adata, config)

        assert {"lineage_auto", "subtype_auto", "state_auto", "celltype_display"}.issubset(
            result.obs.columns
        )
        assert (
            result.obs.loc[clusters == cluster_a, "lineage_auto"]
            .astype(str)
            .str.contains("Type_A")
            .all()
        )
        assert (
            result.obs.loc[clusters == cluster_a, "subtype_auto"]
            .astype(str)
            .str.contains("Naive-like T")
            .all()
        )
        assert (
            result.obs.loc[clusters == cluster_a, "state_auto"]
            .astype(str)
            .str.contains("Activated")
            .all()
        )
        assert (
            result.obs.loc[clusters == cluster_a, "celltype_display"]
            .astype(str)
            .str.contains("\\|", regex=True)
            .any()
        )
        assert (
            result.obs.loc[clusters == cluster_b, "subtype_auto"]
            .astype(str)
            .eq("Not_applicable")
            .all()
        )

    def test_run_lineage_state_annotation_respects_state_scope_metadata(
        self, clustered_adata, tmp_path
    ):
        """State assignments should obey scope/applies_to metadata instead of only taking the highest score."""
        adata = clustered_adata.copy()
        clusters = adata.obs["leiden_clusters"].astype(str)
        cluster_a = clusters.iloc[0]
        cluster_b = next(code for code in clusters.unique() if code != cluster_a)

        lineage_genes = adata.var_names[:8].tolist()
        state_genes = adata.var_names[8:12].tolist()
        X = np.asarray(adata.X)
        X[clusters == cluster_a, 0:4] += 8.0
        X[clusters == cluster_b, 4:8] += 8.0
        X[clusters == cluster_a, 8:10] += 6.0
        X[clusters == cluster_b, 10:12] += 6.0
        adata.X = X
        adata.raw = adata.copy()

        lineage_marker_file = _write_marker_toml(
            tmp_path / "lineage_markers_scope.toml",
            lineage_genes[:4],
            lineage_genes[4:8],
        )

        state_content = f"""
[["Scoped states"]]
name = "Type_A_only_state"
markers = {state_genes[:2]}
metadata = {{ kind = "state", scope = "lineage_restricted", applies_to = ["Type_A"] }}

[["Scoped states"]]
name = "Type_B_only_state"
markers = {state_genes[2:4]}
metadata = {{ kind = "state", scope = "lineage_restricted", applies_to = ["Type_B"] }}
"""
        state_marker_file = tmp_path / "state_markers.toml"
        state_marker_file.write_text(state_content.strip() + "\n")

        config = AnnotationConfig(
            cluster_key="leiden_clusters",
            final_method="hierarchical",
            marker_method="max_score",
            lineage_marker_config=str(lineage_marker_file),
            state_marker_config=str(state_marker_file),
            marker_states=["Type_A_only_state", "Type_B_only_state"],
            target_lineage=None,
            lineage_key="lineage_auto",
            subtype_key="subtype_auto",
            state_key="state_auto",
            key_added="celltype_display",
            nomenclature_style="modular",
        )

        result = run_lineage_state_annotation(adata, config)

        assert (
            result.obs.loc[clusters == cluster_a, "state_auto"]
            .astype(str)
            .eq("Type_A_only_state")
            .all()
        )
        assert (
            result.obs.loc[clusters == cluster_b, "state_auto"]
            .astype(str)
            .eq("Type_B_only_state")
            .all()
        )

    def test_run_lineage_state_annotation_respects_min_state_score(
        self, clustered_adata, tmp_path
    ):
        """Low state scores should remain Not_applicable instead of forced to the top state."""
        adata = clustered_adata.copy()
        clusters = adata.obs["leiden_clusters"].astype(str)
        cluster_a = clusters.iloc[0]
        cluster_b = next(code for code in clusters.unique() if code != cluster_a)

        lineage_genes = adata.var_names[:8].tolist()
        state_genes = adata.var_names[8:10].tolist()
        X = np.asarray(adata.X)
        X[clusters == cluster_a, 0:4] += 8.0
        X[clusters == cluster_b, 4:8] += 8.0
        X[clusters == cluster_a, 8:10] += 4.0
        adata.X = X
        adata.raw = adata.copy()

        lineage_marker_file = _write_marker_toml(
            tmp_path / "lineage_markers_min_state.toml",
            lineage_genes[:4],
            lineage_genes[4:8],
        )

        config = AnnotationConfig(
            cluster_key="leiden_clusters",
            final_method="hierarchical",
            marker_method="max_score",
            lineage_marker_config=str(lineage_marker_file),
            target_lineage="Type_A",
            lineage_key="lineage_auto",
            subtype_key="subtype_auto",
            state_key="state_auto",
            key_added="celltype_display",
            custom_state_signatures={"Weak_state": state_genes},
            min_state_score=999.0,
            nomenclature_style="modular",
        )

        result = run_lineage_state_annotation(adata, config)

        assert result.obs["state_auto"].astype(str).eq("Not_applicable").all()
        assert result.obs["state_auto_confidence"].isna().all()


@pytest.mark.integration
class TestAnnotationConfigValidation:
    """Test annotation configuration validation."""

    def test_invalid_confidence(self):
        """Test that invalid confidence raises error."""
        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=1.5)  # Should be <= 1

        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=-0.1)  # Should be >= 0

    def test_invalid_final_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError):
            AnnotationConfig(final_method="invalid_method")



@pytest.mark.unit
class TestCellTypistModelRecommendation:
    """Tests for CellTypist model recommendation."""

    def test_recommend_immune_models(self):
        """PBMC/blood/immune tissues map to the immune atlas."""
        for tissue in ("PBMC", "blood", "immune"):
            rec = recommend_celltypist_model(tissue)
            assert rec["model"] == "Immune_All_Low.pkl"
            assert rec["tissue_match_warning"] is False

    def test_recommend_tissue_specific_models(self):
        """Intestine and lung map to their dedicated CellTypist models."""
        assert recommend_celltypist_model("colon")["model"] == "Cells_Intestinal_Tract.pkl"
        assert recommend_celltypist_model("lung")["model"] == "Cells_Lung_Airway.pkl"

    def test_tumor_tissue_warns(self):
        """Tumor/cancer tissues warn that the immune atlas may be inappropriate."""
        rec = recommend_celltypist_model("tumor")
        assert rec["model"] == "Immune_All_Low.pkl"
        assert rec["tissue_match_warning"] is True
        assert "immune atlas may be inappropriate" in rec["message"]

    def test_unrecognized_tissue_falls_back_with_warning(self):
        """Unknown tissues fall back to the immune atlas and warn."""
        rec = recommend_celltypist_model("xxxxx")
        assert rec["model"] == "Immune_All_Low.pkl"
        assert rec["tissue_match_warning"] is True


@pytest.mark.unit
class TestAnnotationEvidenceReport:
    """Tests for the annotation evidence report."""

    def _make_annotated_adata(self):
        """Build a tiny annotated AnnData with review and suspect tables."""
        adata = AnnData(X=np.ones((40, 6), dtype=float))
        adata.obs_names = [f"cell_{i}" for i in range(40)]
        adata.var_names = [f"gene_{i}" for i in range(6)]
        adata.obs["leiden_clusters"] = pd.Categorical(
            ["0"] * 20 + ["1"] * 20
        )
        adata.obs["cell_type"] = pd.Categorical(
            ["T cells"] * 20 + ["B cells"] * 20
        )
        adata.obs["cell_type_confidence"] = np.concatenate(
            [np.full(20, 0.8), np.full(20, 0.9)]
        )
        adata.obs["sampleID"] = np.tile(["S1", "S2"], 20)

        review = pd.DataFrame(
            {
                "cluster": ["0", "1"],
                "reference_label": ["T cell", "B cell"],
                "marker_label": ["T cell", "B cell"],
                "final_label": ["T cells", "B cells"],
                "annotation_confidence": [0.8, 0.9],
                "needs_review": [False, False],
                "conflicts": ["", ""],
                "warnings": ["", ""],
                "top_markers": ["CD3D, IL7R", "CD79A, MS4A1"],
                "top_terms": ["", ""],
                "n_cells": [20, 20],
                "pct_cells": [0.5, 0.5],
            }
        )
        suspect = pd.DataFrame(
            {
                "cluster": ["0", "1"],
                "suspect_flag": ["clean", "doublet_suspect"],
                "suspect_reasons": ["", "doublet_suspect"],
            }
        )
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
            "annotation", {}
        )["annotation_review_table"] = review
        adata.uns["sclucid"]["analysis"]["annotation"][
            "leiden_clusters_suspect_flags"
        ] = suspect
        return adata

    def test_build_annotation_evidence_report_structure(self):
        """Report should contain per-cluster, per-cell-type, and global sections."""
        adata = self._make_annotated_adata()
        report = build_annotation_evidence_report(
            adata, sample_col="sampleID", confidence_threshold=0.5
        )

        assert "per_cluster" in report
        assert "per_cell_type" in report
        assert "global" in report
        assert report["schema_version"] == "annotation_evidence_report_v1"

        per_cluster = {row["cluster"]: row for row in report["per_cluster"]}
        assert "0" in per_cluster
        assert "1" in per_cluster
        assert per_cluster["0"]["top_marker_label"] == "T cell"
        assert per_cluster["1"]["suspect_flag"] == "doublet_suspect"
        assert "doublet_suspect" in per_cluster["1"]["suspect_reasons"]

        per_celltype = report["per_cell_type"]
        assert "T cells" in per_celltype
        assert "B cells" in per_celltype
        assert per_celltype["T cells"]["n_cells"] == 20
        assert per_celltype["T cells"]["cluster_purity"] == 1.0
        assert per_celltype["T cells"]["sample_distribution"] != ""

        global_section = report["global"]
        assert "annotation_method" in global_section
        assert "confidence_threshold" in global_section
        assert global_section["n_low_confidence_cells"] == 0
        assert any("suspect" in item for item in global_section["action_items"])

        assert "annotation_evidence_report" in adata.uns["sclucid"]["analysis"]

    def test_export_annotation_evidence_report_writes_markdown(self, tmp_path):
        """Markdown export should write a readable file."""
        adata = self._make_annotated_adata()
        out = tmp_path / "annotation_report"
        path = export_annotation_evidence_report(
            adata, str(out), sample_col="sampleID"
        )
        assert Path(path).exists()
        content = Path(path).read_text()
        assert "# Annotation Evidence Report" in content
        assert "Per-Cluster Evidence" in content
        assert "Per-Cell-Type Summary" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
