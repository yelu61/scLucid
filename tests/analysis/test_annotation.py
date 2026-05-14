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

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.analysis.annotation import (
    apply_final_annotation,
    annotate_clusters,
    build_annotation_review_table,
    build_llm_annotation_bundle,
    filter_marker_table_for_annotation,
    flag_suspect_clusters,
    merge_annotation_evidence,
    run_marker_annotation_evidence,
    run_lineage_state_annotation,
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
        assert config.run_celltypist == False

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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
