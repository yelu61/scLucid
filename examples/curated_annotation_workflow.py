"""Example curated annotation workflow with reusable scLucid utilities.

This script demonstrates a cleaner replacement for long notebook cells that mix:
- marker filtering for dotplots
- marker/enrichment evidence generation
- manual cluster-to-label mapping
- downstream module scoring and composition plotting
"""

from pathlib import Path

import scanpy as sc

import scLucid as scl


def main() -> None:
    data_path = Path("data/clustered_input.h5ad")
    output_dir = Path("results/curated_annotation")
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(data_path)

    # 1) Cluster evidence: markers + enrichment + markdown summary
    de_config = scl.al.DifferentialConfig(
        groupby="leiden_clusters",
        use_raw=adata.raw is not None,
        pval_cutoff=0.05,
    )
    scl.al.find_markers(adata, config=de_config)

    filter_cfg = scl.al.FilterMarkersConfig(
        key="rank_genes_groups",
        key_added="highly_specific_markers_df",
        min_log2fc=0.5,
        max_padj=0.01,
        min_in_group_pct=0.2,
        min_diff_pct=0.1,
        keep_top_n=50,
    )
    highly_specific_markers_df = scl.al.filter_markers(adata, config=filter_cfg)

    enrichment_config = scl.al.EnrichmentConfig(
        de_key="rank_genes_groups_df",
        mode="offline",
        organism="human",
        gene_sets_offline=["go_bp", "reactome"],
    )
    scl.al.run_enrichment(adata, groupby="leiden_clusters", config=enrichment_config)

    scl.al.summarize_markers_and_enrichment(
        adata,
        groupby="leiden_clusters",
        markers_df=highly_specific_markers_df,
        enrichment_key="enrichment",
        summary_file=str(output_dir / "annotation_summary.md"),
    )

    # 2) Optional marker-panel cleanup before custom dotplots
    marker_panels = {
        "T_lineage": ["CD3D", "CD3E", "TRAC", "LCK"],
        "NK_lineage": ["NKG7", "GNLY", "FCGR3A", "KLRD1"],
        "B_lineage": ["MS4A1", "CD79A", "CD74", "HLA-DRA"],
    }
    filtered_panels, missing = scl.ut.filter_marker_dict(
        marker_panels,
        adata.raw.var_names if adata.raw is not None else adata.var_names,
        return_missing=True,
    )
    print("Missing marker genes:", missing)

    scl.al.visualize_markers(
        adata,
        markers=scl.ut.flatten_marker_dict(filtered_panels),
        groupby="leiden_clusters",
        plot_type="dotplot",
        swap_axes=True,
        n_genes_per_group=4,
    )

    # 3) Manual mapping should live in a sidecar table, not a huge notebook cell
    manual_mapping = {
        "0": "Naive T cells",
        "1": "Activated T cells",
        "2": "NK cells",
        "3": "B cells",
    }
    scl.al.apply_annotation_mapping(
        adata,
        cluster_key="leiden_clusters",
        mapping=manual_mapping,
        key_added="cell_type_curated",
    )

    # 4) Annotation review table for manual QC
    mgr = scl.ut.get_marker_manager(species="human", tissue="Blood")
    eval_df = scl.al.evaluate_annotation(
        adata,
        cluster_key="leiden_clusters",
        annotation_key="cell_type_curated",
        marker_config=mgr,
        plot=False,
    )
    eval_df.to_csv(output_dir / "annotation_evaluation.csv", index=False)

    # 5) Downstream state programs with the thin workflow wrapper
    modules = {
        "T_memory": ["CCR7", "LEF1", "SELL", "LTB", "TCF7"],
        "T_activation": ["ICOS", "CD69", "BATF", "FOS", "JUNB"],
        "NK_cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMB", "CTSW"],
    }
    adata, module_results = scl.al.run_module_scoring_workflow(
        adata,
        modules,
        groupby="cell_type_curated",
        sample_col="sampleID",
        condition_col="group",
        use_raw=adata.raw is not None,
    )
    module_results["group_mean_scores"].to_csv(
        output_dir / "module_scores_by_group.csv",
        index=False,
    )

    # 6) Composition plots from pre-aggregated tables
    count_df = (
        adata.obs.groupby(["group", "cell_type_curated"], observed=False)
        .size()
        .reset_index(name="count")
    )
    scl.al.plot_grouped_celltype_counts(
        count_df,
        group_col="group",
        celltype_col="cell_type_curated",
        count_col="count",
        annotate=True,
        out_dir=str(output_dir),
    )

    group_props = (
        adata.obs.groupby(["sampleID", "cell_type_curated"], observed=False)
        .size()
        .unstack(fill_value=0)
    )
    group_props = group_props.div(group_props.sum(axis=1), axis=0)
    sample_to_group = (
        adata.obs[["sampleID", "group"]].drop_duplicates().set_index("sampleID")
    )
    group_props = group_props.join(sample_to_group).groupby("group").mean(numeric_only=True)

    scl.al.plot_grouped_proportion_bar(group_props, out_dir=str(output_dir))
    scl.al.plot_celltype_alluvial(group_props, out_dir=str(output_dir))

    adata.write_h5ad(output_dir / "annotated_curated.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
