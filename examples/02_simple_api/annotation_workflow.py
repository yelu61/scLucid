"""Example curated annotation workflow with reusable scLucid utilities.

This script demonstrates a cleaner replacement for long notebook cells that mix:
- marker filtering for dotplots
- marker/enrichment evidence generation
- manual cluster-to-label mapping
- downstream module scoring and composition plotting

The marker-DE calls in this file are for annotation evidence and exploratory
marker discovery. Formal condition DE should use sample-level pseudobulk via
``scl.al.run_pseudobulk_de``.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl


def main() -> None:
    data_path = Path("data/clustered_input.h5ad")
    output_dir = Path("results/curated_annotation")
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(data_path)

    # 1) Cluster evidence: markers + enrichment + markdown summary
    #    find_markers() is cell-level marker discovery, not publication-grade
    #    condition DE. Its output carries inference_level metadata.
    de_config = scl.al.DifferentialConfig(
        groupby="leiden_clusters",
        use_raw=adata.raw is not None,
        pval_cutoff=0.05,
    )
    marker_df = scl.al.find_markers(adata, config=de_config)
    print(
        "Marker discovery inference level:",
        marker_df.get("inference_level", ["unknown"])[0] if not marker_df.empty else "empty",
    )

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

    # 3) First pass: conservative major-lineage mapping. Subtype/state labels
    #    should be layered after lineage-gated subset refinement.
    manual_mapping = {
        "0": "T cells",
        "1": "T cells",
        "2": "NK cells",
        "3": "B cells",
    }
    scl.al.apply_annotation_mapping(
        adata,
        cluster_key="leiden_clusters",
        mapping=manual_mapping,
        key_added="lineage_curated",
    )
    adata.obs["cell_type_curated"] = adata.obs["lineage_curated"].astype(str)

    # 4) Annotation review table for manual QC
    mgr = scl.ut.get_marker_manager(species="human", tissue="Blood")
    eval_df = scl.al.evaluate_annotation(
        adata,
        cluster_key="leiden_clusters",
        annotation_key="lineage_curated",
        marker_config=mgr,
        plot=False,
    )
    eval_df.to_csv(output_dir / "annotation_evaluation.csv", index=False)

    # 5) Hierarchical annotation: lineage-gated subset refinement.
    #    The subset is extracted from raw/counts when available and reprocessed
    #    independently. Do not write back automatically; review the reconciliation
    #    table first, especially when subset labels conflict with global labels or
    #    when subset review identifies cells to exclude.
    plan = scl.al.build_hierarchical_annotation_plan(
        adata,
        cluster_key="leiden_clusters",
        lineage_key="lineage_curated",
        min_cells_per_lineage=50,
        min_clusters_per_lineage=1,
    )
    plan.to_csv(output_dir / "hierarchical_annotation_plan.csv", index=False)

    subtype_mgr = scl.ut.get_marker_manager(species="human", tissue="Blood", view="subtype_annotation")
    subset_results = scl.al.run_subset_annotation_refinement(
        adata,
        lineage_key="lineage_curated",
        plan=plan,
        counts_layer="counts",
        marker_config=subtype_mgr,
        global_subtype_key="cell_type_curated",
        write_back=False,
        key_added="subset_annotation_refinement",
    )
    reconciliation = scl.al.build_subset_annotation_reconciliation(
        adata,
        subset_results,
        lineage_key="lineage_curated",
        global_subtype_key="cell_type_curated",
        key_added="subset_annotation_refinement",
    )
    reconciliation.to_csv(output_dir / "subset_annotation_reconciliation.csv", index=False)

    # 6) Downstream state programs with the thin workflow wrapper
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

    # 7) Composition plots from pre-aggregated tables
    #    These are visualization summaries. For statistical inference on cell
    #    proportions, prefer sample-level CLR/compositional tests.
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

    if {"sampleID", "group", "cell_type_curated"}.issubset(adata.obs.columns):
        prop_config = scl.al.ProportionConfig(
            celltype_col="cell_type_curated",
            sample_col="sampleID",
            condition_col="group",
            test_method="clr-t-test",
            out_dir=str(output_dir / "proportion_inference"),
        )
        try:
            prop_df, stat_df = scl.al.analyze_celltype_proportion(
                adata,
                method="pseudobulk",
                config=prop_config,
            )
        except Exception as exc:
            print(f"Skipped proportion inference example: {exc}")
        else:
            stat_df.to_csv(output_dir / "celltype_proportion_clr_stats.csv", index=False)
            if "inference_level" in stat_df:
                print(
                    "Proportion inference levels:",
                    sorted(stat_df["inference_level"].dropna().unique()),
                )

    # 8) Optional formal condition DE with sample-level pseudobulk.
    #    Add design_covariates/block_col when batch or paired patient metadata
    #    exist. This is the preferred route for publication-grade condition DE.
    if {"sampleID", "group", "cell_type_curated"}.issubset(adata.obs.columns):
        de_pb_config = scl.al.PseudobulkDEConfig(
            sample_col="sampleID",
            condition_key="group",
            groupby="cell_type_curated",
            contrasts=[("control", "treated")],
            min_cells_per_sample=10,
            method="linear_model_logcpm" if "batch" in adata.obs.columns else "auto",
            design_covariates=["batch"] if "batch" in adata.obs.columns else [],
            block_col="patient_id" if "patient_id" in adata.obs.columns else None,
            key_added="pseudobulk_condition_de",
        )
        try:
            de_pb = scl.al.run_pseudobulk_de(adata, de_pb_config)
        except Exception as exc:
            print(f"Skipped pseudobulk condition DE example: {exc}")
        else:
            de_pb.to_csv(output_dir / "pseudobulk_condition_de.csv", index=False)
            if not de_pb.empty:
                print("Pseudobulk DE inference levels:", sorted(de_pb["inference_level"].unique()))

    adata.write_h5ad(output_dir / "annotated_curated.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
