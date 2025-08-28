"""
Unified Differential Expression and Enrichment Analysis for scRNA-seq.

This module provides:
- Marker gene discovery and DE analysis (find_markers, compare_groups, etc.)
    - Scenario A: Find marker genes for each cluster (one-vs-rest). -find_markers
    - Scenario B: Compare two specific cell types (e.g., CD4+ T vs. CD8+ T). -compare_groups
    - Scenario C: Find markers for a cell type that are conserved across conditions -get_conserved_markers 
        (e.g., find robust T cell markers present in both Control and Treated samples).
    - Scenario D: Compare the same cell type across different conditions -compare_conditions
        (e.g., Treated T-cells vs. Control T-cells).
- DE result filtering (filter_markers)
- Functional enrichment analysis (run_enrichment)
- Combined summary export for AI/manual cell type annotation (summarize_markers_and_enrichment)
- Visualization and cluster characterization

All results and parameters are stored under adata.uns['scrnatk']['analysis']['de'].
"""

import logging
import os
from typing import Dict, List, Literal, Optional, Tuple, Union

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

log = logging.getLogger(__name__)

__all__ = [
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    "run_enrichment",
    "summarize_markers_and_enrichment",
    "characterize_clusters",
    "visualize_markers"
]


# --- Differential Expression Analysis ---

def find_markers(
    adata: AnnData,
    groupby: str,
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon",
    layer: Optional[str] = None,
    key_added: Optional[str] = None,
    use_raw: bool = False,
    min_cells: int = 5,
    groups: Optional[List[str]] = None,
    reference: Optional[str] = "rest",
    fold_change_max: Optional[float] = None,
    pval_cutoff: Optional[float] = None,
    copy: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes for all groups using Scanpy's sc.tl.rank_genes_groups.
    Results are stored in adata.uns['scrnatk']['analysis']['de'].
    Returns a DataFrame of all markers.
    """
    log.info(f"Finding marker genes for groups in '{groupby}' using '{method}' method")
    if copy:
        adata = adata.copy()
    if key_added is None:
        key_added = "rank_genes_groups"
    if groupby not in adata.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in adata.obs")

    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('de', {})
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        layer=layer,
        key_added=key_added,
        use_raw=use_raw,
        pts=True,
        groups=groups,
        reference=reference,
        **kwargs,
    )
    # Get marker results for all groups
    result_dfs = []
    categories = groups if groups is not None else adata.obs[groupby].cat.categories
    for group in categories:
        df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)
        if df.empty:
            continue
        df["group"] = group
        if fold_change_max is not None and fold_change_max > 0:
            df["logfoldchanges"] = df["logfoldchanges"].clip(upper=fold_change_max)
        if pval_cutoff is not None and pval_cutoff > 0:
            df = df[df["pvals_adj"] <= pval_cutoff].copy()
        result_dfs.append(df)
    if not result_dfs:
        log.warning("No valid marker results found for any group")
        return pd.DataFrame()
    full_df = pd.concat(result_dfs, ignore_index=True)
    # Store results
    df_key = f"{key_added}_df"
    adata.uns['scrnatk']['analysis']['de'][key_added] = adata.uns[key_added]
    adata.uns['scrnatk']['analysis']['de'][df_key] = full_df
    log.info(f"Found {len(full_df)} marker genes across {len(result_dfs)} groups. Stored in .uns['scrnatk']['analysis']['de']['{df_key}']")
    return full_df


def filter_markers(
    adata: AnnData,
    key: str = "rank_genes_groups",
    min_log2fc: float = 1.0,
    max_padj: float = 0.05,
    min_in_group_pct: float = 0.25,
    max_out_group_pct: Optional[float] = None,
    min_diff_pct: Optional[float] = None,
    keep_top_n: Optional[int] = None,
    key_added: Optional[str] = None,
) -> pd.DataFrame:
    """
    Filter marker genes by log2FC, adjusted p-value, expression percent, etc.
    Returns filtered DataFrame and stores in adata.uns.
    """
    if key_added is None:
        key_added = f"{key}_filtered_df"
    df_key = f"{key}_df"
    if 'scrnatk' in adata.uns and 'analysis' in adata.uns['scrnatk'] and \
        'de' in adata.uns['scrnatk']['analysis'] and df_key in adata.uns['scrnatk']['analysis']['de']:
        df = adata.uns['scrnatk']['analysis']['de'][df_key].copy()
    elif df_key in adata.uns:
        df = adata.uns[df_key].copy()
    else:
        raise KeyError(f"Results DataFrame not found at `adata.uns['{df_key}']`. Run `find_markers` first.")
    if df.empty:
        return pd.DataFrame()
    # Filters
    filt = (df["logfoldchanges"] >= min_log2fc) & (df["pvals_adj"] <= max_padj) & ((df["pct_nz_group"] / 100) >= min_in_group_pct)
    if max_out_group_pct is not None and "pct_nz_reference" in df.columns:
        filt &= (df["pct_nz_reference"] / 100 <= max_out_group_pct)
    if min_diff_pct is not None and "pct_nz_reference" in df.columns:
        filt &= ((df["pct_nz_group"] - df["pct_nz_reference"]) / 100 >= min_diff_pct)
    filtered_df = df[filt].copy()
    if keep_top_n is not None and keep_top_n > 0:
        top_per_group = []
        for group in filtered_df["group"].unique():
            group_df = filtered_df[filtered_df["group"] == group]
            top_group = group_df.sort_values("logfoldchanges", ascending=False).head(keep_top_n)
            top_per_group.append(top_group)
        filtered_df = pd.concat(top_per_group, ignore_index=True)
    adata.uns['scrnatk']['analysis']['de'][key_added] = filtered_df
    log.info(f"Filtered markers: {len(filtered_df)} genes stored in .uns['scrnatk']['analysis']['de']['{key_added}']")
    return filtered_df


def compare_groups(
    adata: AnnData,
    groupby: str,
    group1: str,
    group2: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    n_genes: int = 50,
    min_log2fc: float = 0.5,
    max_padj: float = 0.05,
    min_in_group_pct: float = 0.1,
    plot: bool = True,
    save_path: Optional[str] = None,
    key_added: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compare two groups (e.g., cell types or conditions) for DE genes.
    """
    if key_added is None:
        key_added = f"compare_{group1}_vs_{group2}"
    temp = adata[adata.obs[groupby].isin([group1, group2])].copy()
    temp.obs["compare_groups"] = temp.obs[groupby].map({group1: "group1", group2: "group2"}).astype("category")
    sc.tl.rank_genes_groups(
        temp,
        groupby="compare_groups",
        groups=["group1", "group2"],
        reference="rest",
        method="wilcoxon",
        layer=layer,
        use_raw=use_raw,
        pts=True,
    )
    group1_vs_group2 = sc.get.rank_genes_groups_df(temp, group="group1")
    group1_vs_group2["comparison"] = f"{group1}_vs_{group2}"
    group1_vs_group2["higher_in"] = group1
    group2_vs_group1 = sc.get.rank_genes_groups_df(temp, group="group2")
    group2_vs_group1["comparison"] = f"{group2}_vs_{group1}"
    group2_vs_group1["higher_in"] = group2
    all_results = pd.concat([group1_vs_group2, group2_vs_group1], ignore_index=True)
    filtered_results = all_results[
        (all_results["logfoldchanges"].abs() >= min_log2fc) &
        (all_results["pvals_adj"] <= max_padj) &
        (all_results["pct_nz_group"] / 100 >= min_in_group_pct)
    ].copy()
    if n_genes > 0:
        top_results = []
        for group in [group1, group2]:
            group_results = filtered_results[filtered_results["higher_in"] == group]
            top_results.append(group_results.head(n_genes))
        filtered_results = pd.concat(top_results, ignore_index=True)
    adata.uns['scrnatk']['analysis']['de'][key_added] = filtered_results
    # Volcano plot
    if plot:
        try:
            plt.figure(figsize=(12, 8))
            plt.scatter(
                all_results["logfoldchanges"],
                -np.log10(all_results["pvals_adj"] + 1e-10),
                alpha=0.3,
                s=10,
                color="grey",
                label="Not significant",
            )
            for group, color in zip([group1, group2], ["red", "blue"]):
                sig = filtered_results[filtered_results["higher_in"] == group]
                if not sig.empty:
                    plt.scatter(
                        sig["logfoldchanges"],
                        -np.log10(sig["pvals_adj"] + 1e-10),
                        alpha=0.7, s=30, color=color, label=f"Higher in {group}"
                    )
            plt.axvline(x=min_log2fc, color="grey", linestyle="--", alpha=0.5)
            plt.axvline(x=-min_log2fc, color="grey", linestyle="--", alpha=0.5)
            plt.axhline(y=-np.log10(max_padj), color="grey", linestyle="--", alpha=0.5)
            plt.xlabel("Log2 Fold Change")
            plt.ylabel("-log10(Adjusted p-value)")
            plt.title(f"Differential Expression: {group1} vs {group2}")
            plt.legend()
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300)
            else:
                plt.show()
        except Exception as e:
            log.warning(f"Error creating volcano plot: {str(e)}")
    return filtered_results


def compare_conditions(
    adata: AnnData,
    groupby: str,
    group_name: str,
    condition_key: str,
    condition1: str,
    condition2: str,
    key_added: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Compare two conditions within a specific group (e.g., cell type).
    """
    if group_name not in adata.obs[groupby].unique():
        raise ValueError(f"Group '{group_name}' not found in adata.obs['{groupby}']")
    adata_subset = adata[adata.obs[groupby] == group_name].copy()
    if key_added is None:
        key_added = f"compare_{condition1}_vs_{condition2}_in_{group_name.replace(' ', '_')}"
    return compare_groups(
        adata_subset,
        groupby=condition_key,
        group1=condition1,
        group2=condition2,
        key_added=key_added,
        **kwargs
    )


def get_conserved_markers(
    adata: AnnData,
    groupby: str,
    condition_key: str,
    method: str = "wilcoxon",
    min_cells: int = 10,
    min_conditions: Optional[int] = None,
    min_log2fc: float = 0.5,
    max_padj: float = 0.05,
    min_in_group_pct: float = 0.25,
    layer: Optional[str] = None,
    use_raw: bool = False,
    key_added: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Find marker genes for a group conserved across multiple conditions.
    Returns dict: group -> DataFrame
    """
    if key_added is None:
        key_added = f"conserved_markers_{groupby}_{condition_key}"
    conditions = adata.obs[condition_key].cat.categories
    groups = adata.obs[groupby].cat.categories
    if min_conditions is None:
        min_conditions = max(1, len(conditions) - 1)
    conserved_markers = {}
    for group in groups:
        markers_per_condition = []
        for condition in conditions:
            subset = adata[(adata.obs[groupby] == group) & (adata.obs[condition_key] == condition)]
            if subset.n_obs < min_cells:
                continue
            temp_adata = adata[adata.obs[condition_key] == condition].copy()
            if group not in temp_adata.obs[groupby].cat.categories:
                continue
            sc.tl.rank_genes_groups(
                temp_adata,
                groupby=groupby,
                groups=[group],
                reference="rest",
                method=method,
                layer=layer,
                use_raw=use_raw,
                pts=True,
            )
            df = sc.get.rank_genes_groups_df(temp_adata, group=group)
            df = df[
                (df["logfoldchanges"] >= min_log2fc) &
                (df["pvals_adj"] <= max_padj) &
                (df["pct_nz_group"] / 100 >= min_in_group_pct)
            ]
            if df.empty:
                continue
            df["condition"] = condition
            markers_per_condition.append(df)
        if len(markers_per_condition) < min_conditions:
            continue
        full_df = pd.concat(markers_per_condition)
        gene_counts = full_df.groupby("names").size()
        conserved_genes = gene_counts[gene_counts >= min_conditions].index.tolist()
        conserved_df = full_df[full_df["names"].isin(conserved_genes)]
        agg_df = (
            conserved_df.groupby("names")
            .agg(
                mean_log2fc=("logfoldchanges", "mean"),
                min_log2fc=("logfoldchanges", "min"),
                max_log2fc=("logfoldchanges", "max"),
                min_pval_adj=("pvals_adj", "min"),
                mean_pval_adj=("pvals_adj", "mean"),
                mean_pct_in_group=("pct_nz_group", "mean"),
                n_conditions=("condition", "nunique"),
            )
            .sort_values("mean_log2fc", ascending=False)
        )
        conserved_markers[group] = agg_df
    adata.uns['scrnatk']['analysis']['de'][key_added] = conserved_markers
    return conserved_markers

# --- Enrichment Analysis ---

def run_enrichment(
    adata: AnnData,
    groupby: str,
    de_key: str = "rank_genes_groups",
    organism: str = "Human",
    gene_sets: List[str] = ["GO_Biological_Process_2023"],
    n_top_genes: int = 100,
    key_added: str = "enrichment",
    min_genes: int = 10,
    max_genes: int = 500,
    min_enrichment_score: float = 0.0,
    max_padj: float = 0.05,
    background_genes: Optional[List[str]] = None,
    plot: bool = False,
    save_path: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Enrichment analysis for each group using GSEApy (Enrichr). Returns dict: group -> enrichment DataFrame.
    Results stored in adata.uns['scrnatk']['analysis']['de'][key_added]
    """
    log.info(f"Running enrichment analysis for {groupby} groups")
    # Get marker_df
    de_results_key = f"{de_key}_df"
    if 'scrnatk' in adata.uns and 'analysis' in adata.uns['scrnatk'] and \
        'de' in adata.uns['scrnatk']['analysis'] and de_results_key in adata.uns['scrnatk']['analysis']['de']:
        marker_df = adata.uns['scrnatk']['analysis']['de'][de_results_key]
    elif de_results_key in adata.uns:
        marker_df = adata.uns[de_results_key]
    else:
        raise KeyError(f"DE results not found. Run `find_markers` first.")
    clusters = marker_df["group"].unique()
    if background_genes is None:
        background_genes = list(adata.var_names)
    if plot and save_path is not None:
        os.makedirs(save_path, exist_ok=True)
    enrichment_results = {}
    for cluster in clusters:
        try:
            cluster_df = marker_df[marker_df["group"] == cluster]
            gene_list = cluster_df.sort_values("logfoldchanges", ascending=False)["names"].head(n_top_genes).tolist()
            if len(gene_list) < min_genes:
                continue
            if len(gene_list) > max_genes:
                gene_list = gene_list[:max_genes]
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets,
                organism=organism,
                background=background_genes,
                outdir=None,
                cutoff=max_padj,
            )
            results = enr.results
            if results.empty:
                enrichment_results[cluster] = pd.DataFrame()
                continue
            if min_enrichment_score > 0:
                results = results[results["Combined Score"] >= min_enrichment_score]
            enrichment_results[cluster] = results
            if plot:
                top_pathways = results.head(20)
                if not top_pathways.empty:
                    plt.figure(figsize=(12, min(10, max(4, len(top_pathways) * 0.3))))
                    plt.barh(
                        top_pathways["Term"].str.split(" ").str[:5].str.join(" "),
                        -np.log10(top_pathways["Adjusted P-value"]),
                        color="skyblue",
                    )
                    plt.xlabel("-log10(Adjusted P-value)")
                    plt.title(f"Top Enriched Pathways for {cluster}")
                    plt.tight_layout()
                    if save_path is not None:
                        safe_cluster = str(cluster).replace("/", "_").replace(" ", "_")
                        plt.savefig(f"{save_path}/{safe_cluster}_enrichment.png", dpi=300)
                        plt.close()
                    else:
                        plt.show()
        except Exception as e:
            log.warning(f"Error analyzing group '{cluster}': {str(e)}")
            enrichment_results[cluster] = pd.DataFrame()
    adata.uns['scrnatk']['analysis']['de'][key_added] = {
        'results': enrichment_results,
        'params': {
            'groupby': groupby,
            'de_key': de_key,
            'gene_sets': gene_sets,
            'organism': organism,
        }
    }
    return enrichment_results

# --- Marker + Enrichment Summary for AI/manual annotation ---

def summarize_markers_and_enrichment(
    adata: AnnData,
    markers_df: pd.DataFrame,
    enrichment_dict: Dict[str, pd.DataFrame],
    groupby: str = "leiden",
    n_markers: int = 10,
    n_terms: int = 5,
    summary_file: Optional[str] = None,
) -> Dict[str, str]:
    """
    Export each group's top marker and enrichment summary (markdown), for AI/manual annotation.
    Returns dict: {cluster: markdown_str}
    """
    summary = {}
    groups = markers_df["group"].unique()
    lines = []
    for g in groups:
        top_genes = markers_df[markers_df["group"] == g].sort_values("logfoldchanges", ascending=False).head(n_markers)["names"].tolist()
        top_terms = (
            enrichment_dict.get(g, pd.DataFrame())
            .sort_values("Adjusted P-value")
            .head(n_terms)["Term"].tolist()
            if g in enrichment_dict and not enrichment_dict[g].empty else []
        )
        s = f"### Cluster {g}\nTop markers: {', '.join(top_genes)}\nTop pathways: {', '.join(top_terms)}"
        summary[g] = s
        lines.append(s)
    if summary_file:
        with open(summary_file, "w") as f:
            f.write("\n\n".join(lines))
    log.info(f"Marker + enrichment summaries exported to {summary_file or '[not saved]'}")
    return summary

# --- Batch Visualization and Cluster Characterization ---

def visualize_markers(
    adata: AnnData,
    markers: Union[pd.DataFrame, Dict[str, List[str]], List[str]],
    groupby: Optional[str] = None,
    n_genes_per_group: int = 5,
    plot_type: Literal[
        "dotplot", "heatmap", "stacked_violin", "violin", "matrixplot"
    ] = "dotplot",
    dendrogram: bool = False,
    standard_scale: Optional[Literal["var", "group"]] = None,
    swap_axes: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    save_path: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """
    Visualize marker genes across groups using various plot types.
    """
    gene_list = []
    if isinstance(markers, pd.DataFrame):
        if "group" not in markers.columns or "names" not in markers.columns:
            raise ValueError("DataFrame must contain 'group' and 'names' columns")
        if n_genes_per_group > 0:
            for group in markers["group"].unique():
                group_markers = markers[markers["group"] == group]
                if "logfoldchanges" in group_markers.columns:
                    group_markers = group_markers.sort_values("logfoldchanges", ascending=False)
                top_genes = group_markers["names"].head(n_genes_per_group).tolist()
                gene_list.extend(top_genes)
        else:
            gene_list = markers["names"].tolist()
    elif isinstance(markers, dict):
        if n_genes_per_group > 0:
            for genes in markers.values():
                gene_list.extend(genes[:n_genes_per_group])
        else:
            for genes in markers.values():
                gene_list.extend(genes)
    elif isinstance(markers, (list, tuple)):
        gene_list = list(markers)
        if groupby is None:
            raise ValueError("groupby must be specified when markers is a list")
    else:
        raise TypeError("markers must be a DataFrame, dictionary, or list")
    gene_list = list(dict.fromkeys(gene_list))  # Remove duplicates
    # Remove genes not in dataset
    gene_list = [gene for gene in gene_list if gene in adata.var_names]
    if not gene_list:
        raise ValueError("No valid genes found for visualization")
    # Plot
    if figsize is None:
        n_groups = len(adata.obs[groupby].cat.categories) if groupby and groupby in adata.obs else 1
        n_genes = len(gene_list)
        figsize = (max(6, min(12, n_genes * 0.5)), max(4, min(10, n_groups * 0.3)))
    if plot_type == "dotplot":
        sc.pl.dotplot(adata, var_names=gene_list, groupby=groupby, dendrogram=dendrogram,
                      standard_scale=standard_scale, swap_axes=swap_axes, use_raw=use_raw,
                      layer=layer, figsize=figsize, **kwargs)
    elif plot_type == "heatmap":
        sc.pl.heatmap(adata, var_names=gene_list, groupby=groupby, dendrogram=dendrogram,
                      standard_scale=standard_scale, swap_axes=swap_axes, use_raw=use_raw,
                      layer=layer, figsize=figsize, **kwargs)
    elif plot_type == "stacked_violin":
        sc.pl.stacked_violin(adata, var_names=gene_list, groupby=groupby, dendrogram=dendrogram,
                             standard_scale=standard_scale, swap_axes=swap_axes, use_raw=use_raw,
                             layer=layer, figsize=figsize, **kwargs)
    elif plot_type == "violin":
        sc.pl.violin(adata, keys=gene_list, groupby=groupby, use_raw=use_raw, layer=layer, figsize=figsize, **kwargs)
    elif plot_type == "matrixplot":
        sc.pl.matrixplot(adata, var_names=gene_list, groupby=groupby, dendrogram=dendrogram,
                         standard_scale=standard_scale, swap_axes=swap_axes, use_raw=use_raw,
                         layer=layer, figsize=figsize, **kwargs)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")

def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_method: str = "wilcoxon",
    n_de_genes: int = 25,
    enrich_organism: str = "Human",
    enrich_gene_sets: List[str] = ["GO_Biological_Process_2023"],
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Run DE and enrichment analysis for each cluster to gather evidence for annotation.
    Stores results under adata.uns[key_added].
    """
    log.info(f"Characterizing clusters in '{groupby}'...")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=de_method, use_raw=True, n_genes=n_de_genes)
    clusters = adata.obs[groupby].cat.categories
    characterization_results = {}
    for cluster in clusters:
        de_df = sc.get.rank_genes_groups_df(adata, group=cluster)
        top_genes = de_df.head(n_de_genes)['names'].tolist()
        try:
            enr = gp.enrichr(
                gene_list=top_genes,
                gene_sets=enrich_gene_sets,
                organism=enrich_organism,
                outdir=None,
            )
            enrich_df = enr.results
        except Exception as e:
            log.warning(f"Enrichment analysis failed for cluster {cluster}: {e}")
            enrich_df = pd.DataFrame()
        characterization_results[cluster] = {
            "top_de_genes": de_df.head(n_de_genes),
            "enrichment": enrich_df,
        }
    adata.uns[key_added] = characterization_results
    log.info(f"Cluster characterization complete. Results stored in adata.uns['{key_added}']")
    return adata