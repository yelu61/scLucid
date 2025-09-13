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

All results and parameters are stored under adata.uns['sclucid']['analysis']['de'].
"""

import dataclasses
import logging
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    DifferentialConfig,
    EnrichmentConfig,
    FilterMarkersConfig,
)

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
    "visualize_markers",
]


# --- Differential Expression Analysis ---
def find_markers(
    adata: AnnData,
    config: Optional[DifferentialConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes using a flexible, config-driven workflow.

    This function performs a one-vs-rest differential expression analysis and stores
    both the raw Scanpy output and a processed DataFrame in adata.uns.
    The entire process is controlled by the DifferentialConfig object.

    Args:
        adata: The AnnData object.
        config: A DifferentialConfig object containing all parameters.
        **kwargs: Allows overriding config parameters for interactive use.

    Returns:
        A pandas DataFrame containing the formatted marker genes for all groups.
    """
    # --- 1. Finalize Configuration ---
    if config is None:
        active_config = DifferentialConfig(**kwargs)
    else:
        # Use a copy to avoid modifying the original config object
        active_config = dataclasses.replace(config)

    # Allow kwargs to override parameters in the provided config object
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    # Unpack parameters for use
    groupby = active_config.groupby
    key_added = active_config.key_added or "rank_genes_groups"

    log.info(f"Finding markers for '{groupby}' using method '{active_config.method}'.")

    # --- 2. Run Scanpy's rank_genes_groups ---
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=active_config.method,
        layer=active_config.layer,
        key_added=key_added,
        use_raw=active_config.use_raw,
        pts=True,
        reference=active_config.reference,
        groups=active_config.groups,
    )

    # --- 3. Format, Filter, and Store Results ---
    result_dfs = []
    # Determine the groups that were actually tested
    groups_tested = adata.uns[key_added]["names"].dtype.names

    for group in groups_tested:
        df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)
        if df.empty:
            continue
        df["group"] = group

        if active_config.pval_cutoff is not None:
            df = df[df["pvals_adj"] <= active_config.pval_cutoff].copy()
            if df.empty:
                continue  # Skip if no genes pass the p-value filter

        if active_config.fold_change_max is not None:
            df["logfoldchanges"] = df["logfoldchanges"].clip(
                upper=active_config.fold_change_max
            )

        result_dfs.append(df)

    if not result_dfs:
        log.warning("No valid marker results found for any group after filtering.")
        return pd.DataFrame()

    full_df = pd.concat(result_dfs, ignore_index=True)

    # Store results in the sclucid namespace
    df_key = f"{key_added}_df"
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    adata.uns["sclucid"]["analysis"]["de"][key_added] = adata.uns[
        key_added
    ]  # Store raw scanpy result
    adata.uns["sclucid"]["analysis"]["de"][df_key] = (
        full_df  # Store processed DataFrame
    )

    log.info(
        f"Found {len(full_df)} total marker gene entries across {len(groups_tested)} groups."
    )
    log.info(f"Results stored in .uns['sclucid']['analysis']['de']['{df_key}']")

    return full_df


def filter_markers(
    adata: AnnData,
    config: FilterMarkersConfig,
) -> pd.DataFrame:
    """
    Filter marker genes based on a comprehensive configuration object.

    This function applies a series of filters to differential expression results
    to identify high-confidence marker genes. Results are stored back into adata.uns.

    Args:
        adata: The AnnData object.
        config: A FilterMarkersConfig object containing all filtering parameters.

    Returns:
        A pandas DataFrame containing the filtered marker genes.
    """
    # --- 1. Unpack parameters from config object ---
    key = config.key
    min_log2fc = config.min_log2fc
    max_padj = config.max_padj
    min_in_group_pct = config.min_in_group_pct
    max_out_group_pct = config.max_out_group_pct
    min_diff_pct = config.min_diff_pct
    keep_top_n = config.keep_top_n

    key_added = config.key_added or f"{key}_filtered_df"
    df_key = f"{key}_df"

    log.info(f"Filtering marker genes from '.uns[...][{df_key}]'")

    # --- 2. Retrieve source DataFrame ---
    if (
        "sclucid" in adata.uns
        and "analysis" in adata.uns["sclucid"]
        and "de" in adata.uns["sclucid"]["analysis"]
        and df_key in adata.uns["sclucid"]["analysis"]["de"]
    ):
        df = adata.uns["sclucid"]["analysis"]["de"][df_key].copy()
    else:
        raise KeyError(
            f"Source DataFrame not found at `adata.uns['sclucid']['analysis']['de']['{df_key}']`. Run `find_markers` first."
        )

    if df.empty:
        log.warning("Source marker DataFrame is empty. Returning empty DataFrame.")
        return pd.DataFrame()

    # --- 3. Apply Filters ---
    # Start with a mask of all True
    filt = pd.Series(True, index=df.index)

    # Core filters
    if min_log2fc is not None:
        filt &= df["logfoldchanges"] >= min_log2fc
    if max_padj is not None:
        filt &= df["pvals_adj"] <= max_padj
    if min_in_group_pct is not None:
        # Note: Scanpy's pct_nz_group is already in percentage, so we don't divide by 100
        filt &= df["pct_nz_group"] >= (min_in_group_pct * 100)

    # Specificity filters
    if "pct_nz_reference" in df.columns:
        if max_out_group_pct is not None:
            filt &= df["pct_nz_reference"] <= (max_out_group_pct * 100)
        if min_diff_pct is not None:
            filt &= (df["pct_nz_group"] - df["pct_nz_reference"]) >= (
                min_diff_pct * 100
            )
    else:
        if max_out_group_pct is not None or min_diff_pct is not None:
            log.warning(
                "'pct_nz_reference' not found in marker DataFrame. Skipping specificity filters ('max_out_group_pct', 'min_diff_pct')."
            )

    filtered_df = df[filt].copy()
    log.info(
        f"Retained {len(filtered_df)} genes after applying statistical and specificity filters."
    )

    # --- 4. Keep Top N genes per group ---
    if keep_top_n is not None and keep_top_n > 0:
        # We sort by 'scores' if available (like from logreg), otherwise fall back to logfoldchanges.
        sort_by_col = "scores" if "scores" in filtered_df.columns else "logfoldchanges"

        log.info(
            f"Selecting top {keep_top_n} genes per group, sorted by '{sort_by_col}'."
        )
        top_per_group = []
        for group in filtered_df["group"].unique():
            group_df = filtered_df[filtered_df["group"] == group]
            top_group = group_df.sort_values(sort_by_col, ascending=False).head(
                keep_top_n
            )
            top_per_group.append(top_group)

        if top_per_group:
            filtered_df = pd.concat(top_per_group, ignore_index=True)
        else:
            filtered_df = pd.DataFrame()  # Handle case where no groups have genes left

    # --- 5. Store Results ---
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})[
        key_added
    ] = filtered_df
    log.info(
        f"Final filtered markers: {len(filtered_df)} genes stored in .uns['sclucid']['analysis']['de']['{key_added}']"
    )

    return filtered_df


def compare_groups(adata: AnnData, config: CompareGroupsConfig) -> pd.DataFrame:
    """
    Compare two groups (e.g., cell types or conditions) for DE genes using a config object.
    """
    # --- 1. Unpack parameters from config ---
    groupby = config.groupby
    group1 = config.group1
    group2 = config.group2
    key_added = config.key_added or f"compare_{group1}_vs_{group2}".replace(" ", "_")

    log.info(
        f"Comparing DE genes between '{group1}' and '{group2}' from column '{groupby}'."
    )

    # --- 2. Subset data and run DE ---
    # Use a view and a temporary column to avoid large data copies
    adata_view = adata[adata.obs[groupby].isin([group1, group2])]

    # Create a temporary AnnData for DE to avoid modifying the view's .obs
    temp_adata = sc.AnnData(
        X=adata_view.X, obs=adata_view.obs.copy(), var=adata_view.var
    )
    temp_adata.obs["_compare_groups"] = (
        temp_adata.obs[groupby]
        .map({group1: "group1", group2: "group2"})
        .astype("category")
    )

    sc.tl.rank_genes_groups(
        temp_adata,
        groupby="_compare_groups",
        groups=["group1"],  # Compare group1 against group2 as reference
        reference="group2",
        method=config.method,
        layer=config.layer,
        use_raw=config.use_raw,
        pts=True,
    )

    # --- 3. Format and Filter Results ---
    results_df = sc.get.rank_genes_groups_df(temp_adata, group="group1")
    results_df.rename(columns={"names": "gene"}, inplace=True)

    # Apply filters
    filt = (
        (results_df["logfoldchanges"].abs() >= config.min_log2fc)
        & (results_df["pvals_adj"] <= config.max_padj)
        & (results_df["pct_nz_group"] / 100 >= config.min_in_group_pct)
    )

    filtered_results = results_df[filt].copy()

    # Sort and take top N
    filtered_results.sort_values("logfoldchanges", ascending=False, inplace=True)
    upregulated = filtered_results[filtered_results["logfoldchanges"] > 0].head(
        config.n_top_genes
    )
    downregulated = filtered_results[filtered_results["logfoldchanges"] < 0].tail(
        config.n_top_genes
    )
    final_results = pd.concat([upregulated, downregulated], ignore_index=True)

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})[
        key_added
    ] = final_results
    log.info(
        f"Found {len(final_results)} DE genes. Stored in .uns['...']['{key_added}']"
    )

    # --- 4. Volcano Plot ---
    if config.plot:
        # ... (Volcano plot logic remains largely the same, but uses the config) ...
        plt.figure(figsize=(10, 7))
        plt.scatter(
            results_df["logfoldchanges"],
            -np.log10(results_df["pvals_adj"].clip(1e-300)),
            alpha=0.3,
            s=10,
            color="grey",
            label="Not significant",
            rasterized=True,
        )
        plt.scatter(
            final_results[final_results["logfoldchanges"] > 0]["logfoldchanges"],
            -np.log10(
                final_results[final_results["logfoldchanges"] > 0]["pvals_adj"].clip(
                    1e-300
                )
            ),
            alpha=0.7,
            s=30,
            color="red",
            label=f"Higher in {group1}",
        )
        plt.scatter(
            final_results[final_results["logfoldchanges"] < 0]["logfoldchanges"],
            -np.log10(
                final_results[final_results["logfoldchanges"] < 0]["pvals_adj"].clip(
                    1e-300
                )
            ),
            alpha=0.7,
            s=30,
            color="blue",
            label=f"Higher in {group2}",
        )
        plt.axvline(x=config.min_log2fc, color="grey", linestyle="--", alpha=0.5)
        plt.axvline(x=-config.min_log2fc, color="grey", linestyle="--", alpha=0.5)
        plt.axhline(
            y=-np.log10(config.max_padj), color="grey", linestyle="--", alpha=0.5
        )
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-log10(Adjusted p-value)")
        plt.title(f"Differential Expression: {group1} vs {group2}")
        plt.legend()
        plt.tight_layout()
        if config.save_dir:
            Path(config.save_dir).mkdir(parents=True, exist_ok=True)
            plt.savefig(Path(config.save_dir) / f"{key_added}_volcano.png", dpi=300)
        plt.show()

    return final_results


def compare_conditions(
    adata: AnnData,
    config: CompareConditionsConfig,
) -> pd.DataFrame:
    """
    Compare two conditions within a specific group using a config object.
    """
    # Unpack parameters
    groupby = config.groupby
    group_name = config.group_name
    condition_key = config.condition_key

    log.info(
        f"Comparing conditions '{config.condition1}' vs '{config.condition2}' within cell group '{group_name}'."
    )

    # --- 1. Subset Data ---
    if group_name not in adata.obs[groupby].unique():
        raise ValueError(f"Group '{group_name}' not found in adata.obs['{groupby}']")
    adata_subset = adata[adata.obs[groupby] == group_name].copy()

    # --- 2. Prepare and run compare_groups ---
    # Use the nested config object for the comparison
    comp_config = config.comparison_params
    comp_config.groupby = condition_key  # Set the correct groupby for the subset
    comp_config.group1 = config.condition1
    comp_config.group2 = config.condition2

    # Set a descriptive key_added if not provided
    if config.key_added is None:
        safe_group = group_name.replace(" ", "_").replace("/", "_")
        comp_config.key_added = (
            f"compare_{config.condition1}_vs_{config.condition2}_in_{safe_group}"
        )
    else:
        comp_config.key_added = config.key_added

    results_df = compare_groups(adata_subset, config=comp_config)

    # Store results back into the original adata object for global context
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})[
        comp_config.key_added
    ] = results_df
    log.info(
        f"Stored condition comparison results in original adata.uns['...']['{comp_config.key_added}']"
    )

    return results_df


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
            subset = adata[
                (adata.obs[groupby] == group) & (adata.obs[condition_key] == condition)
            ]
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
                (df["logfoldchanges"] >= min_log2fc)
                & (df["pvals_adj"] <= max_padj)
                & (df["pct_nz_group"] / 100 >= min_in_group_pct)
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
    adata.uns["sclucid"]["analysis"]["de"][key_added] = conserved_markers
    return conserved_markers


# --- Enrichment Analysis ---


def run_enrichment(
    adata: AnnData,
    groupby: str,
    config: EnrichmentConfig,
) -> Dict[str, pd.DataFrame]:
    """
    Enrichment analysis for each group using GSEApy, driven by a configuration object.
    Supports online (Enrichr) and offline (local GMT files) modes.
    """
    # Unpack parameters from the config object
    de_key = config.de_key
    mode = config.mode
    organism = config.organism
    if mode == "online":
        gene_sets_to_use = config.gene_sets_online
    else:  # offline mode
        gene_sets_to_use = config.gene_sets_offline
    gmt_version = config.gmt_version
    custom_gene_sets = config.custom_gene_sets
    n_top_genes = config.n_top_genes
    key_added = config.key_added
    min_genes_for_enrichment = config.min_genes_for_enrichment
    max_padj = config.max_padj
    plot = config.plot
    save_dir = config.save_dir
    n_plot_terms = config.n_plot_terms

    log.info(f"Running enrichment analysis for '{groupby}' groups in '{mode}' mode.")

    de_results_key = (
        de_key  # Assumes de_key is the full key, e.g., 'rank_genes_groups_filtered_df'
    )
    if (
        "sclucid" in adata.uns
        and "analysis" in adata.uns["sclucid"]
        and "de" in adata.uns["sclucid"]["analysis"]
        and de_results_key in adata.uns["sclucid"]["analysis"]["de"]
    ):
        marker_df = adata.uns["sclucid"]["analysis"]["de"][de_results_key]
    else:
        raise KeyError(
            f"DE results not found at .uns['sclucid']['analysis']['de']['{de_results_key}']. Run find_markers/filter_markers first."
        )

    clusters = marker_df["group"].unique()
    background_genes = list(adata.var_names)

    if plot and save_dir:
        Path(save_dir).mkdir(parents=True, exist_ok=True)

    enrichment_results = {}

    # --- Prepare gene_sets for analysis ---
    if isinstance(gene_sets_to_use, str):
        gene_sets_to_use = [gene_sets_to_use]

    # --- Load gene set file(s) for offline mode ---
    gmt_files_to_run = {}
    if mode == "offline":
        if custom_gene_sets and Path(custom_gene_sets).is_file():
            gmt_files_to_run = {"custom": custom_gene_sets}
            log.info(f"Using single custom gene set file from: {custom_gene_sets}")
        else:
            gmt_version = config.gmt_version
            for gs_category in gene_sets_to_use:
                try:
                    filename = f"{organism.lower()}_{gs_category}_{gmt_version}.gmt"
                    file_path = resources.files("scLucid").joinpath(
                        "resources", filename
                    )
                    if file_path.is_file():
                        gmt_files_to_run[gs_category] = str(file_path)
                    else:
                        log.warning(f"Resource file not found: {filename}")
                except Exception as e:
                    log.error(f"Error finding resource file for {gs_category}: {e}")
            if not gmt_files_to_run:
                raise FileNotFoundError(
                    "No valid gene set files found for offline mode."
                )
            log.info(
                f"Loaded {len(gmt_files_to_run)} local gene set(s) for offline analysis: {list(gmt_files_to_run.keys())}"
            )

    for cluster in clusters:
        try:
            gene_list = (
                marker_df[marker_df["group"] == cluster]
                .sort_values("logfoldchanges", ascending=False)["names"]
                .head(n_top_genes)
                .tolist()
            )

            if len(gene_list) < min_genes_for_enrichment:
                log.warning(
                    f"Skipping cluster {cluster}: Not enough genes ({len(gene_list)}) for enrichment."
                )
                enrichment_results[cluster] = pd.DataFrame()
                continue

            if mode == "online":
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=gene_sets_to_use,
                    organism=organism,
                    background=len(background_genes),
                    outdir=None,
                    cutoff=max_padj,
                )
                results = enr.results
            else:  # Offline mode
                all_offline_results = []
                for category, gmt_file in gmt_files_to_run.items():
                    enr = gp.enrich(
                        gene_list=gene_list,
                        gene_sets=gmt_file,
                        background=len(background_genes),
                        outdir=None,
                        cutoff=max_padj,
                    )
                    # Add the category to the results for clarity
                    enr.results["Gene_set"] = category
                    all_offline_results.append(enr.results)
                results = pd.concat(all_offline_results, ignore_index=True)

            sorted_results = (
                results.sort_values("Adjusted P-value")
                if not results.empty
                else pd.DataFrame()
            )
            enrichment_results[cluster] = sorted_results

            if plot and not sorted_results.empty:
                top_pathways = sorted_results.head(n_plot_terms)
                plt.figure(figsize=(10, max(4, len(top_pathways) * 0.4)))
                # Truncate long term names for better visualization
                y_labels = [
                    term[:70] + "..." if len(term) > 70 else term
                    for term in top_pathways["Term"]
                ]
                plt.barh(
                    y_labels,
                    -np.log10(top_pathways["Adjusted P-value"]),
                    color="skyblue",
                )
                plt.xlabel("-log10(Adjusted P-value)")
                plt.title(f"Top {len(top_pathways)} Enriched Pathways for {cluster}")
                plt.gca().invert_yaxis()
                plt.tight_layout()
                if save_dir:
                    safe_cluster = str(cluster).replace("/", "_").replace(" ", "_")
                    plt.savefig(
                        Path(save_dir) / f"{safe_cluster}_enrichment.png",
                        dpi=300,
                        bbox_inches="tight",
                    )
                plt.close()

        except Exception as e:
            log.warning(f"Enrichment failed for group '{cluster}': {str(e)}")
            enrichment_results[cluster] = pd.DataFrame()

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})[
        key_added
    ] = {
        "results": enrichment_results,
        "params": config.to_dict(),  # Store the full config for reproducibility
    }
    return enrichment_results


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
                    group_markers = group_markers.sort_values(
                        "logfoldchanges", ascending=False
                    )
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
        n_groups = (
            len(adata.obs[groupby].cat.categories)
            if groupby and groupby in adata.obs
            else 1
        )
        n_genes = len(gene_list)
        figsize = (max(6, min(12, n_genes * 0.5)), max(4, min(10, n_groups * 0.3)))
    if plot_type == "dotplot":
        sc.pl.dotplot(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "heatmap":
        sc.pl.heatmap(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "stacked_violin":
        sc.pl.stacked_violin(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "violin":
        sc.pl.violin(
            adata,
            keys=gene_list,
            groupby=groupby,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "matrixplot":
        sc.pl.matrixplot(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Run DE and enrichment analysis for each cluster to gather evidence for annotation.
    This function is now a high-level wrapper around find_markers and run_enrichment.
    """
    log.info(f"Characterizing clusters in '{groupby}'...")

    # --- 1. Run Differential Expression ---
    if de_config is None:
        de_config = DifferentialConfig(
            groupby=groupby, use_raw=True
        )  # Sensible defaults
    else:
        de_config.groupby = groupby  # Ensure groupby is consistent

    de_key = de_config.key_added or "rank_genes_groups"
    de_df_key = f"{de_key}_df"
    find_markers(adata, config=de_config)

    # --- 2. Run Enrichment Analysis ---
    # This now leverages our robust, offline-capable run_enrichment function
    if enrichment_config is None:
        # Create a default config that points to the DE results we just generated
        enrichment_config = EnrichmentConfig(de_key=de_df_key)
    else:
        enrichment_config.de_key = de_df_key  # Ensure it uses the correct DE results

    enrichment_results = run_enrichment(
        adata, groupby=groupby, config=enrichment_config
    )

    # --- 3. Consolidate Results ---
    clusters = adata.obs[groupby].cat.categories
    characterization_results = {}
    de_df = adata.uns["sclucid"]["analysis"]["de"][de_df_key]

    for cluster in clusters:
        characterization_results[cluster] = {
            "top_de_genes": de_df[de_df["group"] == cluster],
            "enrichment": enrichment_results.get(cluster, pd.DataFrame()),
        }

    adata.uns[key_added] = characterization_results
    log.info(
        f"Cluster characterization complete. Results stored in adata.uns['{key_added}']"
    )
    return adata


# --- Marker + Enrichment Summary for AI/manual annotation ---


def summarize_markers_and_enrichment(
    adata: AnnData,
    groupby: str,
    markers_df: Optional[pd.DataFrame] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    markers_key: str = "rank_genes_groups_df",
    enrichment_key: str = "enrichment",
    n_markers: int = 10,
    n_terms: int = 5,
    summary_file: Optional[str] = None,
) -> Dict[str, str]:
    """
    Export each group's top marker and enrichment summary (markdown), for AI/manual annotation.

    This function can now automatically retrieve marker and enrichment data from adata.uns
    if they are not provided directly.

    Args:
        adata: The AnnData object.
        groupby: The key in adata.obs for which the analysis was run.
        markers_df: (Optional) A DataFrame with marker genes. If None, tries to load from adata.uns.
        enrichment_dict: (Optional) A dictionary with enrichment results. If None, tries to load from adata.uns.
        markers_key: The key for the marker DataFrame in adata.uns['sclucid']['analysis']['de'].
        enrichment_key: The key for the enrichment results in adata.uns['sclucid']['analysis']['de'].
        n_markers: Number of top markers to include in the summary.
        n_terms: Number of top enrichment terms to include.
        summary_file: (Optional) Path to save the markdown summary file.

    Returns:
        A dictionary mapping each cluster/group to its markdown summary string.
    """
    # If markers_df is not provided, try to get it from the standard location in adata.uns.
    if markers_df is None:
        try:
            log.info(
                f"Attempting to auto-retrieve markers from .uns using key: {markers_key}"
            )
            markers_df = adata.uns["sclucid"]["analysis"]["de"][markers_key]
        except KeyError:
            raise KeyError(
                f"Marker DataFrame not found at .uns['sclucid']['analysis']['de']['{markers_key}']. Please run find_markers or provide the DataFrame directly."
            )

    # If enrichment_dict is not provided, do the same.
    if enrichment_dict is None:
        try:
            log.info(
                f"Attempting to auto-retrieve enrichment results from .uns using key: {enrichment_key}"
            )
            enrichment_dict = adata.uns["sclucid"]["analysis"]["de"][enrichment_key][
                "results"
            ]
        except KeyError:
            log.warning(
                f"Enrichment dictionary not found at .uns['sclucid']['analysis']['de']['{enrichment_key}']. Summary will not include pathways."
            )
            enrichment_dict = {}  # Use an empty dict to prevent errors

    summary = {}
    # Ensure groups are taken from the provided groupby column for consistency
    groups = (
        adata.obs[groupby].cat.categories
        if groupby in adata.obs
        and pd.api.types.is_categorical_dtype(adata.obs[groupby])
        else adata.obs[groupby].unique()
    )

    lines = []
    for g in groups:
        # Gracefully handle cases where a group might not be in the marker df (e.g., if filtered out)
        if g not in markers_df["group"].unique():
            continue

        top_genes = (
            markers_df[markers_df["group"] == g]
            .sort_values("logfoldchanges", ascending=False)
            .head(n_markers)["names"]
            .tolist()
        )

        enrichment_df_for_group = enrichment_dict.get(g, pd.DataFrame())
        top_terms = (
            enrichment_df_for_group.sort_values("Adjusted P-value")
            .head(n_terms)["Term"]
            .tolist()
            if not enrichment_df_for_group.empty
            else []
        )

        s = f"### Cluster {g}\n**Top Markers**: {', '.join(top_genes)}\n**Top Pathways**: {', '.join(top_terms)}"
        summary[g] = s
        lines.append(s)

    if summary_file:
        Path(summary_file).parent.mkdir(parents=True, exist_ok=True)
        with open(summary_file, "w") as f:
            f.write("\n\n---\n\n".join(lines))
        log.info(f"Marker + enrichment summaries exported to {summary_file}")

    return summary
