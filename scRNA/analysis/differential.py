"""
Differential expression analysis functions for single-cell RNA-seq data.

This module provides functions for identifying marker genes and
performing enrichment analysis for cell clusters.
"""

from typing import Dict, Literal, Optional

import pandas as pd
import scanpy as sc

# --- Main DE Functions ---


def find_markers(
    adata: sc.AnnData,
    groupby: str,
    method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon",
    layer: Optional[str] = None,
    key_added: Optional[str] = None,
    use_raw: bool = False,
) -> pd.DataFrame:
    """
    Find marker genes for all groups using `sc.tl.rank_genes_groups`.

    Args:
        adata: AnnData object.
        groupby: Key in `adata.obs` for grouping cells.
        method: Statistical method for differential testing.
        layer: Layer to use for expression values.
        key_added: Key in `adata.uns` to store results. Defaults to 'rank_genes_groups'.
        use_raw: Whether to use `adata.raw` for testing.

    Returns:
        A pandas DataFrame with the full, structured results.
    """
    if key_added is None:
        key_added = "rank_genes_groups"

    print(f"Finding markers for groups in '{groupby}' using '{method}' method...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        layer=layer,
        key_added=key_added,
        use_raw=use_raw,
        pts=True,  # Always compute pts for filtering
    )

    # Structure the results into a single, clean DataFrame
    result_dfs = []
    for group in adata.obs[groupby].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)
        df["group"] = group
        result_dfs.append(df)

    full_df = pd.concat(result_dfs, ignore_index=True)
    adata.uns[f"{key_added}_df"] = full_df  # Store for reference

    print(
        f"Marker gene analysis complete. Results stored in `adata.uns['{key_added}_df']`."
    )
    return full_df


def filter_markers(
    adata: sc.AnnData,
    key: str = "rank_genes_groups",
    min_log2fc: float = 1.0,
    max_padj: float = 0.05,
    min_in_group_pct: float = 0.25,
) -> pd.DataFrame:
    """
    Filter the results of `find_markers` based on standard criteria.

    Args:
        adata: AnnData object after running `find_markers`.
        key: Key in `adata.uns` where the results DataFrame is stored.
        min_log2fc: Minimum log2 fold change.
        max_padj: Maximum adjusted p-value.
        min_in_group_pct: Minimum percentage of cells in the group expressing the gene.

    Returns:
        A filtered pandas DataFrame of high-quality marker genes.
    """
    df_key = f"{key}_df"
    if df_key not in adata.uns:
        raise KeyError(
            f"Results DataFrame not found at `adata.uns['{df_key}']`. Run `find_markers` first."
        )

    df = adata.uns[df_key].copy()

    # Apply filters
    pass_logfc = df["logfoldchanges"] >= min_log2fc
    pass_padj = df["pvals_adj"] <= max_padj
    pass_pct = df["pct_nz_group"] / 100 >= min_in_group_pct

    filtered_df = df[pass_logfc & pass_padj & pass_pct].copy()

    print(
        f"Filtered {len(df)} total markers down to {len(filtered_df)} high-quality markers."
    )
    adata.uns[f"{key}_filtered_df"] = filtered_df
    return filtered_df


def get_conserved_markers(
    adata: sc.AnnData,
    groupby: str,
    condition_key: str,
    **kwargs,
) -> Dict[str, pd.DataFrame]:
    """
    Find marker genes for a group that are conserved across multiple conditions.
    """
    conserved_markers = {}
    for group in adata.obs[groupby].cat.categories:
        print(f"Finding conserved markers for group: {group}")

        markers_per_condition = []
        for condition in adata.obs[condition_key].cat.categories:
            subset = adata[
                (adata.obs[groupby] == group) & (adata.obs[condition_key] == condition)
            ]
            if subset.n_obs < 10:
                continue

            # Run DE within the condition, comparing the group to all other cells in that condition
            temp_adata = adata[adata.obs[condition_key] == condition].copy()

            # Check if group exists in this condition
            if group not in temp_adata.obs[groupby].cat.categories:
                continue

            sc.tl.rank_genes_groups(
                temp_adata, groupby=groupby, groups=[group], reference="rest", **kwargs
            )
            df = sc.get.rank_genes_groups_df(temp_adata, group=group)
            df["condition"] = condition
            markers_per_condition.append(df)

        if not markers_per_condition:
            print(f"  No markers found for group '{group}' in any condition.")
            continue

        # Combine and find conserved genes
        full_df = pd.concat(markers_per_condition)

        # A gene is conserved if it's significant in most conditions
        # Here we define "conserved" as being in the top markers across multiple conditions
        # A more complex statistical combination (e.g., MetaDE) could be used here
        conserved_df = (
            full_df[full_df["pvals_adj"] < 0.05]
            .groupby("names")
            .filter(
                lambda x: len(x) >= adata.obs[condition_key].nunique() - 1
            )  # Present in n-1 conditions
        )

        # Aggregate stats for conserved genes
        agg_df = (
            conserved_df.groupby("names")
            .agg(
                mean_log2fc=("logfoldchanges", "mean"),
                min_pval_adj=("pvals_adj", "min"),
                n_conditions=("condition", "nunique"),
            )
            .sort_values("mean_log2fc", ascending=False)
        )

        conserved_markers[group] = agg_df
        print(f"  Found {len(agg_df)} conserved markers for group '{group}'.")

    return conserved_markers
