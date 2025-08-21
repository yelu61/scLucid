"""
Differential expression analysis functions for single-cell RNA-seq data.

This module provides comprehensive tools for identifying marker genes and
performing enrichment analysis for cell clusters. It supports various statistical
methods for differential testing, filtering strategies for high-quality markers,
and functional annotation of gene signatures.

Key functionalities:
- Scenario A: Find marker genes for each cluster (one-vs-rest). -find_markers
- Scenario B: Compare two specific cell types (e.g., CD4+ T vs. CD8+ T). -compare_groups
- Scenario C: Find markers for a cell type that are conserved across conditions -get_conserved_markers 
  (e.g., find robust T cell markers present in both Control and Treated samples).
- Scenario D: Compare the same cell type across different conditions -compare_conditions
  (e.g., Treated T-cells vs. Control T-cells).
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
try:
    import cosg
except ImportError:
    cosg = None

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "find_markers",
    "compare_groups",
    "get_conserved_markers",
    "compare_conditions",
    "filter_markers",
    "visualize_markers",
    "characterize_clusters"
]


# --- Helper Functions ---
def _create_pseudobulk(
    adata: AnnData,
    sample_key: str,
    cluster_key: str,
    layer: Optional[str] = None,
    min_cells: int = 10,
) -> AnnData:
    """
    Create a pseudobulk AnnData object from single-cell data.

    Args:
        adata: AnnData object.
        sample_key: Key in adata.obs for sample information.
        cluster_key: Key in adata.obs for cluster/cell type information.
        layer: Layer to use for counts (e.g., 'counts'). If None, uses adata.X.
        min_cells: Minimum number of cells to form a pseudobulk sample.

    Returns:
        A new AnnData object where each observation is a pseudobulk sample.
    """
    if layer:
        count_matrix = adata.layers[layer]
    else:
        count_matrix = adata.X

    groups = adata.obs.groupby([sample_key, cluster_key])
    
    pseudobulk_counts = []
    pseudobulk_obs = []

    for (sample, cluster), idx in groups.indices.items():
        if len(idx) >= min_cells:
            # Sum counts for all cells in the group
            pb_counts = count_matrix[idx, :].sum(axis=0)
            
            # Ensure it's a 1D array
            if not isinstance(pb_counts, np.ndarray):
                pb_counts = pb_counts.A.flatten()
            
            pseudobulk_counts.append(pb_counts)
            pseudobulk_obs.append({"sample": sample, "cluster": cluster, "n_cells": len(idx)})

    if not pseudobulk_counts:
        log.warning("No pseudobulk samples were created (check min_cells threshold).")
        return AnnData()
        
    # Create the new AnnData object
    pb_adata = AnnData(
        X=np.vstack(pseudobulk_counts),
        obs=pd.DataFrame(pseudobulk_obs),
        var=adata.var.copy()
    )
    pb_adata.obs.set_index(pb_adata.obs['sample'] + "_" + pb_adata.obs['cluster'], inplace=True)
    
    log.info(f"Created pseudobulk AnnData with {pb_adata.n_obs} samples and {pb_adata.n_vars} genes.")
    return pb_adata

# --- Main Functions ---
def find_markers(
    adata: AnnData,
    groupby: str,
    method: Literal[
        "wilcoxon", "t-test", "logreg", "cosg",
    ] = "wilcoxon",
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
    Find marker genes for all groups using `sc.tl.rank_genes_groups`.

    This function identifies genes that are differentially expressed between
    groups of cells, using various statistical methods.

    Args:
        adata: AnnData object containing gene expression and metadata
        groupby: Key in `adata.obs` for grouping cells
        method: Statistical method for differential testing:
            - 'wilcoxon': Wilcoxon rank-sum test (non-parametric, default)
            - 't-test': t-test with equal variance assumption
            - 'logreg': Logistic regression with L1 penalty
            - 'cosg': Conserved Overlapping Signature Gene (requires `cosg` package)
        layer: Layer in adata.layers to use for expression values
        key_added: Key in `adata.uns` to store results (default: 'rank_genes_groups')
        use_raw: Whether to use `adata.raw` for testing
        min_cells: Minimum number of cells required per group
        groups: Specific groups to test (default: all groups in groupby)
        reference: Group to compare against. If 'rest', compare to average of all other groups
        fold_change_max: Ceiling for fold change (to avoid extreme values)
        pval_cutoff: If provided, filter results by this p-value threshold
        copy: Return a copy of the AnnData object
        **kwargs: Additional parameters passed to `sc.tl.rank_genes_groups`

    Returns:
        A pandas DataFrame with the full, structured results

    Examples:
        >>> # Find markers for all clusters using Wilcoxon test
        >>> marker_df = find_markers(adata, groupby='leiden')
        >>>
        >>> # Find markers for specific clusters using t-test and a specific layer
        >>> marker_df = find_markers(
        ...     adata,
        ...     groupby='cell_types',
        ...     method='t-test',
        ...     groups=['T cells', 'B cells'],
        ...     layer='log1p_norm'
        ... )
    """
    log.info(f"Finding marker genes for groups in '{groupby}' using '{method}' method")

    if copy:
        adata = adata.copy()

    # Set default key if not provided
    if key_added is None:
        key_added = "rank_genes_groups"

    # Check that the groupby column exists
    if groupby not in adata.obs.columns:
        log.error(f"Column '{groupby}' not found in adata.obs")
        raise ValueError(f"Column '{groupby}' not found in adata.obs")

    # Verify that groups exist if specified
    if groups is not None:
        missing_groups = [
            g for g in groups if g not in adata.obs[groupby].cat.categories
        ]
        if missing_groups:
            log.error(f"Groups not found in {groupby}: {', '.join(missing_groups)}")
            raise ValueError(
                f"Groups not found in {groupby}: {', '.join(missing_groups)}"
            )

    # Check for sufficient cells in each group
    if min_cells > 0:
        group_counts = adata.obs[groupby].value_counts()
        small_groups = group_counts[group_counts < min_cells].index.tolist()

        if groups is None and small_groups:
            log.warning(
                f"Groups with fewer than {min_cells} cells will be skipped: {', '.join(small_groups)}"
            )
        elif groups is not None:
            small_requested = [g for g in groups if g in small_groups]
            if small_requested:
                log.warning(
                    f"Requested groups with fewer than {min_cells} cells: {', '.join(small_requested)}"
                )

    # Ensure the structured .uns path exists for storing results
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('de', {})
    
    if method in ["wilcoxon", "t-test", "logreg"]:
        # --- Standard Scanpy methods ---
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method=method,
            layer=layer,
            key_added=key_added,
            use_raw=use_raw,
            pts=True,  # Always compute percentage of cells expressing the gene
            groups=groups,
            reference=reference,
            **kwargs,
        )

        log.info("Differential expression analysis complete")

        # Structure the results into a single, clean DataFrame
        result_dfs = []

        # Get categories - either specified groups or all categories
        if groups is not None:
            categories = groups
        else:
            categories = adata.obs[groupby].cat.categories

        # Process each group
        for group in categories:
            try:
                df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)

                # Skip if no results for this group
                if df.empty:
                    log.warning(f"No results found for group '{group}'")
                    continue

                df["group"] = group

                # Apply fold change ceiling if requested
                if fold_change_max is not None and fold_change_max > 0:
                    df["logfoldchanges"] = df["logfoldchanges"].clip(
                        upper=fold_change_max
                    )

                # Apply p-value cutoff if requested
                if pval_cutoff is not None and pval_cutoff > 0:
                    df = df[df["pvals_adj"] <= pval_cutoff].copy()

                result_dfs.append(df)

            except Exception as e:
                log.error(f"Error processing group '{group}': {str(e)}")

        if not result_dfs:
            log.warning("No valid results found for any group")
            return pd.DataFrame()

        # Combine results from all groups
        full_df = pd.concat(result_dfs, ignore_index=True)

        # Store the results in the AnnData object
        adata.uns['scrnatk']['analysis']['de'][key_added] = adata.uns[key_added]
        adata.unss['scrnatk']['analysis']['de'][f"{key_added}_df"] = full_df

        log.info(
            f"Found {len(full_df)} differentially expressed genes across {len(result_dfs)} groups"
        )
        log.info(f"Results stored in `adata.uns['{key_added}_df']`")

        return full_df

    elif method == "cosg":
        # --- COSG method ---
        if cosg is None:
            raise ImportError("Please install cosg: `pip install cosg`")
        
        log.info("Running COSG for fast marker gene detection...")
        cosg.cosg(
            adata,
            key_added=key_added,
            groupby=groupby,
            use_raw=kwargs.get('use_raw', True),
            n_genes_user=kwargs.get('n_genes', 100)
        )
        
        # --- IMPORTANT: Reformat COSG output to match Scanpy's structure ---
        log.info("Reformatting COSG output to be compatible with Scanpy.")
        cosg_results = adata.uns[key_added]
        
        # Create an empty structured array matching Scanpy's format
        n_genes = cosg_results['names'].shape[0]
        n_groups = cosg_results['names'].shape[1]
        
        scanpy_recarray = np.rec.fromarrays(
            [
                np.empty((n_genes, n_groups), dtype='<U25'), # names
                np.full((n_genes, n_groups), np.nan, dtype='<f4'), # scores
                np.full((n_genes, n_groups), np.nan, dtype='<f4'), # logfoldchanges
                np.full((n_genes, n_groups), 1.0, dtype='<f4'), # pvals
                np.full((n_genes, n_groups), 1.0, dtype='<f4'), # pvals_adj
            ],
            dtype=[('names', '<U25'), ('scores', '<f4'), ('logfoldchanges', '<f4'), ('pvals', '<f4'), ('pvals_adj', '<f4')]
        )
        
        # Populate the structured array with COSG results
        scanpy_recarray['names'] = cosg_results['names']
        scanpy_recarray['scores'] = cosg_results['scores']
        # COSG provides log-normalized expression, not log-fold-change, but we can store it here.
        # For a true LFC, a separate calculation would be needed. This is a reasonable approximation.
        scanpy_recarray['logfoldchanges'] = cosg_results['log_expr']
        
        # Store the reformatted results in our namespace
        adata.uns['scrnatk']['analysis']['de'][key_added] = {
            'params': {'groupby': groupby, 'method': 'cosg'},
            'names': scanpy_recarray['names'].T,
            'scores': scanpy_recarray['scores'].T,
            'logfoldchanges': scanpy_recarray['logfoldchanges'].T,
            'pvals': scanpy_recarray['pvals'].T,
            'pvals_adj': scanpy_recarray['pvals_adj'].T
        }
        # Delete the original COSG output
        del adata.uns[key_added]

    else:
        raise ValueError(f"Unsupported DE method: {method}")

    # --- Create and store the convenient DataFrame version ---
    df_key = f"{key_added}_df"
    # We now read from our namespaced location
    result_df = sc.get.rank_genes_groups_df(adata, key=key_added, group=None, uns_key=['scrnatk', 'analysis', 'de'])
    adata.uns['scrnatk']['analysis']['de'][df_key] = result_df
    
    log.info(f"DE analysis complete. Results stored under .uns['scrnatk']['analysis']['de']['{key_added}']")
    
    return adata
    

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
    Filter the results of `find_markers` based on standard criteria.

    This function applies multiple filters to identify high-quality marker genes
    that are both statistically significant and biologically meaningful.

    Args:
        adata: AnnData object after running `find_markers`
        key: Key in `adata.uns` where the results DataFrame is stored
        min_log2fc: Minimum log2 fold change (effect size)
        max_padj: Maximum adjusted p-value (statistical significance)
        min_in_group_pct: Minimum percentage of cells in the group expressing the gene
        max_out_group_pct: Maximum percentage of cells outside the group expressing the gene
        min_diff_pct: Minimum percentage difference between in-group and out-group expression
        keep_top_n: If provided, keep only the top N genes per group after filtering
        key_added: Key to store filtered results (default: f"{key}_filtered_df")

    Returns:
        A filtered pandas DataFrame of high-quality marker genes

    Examples:
        >>> # Basic filtering with default parameters
        >>> filtered_df = filter_markers(adata)
        >>>
        >>> # More stringent filtering with custom parameters
        >>> filtered_df = filter_markers(
        ...     adata,
        ...     min_log2fc=2.0,
        ...     max_padj=0.01,
        ...     min_in_group_pct=0.5,
        ...     max_out_group_pct=0.1,
        ...     keep_top_n=50
        ... )
    """
    # Determine key for filtered results
    if key_added is None:
        key_added = f"{key}_filtered_df"

    # Input validation
    df_key = f"{key}_df"
    if df_key not in adata.uns:
        log.error(f"Results DataFrame not found at `adata.uns['{df_key}']`")
        raise KeyError(
            f"Results DataFrame not found at `adata.uns['{df_key}']`. Run `find_markers` first."
        )

    log.info("Filtering marker genes with the following criteria:")
    log.info(f"  - Minimum log2 fold change: {min_log2fc}")
    log.info(f"  - Maximum adjusted p-value: {max_padj}")
    log.info(f"  - Minimum in-group expression: {min_in_group_pct * 100:.1f}%")

    if max_out_group_pct is not None:
        log.info(f"  - Maximum out-group expression: {max_out_group_pct * 100:.1f}%")

    if min_diff_pct is not None:
        log.info(f"  - Minimum expression difference: {min_diff_pct * 100:.1f}%")

    if keep_top_n is not None:
        log.info(f"  - Keeping top {keep_top_n} genes per group after filtering")

    # Get the original results DataFrame
    df = adata.uns[df_key].copy()

    if df.empty:
        log.warning("Input DataFrame is empty, returning empty result")
        adata.uns[key_added] = pd.DataFrame()
        return pd.DataFrame()

    # Check required columns
    required_cols = ["logfoldchanges", "pvals_adj", "pct_nz_group"]
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        log.error(f"Required columns missing from DataFrame: {', '.join(missing_cols)}")
        raise ValueError(f"Required columns missing: {', '.join(missing_cols)}")

    # Apply filters
    original_count = len(df)

    # Basic filters (always applied)
    pass_logfc = df["logfoldchanges"] >= min_log2fc
    pass_padj = df["pvals_adj"] <= max_padj
    pass_pct = df["pct_nz_group"] / 100 >= min_in_group_pct

    # Additional filters (if requested)
    if "pct_nz_reference" in df.columns and max_out_group_pct is not None:
        pass_out_pct = df["pct_nz_reference"] / 100 <= max_out_group_pct
    else:
        pass_out_pct = pd.Series(True, index=df.index)

    if "pct_nz_reference" in df.columns and min_diff_pct is not None:
        pct_diff = (df["pct_nz_group"] - df["pct_nz_reference"]) / 100
        pass_diff_pct = pct_diff >= min_diff_pct
    else:
        pass_diff_pct = pd.Series(True, index=df.index)

    # Combine all filters
    filtered_df = df[
        pass_logfc & pass_padj & pass_pct & pass_out_pct & pass_diff_pct
    ].copy()

    # Keep top N genes per group if requested
    if keep_top_n is not None and keep_top_n > 0:
        top_per_group = []
        for group in filtered_df["group"].unique():
            group_df = filtered_df[filtered_df["group"] == group]
            top_group = group_df.sort_values("logfoldchanges", ascending=False).head(
                keep_top_n
            )
            top_per_group.append(top_group)

        filtered_df = pd.concat(top_per_group, ignore_index=True)

    # Log filtering results
    if original_count > 0:
        pct_remaining = len(filtered_df) / original_count * 100
        log.info(
            f"Filtered {original_count} total markers down to {len(filtered_df)} "
            f"high-quality markers ({pct_remaining:.1f}% retained)"
        )
    else:
        log.warning("No markers to filter")

    # Store filtered results
    adata.uns[key_added] = filtered_df

    # Return filtered DataFrame
    return filtered_df


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
    Find marker genes for a group that are conserved across multiple conditions.

    This function identifies genes that are consistently differentially expressed
    in a particular cell type/group across different experimental conditions,
    batches, or datasets.

    Args:
        adata: AnnData object containing gene expression and metadata
        groupby: Key in adata.obs for cell type/cluster annotations
        condition_key: Key in adata.obs containing condition/batch information
        method: Statistical method for differential testing
        min_cells: Minimum number of cells required in a group within a condition
        min_conditions: Minimum number of conditions where a gene must be significant
        min_log2fc: Minimum log2 fold change required in each condition
        max_padj: Maximum adjusted p-value allowed in each condition
        min_in_group_pct: Minimum percentage of cells expressing the gene in the group
        layer: Layer to use for expression values
        use_raw: Whether to use adata.raw for testing
        key_added: Key in adata.uns to store results

    Returns:
        Dictionary mapping group names to DataFrames of conserved markers

    Examples:
        >>> # Find T cell markers conserved across donors
        >>> conserved_markers = get_conserved_markers(
        ...     adata,
        ...     groupby='cell_type',
        ...     condition_key='donor'
        ... )
        >>>
        >>> # Find cluster markers conserved across experimental conditions
        >>> conserved_markers = get_conserved_markers(
        ...     adata,
        ...     groupby='leiden',
        ...     condition_key='condition',
        ...     min_log2fc=1.0,
        ...     max_padj=0.01
        ... )
    """
    log.info("Finding conserved markers across conditions")
    log.info(f"Group variable: '{groupby}', Condition variable: '{condition_key}'")

    # Validate input columns
    for col in [groupby, condition_key]:
        if col not in adata.obs.columns:
            log.error(f"Column '{col}' not found in adata.obs")
            raise ValueError(f"Column '{col}' not found in adata.obs")

    # Set key for storing results
    if key_added is None:
        key_added = f"conserved_markers_{groupby}_{condition_key}"

    # Get unique conditions and groups
    conditions = adata.obs[condition_key].cat.categories
    groups = adata.obs[groupby].cat.categories

    log.info(f"Analyzing {len(groups)} groups across {len(conditions)} conditions")

    # Set minimum conditions if not specified
    if min_conditions is None:
        # Default to all conditions minus 1 (allow one condition to be missing)
        min_conditions = max(1, len(conditions) - 1)
        log.info(f"Setting min_conditions to {min_conditions} (total conditions - 1)")

    # Initialize results dictionary
    conserved_markers = {}

    # Analyze each group
    for group in groups:
        log.info(f"Finding conserved markers for group: {group}")

        markers_per_condition = []
        conditions_analyzed = []

        # Analyze each condition separately
        for condition in conditions:
            # Extract cells for this group and condition
            subset = adata[
                (adata.obs[groupby] == group) & (adata.obs[condition_key] == condition)
            ]

            # Skip if too few cells
            if subset.n_obs < min_cells:
                log.debug(
                    f"  Skipping condition '{condition}': only {subset.n_obs} cells (min: {min_cells})"
                )
                continue

            # Create temporary dataset with only cells from this condition
            try:
                temp_adata = adata[adata.obs[condition_key] == condition].copy()

                # Check if group exists in this condition
                if group not in temp_adata.obs[groupby].cat.categories:
                    log.debug(f"  Group '{group}' not found in condition '{condition}'")
                    continue

                # Run differential expression within this condition
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

                # Extract results for this group
                df = sc.get.rank_genes_groups_df(temp_adata, group=group)

                # Apply basic filtering
                df = df[
                    (df["logfoldchanges"] >= min_log2fc)
                    & (df["pvals_adj"] <= max_padj)
                    & (df["pct_nz_group"] / 100 >= min_in_group_pct)
                ]

                # Skip if no genes pass filters
                if df.empty:
                    log.debug(
                        f"  No significant markers found for '{group}' in '{condition}'"
                    )
                    continue

                # Add condition information
                df["condition"] = condition
                markers_per_condition.append(df)
                conditions_analyzed.append(condition)

                log.debug(f"  Found {len(df)} markers for '{group}' in '{condition}'")

            except Exception as e:
                log.error(f"  Error analyzing condition '{condition}': {str(e)}")

        # Check if we have enough conditions with results
        if len(markers_per_condition) < min_conditions:
            log.warning(
                f"  Insufficient conditions for group '{group}': "
                f"found {len(markers_per_condition)}, need {min_conditions}"
            )
            continue

        # Combine markers from all conditions
        if not markers_per_condition:
            log.warning(f"  No markers found for group '{group}' in any condition")
            continue

        log.info(
            f"  Found markers for '{group}' in {len(conditions_analyzed)}/{len(conditions)} conditions"
        )

        # Combine and find conserved genes
        full_df = pd.concat(markers_per_condition)

        # A gene is conserved if it's significant in at least min_conditions
        gene_counts = full_df.groupby("names").size()
        conserved_genes = gene_counts[gene_counts >= min_conditions].index.tolist()

        log.info(
            f"  {len(conserved_genes)} genes are conserved across at least {min_conditions} conditions"
        )

        if not conserved_genes:
            log.warning(f"  No conserved markers found for group '{group}'")
            continue

        # Filter for only conserved genes
        conserved_df = full_df[full_df["names"].isin(conserved_genes)]

        # Aggregate statistics for conserved genes
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

        # Store results
        conserved_markers[group] = agg_df

        log.info(
            f"  Found {len(agg_df)} conserved markers for group '{group}' "
            f"(mean log2FC: {agg_df['mean_log2fc'].mean():.2f})"
        )

    # Store in AnnData object
    adata.uns[key_added] = conserved_markers

    log.info(
        f"Analysis complete. Found conserved markers for {len(conserved_markers)}/{len(groups)} groups"
    )
    log.info(f"Results stored in adata.uns['{key_added}']")

    return conserved_markers


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
    Directly compare two specific groups to find differentially expressed genes.

    This function performs a head-to-head comparison between two groups, which is
    often more sensitive than the standard one-vs-rest approach for finding markers.

    Args:
        adata: AnnData object containing expression data and metadata
        groupby: Key in adata.obs for group assignments
        group1: Name of the first group to compare
        group2: Name of the second group to compare
        layer: Layer to use for expression values
        use_raw: Whether to use adata.raw for testing
        n_genes: Number of top genes to return for each direction
        min_log2fc: Minimum log2 fold change required
        max_padj: Maximum adjusted p-value allowed
        min_in_group_pct: Minimum percentage of cells expressing the gene
        plot: Whether to create a volcano plot visualization
        save_path: Path to save the volcano plot
        key_added: Key in adata.uns to store results

    Returns:
        DataFrame containing genes differentially expressed between the two groups

    Examples:
        >>> # Compare CD4+ vs CD8+ T cells
        >>> de_df = compare_groups(
        ...     adata,
        ...     groupby='cell_type',
        ...     group1='CD4+ T',
        ...     group2='CD8+ T',
        ...     plot=True
        ... )
        >>>
        >>> # Compare treatment vs control in a specific cell type
        >>> de_df = compare_groups(
        ...     adata[adata.obs['cell_type'] == 'Macrophages'],
        ...     groupby='condition',
        ...     group1='treatment',
        ...     group2='control',
        ...     min_log2fc=1.0
        ... )
    """
    log.info(f"Comparing groups: '{group1}' vs '{group2}'")

    # Set default key_added if not provided
    if key_added is None:
        key_added = f"compare_{group1}_vs_{group2}"

    # Check if groups exist
    if groupby not in adata.obs.columns:
        log.error(f"Column '{groupby}' not found in adata.obs")
        raise ValueError(f"Column '{groupby}' not found in adata.obs")

    # Check if groups exist in the data
    unique_groups = adata.obs[groupby].unique()

    if group1 not in unique_groups:
        log.error(f"Group '{group1}' not found in {groupby}")
        raise ValueError(f"Group '{group1}' not found in {groupby}")

    if group2 not in unique_groups:
        log.error(f"Group '{group2}' not found in {groupby}")
        raise ValueError(f"Group '{group2}' not found in {groupby}")

    # Create a temporary copy with only the two groups
    temp = adata[adata.obs[groupby].isin([group1, group2])].copy()

    # Check if we have enough cells
    group1_cells = np.sum(temp.obs[groupby] == group1)
    group2_cells = np.sum(temp.obs[groupby] == group2)

    log.info(f"Group '{group1}': {group1_cells} cells")
    log.info(f"Group '{group2}': {group2_cells} cells")

    if group1_cells < 5 or group2_cells < 5:
        log.warning(
            f"Very few cells in one or both groups: {group1}={group1_cells}, {group2}={group2_cells}"
        )

    # Create a new binary grouping for the comparison
    temp.obs["compare_groups"] = temp.obs[groupby].map(
        {group1: "group1", group2: "group2"}
    )
    temp.obs["compare_groups"] = temp.obs["compare_groups"].astype("category")

    # Run differential expression analysis
    try:
        sc.tl.rank_genes_groups(
            temp,
            groupby="compare_groups",
            groups=["group1", "group2"],
            reference="rest",  # Each group compared to the other
            method="wilcoxon",
            layer=layer,
            use_raw=use_raw,
            pts=True,
        )

        # Get results for both directions
        group1_vs_group2 = sc.get.rank_genes_groups_df(temp, group="group1")
        group1_vs_group2["comparison"] = f"{group1}_vs_{group2}"
        group1_vs_group2["higher_in"] = group1

        group2_vs_group1 = sc.get.rank_genes_groups_df(temp, group="group2")
        group2_vs_group1["comparison"] = f"{group2}_vs_{group1}"
        group2_vs_group1["higher_in"] = group2

        # Combine results
        all_results = pd.concat([group1_vs_group2, group2_vs_group1], ignore_index=True)

        # Apply filtering
        filtered_results = all_results[
            (all_results["logfoldchanges"].abs() >= min_log2fc)
            & (all_results["pvals_adj"] <= max_padj)
            & (all_results["pct_nz_group"] / 100 >= min_in_group_pct)
        ].copy()

        # Sort by absolute fold change
        filtered_results["abs_log2fc"] = filtered_results["logfoldchanges"].abs()
        filtered_results = filtered_results.sort_values("abs_log2fc", ascending=False)

        # Limit to top N genes in each direction if requested
        if n_genes > 0:
            top_results = []
            for group in [group1, group2]:
                group_results = filtered_results[filtered_results["higher_in"] == group]
                top_results.append(group_results.head(n_genes))

            filtered_results = pd.concat(top_results, ignore_index=True)

        # Store results
        adata.uns[key_added] = filtered_results

        # Create a volcano plot if requested
        if plot:
            try:
                # Get all results for plotting
                all_results_g1 = sc.get.rank_genes_groups_df(temp, group="group1")
                all_results_g1["higher_in"] = group1

                # Set up the plot
                plt.figure(figsize=(12, 8))

                # Add points for genes that don't pass thresholds
                non_sig = all_results_g1[
                    (all_results_g1["logfoldchanges"].abs() < min_log2fc)
                    | (all_results_g1["pvals_adj"] > max_padj)
                ]

                plt.scatter(
                    non_sig["logfoldchanges"],
                    -np.log10(non_sig["pvals_adj"] + 1e-10),
                    alpha=0.3,
                    s=10,
                    color="grey",
                    label="Not significant",
                )

                # Add points for genes higher in group1
                sig_g1 = filtered_results[filtered_results["higher_in"] == group1]
                if not sig_g1.empty:
                    plt.scatter(
                        sig_g1["logfoldchanges"],
                        -np.log10(sig_g1["pvals_adj"] + 1e-10),
                        alpha=0.7,
                        s=30,
                        color="red",
                        label=f"Higher in {group1}",
                    )

                # Add points for genes higher in group2
                sig_g2 = filtered_results[filtered_results["higher_in"] == group2]
                if not sig_g2.empty:
                    plt.scatter(
                        sig_g2["logfoldchanges"],
                        -np.log10(sig_g2["pvals_adj"] + 1e-10),
                        alpha=0.7,
                        s=30,
                        color="blue",
                        label=f"Higher in {group2}",
                    )

                # Add labels for top genes
                top_genes = pd.concat(
                    [
                        sig_g1.head(10) if not sig_g1.empty else pd.DataFrame(),
                        sig_g2.head(10) if not sig_g2.empty else pd.DataFrame(),
                    ]
                )

                for _, gene in top_genes.iterrows():
                    plt.text(
                        gene["logfoldchanges"],
                        -np.log10(gene["pvals_adj"] + 1e-10),
                        gene["names"],
                        fontsize=9,
                        ha="center",
                        va="bottom",
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.7),
                    )

                # Add threshold lines
                plt.axvline(x=min_log2fc, color="grey", linestyle="--", alpha=0.5)
                plt.axvline(x=-min_log2fc, color="grey", linestyle="--", alpha=0.5)
                plt.axhline(
                    y=-np.log10(max_padj), color="grey", linestyle="--", alpha=0.5
                )

                # Add labels and title
                plt.xlabel("Log2 Fold Change")
                plt.ylabel("-log10(Adjusted p-value)")
                plt.title(f"Differential Expression: {group1} vs {group2}")
                plt.legend()
                plt.grid(alpha=0.2)

                # Save or show the plot
                if save_path is not None:
                    plt.tight_layout()
                    plt.savefig(save_path, dpi=300)
                    plt.close()
                    log.info(f"Saved volcano plot to {save_path}")
                else:
                    plt.tight_layout()
                    plt.show()

            except Exception as e:
                log.warning(f"Error creating volcano plot: {str(e)}")

        # Log results
        log.info(f"Found {len(filtered_results)} differentially expressed genes")
        log.info(
            f"  {len(sig_g1) if 'sig_g1' in locals() else 0} genes higher in {group1}"
        )
        log.info(
            f"  {len(sig_g2) if 'sig_g2' in locals() else 0} genes higher in {group2}"
        )
        log.info(f"Results stored in adata.uns['{key_added}']")

        return filtered_results

    except Exception as e:
        log.error(f"Error in differential expression analysis: {str(e)}")
        raise RuntimeError(f"Comparison failed: {str(e)}")


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
    Compare two conditions within a specific cell type or group.

    This function simplifies the common task of finding differentially expressed
    genes between conditions (e.g., 'Treated' vs 'Control') for a single
    cell population (e.g., 'T cells').

    Args:
        adata: AnnData object.
        groupby: The key in .obs for the cell type/group to isolate (e.g., 'cell_type').
        group_name: The specific cell type/group to analyze (e.g., 'T cells').
        condition_key: The key in .obs that defines the conditions to compare (e.g., 'treatment_status').
        condition1: The first condition name (e.g., 'Treated').
        condition2: The second condition name (e.g., 'Control').
        key_added: Key in .uns to store the results.
        **kwargs: Additional arguments passed to the underlying `compare_groups` function.

    Returns:
        A DataFrame of differentially expressed genes between the two conditions for the specified cell type.
    """
    log.info(f"Initiating DE analysis for '{group_name}' between conditions '{condition1}' and '{condition2}'")
    
    # 1. Subset the data to the cell type of interest
    if group_name not in adata.obs[groupby].unique():
        raise ValueError(f"Group '{group_name}' not found in adata.obs['{groupby}']")
        
    adata_subset = adata[adata.obs[groupby] == group_name].copy()
    log.info(f"Subsetted data to {adata_subset.n_obs} cells of type '{group_name}'.")

    # 2. Use the existing `compare_groups` function on the subset
    # This is a great example of code reuse and building high-level functions
    # from existing components.
    if key_added is None:
        key_added = f"compare_{condition1}_vs_{condition2}_in_{group_name.replace(' ', '_')}"

    # The compare_groups function returns the final DataFrame and also modifies adata.uns
    # We will run it on the subset and then add the results back to the original adata object.
    de_df = compare_groups(
        adata_subset,
        groupby=condition_key,
        group1=condition1,
        group2=condition2,
        key_added=key_added,
        **kwargs
    )
    
    # 3. Store the results back in the original AnnData object for provenance
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('de', {})
    adata.uns['scrnatk']['analysis']['de'][key_added] = adata_subset.uns['scrnatk']['analysis']['de'][key_added]

    log.info(f"Comparison complete. Results stored in .uns['scrnatk']['analysis']['de']['{key_added}']")
    
    return de_df

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

    This function creates visualizations for marker genes, showing their expression
    patterns across different cell groups. It supports various input formats for
    the markers parameter.

    Args:
        adata: AnnData object containing expression data
        markers: Marker genes to visualize, in one of these formats:
            - DataFrame from find_markers() or filter_markers()
            - Dictionary mapping groups to lists of marker genes
            - List of marker gene names
        groupby: Key in adata.obs for grouping cells (required if markers is a list)
        n_genes_per_group: Number of marker genes to show per group
        plot_type: Type of visualization to create
        dendrogram: Whether to include a dendrogram in the plot
        standard_scale: How to scale the data ('var', 'group', or None)
        swap_axes: Whether to swap the x and y axes
        layer: Layer to use for expression values
        use_raw: Whether to use adata.raw for visualization
        save_path: Path to save the visualization
        figsize: Figure size (width, height) in inches
        **kwargs: Additional parameters passed to the plotting function

    Examples:
        >>> # Visualize markers from find_markers() result
        >>> marker_df = find_markers(adata, groupby='leiden')
        >>> visualize_markers(adata, marker_df, n_genes_per_group=10, plot_type='dotplot')
        >>>
        >>> # Visualize specific genes across cell types
        >>> genes = ['CD3D', 'CD4', 'CD8A', 'MS4A1', 'CD19', 'FCGR3A', 'CD14']
        >>> visualize_markers(adata, genes, groupby='cell_type', plot_type='stacked_violin')
        >>>
        >>> # Visualize markers from a dictionary
        >>> marker_dict = {
        ...     'T cells': ['CD3D', 'CD3E', 'CD2'],
        ...     'B cells': ['MS4A1', 'CD79A', 'CD19'],
        ...     'NK cells': ['NCAM1', 'NKG7', 'KLRF1']
        ... }
        >>> visualize_markers(adata, marker_dict, plot_type='heatmap')
    """
    log.info(f"Creating {plot_type} visualization for marker genes")

    # Process input markers based on type
    gene_list = []

    if isinstance(markers, pd.DataFrame):
        # Case 1: DataFrame from find_markers or filter_markers
        if "group" not in markers.columns or "names" not in markers.columns:
            log.error("DataFrame must contain 'group' and 'names' columns")
            raise ValueError("DataFrame must contain 'group' and 'names' columns")

        # Extract top genes per group
        if n_genes_per_group > 0:
            for group in markers["group"].unique():
                group_markers = markers[markers["group"] == group]
                if not group_markers.empty:
                    # Sort by logfoldchanges if available, otherwise by names
                    if "logfoldchanges" in group_markers.columns:
                        group_markers = group_markers.sort_values(
                            "logfoldchanges", ascending=False
                        )

                    # Get top n genes
                    top_genes = group_markers["names"].head(n_genes_per_group).tolist()
                    gene_list.extend(top_genes)
        else:
            # Use all genes
            gene_list = markers["names"].tolist()

    elif isinstance(markers, dict):
        # Case 2: Dictionary mapping groups to gene lists
        if n_genes_per_group > 0:
            for group, genes in markers.items():
                gene_list.extend(genes[:n_genes_per_group])
        else:
            # Use all genes
            for genes in markers.values():
                gene_list.extend(genes)

    elif isinstance(markers, (list, tuple)):
        # Case 3: Simple list of genes
        gene_list = list(markers)

        # Check if groupby is provided
        if groupby is None:
            log.error("groupby must be specified when markers is a list")
            raise ValueError("groupby must be specified when markers is a list")
    else:
        log.error("markers must be a DataFrame, dictionary, or list")
        raise TypeError("markers must be a DataFrame, dictionary, or list")

    # Remove duplicates while preserving order
    gene_list = list(dict.fromkeys(gene_list))

    # Check that genes exist in the dataset
    missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
    if missing_genes:
        log.warning(f"{len(missing_genes)} genes not found in dataset")
        if len(missing_genes) < 10:
            log.warning(f"Missing genes: {', '.join(missing_genes)}")
        else:
            log.warning(f"First 10 missing genes: {', '.join(missing_genes[:10])}...")

        # Remove missing genes
        gene_list = [gene for gene in gene_list if gene in adata.var_names]

        if not gene_list:
            log.error("No valid genes found for visualization")
            raise ValueError("No valid genes found for visualization")

    log.info(
        f"Visualizing {len(gene_list)} marker genes across {groupby or 'unknown'} groups"
    )

    # Create the plot
    try:
        # Set default figure size if not provided
        if figsize is None:
            if plot_type in ["dotplot", "matrixplot"]:
                n_groups = len(adata.obs[groupby].cat.categories) if groupby else 1
                n_genes = len(gene_list)

                if swap_axes:
                    figsize = (
                        max(6, min(12, n_groups * 0.5)),
                        max(4, min(10, n_genes * 0.3)),
                    )
                else:
                    figsize = (
                        max(6, min(12, n_genes * 0.5)),
                        max(4, min(10, n_groups * 0.3)),
                    )
            else:
                figsize = (12, 8)

        # Create appropriate plot
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
            log.error(f"Unknown plot type: {plot_type}")
            raise ValueError(f"Unknown plot type: {plot_type}")

        # Save figure if requested
        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            log.info(f"Saved visualization to {save_path}")

        log.info(f"Created {plot_type} visualization for {len(gene_list)} genes")

    except Exception as e:
        log.error(f"Error creating visualization: {str(e)}")
        raise RuntimeError(f"Visualization failed: {str(e)}")


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

    This function automates the process of finding top marker genes and enriched
    pathways for each cluster, storing the results in a structured format within
    adata.uns for easy inspection.

    Args:
        adata: AnnData object.
        groupby: Key in adata.obs for grouping cells (e.g., 'leiden').
        de_method: Method for differential expression.
        n_de_genes: Number of top DE genes to store and use for enrichment.
        enrich_organism: Organism for GSEApy ('Human' or 'Mouse').
        enrich_gene_sets: Gene sets to use for enrichment.
        key_added: Key in adata.uns to store the results.

    Returns:
        AnnData object with characterization results in adata.uns[key_added].
    """
    log.info(f"Characterizing clusters in '{groupby}'...")
    
    # 1. Run Differential Expression
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=de_method, use_raw=True, n_genes=n_de_genes)
    
    clusters = adata.obs[groupby].cat.categories
    characterization_results = {}

    for cluster in clusters:
        log.info(f"  - Analyzing cluster: {cluster}")
        de_df = sc.get.rank_genes_groups_df(adata, group=cluster)
        top_genes = de_df.head(n_de_genes)['names'].tolist()
        
        # 2. Run Enrichment Analysis
        try:
            enr = gp.enrichr(
                gene_list=top_genes,
                gene_sets=enrich_gene_sets,
                organism=enrich_organism,
                outdir=None,
            )
            enrich_df = enr.results
        except Exception as e:
            log.warning(f"    Enrichment analysis failed for cluster {cluster}: {e}")
            enrich_df = pd.DataFrame()
            
        characterization_results[cluster] = {
            "top_de_genes": de_df.head(n_de_genes),
            "enrichment": enrich_df,
        }
        
    adata.uns[key_added] = characterization_results
    log.info(f"Cluster characterization complete. Results stored in adata.uns['{key_added}']")
    
    return adata


