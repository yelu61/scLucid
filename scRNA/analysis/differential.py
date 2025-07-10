"""
Differential expression analysis functions for single-cell RNA-seq data.

This module provides functions for identifying marker genes and
performing enrichment analysis for cell clusters.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
from typing import Optional, List, Dict, Tuple, Literal, Union
import matplotlib.pyplot as plt

from .manager import Manager


def find_markers(
    adata,
    groupby: str,
    groups: Optional[Union[str, List[str]]] = None,
    reference: Optional[str] = 'rest',
    method: Literal['wilcoxon', 't-test', 'logreg'] = 'wilcoxon',
    layer: Optional[str] = None,
    key_added: Optional[str] = None,
    pts: bool = True,
    min_fold_change: float = 1.5,
    max_pval: float = 0.05,
    max_adj_pval: float = 0.1,
    min_in_group_fraction: float = 0.25,
    max_out_group_fraction: float = 0.5,
    n_genes: int = 100,
    filter_genes: bool = True,
    plot: bool = True,
    use_raw: bool = False
) -> sc.AnnData:
    """
    Find marker genes for groups using various statistical tests.
    
    Identifies genes that are differentially expressed between groups,
    with options for filtering and visualization.
    
    Args:
        adata: AnnData object
        groupby: Key in adata.obs for grouping cells
        groups: Group(s) for which to find marker genes. If None, find for all groups
        reference: Reference group or 'rest' for all other groups
        method: Method for differential testing
        layer: Layer in adata.layers to use for expression values
        key_added: Key under which to add the results (defaults to "rank_genes_groups")
        pts: Whether to compute and plot pct_change
        min_fold_change: Minimum fold change for a gene to be considered
        max_pval: Maximum p-value for a gene to be considered
        max_adj_pval: Maximum adjusted p-value for a gene to be considered
        min_in_group_fraction: Minimum fraction of cells in group expressing the gene
        max_out_group_fraction: Maximum fraction of cells outside group expressing the gene
        n_genes: Number of top genes to return
        filter_genes: Whether to filter the results based on thresholds
        plot: Whether to plot heatmap of top marker genes
        use_raw: Whether to use raw data for differential testing
        
    Returns:
        AnnData with results stored in uns[key_added]
    """
    # Add protection against missing values
    if use_raw and adata.raw is None:
        print("Warning: use_raw=True but adata.raw is None. Using adata.X instead.")
        use_raw = False
    
    # Ensure key_added is defined
    if key_added is None:
        key_added = "rank_genes_groups"
    
    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        method=method,
        layer=layer,
        key_added=key_added,
        pts=pts,
        use_raw=use_raw
    )
    
    # Get results as DataFrame
    results = sc.get.rank_genes_groups_df(adata, group=None, key=key_added)
    
    # Add additional filter columns if pts was computed
    if pts:
        if 'pct_nz_group' in results.columns and 'pct_nz_reference' in results.columns:
            results['in_group_fraction'] = results['pct_nz_group'] / 100
            results['out_group_fraction'] = results['pct_nz_reference'] / 100
        
        # Calculate fold change (if logfoldchanges not already present)
        if 'logfoldchanges' not in results.columns:
            # Approximate fold change from score
            results['logfoldchanges'] = results['scores']
        
        # Convert log fold change to regular fold change
        results['fold_change'] = np.exp(results['logfoldchanges'])
        
        # Add pass/fail filters
        if filter_genes:
            results['pass_fc'] = results['fold_change'] >= min_fold_change
            results['pass_pval'] = results['pvals'] <= max_pval
            
            # Calculate adjusted p-values if not present
            if 'pvals_adj' not in results.columns:
                # Group by group and adjust p-values within each group
                groups = results['group'].unique()
                adj_pvals = []
                
                for group in groups:
                    group_mask = results['group'] == group
                    group_pvals = results.loc[group_mask, 'pvals'].values
                    group_adj_pvals = multipletests(group_pvals, method='fdr_bh')[1]
                    adj_pvals.extend(group_adj_pvals)
                
                results['pvals_adj'] = adj_pvals
            
            results['pass_adj_pval'] = results['pvals_adj'] <= max_adj_pval
            
            if 'in_group_fraction' in results.columns and 'out_group_fraction' in results.columns:
                results['pass_in_fraction'] = results['in_group_fraction'] >= min_in_group_fraction
                results['pass_out_fraction'] = results['out_group_fraction'] <= max_out_group_fraction
                results['pass_all'] = (
                    results['pass_fc'] & 
                    results['pass_pval'] & 
                    results['pass_adj_pval'] & 
                    results['pass_in_fraction'] & 
                    results['pass_out_fraction']
                )
            else:
                results['pass_all'] = results['pass_fc'] & results['pass_pval'] & results['pass_adj_pval']
    
    # Store full results in adata
    adata.uns[f"{key_added}_dataframe"] = results
    
    # Get filtered results if requested
    if filter_genes and 'pass_all' in results.columns:
        filtered_results = results[results['pass_all']].sort_values(['group', 'pvals'])
        adata.uns[f"{key_added}_filtered"] = filtered_results
    
    # Plot heatmap of top marker genes
    if plot:
        if groups is None:
            # Get unique groups
            plot_groups = adata.obs[groupby].cat.categories.tolist()
        elif isinstance(groups, str):
            plot_groups = [groups]
        else:
            plot_groups = list(groups)
        
        # Limit to a reasonable number of groups for visualization
        if len(plot_groups) > 10:
            print(f"Too many groups ({len(plot_groups)}) for heatmap. Showing first 10.")
            plot_groups = plot_groups[:10]
        
        # Get top genes for each group
        plot_genes = []
        for group in plot_groups:
            if filter_genes and 'pass_all' in results.columns:
                # Get filtered top genes
                group_genes = (
                    results[(results['group'] == group) & results['pass_all']]
                    .sort_values('pvals')
                    .head(n_genes)['names']
                    .tolist()
                )
            else:
                # Get top genes by p-value
                group_genes = (
                    results[results['group'] == group]
                    .sort_values('pvals')
                    .head(n_genes)['names']
                    .tolist()
                )
            
            # Add top N genes
            plot_genes.extend(group_genes[:min(10, len(group_genes))])
        
        # Deduplicate genes
        plot_genes = list(dict.fromkeys(plot_genes))
        
        # Plot heatmap
        if plot_genes:
            sc.pl.heatmap(
                adata, plot_genes, groupby=groupby, 
                dendrogram=True, standard_scale='var',
                swap_axes=True, show_gene_labels=True,
                use_raw=use_raw, layer=layer,
                figsize=(12, min(14, len(plot_genes)/2))
            )
        else:
            print("No marker genes found that pass filters.")
    
    return adata


def marker_enrichment_analysis(
    adata,
    de_genes: List[str],
    marker_config: str,
    method: Literal['hypergeometric', 'fisher', 'binomial'] = 'hypergeometric',
    universe: Optional[List[str]] = None,
    min_genes: int = 3,
    pval_cutoff: float = 0.05,
    plot: bool = True,
    n_top: int = 20
) -> pd.DataFrame:
    """
    Perform enrichment analysis of DE genes against known markers.
    
    Args:
        adata: AnnData object
        de_genes: List of differentially expressed genes
        marker_config: Path to marker configuration file
        method: Statistical method for enrichment analysis
        universe: Background gene set (defaults to all genes in adata)
        min_genes: Minimum number of marker genes required for analysis
        pval_cutoff: P-value cutoff for significance
        plot: Whether to plot enrichment results
        n_top: Number of top enriched cell types to plot
        
    Returns:
        DataFrame with enrichment results
    """
    # Load markers
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Set background gene universe
    if universe is None:
        universe = adata.var_names.tolist()
    
    # Get total number of genes in universe
    N = len(universe)
    
    # Get total number of DE genes
    n = len(de_genes)
    
    # Prepare results
    results = []
    
    # Calculate enrichment for each cell type
    for cell_type, cell in mgr.CELLS.items():
        markers = cell.markers
        
        # Skip if not enough markers
        if len(markers) < min_genes:
            continue
        
        # Get number of markers in universe
        K = len(set(markers).intersection(universe))
        
        # Skip if no markers in universe
        if K == 0:
            continue
        
        # Get number of markers in DE genes
        k = len(set(markers).intersection(de_genes))
        
        # Skip if no markers in DE genes
        if k == 0:
            continue
        
        # Calculate enrichment
        if method == 'hypergeometric':
            # P(X >= k) using hypergeometric distribution
            pval = stats.hypergeom.sf(k-1, N, K, n)
        elif method == 'fisher':
            # Fisher's exact test
            contingency = np.array([
                [k, K-k],
                [n-k, N-K-(n-k)]
            ])
            _, pval = stats.fisher_exact(contingency, alternative='greater')
        elif method == 'binomial':
            # Binomial test with probability p = K/N
            pval = stats.binom_test(k, n, K/N, alternative='greater')
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Calculate enrichment score
        # Odds ratio: (k/n) / (K/N) = k*N / (n*K)
        enrichment_score = (k * N) / (n * K) if n * K > 0 else 0
        
        # Calculate adjusted p-value (to be filled later)
        results.append({
            'cell_type': cell_type,
            'n_markers': len(markers),
            'markers_in_universe': K,
            'markers_in_de': k,
            'p_value': pval,
            'enrichment_score': enrichment_score,
            'level': cell.level,
            'color': cell.color
        })
    
    # Create DataFrame
    result_df = pd.DataFrame(results)
    
    # If empty, return empty DataFrame
    if result_df.empty:
        print("No enrichment results found.")
        return result_df
    
    # Calculate adjusted p-values
    result_df['p_value_adj'] = multipletests(result_df['p_value'], method='fdr_bh')[1]
    
    # Add significance flag
    result_df['significant'] = result_df['p_value_adj'] < pval_cutoff
    
    # Sort by p-value
    result_df = result_df.sort_values('p_value')
    
    # Plot results
    if plot and not result_df.empty:
        # Get top enriched cell types
        plot_df = result_df.head(n_top)
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.3)))
        
        # Create bar plot
        bars = ax.barh(
            plot_df['cell_type'],
            plot_df['enrichment_score'],
            color=[c if c else '#1f77b4' for c in plot_df['color']],
            alpha=0.7
        )
        
        # Add p-value stars
        for i, (_, row) in enumerate(plot_df.iterrows()):
            stars = ''
            if row['p_value_adj'] < 0.001:
                stars = '***'
            elif row['p_value_adj'] < 0.01:
                stars = '**'
            elif row['p_value_adj'] < 0.05:
                stars = '*'
            
            if stars:
                ax.text(
                    row['enrichment_score'] + 0.1,
                    i,
                    stars,
                    va='center'
                )
        
        # Customize plot
        ax.set_xlabel('Enrichment Score')
        ax.set_ylabel('Cell Type')
        ax.set_title('Cell Type Marker Enrichment')
        
        # Add count annotations
        for i, (_, row) in enumerate(plot_df.iterrows()):
            ax.text(
                row['enrichment_score'] / 2,
                i,
                f"{row['markers_in_de']}/{row['markers_in_universe']}",
                color='white',
                va='center',
                ha='center',
                fontweight='bold'
            )
        
        # Add legend for significance
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='*', color='w', markerfacecolor='black', markersize=10, label='p < 0.05'),
            Line2D([0], [0], marker='**', color='w', markerfacecolor='black', markersize=10, label='p < 0.01'),
            Line2D([0], [0], marker='***', color='w', markerfacecolor='black', markersize=10, label='p < 0.001')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.show()
    
    return result_df


def get_conserved_markers(
    adata,
    groupby: str,
    groups: Optional[List[str]] = None,
    condition_key: str = "condition",
    conditions: Optional[List[str]] = None,
    method: Literal['wilcoxon', 't-test'] = 'wilcoxon',
    min_fold_change: float = 1.5,
    max_pval: float = 0.05,
    min_significant_conditions: int = 1,
    pts: bool = True,
    n_genes: int = 50,
    layer: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Find conserved marker genes across conditions.
    
    Identifies genes that are consistently differentially expressed
    in a group across multiple conditions.
    
    Args:
        adata: AnnData object
        groupby: Key in adata.obs for grouping cells
        groups: Groups for which to find conserved markers (if None, all groups)
        condition_key: Key in adata.obs for condition information
        conditions: Conditions to compare (if None, all conditions)
        method: Method for differential testing
        min_fold_change: Minimum fold change for a gene to be considered
        max_pval: Maximum p-value for a gene to be considered
        min_significant_conditions: Minimum number of conditions where a gene must be significant
        pts: Whether to compute percentage of cells expressing the gene
        n_genes: Number of top genes to return for each group
        layer: Layer to use for differential testing
        
    Returns:
        Dictionary mapping groups to DataFrames of conserved markers
    """
    # Check if condition_key exists
    if condition_key not in adata.obs.columns:
        raise ValueError(f"Condition key '{condition_key}' not found in adata.obs")
    
    # Get groups if not specified
    if groups is None:
        groups = adata.obs[groupby].cat.categories.tolist()
    elif isinstance(groups, str):
        groups = [groups]
    
    # Get conditions if not specified
    if conditions is None:
        conditions = adata.obs[condition_key].cat.categories.tolist()
    elif isinstance(conditions, str):
        conditions = [conditions]
    
    # Check if we have enough conditions
    if len(conditions) < 2:
        raise ValueError("At least 2 conditions are required for conserved marker analysis")
    
    # Dictionary to store results
    conserved_markers = {}
    
    # For each group, find markers in each condition
    for group in groups:
        print(f"Finding conserved markers for group '{group}'...")
        
        # Dictionary to store markers for this group across conditions
        group_markers = {}
        
        # For each condition, find markers
        for condition in conditions:
            # Subset data for this condition
            condition_mask = adata.obs[condition_key] == condition
            if sum(condition_mask) == 0:
                print(f"  No cells found for condition '{condition}', skipping")
                continue
                
            # Create temporary AnnData with only this condition
            temp_adata = adata[condition_mask].copy()
            
            # Check if the group exists in this condition
            if group not in temp_adata.obs[groupby].cat.categories:
                print(f"  Group '{group}' not found in condition '{condition}', skipping")
                continue
            
            # Run find_markers
            try:
                sc.tl.rank_genes_groups(
                    temp_adata,
                    groupby=groupby,
                    groups=group,
                    reference='rest',
                    method=method,
                    pts=pts,
                    layer=layer
                )
                
                # Get results as DataFrame
                result = sc.get.rank_genes_groups_df(temp_adata, group=group)
                
                # Filter by fold change and p-value
                if 'logfoldchanges' in result.columns:
                    result['fold_change'] = np.exp(result['logfoldchanges'])
                    result = result[
                        (result['fold_change'] >= min_fold_change) &
                        (result['pvals'] <= max_pval)
                    ]
                else:
                    # No fold change available, filter by p-value only
                    result = result[result['pvals'] <= max_pval]
                
                # Store results
                if not result.empty:
                    group_markers[condition] = result
                else:
                    print(f"  No significant markers found for group '{group}' in condition '{condition}'")
            
            except Exception as e:
                print(f"  Error finding markers for group '{group}' in condition '{condition}': {str(e)}")
        
        # Find conserved markers across conditions
        if group_markers:
            # Get all marker genes
            all_markers = set()
            for df in group_markers.values():
                all_markers.update(df['names'])
            
            # Count in how many conditions each gene is significant
            gene_counts = {gene: 0 for gene in all_markers}
            gene_scores = {gene: 0.0 for gene in all_markers}
            gene_fold_changes = {gene: [] for gene in all_markers}
            
            for condition, df in group_markers.items():
                for _, row in df.iterrows():
                    gene = row['names']
                    gene_counts[gene] += 1
                    gene_scores[gene] += -np.log10(row['pvals'])
                    if 'fold_change' in row:
                        gene_fold_changes[gene].append(row['fold_change'])
            
            # Filter genes by minimum number of conditions
            conserved = [
                gene for gene in all_markers 
                if gene_counts[gene] >= min_significant_conditions
            ]
            
            # Create DataFrame with conserved markers
            if conserved:
                conserved_df = pd.DataFrame({
                    'gene': conserved,
                    'n_conditions': [gene_counts[gene] for gene in conserved],
                    'avg_log10_pval': [gene_scores[gene] / gene_counts[gene] for gene in conserved],
                    'avg_fold_change': [
                        np.mean(gene_fold_changes[gene]) if gene_fold_changes[gene] else np.nan
                        for gene in conserved
                    ]
                })
                
                # Sort by number of conditions and then by average score
                conserved_df = conserved_df.sort_values(
                    ['n_conditions', 'avg_log10_pval'], 
                    ascending=[False, False]
                )
                
                # Keep top n_genes
                conserved_df = conserved_df.head(n_genes)
                
                # Store results
                conserved_markers[group] = conserved_df
                
                print(f"  Found {len(conserved_df)} conserved markers for group '{group}'")
            else:
                print(f"  No conserved markers found for group '{group}'")
        else:
            print(f"  No markers found for group '{group}' in any condition")
    
    return conserved_markers