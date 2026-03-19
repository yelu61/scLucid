"""
Differential expression analysis for pyMonocle3 (R-free)
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ranksums, mannwhitneyu, ttest_ind, chi2_contingency
from typing import Optional, List, Union
import logging

from .core import CellDataSet

log = logging.getLogger(__name__)


def top_markers(
    cds: CellDataSet,
    group_cells_by: str = "cluster",
    genes_to_test_per_group: Optional[int] = None,
    reduction_method: str = "PCA",
    marker_sig_test: bool = True,
    reference_cells: Optional[pd.Series] = None,
    speedglm_maxiter: int = 25,
    cores: int = 1,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Find top marker genes for each group

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    group_cells_by : str
        Column to group cells by
    genes_to_test_per_group : int, optional
        Number of top genes to test per group
    reduction_method : str
        Reduction method for residuals
    marker_sig_test : bool
        Perform significance testing
    reference_cells : pd.Series, optional
        Reference cells for comparison
    speedglm_maxiter : int
        Maximum iterations for GLM
    cores : int
        Number of cores
    verbose : bool
        Verbose output

    Returns
    -------
    pd.DataFrame
        Top markers with statistics
    """
    if group_cells_by not in cds.cell_metadata.columns:
        raise ValueError(f"Column '{group_cells_by}' not found in cell metadata")

    groups = cds.cell_metadata[group_cells_by].unique()
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    markers = []

    for group in groups:
        group_cells = cds.cell_metadata[group_cells_by] == group
        other_cells = ~group_cells

        group_expr = expr[:, group_cells]
        other_expr = expr[:, other_cells]

        for i, gene in enumerate(cds.gene_metadata.index):
            g_expr = group_expr[i, :]
            o_expr = other_expr[i, :]

            # Skip if not expressed
            if np.sum(g_expr) == 0:
                continue

            # Calculate statistics
            group_mean = np.mean(g_expr)
            other_mean = np.mean(o_expr)
            group_pct = np.mean(g_expr > 0)
            other_pct = np.mean(o_expr > 0)

            # Fold change
            fc = group_mean / (other_mean + 1e-10)
            log2fc = np.log2(fc + 1e-10)

            # Statistical test
            if marker_sig_test:
                try:
                    stat, pval = ranksums(g_expr, o_expr)
                except:
                    pval = 1.0
            else:
                pval = 1.0

            markers.append({
                'gene': gene,
                'cell_group': group,
                'mean_expr': group_mean,
                'mean_other': other_mean,
                'pct_expr': group_pct,
                'pct_other': other_pct,
                'log2fc': log2fc,
                'pval': pval,
                'qval': pval,  # Placeholder, should do FDR correction
            })

    markers_df = pd.DataFrame(markers)

    if len(markers_df) > 0:
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests

        markers_df['qval'] = multipletests(
            markers_df['pval'].fillna(1),
            method='fdr_bh'
        )[1]

        markers_df = markers_df.sort_values(['cell_group', 'qval'])

    log.info(f"Found {len(markers_df)} marker gene entries for {len(groups)} groups")

    return markers_df


def aggregate_gene_expression(
    cds: CellDataSet,
    gene_group_df: pd.DataFrame,
    cell_group_df: pd.DataFrame,
    scale_agg_values: bool = True,
    max_agg_value: float = 3,
) -> pd.DataFrame:
    """
    Aggregate gene expression by groups

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    gene_group_df : pd.DataFrame
        Gene groups
    cell_group_df : pd.DataFrame
        Cell groups
    scale_agg_values : bool
        Scale aggregated values
    max_agg_value : float
        Maximum value after scaling

    Returns
    -------
    pd.DataFrame
        Aggregated expression matrix
    """
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    # Convert to DataFrame
    expr_df = pd.DataFrame(
        expr,
        index=cds.gene_metadata.index,
        columns=cds.cell_metadata.index
    )

    # Aggregate by cell groups
    cell_groups = cell_group_df.groupby(cell_group_df.columns[0])

    agg_expr = []
    for cell_group_name, cell_group_df in cell_groups:
        cells_in_group = cell_group_df.index
        group_expr = expr_df[cells_in_group].mean(axis=1)
        agg_expr.append(group_expr)

    result = pd.DataFrame(
        agg_expr,
        index=[name for name, _ in cell_groups],
        columns=cds.gene_metadata.index
    )

    # Scale if requested
    if scale_agg_values:
        result = result.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)
        result = result.clip(-max_agg_value, max_agg_value)

    return result.T


def compare_genes(
    cds: CellDataSet,
    group1_cells: List[str],
    group2_cells: List[str],
    test_method: str = "wilcoxon",
    min_pct: float = 0.1,
) -> pd.DataFrame:
    """
    Compare gene expression between two cell groups

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    group1_cells : list
        Cell IDs for group 1
    group2_cells : list
        Cell IDs for group 2
    test_method : str
        Statistical test method
    min_pct : float
        Minimum percentage of cells expressing the gene

    Returns
    -------
    pd.DataFrame
        Differential expression results
    """
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    # Get indices
    cell_idx_map = {c: i for i, c in enumerate(cds.cell_metadata.index)}
    group1_idx = [cell_idx_map[c] for c in group1_cells if c in cell_idx_map]
    group2_idx = [cell_idx_map[c] for c in group2_cells if c in cell_idx_map]

    results = []

    for i, gene in enumerate(cds.gene_metadata.index):
        g1_expr = expr[i, group1_idx]
        g2_expr = expr[i, group2_idx]

        # Percentage expressed
        g1_pct = np.mean(g1_expr > 0)
        g2_pct = np.mean(g2_expr > 0)

        if max(g1_pct, g2_pct) < min_pct:
            continue

        # Calculate statistics
        g1_mean = np.mean(g1_expr)
        g2_mean = np.mean(g2_expr)

        log2fc = np.log2((g1_mean + 1e-10) / (g2_mean + 1e-10))

        # Statistical test
        if test_method == "wilcoxon":
            try:
                stat, pval = ranksums(g1_expr, g2_expr)
            except:
                pval = 1.0
        elif test_method == "t-test":
            try:
                stat, pval = ttest_ind(g1_expr, g2_expr)
            except:
                pval = 1.0
        else:
            pval = 1.0

        results.append({
            'gene': gene,
            'group1_mean': g1_mean,
            'group2_mean': g2_mean,
            'group1_pct': g1_pct,
            'group2_pct': g2_pct,
            'log2fc': log2fc,
            'pval': pval,
        })

    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        results_df['qval'] = multipletests(results_df['pval'].fillna(1), method='fdr_bh')[1]
        results_df = results_df.sort_values('qval')

    return results_df


def pseudotime_de(
    cds: CellDataSet,
    fullModelFormulaStr: str = "~ splines::ns(pseudotime, df=3)",
    reducedModelFormulaStr: str = "~ 1",
    cores: int = 1,
) -> pd.DataFrame:
    """
    Differential expression along pseudotime

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet with pseudotime
    fullModelFormulaStr : str
        Formula for full model
    reducedModelFormulaStr : str
        Formula for reduced model
    cores : int
        Number of cores

    Returns
    -------
    pd.DataFrame
        Pseudotime DE results
    """
    if 'pseudotime' not in cds.cell_metadata.columns:
        raise ValueError("No pseudotime found. Run order_cells first.")

    # For now, simplified version using correlation
    # Full implementation would use regression models

    pseudotime = cds.cell_metadata['pseudotime'].values
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    results = []

    for i, gene in enumerate(cds.gene_metadata.index):
        gene_expr = expr[i, :]

        # Calculate correlation with pseudotime
        correlation = np.corrcoef(gene_expr, pseudotime)[0, 1]

        # Fit simple polynomial
        coeffs = np.polyfit(pseudotime, gene_expr, 3)
        fitted = np.polyval(coeffs, pseudotime)

        # R-squared
        ss_res = np.sum((gene_expr - fitted) ** 2)
        ss_tot = np.sum((gene_expr - np.mean(gene_expr)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        results.append({
            'gene': gene,
            'correlation': correlation,
            'r_squared': r_squared,
            'trend': 'up' if correlation > 0 else 'down',
        })

    return pd.DataFrame(results).sort_values('r_squared', ascending=False)


def calculate_gene_modules(
    cds: CellDataSet,
    genes: Optional[List[str]] = None,
    resolution: float = 1.0,
    k: int = 10,
) -> pd.DataFrame:
    """
    Group genes into modules based on expression correlation

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    genes : list, optional
        Genes to include
    resolution : float
        Resolution for clustering
    k : int
        Number of neighbors

    Returns
    -------
    pd.DataFrame
        Gene module assignments
    """
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    # Subset genes if specified
    if genes is not None:
        gene_mask = cds.gene_metadata.index.isin(genes)
        expr = expr[gene_mask, :]
        gene_names = cds.gene_metadata.index[gene_mask]
    else:
        gene_names = cds.gene_metadata.index

    # Calculate correlation matrix
    # Use a subset of cells for efficiency
    if expr.shape[1] > 5000:
        cell_subset = np.random.choice(expr.shape[1], 5000, replace=False)
        expr_subset = expr[:, cell_subset]
    else:
        expr_subset = expr

    # Calculate correlations
    corr_matrix = np.corrcoef(expr_subset)
    corr_matrix = np.nan_to_num(corr_matrix)

    # Convert to distance
    distance_matrix = 1 - np.abs(corr_matrix)

    # Cluster genes using hierarchical clustering
    from scipy.cluster.hierarchy import linkage, fcluster

    linkage_matrix = linkage(distance_matrix, method='ward')
    clusters = fcluster(linkage_matrix, t=resolution, criterion='maxclust')

    results = pd.DataFrame({
        'gene': gene_names,
        'module': clusters,
    })

    log.info(f"Identified {len(np.unique(clusters))} gene modules")

    return results
