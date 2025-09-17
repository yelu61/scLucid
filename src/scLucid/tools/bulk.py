"""
A bridge for running R-based single-cell analysis tools.

This module uses rpy2 to create a robust interface between Python/AnnData and
popular R/Bioconductor packages like CopyKAT, Monocle3, and CellChat.
"""

import logging
from typing import Literal, Optional

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu, pearsonr, spearmanr
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

log = logging.getLogger(__name__)

# Activate R-to-Python data conversions
pandas2ri.activate()
numpy2ri.activate()

def deconvolve_bulk(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    sample_key: str, # Required for BayesPrism
    method: Literal["DWLS", "BayesPrism"] = "BayesPrism",
    r_tools_instance: Optional[RTools] = None,
    key_added: str = "bulk_deconvolution"
) -> AnnData:
    """
    Estimate cell type proportions in bulk RNA-seq data.

    Args:
        ... (docstring) ...
        method: Deconvolution method to use.
        ...
    """
    if r_tools_instance is None:
        r_tools_instance = RTools()

    common_genes = adata_ref.var_names.intersection(bulk_data.index)
    log.info(f"Found {len(common_genes)} common genes for deconvolution.")
    
    adata_ref_sub = adata_ref[:, common_genes].copy()
    bulk_data_sub = bulk_data.loc[common_genes]

    # --- Method Dispatch ---
    if method == "Bisque":
        proportions_df = r_tools_instance.run_bisque_deconvolution(
            adata_ref=adata_ref_sub, bulk_data=bulk_data_sub,
            cell_type_key=cell_type_key, sample_key=sample_key
        )
    elif method == "DWLS":
        proportions_df = r_tools_instance.run_dwls_deconvolution(
            adata_ref=adata_ref_sub, bulk_data=bulk_data_sub, cell_type_key=cell_type_key
        )
    elif method == "BayesPrism":
        proportions_df = r_tools_instance.run_bayesprism_deconvolution(
            adata_ref=adata_ref_sub, bulk_data=bulk_data_sub,
            cell_type_key=cell_type_key, sample_key=sample_key
        )
    else:
        raise ValueError(f"Unsupported deconvolution method: {method}")

    # Store results
    adata_ref.uns.setdefault('sclucid', {}).setdefault('tools', {})
    adata_ref.uns['sclucid']['tools'][key_added] = {
        'proportions': proportions_df,
        'params': {'method': method, 
                   }
    }
    log.info(f"Deconvolution results stored in .uns.")
    
    return adata_ref


def differential_abundance(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    group1: str,
    group2: str,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon"
) -> pd.DataFrame:
    """
    Perform differential abundance analysis on deconvolution results.

    Args:
        proportions_df: DataFrame of cell type proportions from deconvolution.
        metadata_df: DataFrame with clinical/experimental metadata for bulk samples.
                     Must be indexed by sample ID.
        group_col: Column in metadata_df defining the groups to compare.
        group1: First group for comparison.
        group2: Second group for comparison.
        method: Statistical test to use.

    Returns:
        A DataFrame with differential abundance results for each cell type.
    """
    # Align data
    data = proportions_df.join(metadata_df)
    
    group1_samples = data[data[group_col] == group1].index
    group2_samples = data[data[group_col] == group2].index
    
    results = []
    for cell_type in proportions_df.columns:
        scores1 = data.loc[group1_samples, cell_type].dropna()
        scores2 = data.loc[group2_samples, cell_type].dropna()

        if len(scores1) < 2 or len(scores2) < 2:
            continue
            
        if method == "wilcoxon":
            stat, pval = mannwhitneyu(scores1, scores2)
        else:
            stat, pval = ttest_ind(scores1, scores2)
            
        results.append({
            "cell_type": cell_type,
            "statistic": stat,
            "pvalue": pval,
            "mean_abundance_group1": scores1.mean(),
            "mean_abundance_group2": scores2.mean(),
            "log2fc_abundance": np.log2(scores1.mean() / scores2.mean())
        })
        
    results_df = pd.DataFrame(results).sort_values("pvalue")
    return results_df

def correlate_abundance_with_clinical(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    clinical_variable: str,
    method: Literal["pearson", "spearman"] = "spearman"
) -> pd.DataFrame:
    """
    Correlate cell type abundance with a continuous clinical variable.

    Args:
        proportions_df: DataFrame of cell type proportions.
        metadata_df: DataFrame with clinical metadata, indexed by sample ID.
        clinical_variable: A continuous variable in metadata_df to correlate with.
        method: Correlation method to use.

    Returns:
        A DataFrame with correlation results for each cell type.
    """
    data = proportions_df.join(metadata_df)
    
    results = []
    for cell_type in proportions_df.columns:
        subset = data[[cell_type, clinical_variable]].dropna()
        if len(subset) < 5: # Need a few points to correlate
            continue
        
        if method == "pearson":
            corr, pval = pearsonr(subset[cell_type], subset[clinical_variable])
        else:
            corr, pval = spearmanr(subset[cell_type], subset[clinical_variable])
            
        results.append({
            "cell_type": cell_type,
            "clinical_variable": clinical_variable,
            "correlation_coefficient": corr,
            "pvalue": pval
        })
        
    return pd.DataFrame(results).sort_values("pvalue")