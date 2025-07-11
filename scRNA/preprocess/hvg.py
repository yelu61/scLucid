"""
Functions for identifying highly variable genes in single-cell RNA-seq data.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import itertools
from typing import Optional, List, Literal
from .utils.anndata_helpers import use_layer_as_X

# --- Helper Function for the 'custom' method ---

def _get_sample_specific_genes(
    adata: sc.AnnData,
    sample_key: str,
    n_specific_genes: int,
) -> List[str]:
    """
    Identifies top genes specific to each sample using rank_genes_groups.
    
    This is a helper for the 'custom' HVG method to find genes that might
    drive batch effects rather than shared biological variation.
    """
    # Ensure there are multiple samples to compare
    if adata.obs[sample_key].nunique() <= 1:
        print("Info: Only one sample group found. Skipping sample-specific gene identification.")
        return []

    try:
        # Use a copy to avoid modifying the original object's .uns field
        temp_adata = adata.copy()
        
        # rank_genes_groups requires a 'log1p' normalized layer in .X
        if 'log1p' not in str(temp_adata.X[:5]): # Heuristic check
             sc.pp.log1p(temp_adata)

        sc.tl.rank_genes_groups(
            temp_adata,
            groupby=sample_key,
            method="t-test",
            n_genes=n_specific_genes,
        )
        
        # Extract gene names from the result
        specific_genes_df = pd.DataFrame(temp_adata.uns["rank_genes_groups"]["names"])
        return list(np.unique(specific_genes_df.values.flatten()))

    except Exception as e:
        print(f"Warning: Could not identify sample-specific genes due to an error: {e}")
        return []

# --- Main Functions ---

def annotate_hvg(
    adata: sc.AnnData,
    method: Literal['scanpy', 'custom'] = 'scanpy',
    layer: str = "log1p_norm",
    n_top_genes: Optional[int] = 2000,
    # Scanpy-specific args
    flavor: Literal['seurat', 'seurat_v3', 'pearson_residuals'] = 'seurat',
    batch_key: Optional[str] = None,
    # Custom-specific args
    sample_key: str = "sampleID",
    min_n_samples: int = 2,
    n_highly_expressed_genes: int = 100,
    n_specific_genes: int = 20,
    exclude_mt_ribo: bool = True,
) -> sc.AnnData:
    """
    Annotates highly variable genes using one of two methods.

    This function adds a boolean column to `adata.var` indicating HVGs.

    Args:
        adata: AnnData object.
        method: Method for HVG selection ('scanpy' or 'custom').
        layer: Layer to use for HVG detection.
        n_top_genes: Number of top HVGs to select. Behavior depends on the method.
        flavor: For 'scanpy' method, the flavor of HVG selection.
        batch_key: For 'scanpy' method, batch key for batch-aware HVG selection.
        sample_key: For 'custom' method, the key in adata.obs for per-sample HVG selection.
        min_n_samples: For 'custom' method, min number of samples a gene must be HVG in.
        n_highly_expressed_genes: For 'custom' method, number of top expressed genes to exclude.
        n_specific_genes: For 'custom' method, number of sample-specific genes to exclude.
        exclude_mt_ribo: Whether to exclude mitochondrial and ribosomal genes.

    Returns:
        AnnData object with HVG annotations in `adata.var`.
    """
    output_key = f"highly_variable_{method}"
    
    with use_layer_as_X(adata, layer):
        if method == 'scanpy':
            print(f"Running HVG selection with 'scanpy' method (flavor: {flavor}).")
            if flavor in ['seurat_v3', 'pearson_residuals'] and n_top_genes is None:
                raise ValueError(f"flavor='{flavor}' requires `n_top_genes` to be set.")
            
            sc.pp.highly_variable_genes(
                adata,
                flavor=flavor,
                n_top_genes=n_top_genes,
                batch_key=batch_key,
                inplace=True
            )
            # Rename the default output for clarity
            adata.var[output_key] = adata.var['highly_variable']
            del adata.var['highly_variable']

        elif method == 'custom':
            print("Running HVG selection with 'custom' multi-criteria method.")
            if sample_key not in adata.obs.columns:
                raise KeyError(f"Sample key '{sample_key}' not found in adata.obs")

            # 1. Find HVGs per sample
            hvg_masks = []
            for sample in adata.obs[sample_key].unique():
                sample_adata = adata[adata.obs[sample_key] == sample].copy()
                if sample_adata.n_obs > 10: # Min cells for seurat flavor
                    sc.pp.highly_variable_genes(sample_adata, n_top_genes=n_top_genes, inplace=True)
                    hvg_masks.append(sample_adata.var['highly_variable'])
            
            if not hvg_masks:
                raise ValueError("No samples had enough cells to compute HVGs.")

            # 2. Combine HVGs across samples
            combined_hvgs = pd.concat(hvg_masks, axis=1).sum(axis=1) >= min_n_samples
            print(f"Found {combined_hvgs.sum()} genes considered HVG in at least {min_n_samples} samples.")

            # 3. Identify genes to exclude
            # Highly expressed genes
            gene_expr = np.array(adata.X.sum(axis=0)).flatten()
            top_expr_genes = adata.var_names[np.argsort(-gene_expr)[:n_highly_expressed_genes]]
            
            # Sample-specific genes
            specific_genes = _get_sample_specific_genes(adata, sample_key, n_specific_genes)
            
            exclude_genes = set(top_expr_genes) | set(specific_genes)
            print(f"Identified {len(exclude_genes)} highly expressed or sample-specific genes to exclude.")

            # 4. Final HVG list
            adata.var[output_key] = combined_hvgs & ~adata.var_names.isin(exclude_genes)

        else:
            raise ValueError(f"Unknown method '{method}'. Choose 'scanpy' or 'custom'.")

    # Common exclusion step for both methods
    if exclude_mt_ribo:
        mt_mask = adata.var_names.str.contains(r'^(MT|mt)-', regex=True)
        ribo_mask = adata.var_names.str.contains(r'^(RP[SL]|Rp[sl])', regex=True)
        exclude_mask = mt_mask | ribo_mask
        adata.var[output_key] &= ~exclude_mask
        print(f"Excluded {exclude_mask.sum()} mitochondrial/ribosomal genes.")

    print(f"Final number of highly variable genes ({method} method): {adata.var[output_key].sum()}")
    return adata


def select_hvg(
    adata: sc.AnnData,
    hvg_key: str,
    n_top_genes: Optional[int] = None,
    subset: bool = True,
) -> Optional[sc.AnnData]:
    """
    Selects a final set of HVGs and optionally subsets the AnnData object.

    Args:
        adata: AnnData object with HVG annotations.
        hvg_key: The column in `adata.var` to use (e.g., 'highly_variable_custom').
        n_top_genes: If provided, select the top N genes based on normalized dispersion.
        subset: If True, subsets `adata` to the selected HVGs and returns it.
                If False, adds a new boolean mask to `adata.var` and returns None.

    Returns:
        A subsetted AnnData object if `subset=True`, otherwise None.
    """
    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in `adata.var`. Please run `annotate_hvg` first.")

    final_hvg_mask = adata.var[hvg_key].copy()
    
    if n_top_genes is not None:
        if "dispersions_norm" not in adata.var:
            raise KeyError("`dispersions_norm` not found in `adata.var`. Needed to select top N genes.")
        
        # Sort genes by dispersion, but only consider those already marked as HVG
        dispersions = adata.var.loc[final_hvg_mask, 'dispersions_norm']
        top_genes = dispersions.nlargest(n_top_genes).index
        
        # Create a new mask with only the top N genes
        final_hvg_mask = adata.var_names.isin(top_genes)
        print(f"Selected top {final_hvg_mask.sum()} HVGs based on dispersion.")

    # Update adata.var['highly_variable'] which is the default for downstream tools
    adata.var['highly_variable'] = final_hvg_mask
    print(f"Adata.var['highly_variable'] updated with {final_hvg_mask.sum()} genes.")

    if subset:
        print("Subsetting AnnData object to selected HVGs.")
        return adata[:, adata.var['highly_variable']].copy()
    
    return None # Return None if inplace
