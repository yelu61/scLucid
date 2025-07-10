"""
Functions for identifying highly variable genes in single-cell RNA-seq data.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import itertools
from typing import Optional, List, Literal

def annotate_hvg(
    adata: sc.AnnData,
    method: str = 'scanpy',
    flavor: str = 'seurat',
    n_top_genes_scanpy: int = None,
    sample_key: str = "sampleID",
    min_cells_per_sample: int = 10,
    n_top_genes_custom: int = None,
    min_n_samples: int = 2,
    n_highly_expressed_genes: int = 100,
    n_specific_genes: int = 20,
    exclude_mt: bool = False,
    exclude_ribo: bool = False,
    batch_key: str = None,
    layer: str = "log1p_norm",
):
    """
    Annotate highly variable genes.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        method (str, optional): Method to use for identifying highly variable genes.
            Can be 'scanpy' or 'custom'. Defaults to 'scanpy'.
        flavor (str, optional): Flavor to use for Scanpy's built-in method.
            Can be 'seurat', 'seurat_v3' or 'pearson_residuals'. Defaults to 'seurat'.
        n_top_genes_scanpy (int or None, optional): For Scanpy's built-in method, the number of highly variable genes to keep. Defaults to 'None'.
            If None and flavor='seurat', Scanpy will automatically determine the number.
            If flavor='seurat_v3' or flavor='pearson_residuals', this argument is mandatory and must be provided.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        min_cells_per_sample (int, optional): Minimum number of cells per sample. Defaults to 10.
        n_top_genes_custom (int or None, optional): For the custom multi-criteria method, the number of highly variable genes to keep per sample. Defaults to None.
        min_n_samples (int, optional): For the custom multi-criteria method, minimum number of samples where a gene must be highly variable. Defaults to 2.
        n_highly_expressed_genes (int, optional): For the custom multi-criteria method, the number of top expressed genes to exclude from the highly variable gene set. Defaults to 100.
        n_specific_genes (int, optional): For the custom multi-criteria method, number of sample-specific genes to exclude from the highly variable gene set. Defaults to 20.
        exclude_mt (bool, optional): Whether to exclude mitochondrial genes from the highly variable gene set. Defaults to False.
        exclude_ribo (bool, optional): Whether to exclude ribosomal genes from the highly variable gene set. Defaults to False.
        batch_key (str, optional): The key in adata.obs to identify different batches for batch correction. Defaults to None.
        layer (str, optional): Layer to use for HVG detection. Defaults to "log1p_norm".

    Returns:
        adata (AnnData): AnnData object with highly variable genes identified.
    """
    # Backup X if we're using a layer
    if layer is not None and layer in adata.layers:
        X_backup = adata.X.copy()
        adata.X = adata.layers[layer].copy()
        restore_X = True
    else:
        restore_X = False

    if method == 'scanpy':
        # Input parameter validation for Scanpy method
        allowed_flavors = ['seurat', 'seurat_v3', 'pearson_residuals']
        if flavor not in allowed_flavors:
            raise ValueError(f"Invalid flavor '{flavor}'. Allowed values are: {', '.join(allowed_flavors)}")

        if flavor in ['seurat_v3', 'pearson_residuals'] and n_top_genes_scanpy is None:
            raise ValueError(f"`pp.highly_variable_genes` requires the argument `n_top_genes` for `flavor='{flavor}'`")

        if batch_key is not None and batch_key not in adata.obs.columns:
            raise KeyError(f"Key '{batch_key}' not found in adata.obs")
        
        if flavor == 'pearson_residuals':
            sc.experimental.pp.highly_variable_genes(
                adata,
                flavor=flavor,
                n_top_genes=n_top_genes_scanpy,
                batch_key=batch_key,
            )

        elif flavor in ['seurat', 'seurat_v3']:
            sc.pp.highly_variable_genes(
                adata,
                flavor=flavor,
                n_top_genes=n_top_genes_scanpy,
                batch_key=batch_key,
            )

        adata.var['highly_variable_scanpy'] = adata.var['highly_variable'].astype(bool)

        # Exclude mitochondrial and ribosomal genes
        if exclude_mt or exclude_ribo:
            mt_patterns = [r"^MT-", r"^mt-", r"^mt:", r"^MT:"] if exclude_mt else []
            rp_patterns = [r"^RP[SL]", r"^Rp[sl]"] if exclude_ribo else []
            exclude_patterns = mt_patterns + rp_patterns

            for pattern in exclude_patterns:
                adata.var['highly_variable_scanpy'] = adata.var['highly_variable_scanpy'] & ~adata.var_names.str.contains(pattern)

        print(f"Number of highly variable genes (Scanpy method): {sum(adata.var['highly_variable_scanpy'])}")

    elif method == 'custom':
        # Input parameter validation for custom method
        if sample_key not in adata.obs.columns:
            raise KeyError(f"Key '{sample_key}' not found in adata.obs")

        highly_variable_masks = []
        valid_samples = []

        for sample in adata.obs[sample_key].unique():
            sample_adata = adata[adata.obs[sample_key] == sample, :]
            
            if sample_adata.n_obs < min_cells_per_sample:
                print(f"Warning: Sample '{sample}' has fewer than {min_cells_per_sample} cells. Skipping this sample.")
                continue

            try:
                highly_variable_mask = sc.pp.highly_variable_genes(
                    sample_adata,
                    flavor='seurat',
                    n_top_genes=n_top_genes_custom,
                    inplace=False,
                )["highly_variable"]
                highly_variable_masks.append(highly_variable_mask)
                valid_samples.append(sample)
            except Exception as e:
                print(f"Error processing sample '{sample}': {str(e)}")

        if not highly_variable_masks:
            raise ValueError("No valid samples to process. Check your data and min_cells_per_sample threshold.")

        highly_variable = np.sum(highly_variable_masks, axis=0) >= min(min_n_samples, len(valid_samples))
        adata.var["highly_variable_custom"] = pd.Series(highly_variable, index=adata.var_names)
        print(f"Number of highly variable genes across samples (before filtering): {sum(adata.var['highly_variable_custom'])}")

        # Filter out highly expressed genes
        if sparse.issparse(adata.X):
            gene_expression_sum = np.array(adata.X.sum(axis=0)).flatten()
        else:
            gene_expression_sum = np.sum(adata.X, axis=0)
        gene_names = adata.var_names
        gene_expression_df = pd.DataFrame({'gene': gene_names, 'expression_sum': gene_expression_sum})
        top_genes = gene_expression_df.nlargest(n_highly_expressed_genes, 'expression_sum')['gene'].tolist()
        print(f"Length of top_genes: {len(top_genes)}")
        print(f"Top expressed genes: {', '.join(top_genes[:5])}...")

        # Identify sample-specific genes
        specific_genes = []
        if len(valid_samples) > 1:  # 只在有多个有效样本时执行
            try:
                valid_adata = adata[adata.obs[sample_key].isin(valid_samples)].copy()
            
                # 检查每个样本是否有足够的变异
                if sparse.issparse(valid_adata.X):
                    sample_vars = valid_adata.X.power(2).mean(axis=0).A1 - np.power(valid_adata.X.mean(axis=0).A1, 2)
                else:
                    sample_vars = np.var(valid_adata.X, axis=0)
            
                if np.any(sample_vars == 0):
                    print("Warning: Some genes have zero variance. Removing these genes for rank_genes_groups.")
                    valid_genes = sample_vars > 0
                    valid_adata = valid_adata[:, valid_genes]
            
                sc.tl.rank_genes_groups(
                    valid_adata,
                    groupby=sample_key, 
                    method="t-test", 
                    n_genes=n_specific_genes,
                )

                # 检查 rank_genes_groups 的结果是否按预期存储
                if 'rank_genes_groups' in valid_adata.uns and 'names' in valid_adata.uns['rank_genes_groups']:
                    specific_genes = list(
                        set(itertools.chain.from_iterable(valid_adata.uns["rank_genes_groups"]["names"]))
                    )
                else:
                    print("Warning: rank_genes_groups results not found in the expected format.")
                    print("Available keys in valid_adata.uns['rank_genes_groups']:", valid_adata.uns['rank_genes_groups'].keys())
            except Exception as e:
                print(f"Error in rank_genes_groups: {str(e)}")
    
            if not specific_genes:
                print("No sample-specific genes identified. Proceeding without sample-specific genes.")
            
        else:
            print("Not enough valid samples to perform rank_genes_groups. Skipping this step.")

        print(f"Length of specific_genes: {len(specific_genes)}")
        if specific_genes:
            print(f"Sample-specific genes: {', '.join(specific_genes[:5])}...")

        combined_genes = list(set(top_genes + specific_genes))
        print(f"Number of highest expressed or highly sample-specific genes: {len(combined_genes)}")

        final_highly_variable = adata.var["highly_variable_custom"] & ~adata.var_names.isin(combined_genes)

        # Exclude mitochondrial and ribosomal genes
        if exclude_mt or exclude_ribo:
            mt_patterns = [r"^MT-", r"^mt-", r"^mt:", r"^MT:"] if exclude_mt else []
            rp_patterns = [r"^RP[SL]", r"^Rp[sl]"] if exclude_ribo else []
            exclude_patterns = mt_patterns + rp_patterns
            for pattern in exclude_patterns:
                final_highly_variable = final_highly_variable & ~adata.var_names.str.contains(pattern)

        adata.var["highly_variable_custom"] = final_highly_variable
        print(f"Number of final highly variable genes (custom method): {sum(adata.var['highly_variable_custom'])}")

    # Restore original X if we used a layer
    if restore_X:
        adata.X = X_backup

    return adata


def select_hvg(
    adata: sc.AnnData,
    method: Literal["scanpy", "custom", "combined", "intersection"] = "scanpy",
    subset: bool = True,
    key: Optional[str] = None,
    n_top_genes: Optional[int] = None,
) -> sc.AnnData:
    """
    Select highly variable genes (HVGs) based on previous annotations.
    
    Args:
        adata: AnnData object with HVG annotations (run annotate_hvg first)
        method: Which method's results to use
        subset: Whether to subset adata to only include HVGs
        key: If provided, use this key in adata.var instead of the method
        n_top_genes: If provided, select only this many top genes
        
    Returns:
        AnnData with selected HVGs
    """
    # Determine which key to use
    if key is not None:
        if key not in adata.var:
            raise KeyError(f"Key '{key}' not found in adata.var")
        hvg_key = key
    else:
        if method == "scanpy":
            if "highly_variable_scanpy" not in adata.var:
                raise KeyError("Run annotate_hvg with method='scanpy' first")
            hvg_key = "highly_variable_scanpy"
        elif method == "custom":
            if "highly_variable_custom" not in adata.var:
                raise KeyError("Run annotate_hvg with method='custom' first")
            hvg_key = "highly_variable_custom"
        elif method == "combined":
            # Union of scanpy and custom HVGs
            if "highly_variable_scanpy" not in adata.var or "highly_variable_custom" not in adata.var:
                raise KeyError("Run annotate_hvg with both method='scanpy' and method='custom' first")
            adata.var["highly_variable_combined"] = (
                adata.var["highly_variable_scanpy"] | adata.var["highly_variable_custom"]
            )
            hvg_key = "highly_variable_combined"
        elif method == "intersection":
            # Intersection of scanpy and custom HVGs
            if "highly_variable_scanpy" not in adata.var or "highly_variable_custom" not in adata.var:
                raise KeyError("Run annotate_hvg with both method='scanpy' and method='custom' first")
            adata.var["highly_variable_intersection"] = (
                adata.var["highly_variable_scanpy"] & adata.var["highly_variable_custom"]
            )
            hvg_key = "highly_variable_intersection"
        else:
            raise ValueError(f"Unknown method: {method}")
    
    # If n_top_genes is provided, select top genes
    if n_top_genes is not None:
        if "dispersions_norm" in adata.var:
            dispersion_key = "dispersions_norm"
        elif "dispersions" in adata.var:
            dispersion_key = "dispersions"
        else:
            raise KeyError("No dispersion values found in adata.var")
        
        # Create mask of HVGs sorted by dispersion
        hvg_mask = adata.var[hvg_key]
        dispersions = adata.var[dispersion_key]
        top_genes = dispersions[hvg_mask].sort_values(ascending=False).index[:n_top_genes]
        
        # Create new mask with only top genes
        new_mask = pd.Series(False, index=adata.var_names)
        new_mask.loc[top_genes] = True
        
        # Store the result
        new_key = f"{hvg_key}_top{n_top_genes}"
        adata.var[new_key] = new_mask
        hvg_key = new_key
    
    # Print stats
    n_hvgs = sum(adata.var[hvg_key])
    print(f"Selected {n_hvgs} highly variable genes using {hvg_key}")
    
    # Subset if requested
    if subset:
        print(f"Subsetting adata to {n_hvgs} highly variable genes")
        adata = adata[:, adata.var[hvg_key]].copy()
    
    return adata