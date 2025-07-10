"""
Normalization functions for single-cell RNA-seq data.
"""

import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import numpy as np
from scipy import sparse
from typing import Optional, List, Union, Literal

def normalize_data(
    adata: sc.AnnData,
    method: Literal["standard", "scran", "pearson_residuals"] = "standard",
    layer: Optional[str] = None,
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    log_transform: bool = True,
    output_layer: str = "log1p_norm",
    plot: bool = True,
    save_dir: Optional[str] = None,
):
    """
    Normalize and log-transform the single-cell data.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        method (str): Normalization method to use. Options: "standard", "scran", "pearson_residuals"
        layer (str, optional): Name of the layer in adata.layers containing the raw count data. 
                               If None, use adata.X. Defaults to None.
        target_sum (float, optional): Total count to which each cell should be normalized. Defaults to 1e4.
            If None, after normalization, each observation (cell) has a total count equal to the median of total counts.
            If target_sum=1e6, this is CPM normalization.
        exclude_highly_expressed (bool, optional): Whether to exclude very highly expressed genes for 
                                                 the normalization factor computation. Defaults to False.
        max_fraction (float, optional): If exclude_highly_expressed=True, consider cells as highly expressed 
                                       that have more counts than max_fraction of the original total counts. Defaults to 0.05.
        log_transform (bool): Whether to log-transform the data after normalization. Defaults to True.
        output_layer (str): Name of the output layer. Defaults to "log1p_norm".
        plot (bool, optional): Whether to plot the distributions before and after normalization. Defaults to True.
        save_dir (str, optional): Directory to save plots. If None, plots are not saved. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with normalized and log-transformed data.

    Raises:
        ValueError: If `target_sum` is negative or `max_fraction` is not between 0 and 1.
    """
    from scipy import sparse
    
    # Check input parameter validity
    if target_sum is not None and target_sum <= 0:
        raise ValueError("target_sum must be a positive number.")
    if not 0 < max_fraction < 1:
        raise ValueError("max_fraction must be between 0 and 1 (exclusive).")
    
    # Check normalization method
    valid_methods = ["standard", "scran", "pearson_residuals"]
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}")
    
    # Get raw data
    if layer is None:
        X_raw = adata.X.copy()
    else:
        if layer not in adata.layers:
            print(f"Layer '{layer}' not found in adata.layers. Creating a new layer.")
            adata.layers[layer] = adata.X.copy()
        X_raw = adata.layers[layer].copy()
    
    # Save raw counts if 'raw_counts' layer doesn't exist
    if 'raw_counts' not in adata.layers:
        adata.layers['raw_counts'] = X_raw.copy()
        print("Raw counts saved in adata.layers['raw_counts']")

    print(f"Begin data normalization using {method} method.")

    # Normalize data using the selected method
    if method == "standard":
        # Standard normalization
        if layer is None:
            X_norm = sc.pp.normalize_total(
                adata, 
                target_sum=target_sum, 
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                inplace=False, 
            )['X']
        else:
            X_norm = sc.pp.normalize_total(
                adata, 
                layer=layer,
                target_sum=target_sum, 
                exclude_highly_expressed=exclude_highly_expressed,
                max_fraction=max_fraction,
                inplace=False, 
            )['X']
    
    elif method == "scran":
        try:
            import rpy2
            import rpy2.robjects as ro
            from rpy2.robjects.packages import importr
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
            
            # Use scran from R via rpy2
            importr('scran')
            # Simplified example - in practice this would use rpy2 to call R's scran
            print("Using scran for normalization...")
            # This is a placeholder - actual implementation would involve R code
            X_norm = sc.pp.normalize_total(adata, layer=layer, inplace=False)['X']
            print("Note: This is currently a placeholder. Full scran implementation requires R integration.")
        except ImportError:
            raise ValueError("scran method requires rpy2. Please install with 'pip install rpy2'")
    
    elif method == "pearson_residuals":
        try:
            from scipy import sparse
            
            # Calculate means and variances
            if sparse.issparse(X_raw):
                means = np.array(X_raw.mean(axis=0)).flatten()
                X_raw_dense = X_raw.toarray()
            else:
                means = np.mean(X_raw, axis=0)
                X_raw_dense = X_raw
            
            # Calculate Pearson residuals: (x - mean) / sqrt(mean)
            residuals = (X_raw_dense - means) / np.sqrt(means + 0.1)  # Add 0.1 to avoid division by zero
            
            # Cap values to avoid extreme residuals
            residuals = np.clip(residuals, -10, 10)
            
            if sparse.issparse(X_raw):
                X_norm = sparse.csr_matrix(residuals)
            else:
                X_norm = residuals
                
        except Exception as e:
            raise ValueError(f"Error computing Pearson residuals: {str(e)}")
    
    # Store the normalized result
    adata.layers["normalized"] = X_norm.copy()
    
    # Log transform if requested
    if log_transform and method != "pearson_residuals":  # No need to log Pearson residuals
        print("Applying log1p transformation.")
        if sparse.issparse(X_norm):
            X_log = X_norm.copy()
            X_log.data = np.log1p(X_log.data)
        else:
            X_log = np.log1p(X_norm)
        adata.layers[output_layer] = X_log
    else:
        adata.layers[output_layer] = X_norm
    
    print("Normalization complete.")

    if plot:
        # Visualize the distributions
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        fig.suptitle(f"Data Distributions Before and After {method.capitalize()} Normalization", fontsize=16)

        # Plot total counts before normalization
        if sparse.issparse(X_raw):
            raw_sums = np.array(X_raw.sum(axis=1)).flatten()
        else:
            raw_sums = np.sum(X_raw, axis=1)
            
        sns.histplot(
            raw_sums,
            bins=100,
            kde=True,
            ax=axes[0],
            color="navy",
        )
        axes[0].set_title("Total Counts (Before Normalization)", fontsize=14)
        axes[0].set_xlabel("Total Counts", fontsize=12)

        # Plot values after normalization
        if sparse.issparse(X_norm):
            norm_sums = np.array(adata.layers[output_layer].sum(axis=1)).flatten()
        else:
            norm_sums = np.sum(adata.layers[output_layer], axis=1)
            
        sns.histplot(
            norm_sums,
            bins=100,
            kde=True,
            ax=axes[1],
            color="crimson",
        )
        title_suffix = " (Log-Transformed)" if log_transform and method != "pearson_residuals" else ""
        axes[1].set_title(f"After {method.capitalize()} Normalization{title_suffix}", fontsize=14)
        axes[1].set_xlabel("Values", fontsize=12)

        plt.tight_layout()
        
        if save_dir:
            import os
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f"normalization_{method}.png"), dpi=300)
        
        plt.show()

    return adata

def regress_out(
    adata: sc.AnnData,
    keys: list,
    layer: Optional[str] = "log1p_norm",
    n_jobs: int = None,
    output_layer: str = "regressed_out",
) -> sc.AnnData:
    """
    Regress out unwanted sources of variation.
    
    Args:
        adata: AnnData object
        keys: Variables to regress out (must be in adata.obs)
        layer: Layer to use as input
        n_jobs: Number of parallel jobs
        output_layer: Layer to store result
        
    Returns:
        AnnData with regressed out variables
    """
    # Check if keys exist
    for key in keys:
        if key not in adata.obs:
            raise ValueError(f"Key '{key}' not found in adata.obs")
    
    # Determine data source
    if layer is not None and layer in adata.layers:
        print(f"Using layer '{layer}' for regression")
        X_backup = adata.X.copy()
        adata.X = adata.layers[layer].copy()
    
    # Regress out the specified variables
    print(f"Regressing out: {', '.join(keys)}")
    sc.pp.regress_out(adata, keys=keys, n_jobs=n_jobs)
    
    # Store result
    if output_layer is not None:
        adata.layers[output_layer] = adata.X.copy()
        print(f"Regressed data stored in adata.layers['{output_layer}']")
    
    # Restore original data if needed
    if layer is not None and layer in adata.layers:
        adata.X = X_backup
    
    return adata