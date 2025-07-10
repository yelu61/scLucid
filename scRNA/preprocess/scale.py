"""
Scaling and cell cycle analysis for single-cell RNA-seq data.

This module provides functions for scaling data and scoring cell cycle phases.
"""

import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from scipy import sparse
from typing import Optional, List, Union, Literal

def scale_data(
    adata: sc.AnnData,
    max_value: float = 10.0,
    vars_to_regress: Optional[List[str]] = None,
    layer: Optional[str] = "log1p_norm",
    output_layer: str = "scaled",
    zero_center: bool = True,
) -> sc.AnnData:
    """
    Scale the data to unit variance and zero mean (if zero_center=True).
    
    Args:
        adata: AnnData object
        max_value: Truncate values exceeding this value
        vars_to_regress: Variables to regress out before scaling
        layer: Layer to scale. If None, use adata.X
        output_layer: Layer to store scaled data
        zero_center: Whether to center data to zero mean
        
    Returns:
        AnnData with scaled data
        
    Example:
        >>> adata = pp.normalize_data(adata)
        >>> adata = pp.scale_data(adata, layer="log1p_norm", output_layer="scaled")
    """
    # Store original data
    X_backup = None
    if layer is not None and layer in adata.layers:
        X_backup = adata.X.copy()
        adata.X = adata.layers[layer].copy()
    
    # Optionally regress out variables
    if vars_to_regress is not None:
        for var in vars_to_regress:
            if var not in adata.obs:
                raise ValueError(f"Variable '{var}' not found in adata.obs")
        
        print(f"Regressing out variables: {', '.join(vars_to_regress)}")
        sc.pp.regress_out(adata, vars_to_regress)
    
    # Scale data
    print("Scaling data...")
    sc.pp.scale(adata, max_value=max_value, zero_center=zero_center)
    
    # Store results in specified layer
    adata.layers[output_layer] = adata.X.copy()
    print(f"Scaled data stored in layer '{output_layer}'")
    
    # Restore original data
    if X_backup is not None:
        adata.X = X_backup
    
    return adata