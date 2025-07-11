"""
Data scaling function for single-cell RNA-seq data.
"""

import scanpy as sc
import numpy as np
from typing import Optional
from .utils.anndata_helpers import use_layer_as_X

__all__ = [
    "scale_data",
]

def scale_data(
    adata: sc.AnnData,
    layer: str = "log1p_norm",
    output_layer: str = "scaled",
    max_value: Optional[float] = 10.0,
    zero_center: bool = True,
) -> sc.AnnData:
    """
    Scales data to unit variance and optionally zero mean.

    This function is a wrapper around `scanpy.pp.scale`. It operates on a
    specified input layer and saves the result to a new output layer. It is
    recommended to run this on log-normalized data.

    Note: This function modifies the AnnData object in place by adding a new layer.

    Args:
        adata: AnnData object.
        layer: Layer to use for scaling. Typically log-normalized data.
        output_layer: Layer to store the scaled data.
        max_value: Clip (truncate) values exceeding this value. If None, no clipping is performed.
        zero_center: If True, center the data to zero mean. `scanpy.pp.scale` default is True.

    Returns:
        The modified AnnData object with the new scaled layer.
        
    Example:
        >>> # Standard workflow: normalize, then scale
        >>> adata = pp.normalize_data(adata)
        >>> adata = pp.scale_data(adata, layer="log1p_norm", output_layer="scaled")
    """
    print(f"Scaling data from layer '{layer}' and saving to '{output_layer}'.")

    with use_layer_as_X(adata, layer):
        # Heuristic check: Warn user if data looks like raw counts
        # np.max works efficiently on both sparse and dense matrices
        if np.max(adata.X) > 100:
            print(f"Warning: Max value in layer '{layer}' is > 100. "
                  "Scaling is typically performed on log-normalized data, not raw counts.")

        # sc.pp.scale modifies adata.X in place
        sc.pp.scale(adata, max_value=max_value, zero_center=zero_center)
        
        # Store the result from the modified adata.X into the output layer
        adata.layers[output_layer] = adata.X.copy()
    
    print("Scaling complete.")
    return adata
