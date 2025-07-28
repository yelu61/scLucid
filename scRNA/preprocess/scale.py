"""
Data scaling functions for single-cell RNA-seq data.

This module provides functions for scaling gene expression data,
including z-score normalization and robust scaling options.
The main function, scale_data, serves as a wrapper around scanpy's scaling
functionality with enhanced error handling and options.
"""

import logging
import os
from typing import Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from ..utils.utils import use_layer_as_X

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "scale_data",
    "plot_scaling_effect",
]


def scale_data(
    adata: AnnData,
    layer: Optional[str] = "log1p_norm",
    output_layer: str = "scaled",
    max_value: Optional[float] = 10.0,
    zero_center: bool = True,
    vars_to_regress: Optional[list] = None,
    scale_method: Literal["zscore", "robust", "minmax"] = "zscore",
    subset_highly_variable: bool = False,
    plot: bool = False,
    save_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Scales gene expression data to unit variance and optionally zero mean.

    This function is a wrapper around `scanpy.pp.scale` with additional options
    and better error handling. It operates on a specified input layer and saves
    the result to a new output layer. It is recommended to run this on
    log-normalized data.

    Args:
        adata: AnnData object.
        layer: Layer to use for scaling. Typically log-normalized data.
              If None, uses adata.X.
        output_layer: Layer to store the scaled data.
        max_value: Clip values exceeding this value. If None, no clipping is performed.
        zero_center: If True, center the data to zero mean.
        vars_to_regress: Variables to regress out during scaling (must be in adata.obs).
                        This integrates scaling and regression in one step.
        scale_method: Scaling method to use:
                     - "zscore": Standard z-score normalization (default)
                     - "robust": Robust scaling using median and MAD
                     - "minmax": Min-max scaling to [0,1] range
        subset_highly_variable: Whether to scale only highly variable genes.
                               If True, only scales genes where adata.var['highly_variable'] is True.
        plot: Whether to generate plots showing the effect of scaling.
        save_dir: Directory to save plots. If None, plots are not saved.
        force: Whether to overwrite existing output_layer if it exists.

    Returns:
        The modified AnnData object with the new scaled layer.

    Raises:
        ValueError: If the input layer doesn't exist or parameters are invalid.
        RuntimeError: If scaling fails due to computational issues.

    Examples:
        >>> # Standard workflow: normalize, then scale
        >>> adata = pp.normalize_data(adata)
        >>> adata = pp.scale_data(adata, layer="log1p_norm", output_layer="scaled")
        >>>
        >>> # Scale only highly variable genes
        >>> adata = pp.scale_data(adata, subset_highly_variable=True)
        >>>
        >>> # Use robust scaling (less sensitive to outliers)
        >>> adata = pp.scale_data(adata, scale_method="robust", max_value=None)
    """
    # Parameter validation
    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    if layer is not None and layer not in adata.layers:
        available_layers = list(adata.layers.keys())
        raise ValueError(
            f"Layer '{layer}' not found in adata.layers. "
            f"Available layers: {available_layers}"
        )

    # Validate vars_to_regress
    if vars_to_regress:
        missing_vars = [var for var in vars_to_regress if var not in adata.obs.columns]
        if missing_vars:
            raise ValueError(
                f"The following variables for regression are missing from adata.obs: "
                f"{', '.join(missing_vars)}"
            )

    # Validate subset_highly_variable
    if subset_highly_variable and "highly_variable" not in adata.var:
        raise ValueError(
            "subset_highly_variable=True but 'highly_variable' not found in adata.var. "
            "Run 'pp.highly_variable_genes()' first."
        )

    # Log the operation
    log.info(
        f"Scaling data from {'adata.X' if layer is None else f"layer '{layer}'"} "
        f"using {scale_method} scaling"
    )
    log.info(
        f"Parameters: max_value={max_value}, zero_center={zero_center}, "
        f"subset_highly_variable={subset_highly_variable}"
    )

    if vars_to_regress:
        log.info(
            f"Will regress out the following variables during scaling: {vars_to_regress}"
        )

    # Create a mask for highly variable genes if needed
    if subset_highly_variable:
        gene_mask = adata.var["highly_variable"].copy()
        n_hvgs = gene_mask.sum()
        log.info(
            f"Scaling only {n_hvgs} highly variable genes ({n_hvgs / adata.n_vars:.1%} of all genes)"
        )
    else:
        gene_mask = None

    # Check the input data for potential issues
    try:
        with use_layer_as_X(adata, layer):
            # Get a small sample to check statistics
            if scipy.sparse.issparse(adata.X):
                max_val = adata.X.max()
                min_val = adata.X.data.min() if adata.X.nnz > 0 else 0
                sparsity = 1 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
            else:
                max_val = np.max(adata.X)
                min_val = np.min(adata.X)
                sparsity = np.mean(adata.X == 0)

            log.info(
                f"Input data statistics: min={min_val:.3f}, max={max_val:.3f}, sparsity={sparsity:.2%}"
            )

            # Heuristic check: Warn user if data looks like raw counts
            if max_val > 100:
                log.warning(
                    f"Max value in layer '{layer}' is > 100. "
                    "Scaling is typically performed on log-normalized data, not raw counts."
                )

            # Additional warning for non-logged data
            if min_val >= 0 and max_val > 30 and sparsity > 0.5:
                log.warning(
                    "Data appears to be raw counts rather than log-normalized values. "
                    "Consider running normalize_data with log_transform=True first."
                )

            # Store original data for plotting if needed
            if plot:
                if scipy.sparse.issparse(adata.X):
                    original_data = adata.X.copy()
                else:
                    original_data = adata.X.copy()

            # Perform scaling based on the method
            if scale_method == "zscore":
                # Standard z-score scaling (scanpy default)
                sc.pp.scale(
                    adata,
                    max_value=max_value,
                    zero_center=zero_center,
                    vars_to_regress=vars_to_regress,
                    mask_genes=gene_mask,
                )

            elif scale_method == "robust":
                # Robust scaling using median and MAD
                _robust_scale(
                    adata,
                    max_value=max_value,
                    zero_center=zero_center,
                    mask_genes=gene_mask,
                )

            elif scale_method == "minmax":
                # Min-max scaling to [0,1]
                _minmax_scale(adata, feature_range=(0, 1), mask_genes=gene_mask)

            else:
                raise ValueError(
                    f"Unknown scale_method '{scale_method}'. "
                    f"Valid options: 'zscore', 'robust', 'minmax'"
                )

            # Store the result from the modified adata.X into the output layer
            adata.layers[output_layer] = adata.X.copy()

            # Restore original .X if using a layer
            if layer is not None:
                adata.X = adata.layers[layer].copy()

    except Exception as e:
        log.error(f"Scaling failed: {str(e)}")
        raise RuntimeError(f"Failed to scale data: {str(e)}")

    # Store metadata about the scaling
    adata.uns["scaling"] = {
        "method": scale_method,
        "input_layer": layer,
        "output_layer": output_layer,
        "zero_center": zero_center,
        "max_value": max_value,
        "vars_regressed": vars_to_regress,
        "subset_highly_variable": subset_highly_variable,
    }

    # Get statistics on the scaled data
    scaled_data = adata.layers[output_layer]
    if scipy.sparse.issparse(scaled_data):
        scaled_max = scaled_data.max()
        scaled_min = scaled_data.data.min() if scaled_data.nnz > 0 else 0
        scaled_mean = scaled_data.mean()
        scaled_std = np.sqrt(scaled_data.power(2).mean() - scaled_mean**2)
    else:
        scaled_max = np.max(scaled_data)
        scaled_min = np.min(scaled_data)
        scaled_mean = np.mean(scaled_data)
        scaled_std = np.std(scaled_data)

    log.info(
        f"Scaling complete. Output statistics: min={scaled_min:.3f}, max={scaled_max:.3f}, "
        f"mean={scaled_mean:.3f}, std={scaled_std:.3f}"
    )
    log.info(f"Scaled data stored in adata.layers['{output_layer}']")

    # Generate plots if requested
    if plot:
        try:
            plot_scaling_effect(
                adata,
                original_data=original_data,
                scaled_layer=output_layer,
                save_dir=save_dir,
            )
        except Exception as e:
            log.warning(f"Failed to generate scaling effect plots: {str(e)}")

    return adata


def _robust_scale(
    adata: AnnData,
    max_value: Optional[float] = None,
    zero_center: bool = True,
    mask_genes: Optional[np.ndarray] = None,
) -> None:
    """
    Scale data using robust scaling (median and MAD).

    This is less sensitive to outliers than standard z-score scaling.

    Args:
        adata: AnnData object
        max_value: Maximum value to clip at
        zero_center: Whether to center to zero median
        mask_genes: Boolean mask for genes to scale

    Note:
        This function modifies adata.X in place.
    """
    # Handle gene mask
    if mask_genes is not None:
        genes_to_scale = np.where(mask_genes)[0]
    else:
        genes_to_scale = np.arange(adata.n_vars)

    # Get the data matrix
    X = adata.X

    # Convert to dense if sparse
    if scipy.sparse.issparse(X):
        # For sparse matrices, we operate only on the genes we need to scale
        if len(genes_to_scale) < adata.n_vars:
            # Only convert the necessary genes to dense
            for gene_idx in genes_to_scale:
                gene_vector = X[:, gene_idx].toarray().flatten()

                # Compute median and MAD
                med = np.median(gene_vector)
                mad = np.median(np.abs(gene_vector - med))

                # Handle zero MAD
                if mad == 0:
                    mad = 1.0

                # Center if requested
                if zero_center:
                    gene_vector = gene_vector - med

                # Scale by MAD
                gene_vector = gene_vector / mad

                # Clip if requested
                if max_value is not None:
                    gene_vector = np.clip(gene_vector, -max_value, max_value)

                # Put back in the matrix
                X[:, gene_idx] = scipy.sparse.csr_matrix(gene_vector).T
        else:
            # Convert the whole matrix
            X_dense = X.toarray()

            # Compute median and MAD for each gene
            gene_medians = np.median(X_dense, axis=0)
            gene_mads = np.median(np.abs(X_dense - gene_medians), axis=0)

            # Replace zero MADs with 1
            gene_mads[gene_mads == 0] = 1.0

            # Center if requested
            if zero_center:
                X_dense = X_dense - gene_medians

            # Scale by MAD
            X_dense = X_dense / gene_mads

            # Clip if requested
            if max_value is not None:
                X_dense = np.clip(X_dense, -max_value, max_value)

            # Convert back to sparse
            adata.X = scipy.sparse.csr_matrix(X_dense)
    else:
        # For dense matrices
        if len(genes_to_scale) < adata.n_vars:
            # Only scale selected genes
            gene_medians = np.zeros(adata.n_vars)
            gene_mads = np.ones(adata.n_vars)

            for gene_idx in genes_to_scale:
                gene_vector = X[:, gene_idx]
                gene_medians[gene_idx] = np.median(gene_vector)
                gene_mads[gene_idx] = np.median(
                    np.abs(gene_vector - gene_medians[gene_idx])
                )
                if gene_mads[gene_idx] == 0:
                    gene_mads[gene_idx] = 1.0

            # Center if requested
            if zero_center:
                X = X - gene_medians

            # Scale by MAD
            X = X / gene_mads

            # Clip if requested
            if max_value is not None:
                X = np.clip(X, -max_value, max_value)

            adata.X = X
        else:
            # Scale all genes
            gene_medians = np.median(X, axis=0)
            gene_mads = np.median(np.abs(X - gene_medians), axis=0)

            # Replace zero MADs with 1
            gene_mads[gene_mads == 0] = 1.0

            # Center if requested
            if zero_center:
                X = X - gene_medians

            # Scale by MAD
            X = X / gene_mads

            # Clip if requested
            if max_value is not None:
                X = np.clip(X, -max_value, max_value)

            adata.X = X


def _minmax_scale(
    adata: AnnData,
    feature_range: tuple = (0, 1),
    mask_genes: Optional[np.ndarray] = None,
) -> None:
    """
    Scale data to a fixed range using min-max scaling.

    Args:
        adata: AnnData object
        feature_range: (min, max) tuple of the desired range
        mask_genes: Boolean mask for genes to scale

    Note:
        This function modifies adata.X in place.
    """
    # Handle gene mask
    if mask_genes is not None:
        genes_to_scale = np.where(mask_genes)[0]
    else:
        genes_to_scale = np.arange(adata.n_vars)

    # Get the data matrix
    X = adata.X

    # Extract range parameters
    min_val, max_val = feature_range

    # Convert to dense if sparse
    if scipy.sparse.issparse(X):
        # For sparse matrices, we operate only on the genes we need to scale
        if len(genes_to_scale) < adata.n_vars:
            # Only convert the necessary genes to dense
            for gene_idx in genes_to_scale:
                gene_vector = X[:, gene_idx].toarray().flatten()

                # Compute min and max
                x_min = np.min(gene_vector)
                x_max = np.max(gene_vector)

                # Handle case when all values are the same
                if x_max == x_min:
                    gene_vector[:] = (min_val + max_val) / 2
                else:
                    # Scale to [0, 1]
                    gene_vector = (gene_vector - x_min) / (x_max - x_min)

                    # Scale to target range
                    gene_vector = gene_vector * (max_val - min_val) + min_val

                # Put back in the matrix
                X[:, gene_idx] = scipy.sparse.csr_matrix(gene_vector).T
        else:
            # Convert the whole matrix
            X_dense = X.toarray()

            # Compute min and max for each gene
            gene_mins = np.min(X_dense, axis=0)
            gene_maxs = np.max(X_dense, axis=0)

            # Handle case when all values are the same for a gene
            equal_genes = gene_maxs == gene_mins
            gene_range = gene_maxs - gene_mins
            gene_range[equal_genes] = 1.0  # Avoid division by zero

            # Scale to [0, 1]
            X_dense = (X_dense - gene_mins) / gene_range

            # Set constant genes to mid-range
            for gene_idx in np.where(equal_genes)[0]:
                X_dense[:, gene_idx] = 0.5

            # Scale to target range
            X_dense = X_dense * (max_val - min_val) + min_val

            # Convert back to sparse
            adata.X = scipy.sparse.csr_matrix(X_dense)
    else:
        # For dense matrices
        if len(genes_to_scale) < adata.n_vars:
            # Only scale selected genes
            gene_mins = np.zeros(adata.n_vars)
            gene_maxs = np.zeros(adata.n_vars)

            for gene_idx in genes_to_scale:
                gene_vector = X[:, gene_idx]
                gene_mins[gene_idx] = np.min(gene_vector)
                gene_maxs[gene_idx] = np.max(gene_vector)

            # Handle case when all values are the same for a gene
            equal_genes = gene_maxs == gene_mins
            gene_range = gene_maxs - gene_mins
            gene_range[equal_genes] = 1.0  # Avoid division by zero

            # Scale to [0, 1]
            X_scaled = np.zeros_like(X)
            X_scaled = (X - gene_mins) / gene_range

            # Set constant genes to mid-range
            for gene_idx in np.where(equal_genes)[0]:
                if gene_idx in genes_to_scale:
                    X_scaled[:, gene_idx] = 0.5

            # Scale to target range
            X_scaled = X_scaled * (max_val - min_val) + min_val

            # Only replace the scaled genes
            X = X_scaled

        else:
            # Scale all genes
            gene_mins = np.min(X, axis=0)
            gene_maxs = np.max(X, axis=0)

            # Handle case when all values are the same for a gene
            equal_genes = gene_maxs == gene_mins
            gene_range = gene_maxs - gene_mins
            gene_range[equal_genes] = 1.0  # Avoid division by zero

            # Scale to [0, 1]
            X = (X - gene_mins) / gene_range

            # Set constant genes to mid-range
            for gene_idx in np.where(equal_genes)[0]:
                X[:, gene_idx] = 0.5

            # Scale to target range
            X = X * (max_val - min_val) + min_val

        adata.X = X


def plot_scaling_effect(
    adata: AnnData,
    original_data: Union[np.ndarray, scipy.sparse.spmatrix],
    scaled_layer: str = "scaled",
    n_genes: int = 5,
    gene_subset: Optional[list] = None,
    save_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot the effect of scaling on gene expression distributions.

    This function creates before/after distribution plots to visualize
    how scaling affects gene expression values.

    Args:
        adata: AnnData object with scaled data
        original_data: Original data before scaling
        scaled_layer: Layer containing scaled data
        n_genes: Number of top variable genes to plot
        gene_subset: Specific genes to plot instead of top variable
        save_dir: Directory to save the plot

    Returns:
        matplotlib Figure object
    """
    if scaled_layer not in adata.layers:
        raise ValueError(f"Scaled layer '{scaled_layer}' not found in adata.layers")

    # Determine genes to plot
    if gene_subset is not None:
        genes_to_plot = [g for g in gene_subset if g in adata.var_names]
        if not genes_to_plot:
            raise ValueError("None of the specified genes were found in the data")
    else:
        # Find the most variable genes
        if scipy.sparse.issparse(original_data):
            gene_vars = np.var(original_data.toarray(), axis=0)
        else:
            gene_vars = np.var(original_data, axis=0)

        top_idx = np.argsort(-gene_vars)[:n_genes]
        genes_to_plot = adata.var_names[top_idx].tolist()

    # Create the plot
    fig, axes = plt.subplots(
        len(genes_to_plot), 2, figsize=(12, 3 * len(genes_to_plot))
    )
    fig.suptitle("Effect of Scaling on Gene Expression Distributions", fontsize=16)

    # Ensure axes is 2D even with a single gene
    if len(genes_to_plot) == 1:
        axes = np.array([axes])

    # Plot each gene
    for i, gene in enumerate(genes_to_plot):
        gene_idx = np.where(adata.var_names == gene)[0][0]

        # Get original data for this gene
        if scipy.sparse.issparse(original_data):
            orig_data = original_data[:, gene_idx].toarray().flatten()
        else:
            orig_data = original_data[:, gene_idx].flatten()

        # Get scaled data for this gene
        scaled_data = adata.layers[scaled_layer][:, gene_idx]
        if scipy.sparse.issparse(scaled_data):
            scaled_data = scaled_data.toarray().flatten()
        else:
            scaled_data = scaled_data.flatten()

        # Plot original distribution
        sns.histplot(orig_data, bins=30, kde=True, ax=axes[i, 0])
        axes[i, 0].set_title(f"{gene} - Before Scaling")
        axes[i, 0].text(
            0.05,
            0.95,
            f"Mean: {np.mean(orig_data):.2f}\nStd: {np.std(orig_data):.2f}",
            transform=axes[i, 0].transAxes,
            va="top",
        )

        # Plot scaled distribution
        sns.histplot(scaled_data, bins=30, kde=True, ax=axes[i, 1])
        axes[i, 1].set_title(f"{gene} - After Scaling")
        axes[i, 1].text(
            0.05,
            0.95,
            f"Mean: {np.mean(scaled_data):.2f}\nStd: {np.std(scaled_data):.2f}",
            transform=axes[i, 1].transAxes,
            va="top",
        )

    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, "scaling_effect.png"), dpi=300)

    return fig
