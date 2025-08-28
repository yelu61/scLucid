"""
Data scaling functions for single-cell RNA-seq data.

This module provides functions for scaling gene expression data,
including z-score normalization and robust scaling options.
The main function, scale_data, serves as a wrapper around scanpy's scaling
functionality with enhanced error handling and options.
"""

import logging
import os
from typing import List, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from ..utils.utils import use_layer_as_X
from .config import ScalingConfig

# Configure logging
log = logging.getLogger(__name__)

__all__ = ["scale_data", "plot_scaling_effect", "regress_out"]


# --- Helper Functions ---#
def _robust_scale(
    adata: AnnData,
    max_value: Optional[float] = None,
    zero_center: bool = True,
) -> None:
    """
    Scale data using robust scaling (median and MAD).
    This is less sensitive to outliers than standard z-score scaling.

    Args:
        adata: AnnData object
        max_value: Maximum value to clip at
        zero_center: Whether to center to zero median

    Note:
        This function modifies adata.X in place.
    """
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    gene_medians = np.median(X, axis=0)
    gene_mads = np.median(np.abs(X - gene_medians), axis=0)
    gene_mads[gene_mads == 0] = 1.0
    if zero_center:
        X = X - gene_medians
    X = X / gene_mads
    if max_value is not None:
        X = np.clip(X, -max_value, max_value)
    adata.X = X


def _minmax_scale(
    adata: AnnData,
    feature_range: tuple = (0, 1),
) -> None:
    """
    Scale data to a fixed range using min-max scaling.

    Args:
        adata: AnnData object
        feature_range: (min, max) tuple of the desired range

    Note:
        This function modifies adata.X in place.
    """
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    gene_mins = np.min(X, axis=0)
    gene_maxs = np.max(X, axis=0)
    equal_genes = gene_maxs == gene_mins
    gene_range = gene_maxs - gene_mins
    gene_range[equal_genes] = 1.0
    X = (X - gene_mins) / gene_range
    for gene_idx in np.where(equal_genes)[0]:
        X[:, gene_idx] = 0.5
    min_val, max_val = feature_range
    X = X * (max_val - min_val) + min_val
    adata.X = X


# --- Main Functions ---#
def scale_data(
    adata: AnnData,
    input_layer: Optional[str] = "log1p_norm",
    output_layer: str = "scaled",
    max_value: Optional[float] = 10.0,
    zero_center: bool = True,
    vars_to_regress: Optional[List[str]] = None,
    scale_method: Literal["zscore", "robust", "minmax"] = "zscore",
    subset_highly_variable: bool = False,
    hvg_key: str = "highly_variable",
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
        input_layer: Layer to use for scaling. Typically log-normalized data.
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
                                If True, only scales genes where adata.var[hvg_key] is True.
        hvg_key: The column in adata.var to use as the HVG mask (default "highly_variable").
        plot: Whether to generate plots showing the effect of scaling.
        save_dir: Directory to save plots. If None, plots are not saved.
        force: Whether to overwrite existing output_layer if it exists.

    Returns:
        The modified AnnData object with the new scaled layer.

    Raises:
        ValueError: If the input layer doesn't exist or parameters are invalid.
        RuntimeError: If scaling fails due to computational issues.

    Examples:
        >>> adata = pp.normalize_data(adata)
        >>> adata = pp.scale_data(adata, input_layer="log1p_norm", output_layer="scaled")
        >>> adata = pp.scale_data(adata, subset_highly_variable=True, hvg_key="highly_variable_custom")
        >>> adata = pp.scale_data(adata, scale_method="robust", max_value=None)
    """
    # Parameter validation
    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    if input_layer is not None and input_layer not in adata.layers:
        available_layers = list(adata.layers.keys())
        raise ValueError(
            f"Layer '{input_layer}' not found in adata.layers. "
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
    if subset_highly_variable:
        if hvg_key not in adata.var:
            raise ValueError(
                f"subset_highly_variable=True but '{hvg_key}' not found in adata.var. "
                "Run HVG selection first or specify the correct hvg_key."
            )
        gene_mask = adata.var[hvg_key].values
        n_hvgs = gene_mask.sum()
        if n_hvgs == 0:
            raise ValueError(f"No genes are marked as True in adata.var['{hvg_key}']")
        log.info(
            f"Scaling only {n_hvgs} highly variable genes ({n_hvgs / adata.n_vars:.1%} of all genes) using key '{hvg_key}'"
        )
    else:
        gene_mask = None

    log.info(
        f"Scaling data from {'adata.X' if input_layer is None else f'layer {input_layer}'} "
        f"using {scale_method} scaling"
    )
    log.info(
        f"Parameters: max_value={max_value}, zero_center={zero_center}, "
        f"subset_highly_variable={subset_highly_variable}, hvg_key={hvg_key}"
    )
    if vars_to_regress:
        log.info(
            f"Will regress out the following variables during scaling: {vars_to_regress}"
        )

    try:
        with use_layer_as_X(adata, input_layer):
            # Save original data for plotting if requested
            if plot:
                if scipy.sparse.issparse(adata.X):
                    original_data = adata.X.copy()
                else:
                    original_data = adata.X.copy()

            # Regression (if requested)
            if vars_to_regress:
                if scipy.sparse.issparse(adata.X):
                    adata.X = adata.X.toarray()
                sc.pp.regress_out(adata, keys=vars_to_regress)

            # Scaling
            if scale_method == "zscore":
                if subset_highly_variable:
                    adata_hvg = adata[:, gene_mask].copy()
                    sc.pp.scale(
                        adata_hvg,
                        max_value=max_value,
                        zero_center=zero_center,
                    )
                    X_scaled = (
                        adata.X.toarray()
                        if scipy.sparse.issparse(adata.X)
                        else adata.X.copy()
                    )
                    X_scaled[:, gene_mask] = (
                        adata_hvg.X.toarray()
                        if scipy.sparse.issparse(adata_hvg.X)
                        else adata_hvg.X
                    )
                    if scipy.sparse.issparse(adata.X):
                        adata.layers[output_layer] = scipy.sparse.csr_matrix(X_scaled)
                    else:
                        adata.layers[output_layer] = X_scaled
                else:
                    sc.pp.scale(
                        adata,
                        max_value=max_value,
                        zero_center=zero_center,
                    )
                    adata.layers[output_layer] = adata.X.copy()

            elif scale_method == "robust":
                if subset_highly_variable:
                    adata_hvg = adata[:, gene_mask].copy()
                    _robust_scale(
                        adata_hvg,
                        max_value=max_value,
                        zero_center=zero_center,
                    )
                    X_scaled = (
                        adata.X.toarray()
                        if scipy.sparse.issparse(adata.X)
                        else adata.X.copy()
                    )
                    X_scaled[:, gene_mask] = (
                        adata_hvg.X.toarray()
                        if scipy.sparse.issparse(adata_hvg.X)
                        else adata_hvg.X
                    )
                    if scipy.sparse.issparse(adata.X):
                        adata.layers[output_layer] = scipy.sparse.csr_matrix(X_scaled)
                    else:
                        adata.layers[output_layer] = X_scaled
                else:
                    _robust_scale(
                        adata,
                        max_value=max_value,
                        zero_center=zero_center,
                    )
                    adata.layers[output_layer] = adata.X.copy()

            elif scale_method == "minmax":
                if subset_highly_variable:
                    adata_hvg = adata[:, gene_mask].copy()
                    _minmax_scale(adata_hvg, feature_range=(0, 1))
                    X_scaled = (
                        adata.X.toarray()
                        if scipy.sparse.issparse(adata.X)
                        else adata.X.copy()
                    )
                    X_scaled[:, gene_mask] = (
                        adata_hvg.X.toarray()
                        if scipy.sparse.issparse(adata_hvg.X)
                        else adata_hvg.X
                    )
                    if scipy.sparse.issparse(adata.X):
                        adata.layers[output_layer] = scipy.sparse.csr_matrix(X_scaled)
                    else:
                        adata.layers[output_layer] = X_scaled
                else:
                    _minmax_scale(adata, feature_range=(0, 1))
                    adata.layers[output_layer] = adata.X.copy()
            else:
                raise ValueError(
                    f"Unknown scale_method '{scale_method}'. "
                    f"Valid options: 'zscore', 'robust', 'minmax'"
                )

            # Restore original .X if using a layer
            if input_layer is not None:
                adata.X = adata.layers[input_layer].copy()

    except Exception as e:
        log.error(f"Scaling failed: {str(e)}")
        raise RuntimeError(f"Failed to scale data: {str(e)}")

    # Store method information in uns
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["scaling"] = {
        "method": scale_method,
        "input_layer": input_layer,
        "output_layer": output_layer,
        "zero_center": zero_center,
        "max_value": max_value,
        "vars_regressed": vars_to_regress,
        "subset_highly_variable": subset_highly_variable,
        "hvg_key": hvg_key,
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


def plot_scaling_effect(
    adata: AnnData,
    original_data: Union[np.ndarray, scipy.sparse.spmatrix],
    scaled_layer: str = "scaled",
    n_genes: int = 5,
    gene_subset: Optional[List[str]] = None,
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

    if gene_subset is not None:
        genes_to_plot = [g for g in gene_subset if g in adata.var_names]
        if not genes_to_plot:
            raise ValueError("None of the specified genes were found in the data")
    else:
        # Find top variable genes in original data
        data = (
            original_data.toarray()
            if scipy.sparse.issparse(original_data)
            else original_data
        )
        gene_vars = np.var(data, axis=0)
        top_idx = np.argsort(-gene_vars)[:n_genes]
        genes_to_plot = adata.var_names[top_idx].tolist()

    fig, axes = plt.subplots(
        len(genes_to_plot), 2, figsize=(12, 3 * len(genes_to_plot))
    )
    fig.suptitle("Effect of Scaling on Gene Expression Distributions", fontsize=16)
    if len(genes_to_plot) == 1:
        axes = np.array([axes])

    for i, gene in enumerate(genes_to_plot):
        gene_idx = np.where(adata.var_names == gene)[0][0]
        orig = original_data[:, gene_idx]
        if scipy.sparse.issparse(orig):
            orig = orig.toarray().flatten()
        else:
            orig = orig.flatten()
        scaled = adata.layers[scaled_layer][:, gene_idx]
        if scipy.sparse.issparse(scaled):
            scaled = scaled.toarray().flatten()
        else:
            scaled = scaled.flatten()
        sns.histplot(orig, bins=30, kde=True, ax=axes[i, 0])
        axes[i, 0].set_title(f"{gene} - Before Scaling")
        axes[i, 0].text(
            0.05,
            0.95,
            f"Mean: {np.mean(orig):.2f}\nStd: {np.std(orig):.2f}",
            transform=axes[i, 0].transAxes,
            va="top",
        )
        sns.histplot(scaled, bins=30, kde=True, ax=axes[i, 1])
        axes[i, 1].set_title(f"{gene} - After Scaling")
        axes[i, 1].text(
            0.05,
            0.95,
            f"Mean: {np.mean(scaled):.2f}\nStd: {np.std(scaled):.2f}",
            transform=axes[i, 1].transAxes,
            va="top",
        )

    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, "scaling_effect.png"), dpi=300)

    return fig


def regress_out(
    adata: AnnData,
    config: ScalingConfig,
    input_layer: Optional[str] = "log1p_norm",
    output_layer: Optional[str] = "regressed_out",
    force: bool = False,
    update_X: bool = True,
) -> Optional[AnnData]:
    """
    Regress out unwanted sources of variation from gene expression data.

    Args:
        adata: The AnnData object to process.
        config: A ScalingConfig object containing `vars_to_regress`.
        input_layer: Layer containing the data to be corrected (e.g., 'log1p_norm').
        output_layer: Layer to store the regressed-out data.
        force: If True, overwrite the `output_layer` if it already exists.
        update_X: If True, also update adata.X to the regressed layer.

    Returns:
        The modified AnnData object.

    Raises:
        ValueError: If keys are not found in adata.obs or layers are invalid.
    """
    # Parameter Validation
    if not config.vars_to_regress:
        log.info(
            "No variables specified in `config.vars_to_regress`. Skipping regression."
        )
        if input_layer != output_layer:
            adata.layers[output_layer] = adata.layers[input_layer].copy()
        if update_X:
            adata.X = adata.layers[output_layer].copy()
        return adata

    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    missing_keys = [key for key in config.vars_to_regress if key not in adata.obs]
    if missing_keys:
        raise ValueError(f"Keys not found in adata.obs: {', '.join(missing_keys)}")

    if input_layer not in adata.layers:
        raise ValueError(f"Input layer '{input_layer}' not found in adata.layers.")

    log.info(
        f"Regressing out: {', '.join(config.vars_to_regress)} from layer '{input_layer}'"
    )

    # Use a temporary object for the regression
    temp_adata = AnnData(X=adata.layers[input_layer].copy(), obs=adata.obs.copy())

    try:
        sc.pp.regress_out(temp_adata, keys=config.vars_to_regress)
    except Exception as e:
        log.error(f"Regression failed: {e}")
        raise RuntimeError("Failed to regress out variables.")

    adata.layers[output_layer] = temp_adata.X.copy()
    log.info(f"Regression complete. Results stored in 'adata.layers[{output_layer}]'.")

    if update_X:
        adata.X = adata.layers[output_layer].copy()

    # Store info to .uns for traceability
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["regress_out"] = {
        "vars_to_regress": config.vars_to_regress,
        "input_layer": input_layer,
        "output_layer": output_layer,
        "params": config.__dict__,
    }

    return adata
