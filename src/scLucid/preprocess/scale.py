"""
Data scaling and regression functions for single-cell RNA-seq data.

This module provides flexible, config-driven functions for scaling and
regressing out covariates, ensuring consistency with the scLucid workflow.
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from .config import ScalingConfig, apply_config_overrides
from .utils import validate_matrix_input

log = logging.getLogger(__name__)

__all__ = ["scale_data", "regress_out", "plot_scaling_effect"]


# --- Helper functions for different scaling methods ---
def _robust_scale(X: np.ndarray, max_value: Optional[float]) -> np.ndarray:
    """Robustly scales a dense matrix X."""
    gene_medians = np.median(X, axis=0)
    # MAD calculation corrected for robustness
    gene_mads = np.median(np.abs(X - gene_medians), axis=0)
    gene_mads[gene_mads == 0] = 1e-8  # Avoid division by zero

    X_scaled = (X - gene_medians) / gene_mads

    if max_value is not None:
        X_scaled = np.clip(X_scaled, -max_value, max_value)
    return X_scaled


def _minmax_scale(X: np.ndarray) -> np.ndarray:
    """Scales a dense matrix X to a [0, 1] range."""
    gene_mins = np.min(X, axis=0)
    gene_ranges = np.max(X, axis=0) - gene_mins
    gene_ranges[gene_ranges == 0] = 1e-8  # Avoid division by zero

    return (X - gene_mins) / gene_ranges


def _robust_scale_sparse(
    X: scipy.sparse.spmatrix, max_value: Optional[float]
) -> scipy.sparse.spmatrix:
    """
    Robust scaling for sparse matrices (memory efficient).

    Uses median and MAD computed in a sparse-aware manner.
    """
    import scipy.sparse as sp

    if not sp.issparse(X):
        return _robust_scale(X, max_value)  # Fallback to dense version

    # Convert to CSC format for efficient column operations
    X_csc = X.tocsc()

    n_genes = X_csc.shape[1]
    medians = np.zeros(n_genes)
    mads = np.zeros(n_genes)

    # Compute median and MAD per gene
    for i in range(n_genes):
        col_data = X_csc.getcol(i).data  # Only non-zero values

        if len(col_data) > 0:
            medians[i] = np.median(col_data)

            # MAD calculation
            deviations = np.abs(col_data - medians[i])
            mads[i] = np.median(deviations)
        else:
            medians[i] = 0
            mads[i] = 1e-8

    # Avoid division by zero
    mads[mads == 0] = 1e-8

    # Scale: (X - median) / MAD
    X_scaled = X_csc.copy()

    # Subtract median (broadcast-safe for sparse)
    for i in range(n_genes):
        col = X_scaled.getcol(i)
        col.data = (col.data - medians[i]) / mads[i]
        X_scaled[:, i] = col

    # Clip values if max_value specified
    if max_value is not None:
        X_scaled.data = np.clip(X_scaled.data, -max_value, max_value)

    return X_scaled.tocsr()  # Convert back to CSR


def _minmax_scale_sparse(X: scipy.sparse.spmatrix) -> scipy.sparse.spmatrix:
    """
    MinMax scaling for sparse matrices.
    """
    import scipy.sparse as sp

    if not sp.issparse(X):
        return _minmax_scale(X)

    X_csc = X.tocsc()
    n_genes = X_csc.shape[1]

    mins = np.zeros(n_genes)
    maxs = np.zeros(n_genes)

    for i in range(n_genes):
        col_data = X_csc.getcol(i).data
        if len(col_data) > 0:
            mins[i] = col_data.min()
            maxs[i] = col_data.max()

    ranges = maxs - mins
    ranges[ranges == 0] = 1e-8

    X_scaled = X_csc.copy()
    for i in range(n_genes):
        col = X_scaled.getcol(i)
        col.data = (col.data - mins[i]) / ranges[i]
        X_scaled[:, i] = col

    return X_scaled.tocsr()


# --- Main Functions ---
def scale_data(
    adata: AnnData,
    config: Optional[ScalingConfig] = None,
    output_layer: Optional[str] = "scaled",
    **kwargs,
) -> AnnData:
    """
    Scales gene expression data using the specified method.

    Operates on adata.X and modifies it in place. Assumes adata has been
    subsetted to the desired features (e.g., HVGs).

    Args:
        adata: AnnData object (will be modified in place).
        config: A ScalingConfig object. If None, a default config is used.
        **kwargs: Keyword arguments to override parameters in the config object
                  (e.g., `max_value=15`, `scale_method='robust'`).

    Returns:
        The modified AnnData object with scaled adata.X.
        If `output_layer` is not None, scaled values are also stored in
        `adata.layers[output_layer]` for downstream compatibility.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = ScalingConfig()
    else:
        active_config = apply_config_overrides(config, **kwargs)

    log.info(
        f"Scaling data in .X (shape: {adata.shape}) using '{active_config.scale_method}' method."
    )

    # --- 2. Validate input matrix ---
    validate_matrix_input(adata.X, name="adata.X", allow_negative=True)

    # --- 3. Apply the scaling method ---
    if active_config.regress_in_scale:
        vars_reg = active_config.vars_to_regress_in_scale or active_config.vars_to_regress or []
        vars_reg = [v for v in vars_reg if v]

        if vars_reg:
            missing = [k for k in vars_reg if k not in adata.obs.columns]
            if missing:
                log.warning(f"Vars to regress not found: {missing}")
                vars_reg = [k for k in vars_reg if k in adata.obs.columns]

            if vars_reg:
                # === IMPROVED: Use config setting ===
                input_layer = active_config.input_layer_for_regress

                if input_layer not in adata.layers:
                    raise ValueError(
                        f"Layer '{input_layer}' specified in config.input_layer_for_regress "
                        f"not found in adata.layers. Available: {list(adata.layers.keys())}"
                    )

                X_in = adata.layers[input_layer].copy()
                temp = AnnData(X=X_in, obs=adata.obs.copy(), var=adata.var.copy())

                log.info(f"Regressing {vars_reg} from layer '{input_layer}' before scaling")
                sc.pp.regress_out(temp, keys=vars_reg)

                adata.X = temp.X.copy()

                # Store metadata
                adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
                    "regress_inline"
                ] = {
                    "vars_to_regress": vars_reg,
                    "input_layer": input_layer,
                    "timestamp": pd.Timestamp.now().isoformat(),
                }
            else:
                log.info("No valid variables to regress in scale step; skipping inline regression.")
        else:
            log.info("regress_in_scale=True but no variables provided; skipping inline regression.")

    # --- 3. Scale the data ---
    original_X = adata.X.copy()

    try:
        if active_config.scale_method == "zscore":
            sc.pp.scale(adata, max_value=active_config.max_value, zero_center=True)

        elif active_config.scale_method == "robust":
            if scipy.sparse.issparse(adata.X):
                log.info("Using sparse-aware robust scaling")
                adata.X = _robust_scale_sparse(adata.X, max_value=active_config.max_value)
            else:
                adata.X = _robust_scale(adata.X, max_value=active_config.max_value)

        elif active_config.scale_method == "minmax":
            if scipy.sparse.issparse(adata.X):
                log.info("Using sparse-aware minmax scaling")
                adata.X = _minmax_scale_sparse(adata.X)
            else:
                adata.X = _minmax_scale(adata.X)

        else:
            raise ValueError(
                f"Unknown scale_method '{active_config.scale_method}'. "
                "Expected one of: zscore, robust, minmax."
            )

    except Exception as e:
        raise RuntimeError(f"[preprocess] Scaling failed: {e}. Check input data format.") from e

    # Backward compatibility: persist scaled matrix into a named layer.
    if output_layer:
        adata.layers[output_layer] = adata.X.copy()

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["scaling"] = {
        "params": active_config.to_dict(),  # Pydantic's built-in serialization
        "output_layer": output_layer,
    }

    log.info("Scaling complete. adata.X has been updated.")
    return adata


def regress_out(
    adata: AnnData,
    config: Optional[ScalingConfig] = None,
    input_layer: str = "normalized",
    output_layer: str = "regressed",
    **kwargs,
) -> AnnData:
    """
    Regress out unwanted sources of variation from gene expression data.

    Notes:
    - Expects input_layer to be log-normalized-like.
    - Variables must exist in adata.obs; typical covariates:
      ['total_counts', 'pct_counts_mt', 'S_score', 'G2M_score', 'cc_diff'].
    """
    if config is None:
        active_config = ScalingConfig()
    else:
        active_config = apply_config_overrides(config, **kwargs)

    vars_to_regress = list(active_config.vars_to_regress or [])
    if not vars_to_regress:
        log.info("No variables specified for regression. Skipping.")
        if input_layer != output_layer and input_layer in adata.layers:
            adata.layers[output_layer] = adata.layers[input_layer].copy()
        return adata

    missing_keys = [k for k in vars_to_regress if k not in adata.obs.columns]
    if missing_keys:
        log.warning(
            f"Variables to regress not found in adata.obs: {missing_keys}. "
            "Proceeding with available variables."
        )
        vars_to_regress = [k for k in vars_to_regress if k in adata.obs.columns]
        if not vars_to_regress:
            log.info("No valid variables left to regress. Skipping.")
            if input_layer != output_layer and input_layer in adata.layers:
                adata.layers[output_layer] = adata.layers[input_layer].copy()
            return adata

    if input_layer not in adata.layers:
        raise ValueError(f"Input layer '{input_layer}' not found in adata.layers.")

    log.info(f"Regressing out: {', '.join(vars_to_regress)} from layer '{input_layer}'")
    temp_adata = AnnData(
        X=adata.layers[input_layer].copy(), obs=adata.obs.copy(), var=adata.var.copy()
    )
    try:
        sc.pp.regress_out(temp_adata, keys=vars_to_regress)
    except Exception as e:
        log.error(f"regress_out failed: {e}")
        raise

    adata.layers[output_layer] = temp_adata.X.copy()
    log.info(f"Regression complete. Results stored in adata.layers['{output_layer}'].")

    # Metadata
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["regress"] = {
        "input_layer": input_layer,
        "output_layer": output_layer,
        "vars_to_regress": vars_to_regress,
        "scanpy_version": getattr(sc, "__version__", "unknown"),
    }
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
        data = original_data.toarray() if scipy.sparse.issparse(original_data) else original_data
        gene_vars = np.var(data, axis=0)
        top_idx = np.argsort(-gene_vars)[:n_genes]
        genes_to_plot = adata.var_names[top_idx].tolist()

    fig, axes = plt.subplots(len(genes_to_plot), 2, figsize=(12, 3 * len(genes_to_plot)))
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
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        figure_path = save_path / "scaling_effect.png"
        plt.savefig(figure_path, dpi=300)
        log.info(f"Saved scaling effect plot to {figure_path}")

    return fig
