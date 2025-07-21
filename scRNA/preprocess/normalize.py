"""
Normalization functions for single-cell RNA-seq data.
"""

import sys

sys.path.append("..")

import os
from typing import List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns

from ..utils import use_layer_as_X


def normalize_data(
    adata: sc.AnnData,
    method: Literal["standard", "scran", "pearson_residuals"] = "standard",
    layer: Optional[str] = "counts",
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
    print(f"Begin data normalization using '{method}' method.")

    # Validate parameters
    if target_sum is not None and target_sum <= 0:
        raise ValueError("target_sum must be a positive number.")
    if not 0 < max_fraction < 1:
        raise ValueError("max_fraction must be between 0 and 1 (exclusive).")

    # Ensure a raw counts layer exists for comparison plotting
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        print("Saved current adata.X to adata.layers['counts'] for reference.")

    input_data_source = adata.layers.get(layer, adata.X)

    # --- Normalization ---
    if method == "standard":
        # Use a temporary AnnData object for normalization to avoid inplace modification issues
        temp_adata = sc.AnnData(input_data_source.copy())
        sc.pp.normalize_total(
            temp_adata,
            target_sum=target_sum,
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            inplace=True,
        )
        X_norm = temp_adata.X.copy()

    elif method == "scran":
        print(
            "Warning: The 'scran' method requires a separate installation of R, rpy2, and Bioconductor's scran package."
        )
        print(
            "This is a complex dependency. For now, this is a placeholder and will proceed with 'standard' normalization."
        )
        # Fallback to standard normalization
        temp_adata = sc.AnnData(input_data_source.copy())
        sc.pp.normalize_total(temp_adata, target_sum=target_sum, inplace=True)
        X_norm = temp_adata.X.copy()
        # try:
        #    import rpy2
        #    import rpy2.robjects as ro
        #    from rpy2.robjects.packages import importr
        #    from rpy2.robjects import pandas2ri
        #    pandas2ri.activate()

        #    # Use scran from R via rpy2
        #    importr('scran')
        #    # Simplified example - in practice this would use rpy2 to call R's scran
        #    print("Using scran for normalization...")
        #    # This is a placeholder - actual implementation would involve R code
        #    X_norm = sc.pp.normalize_total(adata, layer=layer, inplace=False)['X']
        #    print("Note: This is currently a placeholder. Full scran implementation requires R integration.")
        # except ImportError:
        #    raise ValueError("scran method requires rpy2. Please install with 'pip install rpy2'")

    elif method == "pearson_residuals":
        # Pearson residuals normalization is an advanced method that replaces standard normalization and scaling.
        print("Applying Pearson residuals normalization...")
        X_norm = sc.experimental.pp.normalize_pearson_residuals(
            sc.AnnData(input_data_source.copy()), inplace=False
        )["X"]
        # Log transform is not needed for Pearson residuals
        log_transform = False

    # --- Store Results and Log Transform ---
    adata.layers["normalized"] = X_norm.copy()

    if log_transform:
        print("Applying log1p transformation.")
        # sc.pp.log1p is robust for both sparse and dense matrices
        adata.layers[output_layer] = sc.pp.log1p(X_norm, copy=True)
    else:
        adata.layers[output_layer] = X_norm.copy()

    print(
        f"Normalization complete. Final data stored in adata.layers['{output_layer}']"
    )

    # --- Plotting ---
    if plot:
        # Visualize the distributions
        print("Generating comparison plots...")
        plt.rcParams.update(
            {
                "figure.facecolor": "white",
                "axes.facecolor": "white",
                "savefig.facecolor": "white",
                "text.color": "black",
                "axes.labelcolor": "black",
                "axes.edgecolor": "black",
                "xtick.color": "black",
                "ytick.color": "black",
            }
        )
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
        fig.suptitle(
            f"Data Distributions Before and After {method.capitalize()} Normalization",
            fontsize=16,
            color="black",
        )

        # Before
        sns.histplot(
            np.array(input_data_source.sum(axis=1)).flatten(),
            bins=100,
            kde=True,
            ax=axes[0],
            color="navy",
        )
        axes[0].set_title("Before Normalization")
        axes[0].set_xlabel("Total Counts per Cell")

        # After
        final_data = adata.layers[output_layer]
        sns.histplot(
            np.array(final_data.sum(axis=1)).flatten(),
            bins=100,
            kde=True,
            ax=axes[1],
            color="crimson",
        )
        title_suffix = " (Log-Transformed)" if log_transform else ""
        axes[1].set_title(f"After Normalization{title_suffix}")
        axes[1].set_xlabel("Sum of Normalized Values per Cell")

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f"normalization_{method}.png"), dpi=300)
        plt.show()
        plt.close(fig)

    return adata


def regress_out(
    adata: sc.AnnData,
    keys: List[str],
    layer: Optional[str] = "log1p_norm",
    n_jobs: Optional[int] = None,
    output_layer: str = "regressed_out",
) -> sc.AnnData:
    """
    Regress out unwanted sources of variation from a specified layer.

    Args:
        adata: AnnData object
        keys: Variables to regress out (must be in adata.obs)
        layer: Layer to use as input for regression.
        n_jobs: Number of parallel jobs to use
        output_layer: Layer to store the regressed-out data.

    Returns:
        AnnData with regressed out variables
    """
    # Check if keys exist in adata.obs
    missing_keys = [key for key in keys if key not in adata.obs]
    if missing_keys:
        raise ValueError(f"Keys not found in adata.obs: {', '.join(missing_keys)}")

    print(f"Regressing out: {', '.join(keys)} from layer '{layer}'")

    # Use the context manager to safely handle layers
    with use_layer_as_X(adata, layer):
        # The regress_out function modifies adata.X in place
        sc.pp.regress_out(adata, keys=keys, n_jobs=n_jobs)

        # Store the result from the modified adata.X into the output layer
        adata.layers[output_layer] = adata.X.copy()

    print(f"Regressed data stored in adata.layers['{output_layer}']")

    return adata
