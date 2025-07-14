"""
Data integration and batch effect correction methods.
"""

import os
from typing import List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

# ==============================================================================
# Low-Level Integration Wrappers
# ==============================================================================


def _integrate_harmony(
    adata: sc.AnnData,
    batch_key: str,
    basis: str = "X_pca",
    embedding_key: str = "X_harmony",
    **kwargs,
):
    """Internal wrapper for Harmony integration."""
    try:
        from scanpy.external.pp import harmony_integrate
    except ImportError:
        raise ImportError("Please install harmonypy: pip install harmonypy")

    if basis not in adata.obsm:
        raise ValueError(f"Basis '{basis}' not found in adata.obsm. Run PCA first.")

    # Harmony modifies the basis in-place and saves to a new key
    harmony_integrate(adata, key=batch_key, basis=basis, **kwargs)

    # Standardize output key
    adata.obsm[embedding_key] = adata.obsm["X_pca_harmony"].copy()
    if "X_pca_harmony" != embedding_key:
        del adata.obsm["X_pca_harmony"]

    return adata


def _integrate_scanorama(
    adata: sc.AnnData,
    batch_key: str,
    hvg: Optional[List[str]] = None,
    dims: int = 50,
    embedding_key: str = "X_scanorama",
    **kwargs,
):
    """Internal wrapper for Scanorama integration."""
    try:
        import scanorama
    except ImportError:
        raise ImportError("Please install Scanorama: pip install scanorama")

    # Split AnnData object by batch
    batches = adata.obs[batch_key].unique()
    adatas_list = [adata[adata.obs[batch_key] == b].copy() for b in batches]

    # Use HVGs if provided
    if hvg is not None:
        print(f"Subsetting to {len(hvg)} highly variable genes for Scanorama.")
        for ad in adatas_list:
            # Ensure we don't error on missing genes, just use what's available
            genes_in_batch = ad.var_names.intersection(hvg)
            ad._inplace_subset_var(genes_in_batch)

    # Run Scanorama integration
    scanorama.integrate_scanpy(adatas_list, dimred=dims, **kwargs)

    # Concatenate the results back into the original adata object
    # The order is preserved from the initial split
    integrated_embedding = np.concatenate(
        [ad.obsm["X_scanorama"] for ad in adatas_list]
    )
    adata.obsm[embedding_key] = integrated_embedding

    return adata


def _integrate_scvi(
    adata: sc.AnnData,
    batch_key: str,
    layer: Optional[str] = "counts",
    n_latent: int = 30,
    embedding_key: str = "X_scVI",
    **kwargs,
):
    """Internal wrapper for scVI integration."""
    try:
        import scvi
    except ImportError:
        raise ImportError("Please install scvi-tools: pip install scvi-tools")

    # scVI requires setup on the AnnData object
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)

    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(**kwargs)  # Pass other training args like max_epochs

    adata.obsm[embedding_key] = model.get_latent_representation()
    return adata


# ==============================================================================
# High-Level Batch Correction Workflow Function
# ==============================================================================


def batch_correction(
    adata: sc.AnnData,
    batch_key: str,
    method: Literal["harmony", "scanorama", "scvi"] = "harmony",
    use_rep: str = "X_pca",
    hvg_key: Optional[str] = "highly_variable",
    layer: Optional[str] = "counts",
    plot: bool = True,
    save_dir: Optional[str] = None,
    **kwargs,
) -> sc.AnnData:
    """
    Performs batch correction using a specified integration method.

    This function serves as a high-level wrapper that calls the appropriate
    integration tool and visualizes the result.

    Args:
        adata: AnnData object. Must contain a dimensionality reduction (e.g., PCA).
        batch_key: Key in `adata.obs` that denotes the batch.
        method: Integration method to use.
        use_rep: The representation to use as input for Harmony. Defaults to 'X_pca'.
        hvg_key: Key in `adata.var` specifying HVGs, used by Scanorama.
        layer: Layer containing raw counts, used by scVI.
        plot: If True, plots UMAPs before and after correction for comparison.
        save_dir: Directory to save plots.
        **kwargs: Additional arguments to pass to the integration method's train/run function.

    Returns:
        The AnnData object with a new batch-corrected embedding in `adata.obsm`.
    """
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in `adata.obs`.")

    output_key = f"X_{method}"

    # --- Plotting: Before Correction ---
    if plot:
        print("Generating UMAP before batch correction...")
        adata_before = adata.copy()
        sc.pp.neighbors(adata_before, use_rep=use_rep)
        sc.tl.umap(adata_before)

    # --- Integration ---
    print(f"Performing batch correction using '{method}'...")
    if method == "harmony":
        _integrate_harmony(
            adata,
            batch_key=batch_key,
            basis=use_rep,
            embedding_key=output_key,
            **kwargs,
        )

    elif method == "scanorama":
        hvg_list = adata.var_names[adata.var[hvg_key]] if hvg_key in adata.var else None
        _integrate_scanorama(
            adata,
            batch_key=batch_key,
            hvg=hvg_list,
            dims=adata.obsm[use_rep].shape[1],
            embedding_key=output_key,
            **kwargs,
        )

    elif method == "scvi":
        _integrate_scvi(
            adata,
            batch_key=batch_key,
            layer=layer,
            n_latent=adata.obsm[use_rep].shape[1],
            embedding_key=output_key,
            **kwargs,
        )

    else:
        raise ValueError(
            f"Unknown integration method: '{method}'. Choose from 'harmony', 'scanorama', 'scvi'."
        )

    print(
        f"Integration complete. Corrected embedding stored in `adata.obsm['{output_key}']`."
    )

    # --- Plotting: After Correction ---
    if plot:
        print("Generating UMAP after batch correction...")
        # Compute neighbors and UMAP on the new, corrected embedding
        sc.pp.neighbors(adata, use_rep=output_key)
        sc.tl.umap(adata)

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        fig.suptitle("Batch Correction Comparison", fontsize=16)

        sc.pl.umap(
            adata_before,
            color=batch_key,
            ax=axes[0],
            show=False,
            title=f"Before Correction ({use_rep})",
        )
        sc.pl.umap(
            adata,
            color=batch_key,
            ax=axes[1],
            show=False,
            title=f"After {method.capitalize()} Correction",
        )

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(
                os.path.join(save_dir, f"batch_correction_{method}.png"), dpi=300
            )
        plt.show()
        plt.close(fig)

    return adata
