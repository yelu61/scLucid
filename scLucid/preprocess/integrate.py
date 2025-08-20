"""
Data integration and batch effect correction methods.

This module provides functions to correct batch effects and integrate data from
multiple experiments or conditions in single-cell RNA-seq data. It implements
wrappers around popular integration methods including Harmony, Scanorama, and
scVI.
"""

import logging
import os
from typing import Dict, List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = ["batch_correction", "evaluate_integration"]

# ==============================================================================
# Low-Level Integration Wrappers
# ==============================================================================


def _integrate_harmony(
    adata: AnnData,
    batch_key: str,
    basis: str = "X_pca",
    embedding_key: str = "X_harmony",
    max_iter_harmony: int = 20,
    theta: float = 2.0,
    lambda_val: float = 1.0,
    sigma: float = 0.1,
    random_state: int = 42,
    plot_convergence: bool = False,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Internal wrapper for Harmony integration.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key for batch/sample annotation in adata.obs.
    basis : str, optional (default: "X_pca")
        The embedding in adata.obsm to use as input for Harmony.
    embedding_key : str, optional (default: "X_harmony")
        Key under which to store the integrated embedding in adata.obsm.
    max_iter_harmony : int, optional (default: 20)
        Maximum number of iterations for Harmony optimization.
    theta : float, optional (default: 2.0)
        Diversity clustering penalty parameter. Larger values result in more
        diverse clusters. Default: 2.
    lambda_val : float, optional (default: 1.0)
        Ridge regression penalty parameter. Default: 1.
    sigma : float, optional (default: 0.1)
        Soft cluster scatter. Default: 0.1.
    random_state : int, optional (default: 42)
        Random seed for reproducibility.
    plot_convergence : bool, optional (default: False)
        Whether to plot the convergence of the Harmony objective function.
    copy : bool, optional (default: False)
        If True, return a copy of the AnnData object.
    **kwargs
        Additional arguments to scanpy.external.pp.harmony_integrate().

    Returns
    -------
    adata : AnnData
        Annotated data matrix with integrated embedding in .obsm[embedding_key].

    Notes
    -----
    This function requires the 'harmonypy' package to be installed:
    pip install harmonypy

    The input embedding (specified by `basis`) should typically be PCA,
    but can be any low-dimensional representation.
    """
    try:
        from scanpy.external.pp import harmony_integrate
    except ImportError:
        try:
            # Try to import harmonypy directly as a fallback
            import harmonypy

            log.error(
                "Found harmonypy, but failed to import through scanpy. "
                "Please update scanpy: pip install --upgrade scanpy"
            )
            raise ImportError(
                "Found harmonypy, but failed to import through scanpy. "
                "Please update scanpy: pip install --upgrade scanpy"
            )
        except ImportError:
            log.error("Please install harmonypy: pip install harmonypy")
            raise ImportError("Please install harmonypy: pip install harmonypy")

    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

    # Check number of batches
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning(
            f"Only {n_batches} batch found. Harmony works best with multiple batches."
        )

    # Check if basis exists
    if basis not in adata.obsm:
        log.error(f"Basis '{basis}' not found in adata.obsm")
        raise ValueError(
            f"Basis '{basis}' not found in adata.obsm. "
            f"Run dimensionality reduction first."
        )

    # Create working copy if needed
    if copy:
        adata = adata.copy()

    log.info(f"Running Harmony integration using '{batch_key}' as batch key")
    log.info(f"Input dimensions: {adata.obsm[basis].shape}")
    log.info(
        f"Parameters: theta={theta}, lambda={lambda_val}, sigma={sigma}, "
        f"max_iter={max_iter_harmony}"
    )

    # Build full parameter set for harmony_integrate
    full_params = {
        "theta": theta,
        "lamb": lambda_val, 
        "sigma": sigma,
        "max_iter_harmony": max_iter_harmony,
        "random_state": random_state,
        "plot_convergence": plot_convergence,
    }

    # Add user-defined parameters
    full_params.update(kwargs)

    # Run Harmony integration
    try:
        harmony_integrate(
            adata,
            key=batch_key,
            basis=basis,
            adjusted_basis=f"{basis}_harmony",
            **full_params,
        )
    except Exception as e:
        log.error(f"Harmony integration failed: {str(e)}")
        raise RuntimeError(f"Harmony integration failed: {str(e)}")

    # Check if Harmony output was successfully generated
    expected_key = f"{basis}_harmony"
    if expected_key not in adata.obsm:
        log.error(f"Harmony failed to generate output in adata.obsm['{expected_key}']")
        raise RuntimeError(
            f"Harmony failed to generate output in adata.obsm['{expected_key}']"
        )

    # Copy results to user-specified key
    adata.obsm[embedding_key] = adata.obsm[expected_key].copy()

    # If user specified a different key, and it's not the default harmony key, delete default key
    if expected_key != embedding_key:
        del adata.obsm[expected_key]

    # Store integration metadata
    if "integration" not in adata.uns:
        adata.uns["integration"] = {}

    adata.uns["integration"]["harmony"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "params": {
            "theta": theta,
            "lambda": lambda_val,
            "sigma": sigma,
            "max_iter_harmony": max_iter_harmony,
            "random_state": random_state,
        },
        "input_dims": adata.obsm[basis].shape[1],
        "output_dims": adata.obsm[embedding_key].shape[1],
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    log.info(
        f"Harmony integration complete. Results stored in adata.obsm['{embedding_key}']"
    )
    log.info(f"Output dimensions: {adata.obsm[embedding_key].shape}")

    return adata


def _integrate_scanorama(
    adata: AnnData,
    batch_key: str,
    hvg: Optional[List[str]] = None,
    dims: int = 50,
    embedding_key: str = "X_scanorama",
    correct_expression: bool = False,
    return_corrected_expression: bool = False,
    layer: Optional[str] = None,
    knn: int = 20,
    sigma: float = 15,
    approx: bool = True,
    alpha: float = 0.1,
    batch_size: int = 5000,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Internal wrapper for Scanorama integration.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key for batch/sample annotation in adata.obs.
    hvg : List[str], optional (default: None)
        List of highly variable gene names to use for integration.
        If None, all genes will be used.
    dims : int, optional (default: 50)
        Dimensionality of the integrated embedding.
    embedding_key : str, optional (default: "X_scanorama")
        Key under which to store the integrated embedding in adata.obsm.
    correct_expression : bool, optional (default: False)
        Whether to correct the expression values in addition to computing embeddings.
    return_corrected_expression : bool, optional (default: False)
        Whether to return batch-corrected expression values.
        If True, corrected values will be stored in adata.layers['scanorama_corrected'].
    layer : str, optional (default: None)
        If provided, use this layer for integration instead of adata.X.
    knn : int, optional (default: 20)
        Number of nearest neighbors to use for matching cells across datasets.
    sigma : float, optional (default: 15)
        Correction smoothing parameter on Gaussian kernel.
    approx : bool, optional (default: True)
        Use approximate nearest neighbors with ANNOY for scalability.
    alpha : float, optional (default: 0.1)
        Alignment score minimum cutoff.
    batch_size : int, optional (default: 5000)
        The batch size used for the matching algorithm.
    copy : bool, optional (default: False)
        If True, return a copy of the AnnData object.
    **kwargs
        Additional arguments to scanorama.integrate_scanpy().

    Returns
    -------
    adata : AnnData
        Annotated data matrix with integrated embedding in .obsm[embedding_key]
        and optionally batch-corrected expression in .layers['scanorama_corrected'].
    """
    try:
        import scanorama
    except ImportError:
        log.error("Please install Scanorama: pip install scanorama")
        raise ImportError("Please install Scanorama: pip install scanorama")

    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

    # Check if batch_key has more than one unique value
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    if n_batches < 2:
        log.error(
            f"Found only {n_batches} batch. Scanorama requires at least 2 batches"
        )
        raise ValueError(
            f"Found only {n_batches} batch. Scanorama requires at least 2 batches for integration."
        )

    # Create a copy if requested
    if copy:
        adata = adata.copy()

    log.info(f"Running Scanorama integration on {n_batches} batches")
    log.info(f"Parameters: knn={knn}, sigma={sigma}, approx={approx}, alpha={alpha}")

    # Save original data
    original_obs = adata.obs.copy()
    original_var = adata.var.copy()

    # Split AnnData object by batch
    adatas_list = []
    genes_union = set()

    log.info("Splitting data by batch and preparing for integration")
    for i, b in enumerate(batches):
        batch_cells = np.sum(adata.obs[batch_key] == b)
        log.info(f"Preparing batch {i + 1}/{n_batches}: '{b}' with {batch_cells} cells")

        batch_data = adata[adata.obs[batch_key] == b].copy()

        # If a layer is specified, use it instead of adata.X
        if layer is not None:
            if layer not in batch_data.layers:
                log.error(f"Layer '{layer}' not found in adata")
                raise ValueError(f"Layer '{layer}' not found in adata.")
            batch_data.X = batch_data.layers[layer].copy()

        adatas_list.append(batch_data)
        genes_union.update(batch_data.var_names)

    # Use HVGs if provided
    if hvg is not None:
        if not isinstance(hvg, list) and not isinstance(hvg, np.ndarray):
            log.error("hvg must be a list or numpy array of gene names")
            raise TypeError("hvg must be a list or numpy array of gene names")

        log.info(f"Subsetting to {len(hvg)} highly variable genes for Scanorama")
        genes_to_use = list(set(hvg).intersection(genes_union))

        if len(genes_to_use) < len(hvg):
            log.warning(
                f"Only {len(genes_to_use)} out of {len(hvg)} HVGs found in the data"
            )

        if len(genes_to_use) == 0:
            log.error("No genes from the provided hvg list were found in the data")
            raise ValueError(
                "No genes from the provided hvg list were found in the data."
            )
    else:
        # If no HVGs provided, use all genes common to all batches
        common_genes = set.intersection(*[set(ad.var_names) for ad in adatas_list])
        genes_to_use = list(common_genes)
        log.info(
            f"Using {len(genes_to_use)} genes common to all batches for integration"
        )

    # Make sure all batches have the same genes
    for i, ad in enumerate(adatas_list):
        # Get genes to use in this batch
        genes_in_batch = [g for g in genes_to_use if g in ad.var_names]

        if len(genes_in_batch) < len(genes_to_use):
            log.warning(
                f"Batch {i} contains {len(genes_in_batch)}/{len(genes_to_use)} integration genes"
            )

        # Subset the data to the genes in this batch
        adatas_list[i] = ad[:, genes_in_batch].copy()

    # Prepare status message
    if correct_expression:
        log.info("Computing integrated embedding and batch-corrected expression")
    else:
        log.info("Computing integrated embedding only")

    # Run Scanorama integration
    try:
        scanorama.integrate_scanpy(
            adatas_list,
            dimred=dims,
            knn=knn,
            sigma=sigma,
            approx=approx,
            alpha=alpha,
            batch_size=batch_size,
            **kwargs,
        )
    except Exception as e:
        log.error(f"Scanorama integration failed: {str(e)}")
        raise RuntimeError(f"Scanorama integration failed: {str(e)}")

    # Reconstruct the integrated embedding
    # Note: Scanorama returns a list of embeddings, one for each batch
    log.info("Assembling integrated embedding from individual batch results")
    cell_indices = np.concatenate(
        [np.where(adata.obs[batch_key] == b)[0] for b in batches]
    )
    integrated_embedding = np.zeros((adata.shape[0], dims))

    counter = 0
    for i, b in enumerate(batches):
        batch_size = np.sum(adata.obs[batch_key] == b)
        if "X_scanorama" not in adatas_list[i].obsm:
            log.error(f"Scanorama failed to generate embedding for batch {i}")
            raise RuntimeError(f"Scanorama failed to generate embedding for batch {i}")

        integrated_embedding[cell_indices[counter : counter + batch_size]] = (
            adatas_list[i].obsm["X_scanorama"]
        )
        counter += batch_size

    adata.obsm[embedding_key] = integrated_embedding
    log.info(f"Integrated embedding stored in adata.obsm['{embedding_key}']")

    # If requested, also return batch-corrected expression
    if return_corrected_expression:
        if not hasattr(adatas_list[0], "X_scanorama"):
            log.warning(
                "Batch-corrected expression not available. Set correct_expression=True to enable this."
            )
        else:
            log.info("Assembling batch-corrected expression matrix")
            corrected_exp = np.zeros(adata.shape)
            genes_idx = {g: i for i, g in enumerate(adata.var_names)}

            counter = 0
            for i, b in enumerate(batches):
                batch_size = np.sum(adata.obs[batch_key] == b)
                # Mapping the corrected expression to the original gene order
                for j, gene in enumerate(adatas_list[i].var_names):
                    if gene in genes_idx:
                        corrected_exp[
                            cell_indices[counter : counter + batch_size],
                            genes_idx[gene],
                        ] = adatas_list[i].X_scanorama[:, j]
                counter += batch_size

            # Save the corrected expression in a new layer
            adata.layers["scanorama_corrected"] = corrected_exp
            log.info(
                "Batch-corrected expression stored in adata.layers['scanorama_corrected']"
            )

    # Make sure to return the original adata object with the same obs and var
    for col in original_obs.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = original_obs[col].values

    for col in original_var.columns:
        if col not in adata.var.columns:
            adata.var[col] = original_var[col].values

    # Store integration metadata
    if "integration" not in adata.uns:
        adata.uns["integration"] = {}

    adata.uns["integration"]["scanorama"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "params": {
            "knn": knn,
            "sigma": sigma,
            "approx": approx,
            "alpha": alpha,
            "batch_size": batch_size,
            "n_hvgs": len(genes_to_use),
            "layer": layer,
        },
        "output_dims": dims,
        "corrected_expression": return_corrected_expression,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    return adata


def _integrate_scvi(
    adata: AnnData,
    batch_key: str,
    layer: Optional[str] = "counts",
    n_layers: int = 2,
    n_latent: int = 30,
    batch_size: int = 256,
    max_epochs: int = 500,
    embedding_key: str = "X_scVI",
    gene_likelihood: str = "nb",
    #use_gpu: Optional[bool] = None,
    save_model: bool = False,
    model_path: Optional[str] = None,
    plan_kwargs: Optional[dict] = None,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Internal wrapper for scVI integration.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key for batch/sample annotation in adata.obs.
    layer : str, optional (default: "counts")
        Layer containing raw count data. If None, uses .X.
    n_layers : int, optional (default: 2)
        Number of hidden layers used for encoder and decoder NNs.
    n_latent : int, optional (default: 30)
        Dimensionality of the latent space.
    batch_size : int, optional (default: 256)
        Minibatch size for training.
    max_epochs : int, optional (default: 500)
        Maximum number of training epochs.
    embedding_key : str, optional (default: "X_scVI")
        Key under which to store the latent representation in adata.obsm.
    gene_likelihood : str, optional (default: "nb")
        One of 'nb' (negative binomial), 'zinb' (zero-inflated negative binomial), or 'poisson'.
    use_gpu : bool, optional (default: None)
        If True, use GPU if available. If None, automatically detect.
    save_model : bool, optional (default: False)
        If True, save the scVI model.
    model_path : str, optional (default: None)
        Path to save the model. Required if save_model is True.
    plan_kwargs : dict, optional (default: None)
        Keyword args for TrainingPlan.
    copy : bool, optional (default: False)
        If True, return a copy of the AnnData object.
    **kwargs
        Additional arguments to scvi.model.SCVI.train().

    Returns
    -------
    adata : AnnData
        Annotated data matrix with latent representation in .obsm[embedding_key].
    """
    try:
        import scvi
    except ImportError:
        log.error("Please install scvi-tools: pip install scvi-tools")
        raise ImportError("Please install scvi-tools: pip install scvi-tools")

    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

    # Check if the batch_key has more than one unique value
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning(
            f"Only {n_batches} batch found. scVI works best with multiple batches."
        )

    # Create a copy if requested
    if copy:
        adata = adata.copy()

    # Check if layer exists if specified
    if layer is not None and layer not in adata.layers:
        log.error(f"Layer '{layer}' not found in adata.layers")
        raise ValueError(f"Layer '{layer}' not found in adata.layers")

    log.info(f"Setting up scVI with {n_layers} layers and {n_latent} latent dimensions")
    log.info(
        f"Parameters: batch_size={batch_size}, max_epochs={max_epochs}, "
        f"gene_likelihood={gene_likelihood}"
    )

    # scVI requires setup on the AnnData object
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)

    # Prepare the model
    if plan_kwargs is None:
        plan_kwargs = {}

    # Create the scVI model
    model = scvi.model.SCVI(
        adata, n_layers=n_layers, n_latent=n_latent, gene_likelihood=gene_likelihood
    )

    # Train the model
    log.info(
        f"Training scVI model with batch_size={batch_size}, max_epochs={max_epochs}"
    )
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
        plan_kwargs=plan_kwargs,
        **kwargs,
    )  # pass any additional arguments to train()

    log.info("Extracting latent representation")
    # Extract the latent representation
    adata.obsm[embedding_key] = model.get_latent_representation()

    # Save the model if requested
    if save_model:
        if model_path is None:
            log.error("model_path must be provided when save_model=True")
            raise ValueError("model_path must be provided when save_model=True")

        os.makedirs(os.path.dirname(os.path.abspath(model_path)), exist_ok=True)
        log.info(f"Saving model to {model_path}")
        model.save(model_path)

    # Store integration metadata
    if "integration" not in adata.uns:
        adata.uns["integration"] = {}

    adata.uns["integration"]["scvi"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "params": {
            "n_layers": n_layers,
            "n_latent": n_latent,
            "batch_size": batch_size,
            "max_epochs": max_epochs,
            "gene_likelihood": gene_likelihood,
            "layer": layer,
        },
        "output_dims": n_latent,
        "model_saved": save_model,
        "model_path": model_path if save_model else None,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    log.info(
        f"scVI integration complete. Latent representation stored in adata.obsm['{embedding_key}']"
    )
    return adata


def _integrate_bbknn(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
    neighbors_within_batch: int = 3,
    n_pcs: Optional[int] = None,
    trim: Optional[int] = None,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Internal wrapper for BBKNN (Batch Balanced KNN) integration.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key for batch/sample annotation in adata.obs.
    use_rep : str, optional (default: "X_pca")
        The representation to use for computing neighbors.
    neighbors_within_batch : int, optional (default: 3)
        How many top neighbors to report for each batch; total number of neighbors
        will be this times the number of batches.
    n_pcs : int, optional (default: None)
        Number of PCs to use. If None, use all PCs in use_rep.
    trim : int, optional (default: None)
        Trim the neighbors to this number (i.e., report a total of `trim` neighbors).
    copy : bool, optional (default: False)
        If True, return a copy of the AnnData object.
    **kwargs
        Additional arguments to scanpy.external.pp.bbknn().

    Returns
    -------
    adata : AnnData
        Annotated data matrix with the batch-balanced nearest neighbor graph computed.
    """
    try:
        from scanpy.external.pp import bbknn
    except ImportError:
        log.error("Please install bbknn: pip install bbknn")
        raise ImportError("Please install bbknn: pip install bbknn")

    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

    # Check if use_rep exists
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"use_rep '{use_rep}' not found in adata.obsm")

    # Check number of batches
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning(
            f"Only {n_batches} batch found. BBKNN works best with multiple batches."
        )

    # Create a copy if requested
    if copy:
        adata = adata.copy()

    # Determine n_pcs if not provided
    if n_pcs is None:
        n_pcs = adata.obsm[use_rep].shape[1]

    log.info(f"Running BBKNN integration with '{batch_key}' as batch key")
    log.info(
        f"Parameters: neighbors_within_batch={neighbors_within_batch}, n_pcs={n_pcs}"
    )

    # Run BBKNN
    try:
        bbknn(
            adata,
            batch_key=batch_key,
            use_rep=use_rep,
            neighbors_within_batch=neighbors_within_batch,
            n_pcs=n_pcs,
            trim=trim,
            **kwargs,
        )
    except Exception as e:
        log.error(f"BBKNN integration failed: {str(e)}")
        raise RuntimeError(f"BBKNN integration failed: {str(e)}")

    # Store integration metadata
    if "integration" not in adata.uns:
        adata.uns["integration"] = {}

    adata.uns["integration"]["bbknn"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "params": {
            "neighbors_within_batch": neighbors_within_batch,
            "n_pcs": n_pcs,
            "trim": trim,
            "use_rep": use_rep,
        },
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    log.info("BBKNN integration complete. Neighborhood graph stored in adata.obsp")
    return adata


def _integrate_combat(
    adata: AnnData,
    batch_key: str,
    layer: Optional[str] = None,
    covariates: Optional[List[str]] = None,
    inplace: bool = True,
    output_layer: str = "combat_corrected",
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Internal wrapper for ComBat batch correction.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key for batch/sample annotation in adata.obs.
    layer : str, optional (default: None)
        If provided, use this layer for correction instead of adata.X.
    covariates : List[str], optional (default: None)
        List of covariates to preserve during batch correction (must be columns in adata.obs).
    inplace : bool, optional (default: True)
        Whether to update adata.X with the corrected values.
    output_layer : str, optional (default: "combat_corrected")
        If provided, store the corrected data in this layer.
    copy : bool, optional (default: False)
        If True, return a copy of the AnnData object.
    **kwargs
        Additional arguments to scanpy.pp.combat().

    Returns
    -------
    adata : AnnData
        Annotated data matrix with batch-corrected data.
    """
    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")

    # Check covariates
    if covariates is not None:
        missing_covariates = [c for c in covariates if c not in adata.obs.columns]
        if missing_covariates:
            log.error(f"Covariates not found in adata.obs: {missing_covariates}")
            raise ValueError(f"Covariates not found in adata.obs: {missing_covariates}")

    # Check layer
    if layer is not None and layer not in adata.layers:
        log.error(f"Layer '{layer}' not found in adata.layers")
        raise ValueError(f"layer '{layer}' not found in adata.layers")

    # Create a copy if requested
    if copy:
        adata = adata.copy()

    # Check number of batches
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning(
            f"Only {n_batches} batch found. ComBat requires at least 2 batches."
        )
        return adata

    log.info(f"Running ComBat batch correction with '{batch_key}' as batch key")

    # Prepare data for ComBat
    if layer is not None:
        X_combat = adata.layers[layer].copy()
    else:
        X_combat = adata.X.copy()

    # Convert to dense if sparse
    import scipy.sparse

    if scipy.sparse.issparse(X_combat):
        log.info("Converting sparse matrix to dense for ComBat")
        X_combat = X_combat.toarray()

    # Run ComBat
    try:
        # If scanpy's implementation is available
        try:
            import scanpy as sc

            log.info("Using scanpy.pp.combat for batch correction")

            # Create temporary AnnData
            import anndata

            temp_adata = anndata.AnnData(X=X_combat, obs=adata.obs.copy())

            # Run ComBat
            sc.pp.combat(
                temp_adata, key=batch_key, covariates=covariates, inplace=True, **kwargs
            )

            # Get corrected data
            X_corrected = temp_adata.X

        # If scanpy's implementation is not available, try using combat directly
        except (ImportError, AttributeError):
            log.info("Using standalone combat implementation")

            try:
                from combat.pycombat import pycombat
            except ImportError:
                try:
                    import combat as cb

                    pycombat = cb.pycombat
                except ImportError:
                    log.error("Please install combat: pip install combat")
                    raise ImportError("Please install combat: pip install combat")

            # Prepare batch vector
            batch_vec = adata.obs[batch_key].values

            # Prepare covariates matrix
            covariates_mat = None
            if covariates is not None:
                covariates_mat = adata.obs[covariates].values

            # Run ComBat
            X_corrected = pycombat(X_combat.T, batch_vec, covariates_mat).T

    except Exception as e:
        log.error(f"ComBat batch correction failed: {str(e)}")
        raise RuntimeError(f"ComBat batch correction failed: {str(e)}")

    # Store corrected data
    if inplace:
        adata.X = X_corrected
        log.info("Updated adata.X with batch-corrected data")

    if output_layer is not None:
        adata.layers[output_layer] = X_corrected
        log.info(f"Stored batch-corrected data in adata.layers['{output_layer}']")

    # Store integration metadata
    if "integration" not in adata.uns:
        adata.uns["integration"] = {}

    adata.uns["integration"]["combat"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "params": {
            "covariates": covariates,
            "layer": layer,
            "inplace": inplace,
            "output_layer": output_layer,
        },
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    log.info("ComBat batch correction complete")
    return adata


# ==============================================================================
# High-Level Batch Correction Workflow Function
# ==============================================================================


def batch_correction(
    adata: AnnData,
    batch_key: str,
    method: Literal["harmony", "scanorama", "scvi", "bbknn", "combat"] = "harmony",
    use_rep: str = "X_pca",
    hvg_key: Optional[str] = "highly_variable",
    layer: Optional[str] = "counts",
    plot: bool = True,
    save_dir: Optional[str] = None,
    copy: bool = False,
    force: bool = False,
    **kwargs,
) -> AnnData:
    """
    Performs batch correction using a specified integration method.

    This function serves as a high-level wrapper that calls the appropriate
    integration tool and visualizes the result.

    Args:
        adata: AnnData object. Must contain a dimensionality reduction (e.g., PCA).
        batch_key: Key in `adata.obs` that denotes the batch.
        method: Integration method to use:
            - "harmony": Fast and efficient linear integration (Korsunsky et al., 2019)
            - "scanorama": Panoramic stitching of single-cell data (Hie et al., 2019)
            - "scvi": Deep generative model (Lopez et al., 2018)
            - "bbknn": Batch-balanced k-nearest neighbors (Polański et al., 2020)
            - "combat": ComBat batch correction (Johnson et al., 2007)
        use_rep: The representation to use as input (for methods that need it).
        hvg_key: Key in `adata.var` specifying HVGs, used by Scanorama.
        layer: Layer containing raw counts, used by scVI.
        plot: If True, plots UMAPs before and after correction for comparison.
        save_dir: Directory to save plots.
        copy: If True, return a copy instead of modifying the original object.
        force: If True, run even if batch correction has already been applied.
        **kwargs: Additional arguments to pass to the integration method.

    Returns:
        The AnnData object with batch correction applied.

    Examples:
        >>> # Apply Harmony batch correction
        >>> adata = batch_correction(adata, batch_key="sample", method="harmony")
        >>>
        >>> # Apply scVI with custom parameters
        >>> adata = batch_correction(
        ...     adata,
        ...     batch_key="sample",
        ...     method="scvi",
        ...     layer="counts",
        ...     n_latent=20,
        ...     max_epochs=300
        ... )
    """
    # Parameter validation
    if batch_key not in adata.obs:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"Batch key '{batch_key}' not found in `adata.obs`.")

    # Handle copy
    if copy:
        adata = adata.copy()

    # Define output key based on method
    output_key = f"X_{method}"

    # Check if batch correction has already been applied
    if output_key in adata.obsm and not force:
        log.info(
            f"Batch correction using {method} already applied. Use force=True to rerun."
        )
        return adata

    # Check number of batches
    n_batches = adata.obs[batch_key].nunique()
    batch_counts = adata.obs[batch_key].value_counts()

    log.info(f"Performing batch correction using '{method}' method")
    log.info(f"Found {n_batches} batches in '{batch_key}':")
    for batch, count in batch_counts.items():
        log.info(f"  - {batch}: {count} cells ({count / adata.n_obs:.1%})")

    # Check if required dimensionality reduction exists
    if method in ["harmony", "bbknn"] and use_rep not in adata.obsm:
        log.error(f"Required representation '{use_rep}' not found in adata.obsm")
        raise ValueError(
            f"Required representation '{use_rep}' not found in adata.obsm."
        )

    # --- Plotting: Before Correction ---
    if plot:
        log.info("Generating visualization before batch correction")
        adata_before = adata.copy()

        # Compute neighbors and UMAP if needed
        if "neighbors" not in adata_before.uns:
            sc.pp.neighbors(adata_before, use_rep=use_rep)
            log.info(f"Computed neighbors using '{use_rep}'")

        if "X_umap" not in adata_before.obsm:
            sc.tl.umap(adata_before)
            log.info("Computed UMAP embedding")

    # --- Integration ---
    if method == "harmony":
        _integrate_harmony(
            adata,
            batch_key=batch_key,
            basis=use_rep,
            embedding_key=output_key,
            **kwargs,
        )

    elif method == "scanorama":
        # Get HVG list if specified
        hvg_list = None
        if hvg_key is not None and hvg_key in adata.var:
            hvg_list = adata.var_names[adata.var[hvg_key]].tolist()
            log.info(f"Using {len(hvg_list)} highly variable genes from '{hvg_key}'")

        _integrate_scanorama(
            adata,
            batch_key=batch_key,
            hvg=hvg_list,
            dims=adata.obsm[use_rep].shape[1] if use_rep in adata.obsm else 50,
            embedding_key=output_key,
            **kwargs,
        )

    elif method == "scvi":
        _integrate_scvi(
            adata,
            batch_key=batch_key,
            layer=layer,
            n_latent=adata.obsm[use_rep].shape[1] if use_rep in adata.obsm else 30,
            embedding_key=output_key,
            **kwargs,
        )

    elif method == "bbknn":
        _integrate_bbknn(
            adata,
            batch_key=batch_key,
            use_rep=use_rep,
            **kwargs,
        )

        # BBKNN creates a neighborhood graph but not an embedding
        # We'll create a dummy embedding to match the interface of other methods
        if "X_pca" in adata.obsm:
            adata.obsm[output_key] = adata.obsm["X_pca"].copy()
            log.info(f"Created placeholder embedding in adata.obsm['{output_key}']")

    elif method == "combat":
        _integrate_combat(
            adata,
            batch_key=batch_key,
            layer=layer,
            output_layer="combat_corrected",
            **kwargs,
        )

        # ComBat corrects the expression matrix but doesn't create an embedding
        # We'll compute PCA on the corrected data to create an embedding
        log.info("Computing PCA on batch-corrected data")
        # Save the original X
        X_original = adata.X.copy()

        # Use corrected data for PCA
        if "combat_corrected" in adata.layers:
            adata.X = adata.layers["combat_corrected"].copy()

        # Compute PCA
        sc.pp.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1))

        # Store the result as the integration embedding
        adata.obsm[output_key] = adata.obsm["X_pca"].copy()

        # Restore original X
        adata.X = X_original

    else:
        log.error(f"Unknown integration method: '{method}'")
        raise ValueError(
            f"Unknown integration method: '{method}'. "
            f"Choose from 'harmony', 'scanorama', 'scvi', 'bbknn', 'combat'."
        )

    log.info(
        f"Integration complete. Corrected embedding stored in `adata.obsm['{output_key}']`"
    )


    # Store method information in uns
    if "scrnatk" not in adata.uns:
        adata.uns["scrnatk"] = {}
    if "preprocess" not in adata.uns["scrnatk"]:
        adata.uns["scrnatk"]["preprocess"] = {}
    adata.uns["scrnatk"]["preprocess"]["integration"] = {
        "method": method,
        "input_layer": layer,
        "use_rep": use_rep,
        "batch_key": batch_key,
        "hvg_key": hvg_key,
    }
    
    # --- Plotting: After Correction ---
    if plot:
        log.info("Generating visualization after batch correction")

        # Compute neighbors and UMAP on the new, corrected embedding
        if method != "bbknn":  # BBKNN already computes the neighborhood graph
            sc.pp.neighbors(adata, use_rep=output_key)
            log.info(f"Computed neighbors using '{output_key}'")

        sc.tl.umap(adata)
        log.info("Computed UMAP embedding")

        # Create comparison figure
        log.info("Creating batch correction comparison plot")
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
            save_path = os.path.join(save_dir, f"batch_correction_{method}.png")
            plt.savefig(save_path, dpi=300)
            log.info(f"Saved comparison plot to {save_path}")

        plt.show()
        plt.close(fig)

    return adata


def evaluate_integration(
    adata: AnnData,
    batch_key: str,
    label_key: Optional[str] = None,
    integration_method: Optional[str] = None,
    use_rep: Optional[str] = None,
    n_neighbors: int = 30,
    metric: str = "euclidean",
    methods: List[str] = ["silhouette", "kbet", "graph_connectivity"],
    plot: bool = True,
    save_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Evaluates the quality of batch integration using multiple metrics.

    This function computes various metrics to assess how well batches have been
    integrated, optionally including biological label preservation.

    Args:
        adata: AnnData object with integrated data
        batch_key: Key in adata.obs identifying batch information
        label_key: Key in adata.obs identifying biological labels (cell types, etc.)
        integration_method: Name of integration method used (for reporting)
        use_rep: Representation to evaluate. If None, uses active .X or specified embedding
        n_neighbors: Number of neighbors to use for graph-based metrics
        metric: Distance metric to use for neighbor calculations
        methods: List of evaluation methods to use
        plot: Whether to create visualization of results
        save_path: Path to save the evaluation plot

    Returns:
        Dictionary containing computed metrics

    Examples:
        >>> # Evaluate harmony integration with cell type labels
        >>> metrics = evaluate_integration(
        ...     adata,
        ...     batch_key="sample",
        ...     label_key="cell_type",
        ...     integration_method="harmony",
        ...     use_rep="X_harmony"
        ... )
        >>>
        >>> # Evaluate without biological labels
        >>> metrics = evaluate_integration(adata, batch_key="batch")
    """
    try:
        import sklearn.metrics as skm
    except ImportError:
        log.error("sklearn is required for integration evaluation")
        raise ImportError("Please install sklearn: pip install scikit-learn")

    # Parameter validation
    if batch_key not in adata.obs:
        log.error(f"Batch key '{batch_key}' not found in adata.obs")
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")

    if label_key is not None and label_key not in adata.obs:
        log.warning(
            f"Label key '{label_key}' not found in adata.obs. Continuing without labels."
        )
        label_key = None

    # Set a default integration method name if not provided
    if integration_method is None:
        # Try to guess from obsm keys
        for key in adata.obsm.keys():
            if key.startswith("X_") and key != "X_pca":
                integration_method = key[2:]  # Remove "X_" prefix
                break

        if integration_method is None:
            integration_method = "unknown"

    # Determine which embedding to use
    if use_rep is None:
        # Look for integrated embeddings
        for key in ["X_harmony", "X_scanorama", "X_scVI", f"X_{integration_method}"]:
            if key in adata.obsm:
                use_rep = key
                break

        if use_rep is None:
            # Fall back to PCA
            if "X_pca" in adata.obsm:
                use_rep = "X_pca"
                log.warning(
                    "No integration embedding found. Using X_pca for evaluation."
                )
            else:
                log.warning(
                    "No dimensionality reduction found. Using adata.X for evaluation."
                )
    elif use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

    # Get the data to evaluate
    if use_rep is not None:
        log.info(f"Evaluating integration using '{use_rep}' representation")
        X = adata.obsm[use_rep]
    else:
        log.info("Evaluating integration using adata.X")
        X = adata.X

        # Convert to dense if sparse
        import scipy.sparse

        if scipy.sparse.issparse(X):
            X = X.toarray()

    # Initialize results dictionary
    results = {
        "method": integration_method,
        "n_batches": adata.obs[batch_key].nunique(),
        "n_cells": adata.n_obs,
    }

    if label_key is not None:
        results["n_labels"] = adata.obs[label_key].nunique()

    # Compute metrics
    log.info("Computing integration quality metrics")

    # 1. Silhouette score (measures batch mixing)
    if "silhouette" in methods:
        try:
            # Compute silhouette score based on batches (lower is better for integration)
            # We want batches to be well-mixed, so a good integration should have a negative score
            batch_silhouette = skm.silhouette_score(
                X,
                adata.obs[batch_key],
                metric=metric,
                sample_size=min(5000, adata.n_obs),  # Use sampling for large datasets
            )

            # Invert so higher is better
            results["batch_silhouette"] = -batch_silhouette

            # If we have labels, also compute silhouette based on biological labels
            # For biological labels, higher silhouette is better (we want label clusters preserved)
            if label_key is not None:
                label_silhouette = skm.silhouette_score(
                    X,
                    adata.obs[label_key],
                    metric=metric,
                    sample_size=min(5000, adata.n_obs),
                )
                results["label_silhouette"] = label_silhouette

                # Calculate a combined score that rewards both batch mixing and label preservation
                results["overall_silhouette"] = (
                    results["batch_silhouette"] + label_silhouette
                ) / 2

            log.info(
                f"Silhouette score (batch mixing): {results['batch_silhouette']:.4f}"
            )
            if label_key is not None:
                log.info(
                    f"Silhouette score (label preservation): {results['label_silhouette']:.4f}"
                )

        except Exception as e:
            log.warning(f"Failed to compute silhouette score: {str(e)}")

    # 2. k-BET (batch effect test - measures batch mixing in local neighborhoods)
    if "kbet" in methods:
        try:
            # We'll implement a simplified version of k-BET
            # Compute nearest neighbors for each cell
            from sklearn.neighbors import NearestNeighbors

            nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric).fit(X)
            distances, indices = nn.kneighbors(X)

            # For each neighborhood, test if batch distribution matches global distribution
            from scipy.stats import chi2_contingency

            # Get global batch distribution
            batch_categories = (
                adata.obs[batch_key].cat.categories
                if hasattr(adata.obs[batch_key], "cat")
                else sorted(adata.obs[batch_key].unique())
            )
            global_batch_counts = (
                adata.obs[batch_key].value_counts().reindex(batch_categories).fillna(0)
            )
            global_batch_freqs = global_batch_counts / global_batch_counts.sum()

            # Test each neighborhood
            rejection_rate = 0
            n_valid_tests = 0

            # We'll test a random subset for very large datasets
            import random

            test_indices = random.sample(range(adata.n_obs), min(1000, adata.n_obs))

            for i in test_indices:
                # Get batch labels of neighbors
                neighbor_batches = adata.obs[batch_key].iloc[indices[i]]
                local_batch_counts = (
                    neighbor_batches.value_counts().reindex(batch_categories).fillna(0)
                )

                # Skip neighborhoods with too few distinct batches for a meaningful test
                if (local_batch_counts > 0).sum() <= 1:
                    continue

                # Compute chi-square test
                expected_counts = global_batch_freqs * len(neighbor_batches)

                # Filter out categories with too few expected counts for a valid chi-square test
                valid_cats = expected_counts >= 5
                if sum(valid_cats) <= 1:
                    continue

                observed = local_batch_counts[valid_cats].values
                expected = expected_counts[valid_cats].values

                _, p_value, _, _ = chi2_contingency(
                    [observed, expected], correction=True
                )

                # A low p-value means we reject the null hypothesis that the local distribution
                # matches the global distribution (bad for integration)
                if p_value < 0.05:
                    rejection_rate += 1

                n_valid_tests += 1

            if n_valid_tests > 0:
                # Compute rejection rate (lower is better)
                rejection_rate = rejection_rate / n_valid_tests

                # Convert to an acceptance rate (higher is better)
                results["kbet_acceptance"] = 1 - rejection_rate

                log.info(f"k-BET acceptance rate: {results['kbet_acceptance']:.4f}")
            else:
                log.warning("Could not compute k-BET: insufficient valid tests")

        except Exception as e:
            log.warning(f"Failed to compute k-BET: {str(e)}")

    # 3. Graph connectivity (measures label preservation)
    if "graph_connectivity" in methods and label_key is not None:
        try:
            # Build a nearest-neighbor graph
            from sklearn.neighbors import kneighbors_graph

            # Create adjacency matrix
            adjacency = kneighbors_graph(
                X,
                n_neighbors=n_neighbors,
                mode="connectivity",
                metric=metric,
                include_self=False,
            )

            # For each label, compute the largest connected component size
            from scipy.sparse.csgraph import connected_components

            label_categories = (
                adata.obs[label_key].cat.categories
                if hasattr(adata.obs[label_key], "cat")
                else sorted(adata.obs[label_key].unique())
            )
            connectivities = []

            for label in label_categories:
                # Get indices for this label
                label_mask = adata.obs[label_key] == label
                label_indices = np.where(label_mask)[0]

                if len(label_indices) <= 1:
                    continue

                # Extract subgraph for this label
                subgraph = adjacency[label_indices][:, label_indices]

                # Find connected components
                n_components, component_labels = connected_components(
                    subgraph, directed=False
                )

                # Compute largest component size
                component_sizes = np.bincount(component_labels)
                largest_component = component_sizes.max()

                # Compute connectivity score (ratio of largest component to total cells in this label)
                connectivity = largest_component / len(label_indices)
                connectivities.append(connectivity)

            if connectivities:
                # Average connectivity across all labels
                results["graph_connectivity"] = np.mean(connectivities)
                log.info(f"Graph connectivity: {results['graph_connectivity']:.4f}")
            else:
                log.warning("Could not compute graph connectivity: no valid labels")

        except Exception as e:
            log.warning(f"Failed to compute graph connectivity: {str(e)}")

    # 4. Batch ASW (Adjusted Silhouette Width) - more sophisticated batch mixing metric
    if "batch_asw" in methods:
        try:
            # Compute average silhouette width per batch
            batch_categories = (
                adata.obs[batch_key].cat.categories
                if hasattr(adata.obs[batch_key], "cat")
                else sorted(adata.obs[batch_key].unique())
            )
            batch_asw_scores = []

            for batch in batch_categories:
                # Get indices for this batch
                batch_mask = adata.obs[batch_key] == batch
                batch_indices = np.where(batch_mask)[0]
                other_indices = np.where(~batch_mask)[0]

                if len(batch_indices) <= 1 or len(other_indices) <= 1:
                    continue

                # Sample for large datasets
                if len(batch_indices) > 1000:
                    batch_indices = np.random.choice(batch_indices, 1000, replace=False)
                if len(other_indices) > 1000:
                    other_indices = np.random.choice(other_indices, 1000, replace=False)

                # Compute distances from batch cells to other batch cells
                from sklearn.metrics import pairwise_distances

                # Within-batch distances
                batch_X = X[batch_indices]
                within_dists = pairwise_distances(batch_X, metric=metric)
                np.fill_diagonal(within_dists, np.inf)  # Exclude self-distances
                avg_within = np.mean(np.min(within_dists, axis=1))

                # Between-batch distances
                other_X = X[other_indices]
                between_dists = pairwise_distances(batch_X, other_X, metric=metric)
                avg_between = np.mean(np.min(between_dists, axis=1))

                # Compute ASW for this batch
                batch_size = len(batch_indices)
                total_cells = adata.n_obs
                weight = batch_size / total_cells

                asw = (avg_between - avg_within) / max(avg_between, avg_within)
                weighted_asw = weight * asw

                batch_asw_scores.append(weighted_asw)

            if batch_asw_scores:
                # Sum weighted ASW scores
                results["batch_asw"] = sum(batch_asw_scores)
                log.info(f"Batch ASW: {results['batch_asw']:.4f}")
            else:
                log.warning("Could not compute batch ASW: insufficient data")

        except Exception as e:
            log.warning(f"Failed to compute batch ASW: {str(e)}")

    # Plot results if requested
    if plot and results:
        try:
            log.info("Generating integration evaluation plot")

            # Determine which metrics to plot
            plot_metrics = [
                k
                for k in results.keys()
                if k not in ["method", "n_batches", "n_cells", "n_labels"]
            ]

            if not plot_metrics:
                log.warning("No metrics available to plot")
                return results

            # Create a figure
            fig, axes = plt.subplots(
                1, len(plot_metrics), figsize=(4 * len(plot_metrics), 5)
            )

            # Ensure axes is always a list
            if len(plot_metrics) == 1:
                axes = [axes]

            # Plot each metric
            for i, metric_name in enumerate(plot_metrics):
                value = results[metric_name]

                # Create a bar plot
                axes[i].bar([0], [value], color="skyblue", width=0.5)
                axes[i].set_title(metric_name.replace("_", " ").title())
                axes[i].set_ylabel("Score (higher is better)")
                axes[i].set_ylim(0, 1)  # All our metrics are normalized to [0,1]
                axes[i].set_xticks([])

                # Add the value as text
                axes[i].text(
                    0,
                    value / 2,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    fontweight="bold",
                )

                # Add a grid
                axes[i].grid(axis="y", linestyle="--", alpha=0.7)

            # Add a title
            method_name = (
                integration_method.capitalize() if integration_method else "Unknown"
            )
            fig.suptitle(
                f"Integration Quality: {method_name}", fontsize=16, fontweight="bold"
            )

            # Add a text box with summary information
            info_text = (
                f"Data: {adata.n_obs} cells, {adata.n_vars} genes\n"
                f"Batches: {results['n_batches']} ({batch_key})"
            )

            if label_key is not None:
                info_text += f"\nLabels: {results.get('n_labels', 'N/A')} ({label_key})"

            if use_rep:
                info_text += f"\nEmbedding: {use_rep}"

            fig.text(
                0.5,
                0.01,
                info_text,
                ha="center",
                va="bottom",
                bbox=dict(facecolor="white", alpha=0.8, boxstyle="round,pad=0.5"),
            )

            plt.tight_layout(rect=[0, 0.05, 1, 0.95])

            # Save if requested
            if save_path:
                os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                log.info(f"Saved evaluation plot to {save_path}")

            plt.show()

        except Exception as e:
            log.warning(f"Failed to create evaluation plot: {str(e)}")

    return results
