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
    n_clusters: Optional[int] = None,
    max_iter_harmony: int = 20,
    theta: float = 2.0,
    lambda_val: float = 1.0,
    sigma: float = 0.1,
    random_state: int = 42,
    plot_convergence: bool = False,
    verbose: bool = False,
    copy: bool = False,
    **kwargs,
):
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
    n_clusters : int, optional (default: None)
        Number of clusters for Harmony. If None, defaults to min(100, d * 20),
        where d is the number of dimensions in the basis.
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
    verbose : bool, optional (default: False)
        Whether to print progress information.
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
            # 尝试直接导入harmonypy作为备选
            import harmonypy
            raise ImportError(
                "Found harmonypy, but failed to import through scanpy. "
                "Please update scanpy: pip install --upgrade scanpy"
            )
        except ImportError:
            raise ImportError(
                "Please install harmonypy: pip install harmonypy"
            )
    
    # 检查batch_key是否在obs中
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
    
    # 检查批次数量
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        print(f"Warning: Only {n_batches} batch found. Harmony works best with multiple batches.")
        
    # 检查basis是否存在
    if basis not in adata.obsm:
        raise ValueError(f"Basis '{basis}' not found in adata.obsm. Run dimensionality reduction first.")
    
    # 计算自动n_clusters
    if n_clusters is None:
        n_dims = adata.obsm[basis].shape[1]
        n_clusters = min(100, n_dims * 20)
        if verbose:
            print(f"Automatically setting n_clusters to {n_clusters}")
    
    # 创建工作副本如果需要
    if copy:
        adata = adata.copy()
    
    if verbose:
        print(f"Running Harmony integration using '{batch_key}' as batch key...")
        print(f"  - Input dimensions: {adata.obsm[basis].shape}")
        print(f"  - Parameters: theta={theta}, lambda={lambda_val}, sigma={sigma}")
        print(f"  - Max iterations: {max_iter_harmony}")
    
    # 这里构建full_params将所有参数明确传递给harmony_integrate
    full_params = {
        'theta': theta,
        'lambda_': lambda_val,  # 注意这里使用lambda_而不是lambda_val
        'sigma': sigma,
        'n_clusters': n_clusters,
        'max_iter_harmony': max_iter_harmony,
        'random_state': random_state,
        'plot_convergence': plot_convergence,
    }
    
    # 添加用户自定义参数
    full_params.update(kwargs)
    
    # 运行Harmony集成
    try:
        harmony_integrate(
            adata,
            key=batch_key,
            basis=basis,
            adjusted_basis=f"{basis}_harmony",
            **full_params
        )
    except Exception as e:
        raise RuntimeError(f"Harmony integration failed: {str(e)}")
    
    # 检查Harmony输出是否成功生成
    expected_key = f"{basis}_harmony"
    if expected_key not in adata.obsm:
        raise RuntimeError(f"Harmony failed to generate output in adata.obsm['{expected_key}']")
    
    # 将结果复制到用户指定的键
    adata.obsm[embedding_key] = adata.obsm[expected_key].copy()
    
    # 如果用户指定了不同的键，并且不是默认的harmony键，则删除默认键
    if expected_key != embedding_key:
        del adata.obsm[expected_key]
    
    if verbose:
        print(f"Harmony integration complete. Results stored in adata.obsm['{embedding_key}']")
        print(f"  - Output dimensions: {adata.obsm[embedding_key].shape}")
    
    return adata


def _integrate_scanorama(
    adata: sc.AnnData,
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
    **kwargs,
):
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
        raise ImportError(
            "Please install Scanorama: pip install scanorama"
        )
        
    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
        
    # Check if batch_key has more than one unique value
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    if n_batches < 2:
        raise ValueError(f"Found only {n_batches} batch. Scanorama requires at least 2 batches for integration.")
    
    print(f"Running Scanorama integration on {n_batches} batches...")
    
    # Save original data
    original_obs = adata.obs.copy()
    original_var = adata.var.copy()
    
    # Split AnnData object by batch
    adatas_list = []
    genes_union = set()
    for b in batches:
        batch_data = adata[adata.obs[batch_key] == b].copy()
        
        # If a layer is specified, use it instead of adata.X
        if layer is not None:
            if layer not in batch_data.layers:
                raise ValueError(f"Layer '{layer}' not found in adata.")
            batch_data.X = batch_data.layers[layer].copy()
        
        adatas_list.append(batch_data)
        genes_union.update(batch_data.var_names)
    
    # Use HVGs if provided
    if hvg is not None:
        if not isinstance(hvg, list) and not isinstance(hvg, np.ndarray):
            raise TypeError("hvg must be a list or numpy array of gene names")
        
        print(f"Subsetting to {len(hvg)} highly variable genes for Scanorama.")
        genes_to_use = list(set(hvg).intersection(genes_union))
        if len(genes_to_use) < len(hvg):
            print(f"Warning: Only {len(genes_to_use)} out of {len(hvg)} HVGs found in the data.")
        if len(genes_to_use) == 0:
            raise ValueError("No genes from the provided hvg list were found in the data.")
    else:
        # If no HVGs provided, use all genes common to all batches
        common_genes = set.intersection(*[set(ad.var_names) for ad in adatas_list])
        genes_to_use = list(common_genes)
        print(f"Using {len(genes_to_use)} genes common to all batches for integration.")
    
    # Make sure all batches have the same genes
    for i, ad in enumerate(adatas_list):
        # Get genes to use in this batch
        genes_in_batch = [g for g in genes_to_use if g in ad.var_names]
        if len(genes_in_batch) < len(genes_to_use):
            print(f"Batch {i} contains {len(genes_in_batch)}/{len(genes_to_use)} integration genes.")
        
        # Subset the data to the genes in this batch
        adatas_list[i] = ad[:, genes_in_batch].copy()
    
    # Preprocess data for Scanorama
    if correct_expression:
        print("Computing integrated embedding and batch-corrected expression...")
    else:
        print("Computing integrated embedding only...")
    
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
            **kwargs
        )
    except Exception as e:
        raise RuntimeError(f"Scanorama integration failed: {str(e)}")
    
    # Reconstruct the integrated embedding
    # Note: Scanorama returns a list of embeddings, one for each batch
    cell_indices = np.concatenate([np.where(adata.obs[batch_key] == b)[0] for b in batches])
    integrated_embedding = np.zeros((adata.shape[0], dims))
    
    counter = 0
    for i, b in enumerate(batches):
        batch_size = np.sum(adata.obs[batch_key] == b)
        integrated_embedding[cell_indices[counter:counter+batch_size]] = adatas_list[i].obsm["X_scanorama"]
        counter += batch_size
    
    adata.obsm[embedding_key] = integrated_embedding
    print(f"Integrated embedding stored in adata.obsm['{embedding_key}']")
    
    # If requested, also return batch-corrected expression
    if return_corrected_expression:
        if not hasattr(adatas_list[0], 'X_scanorama'):
            print("Warning: Batch-corrected expression not available. Set correct_expression=True to enable this.")
        else:
            corrected_exp = np.zeros(adata.shape)
            genes_idx = {g: i for i, g in enumerate(adata.var_names)}
            
            counter = 0
            for i, b in enumerate(batches):
                batch_size = np.sum(adata.obs[batch_key] == b)
                # Mapping the corrected expression to the original gene order
                for j, gene in enumerate(adatas_list[i].var_names):
                    if gene in genes_idx:
                        corrected_exp[cell_indices[counter:counter+batch_size], genes_idx[gene]] = \
                            adatas_list[i].X_scanorama[:, j]
                counter += batch_size
            
            # Save the corrected expression in a new layer
            adata.layers['scanorama_corrected'] = corrected_exp
            print("Batch-corrected expression stored in adata.layers['scanorama_corrected']")
    
    # Make sure to return the original adata object with the same obs and var
    for col in original_obs.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = original_obs[col].values
    
    for col in original_var.columns:
        if col not in adata.var.columns:
            adata.var[col] = original_var[col].values
    
    return adata


def _integrate_scvi(
    adata: sc.AnnData,
    batch_key: str,
    layer: Optional[str] = "counts",
    n_layers: int = 2, 
    n_latent: int = 30,
    batch_size: int = 256,
    max_epochs: int = 500,
    embedding_key: str = "X_scVI",
    gene_likelihood: str = "nb",
    use_gpu: Optional[bool] = None,
    save_model: bool = False,
    model_path: Optional[str] = None,
    plan_kwargs: Optional[dict] = None,
    **kwargs,
):
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
        raise ImportError(
            "Please install scvi-tools: pip install scvi-tools"
        )

    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
        
    # Check if the batch_key has more than one unique value
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        print(f"Warning: Only {n_batches} batch found. scVI works best with multiple batches.")
    
    # GPU detection
    if use_gpu is None:
        use_gpu = True
        try:
            import torch
            use_gpu = torch.cuda.is_available()
            if use_gpu:
                print("GPU detected. Using GPU for scVI.")
            else:
                print("No GPU detected. Using CPU for scVI.")
        except:
            use_gpu = False
            print("Could not detect GPU. Using CPU for scVI.")
    
    print(f"Setting up scVI with {n_layers} layers and {n_latent} latent dimensions...")
    
    # scVI requires setup on the AnnData object
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)

    # Prepare the model
    if plan_kwargs is None:
        plan_kwargs = {}
    
    # Create the scVI model
    model = scvi.model.SCVI(
        adata, 
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood)
    
    # Train the model
    print(f"Training scVI model with batch_size={batch_size}, max_epochs={max_epochs}...")
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
        plan_kwargs=plan_kwargs,
        **kwargs)  # pass any additional arguments to train()
    
    print("Extracting latent representation...")
    # Extract the latent representation
    adata.obsm[embedding_key] = model.get_latent_representation()
    
    # Save the model if requested
    if save_model:
        if model_path is None:
            raise ValueError("model_path must be provided when save_model=True")
        print(f"Saving model to {model_path}")
        model.save(model_path)
    
    print(f"scVI integration complete. Latent representation stored in adata.obsm['{embedding_key}']")
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
