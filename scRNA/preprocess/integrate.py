"""
Data integration methods for single-cell RNA-seq data.
"""

import anndata as ad
import numpy as np
import scanpy as sc
from typing import Optional, Literal, Union, List
import matplotlib.pyplot as plt

def integrate_scanorama(
    adata: sc.AnnData,
    batch_key: str = "sampleID", 
    dims: int = 50,
    hvg_key: Optional[str] = "highly_variable",
    min_batches: int = 2,
    embedding_key: str = "X_scanorama",
):
    """
    Integrate data using Scanorama.
    
    Args:
        adata: AnnData object
        batch_key: Key in adata.obs that denotes the batch
        dims: Number of dimensions for the integrated embedding
        hvg_key: Key in adata.var to use for highly variable genes
        min_batches: Minimum number of batches in which a gene should be variable
        embedding_key: Key to store the integrated embedding
        
    Returns:
        AnnData with Scanorama integration
    """
    try:
        import scanorama
    except ImportError:
        raise ImportError("Please install Scanorama: pip install scanorama")
    
    # Check batch_key
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
    
    # Use HVGs across multiple batches
    if hvg_key is not None and hvg_key in adata.var:
        # If we have information about number of batches where gene is variable
        if "highly_variable_nbatches" in adata.var:
            var_select = adata.var.highly_variable_nbatches >= min_batches
        else:
            var_select = adata.var[hvg_key]
        
        var_genes = var_select.index[var_select]
        print(f"Using {len(var_genes)} genes for integration")
    else:
        # Use all genes if no HVG information
        var_genes = adata.var_names
        print("No highly variable genes specified. Using all genes.")
    
    # Split per batch into new objects
    batches = adata.obs[batch_key].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        alldata[batch] = adata[adata.obs[batch_key] == batch,]

    # Subset the individual dataset to the variable genes
    alldata2 = dict()
    for ds in alldata.keys():
        print(f"Processing batch: {ds}")
        alldata2[ds] = alldata[ds][:,var_genes]

    # Convert to list of AnnData objects
    adatas = list(alldata2.values())

    # Run scanorama.integrate
    print(f"Running Scanorama with {dims} dimensions...")
    scanorama.integrate_scanpy(adatas, dimred=dims) 
    
    # Get all the integrated matrices
    scanorama_int = [ad_.obsm['X_scanorama'] for ad_ in adatas]

    # Make into one matrix
    all_s = np.concatenate(scanorama_int)
    print(f"Integrated matrix shape: {all_s.shape}")

    # Add to the AnnData object
    adata.obsm[embedding_key] = all_s
    
    print(f"Scanorama integration complete. Results stored in adata.obsm['{embedding_key}']")
    
    return adata


def integrate_scvi(
    adata: sc.AnnData, 
    layer: Optional[str] = "counts", 
    batch_key: str = "sampleID", 
    batch_size: int = 256,
    max_epochs: int = 500,
    n_layers: int = 2, 
    n_latent: int = 30, 
    embedding_key: str = "X_scVI",
):
    """
    Integrate data using scVI.
    
    Args:
        adata: AnnData object
        layer: Layer containing count data
        batch_key: Key in adata.obs that denotes the batch
        batch_size: Batch size for training
        max_epochs: Maximum number of training epochs
        n_layers: Number of layers in the model
        n_latent: Number of latent dimensions
        embedding_key: Key to store the integrated embedding
        
    Returns:
        AnnData with scVI integration
    """
    try:
        import scvi
    except ImportError:
        raise ImportError("Please install scVI: pip install scvi-tools")
    
    # Check batch_key
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
    
    # Setup anndata for scVI
    print("Setting up AnnData for scVI...")
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
    
    # Create and train the model
    print(f"Creating scVI model with {n_latent} latent dimensions...")
    model = scvi.model.SCVI(
        adata, 
        n_layers=n_layers, 
        n_latent=n_latent,
        gene_likelihood="nb")
    
    print(f"Training scVI model (max_epochs={max_epochs})...")
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
    )
    
    # Get latent representation
    print("Extracting latent representation...")
    adata.obsm[embedding_key] = model.get_latent_representation()
    
    print(f"scVI integration complete. Results stored in adata.obsm['{embedding_key}']")
 
    return adata


def integrate_harmony(
    adata: sc.AnnData, 
    batch_key: str = "sampleID", 
    basis: str = 'X_pca',
    embedding_key: str = "X_harmony",
):
    """
    Integrate data using Harmony.
    
    Args:
        adata: AnnData object
        batch_key: Key in adata.obs that denotes the batch
        basis: Representation to use for integration
        embedding_key: Key to store the integrated embedding
        
    Returns:
        AnnData with Harmony integration
    """
    try:
        from scanpy.external.pp import harmony_integrate
    except ImportError:
        raise ImportError("Please install harmonypy: pip install harmonypy")
    
    # Check batch_key
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
    
    # Check if basis exists
    if basis not in adata.obsm:
        raise ValueError(f"Basis '{basis}' not found in adata.obsm. Run PCA first.")
    
    print(f"Running Harmony integration using {basis} as input...")
    harmony_integrate(adata, key=batch_key, basis=basis)
    
    # If the key isn't the desired one, rename it
    if f"X_pca_{batch_key}" in adata.obsm and embedding_key != f"X_pca_{batch_key}":
        adata.obsm[embedding_key] = adata.obsm[f"X_pca_{batch_key}"].copy()
    
    print(f"Harmony integration complete. Results stored in adata.obsm['{embedding_key}']")
    
    return adata


def batch_correction(
    adata: sc.AnnData,
    batch_key: str,
    method: Literal["harmony", "scanorama", "scvi", "bbknn"] = "harmony",
    embedding_key: Optional[str] = None,
    n_pcs: int = 50,
    use_rep: str = "X_pca",
    hvg_key: str = "highly_variable",
    layer: Optional[str] = None,
    plot: bool = True,
    save_dir: Optional[str] = None,
    **kwargs
) -> sc.AnnData:
    """
    Perform batch correction on the AnnData object.
    
    Args:
        adata: AnnData object
        batch_key: Key in adata.obs that denotes the batch
        method: Method for batch correction
        embedding_key: Key to store the corrected embedding
        n_pcs: Number of principal components to use
        use_rep: Representation to use for batch correction
        hvg_key: Key in adata.var for highly variable genes
        layer: Layer to use for scVI integration
        plot: Whether to plot UMAP before and after correction
        save_dir: Directory to save plots
        **kwargs: Additional arguments to pass to the integration method
        
    Returns:
        AnnData with batch-corrected representation
    """
    import gc
    
    # Check if PCA has been computed
    if use_rep == "X_pca" and use_rep not in adata.obsm:
        print("Computing PCA...")
        sc.pp.pca(adata, n_comps=n_pcs)
    
    # Check batch_key
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
    
    # Set default embedding key if not provided
    if embedding_key is None:
        embedding_key = f"X_{method}"
    
    # Make a copy for visualization
    if plot:
        adata_copy = adata.copy()
        sc.pp.neighbors(adata_copy, use_rep=use_rep)
        sc.tl.umap(adata_copy)
    
    print(f"Performing batch correction using {method}...")
    
    if method == "harmony":
        integrate_harmony(adata, batch_key=batch_key, basis=use_rep, embedding_key=embedding_key, **kwargs)
    
    elif method == "scanorama":
        integrate_scanorama(adata, batch_key=batch_key, dims=n_pcs, hvg_key=hvg_key, embedding_key=embedding_key, **kwargs)
    
    elif method == "scvi":
        integrate_scvi(adata, layer=layer, batch_key=batch_key, n_latent=n_pcs, embedding_key=embedding_key, **kwargs)
    
    elif method == "bbknn":
        try:
            from scanpy.external.pp import bbknn
            
            # Run BBKNN
            bbknn(adata, batch_key=batch_key, use_rep=use_rep, **kwargs)
            
            # BBKNN directly modifies the neighbor graph
            print("BBKNN batch correction complete.")
            
        except ImportError:
            raise ImportError("Please install BBKNN: pip install bbknn")
    
    else:
        raise ValueError(f"Unknown batch correction method: {method}")
    
    gc.collect()
    
    # Plot results if requested
    if plot:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # Before correction
        sc.pl.umap(adata_copy, color=batch_key, ax=axes[0], show=False, title="Before Batch Correction")
        
        # After correction
        if method != "bbknn":
            # Compute neighbors and UMAP with corrected embedding
            sc.pp.neighbors(adata, use_rep=embedding_key)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=batch_key, ax=axes[1], show=False, title=f"After {method.upper()} Correction")
        
        plt.tight_layout()
        
        if save_dir:
            import os
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f"batch_correction_{method}.png"), dpi=300)
        
        plt.show()
    
    return adata

def adaptive_batch_correction(
    adata,
    batch_key,
    method="auto",
    n_cells_threshold=5000,
    n_batches_threshold=5,
    **kwargs
):
    """自适应批次校正，根据数据规模选择最佳方法"""
    # 验证批次键
    if batch_key not in adata.obs:
        raise ValueError(f"批次键'{batch_key}'不在adata.obs中")
    
    # 计算批次数和细胞数
    n_batches = len(adata.obs[batch_key].unique())
    n_cells = adata.n_obs
    
    if method == "auto":
        # 根据数据规模自动选择方法
        if n_cells > n_cells_threshold and n_batches > n_batches_threshold:
            print(f"检测到大型数据集({n_cells}细胞, {n_batches}批次), 使用scVI进行批次校正")
            method = "scvi"
        elif n_batches <= 2:
            print(f"检测到较少批次({n_batches}), 使用Harmony进行批次校正")
            method = "harmony"
        else:
            print(f"使用Scanorama进行批次校正")
            method = "scanorama"
    
    print(f"执行批次校正, 方法: {method}")
    return batch_correction(adata, batch_key=batch_key, method=method, **kwargs)