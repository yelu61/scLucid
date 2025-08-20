"""
Normalization functions for single-cell RNA-seq data.

This module provides methods for normalizing raw count data in single-cell RNA-seq
analysis, including standard library size normalization, centered log-ratio (CLR),
and advanced methods like Pearson residuals.
"""

import logging
import os
from typing import List, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from ..utils.utils import use_layer_as_X

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = ["normalize_data", "regress_out"]


# --- Helper Functions ---
def _plot_normalization_global(
    adata: AnnData,
    input_data: Union[np.ndarray, scipy.sparse.spmatrix],
    output_data: Union[np.ndarray, scipy.sparse.spmatrix],
    method: str = "standard",
    log_transformed: bool = True,
    save_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Generate global distribution plots comparing before and after normalization.

    Args:
        adata: AnnData object
        input_data: Original data before normalization
        output_data: Normalized data
        method: Normalization method used
        log_transformed: Whether log transformation was applied
        save_dir: Directory to save plots

    Returns:
        matplotlib Figure object
    """
    # Set up the plotting style
    rc_params = {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
        }

    # Create figure
    with plt.rc_context(rc_params):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
        fig.suptitle(
            f"Data Distributions Before and After {method.capitalize()} Normalization",
            fontsize=16,
            color="black",
        )

        # Get cell sums for before and after
        if scipy.sparse.issparse(input_data):
            before_sums = input_data.sum(axis=1).A1
        else:
            before_sums = input_data.sum(axis=1)

        if scipy.sparse.issparse(output_data):
            after_sums = output_data.sum(axis=1).A1
        else:
            after_sums = output_data.sum(axis=1)

        # Plot before normalization
        sns.histplot(
            before_sums,
            bins=100,
            kde=True,
            ax=axes[0],
            color="navy",
        )
        axes[0].set_title("Before Normalization")
        axes[0].set_xlabel("Total Counts per Cell")
        axes[0].set_ylabel("Frequency")
        axes[0].text(
            0.05,
            0.95,
            f"Mean: {before_sums.mean():.1f}\nMedian: {np.median(before_sums):.1f}",
            transform=axes[0].transAxes,
            va="top",
        )

        # Plot after normalization
        sns.histplot(
            after_sums,
            bins=100,
            kde=True,
            ax=axes[1],
            color="crimson",
        )
        title_suffix = " (Log-Transformed)" if log_transformed else ""
        axes[1].set_title(f"After {method.capitalize()} Normalization{title_suffix}")
        axes[1].set_xlabel("Sum of Normalized Values per Cell")
        axes[1].set_ylabel("Frequency")
        axes[1].text(
            0.05,
            0.95,
            f"Mean: {after_sums.mean():.1f}\nMedian: {np.median(after_sums):.1f}",
            transform=axes[1].transAxes,
            va="top",
        )

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(
                os.path.join(save_dir, f"normalization_{method}_global.png"), dpi=300
            )

    return fig


def _plot_normalization_comparison(
    adata: AnnData,
    original_layer: str,
    normalized_layer: str,
    log_transformed: bool = True,
    n_genes: int = 10,
    gene_subset: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Generate comparison plots showing gene expression distributions before and after normalization.

    This function is useful for visualizing how normalization affects the distribution
    of individual genes, particularly for highly expressed genes that can be most
    affected by normalization.

    Args:
        adata: AnnData object with raw and normalized data
        original_layer: Layer containing original data
        normalized_layer: Layer containing normalized data
        log_transformed: Whether the normalized data is log-transformed
        n_genes: Number of top genes to show (by mean expression)
        gene_subset: Specific genes to show instead of top expressing genes
        save_dir: Directory to save the generated figure

    Returns:
        matplotlib Figure object with the comparison plots

    Raises:
        ValueError: If specified layers don't exist

    Examples:
        >>> # Compare raw counts to normalized data
        >>> fig = plot_normalization_comparison(adata, 'counts', 'log1p_norm')
        >>> plt.show()
        >>>
        >>> # Compare specific genes
        >>> fig = plot_normalization_comparison(
        ...     adata, 'counts', 'log1p_norm',
        ...     gene_subset=['MALAT1', 'GAPDH', 'ACTB']
        ... )
    """
    # Validate inputs
    if original_layer not in adata.layers:
        raise ValueError(f"Original layer '{original_layer}' not found in adata.layers")
    if normalized_layer not in adata.layers:
        raise ValueError(
            f"Normalized layer '{normalized_layer}' not found in adata.layers"
        )

    # Import scipy here to avoid importing at the top level
    import scipy.sparse

    # Select genes to plot
    if gene_subset is not None:
        # Use user-specified genes
        genes_to_plot = [g for g in gene_subset if g in adata.var_names]
        if not genes_to_plot:
            raise ValueError("None of the specified genes were found in the data")
        if len(genes_to_plot) < len(gene_subset):
            log.warning(
                f"Only {len(genes_to_plot)}/{len(gene_subset)} specified genes were found"
            )
    else:
        # Select top expressed genes
        mean_expr = np.array(adata.layers[original_layer].mean(axis=0)).flatten()
        top_genes_idx = np.argsort(-mean_expr)[:n_genes]
        genes_to_plot = adata.var_names[top_genes_idx].tolist()

    # Create the figure
    n_to_plot = len(genes_to_plot)
    fig, axes = plt.subplots(n_to_plot, 2, figsize=(14, 3 * n_to_plot))

    # Ensure axes is 2D even with a single gene
    if n_to_plot == 1:
        axes = np.array([axes])

    # Define helper function to safely get data
    def get_data(layer, gene):
        gene_idx = adata.var.index.get_loc(gene) 
        data = adata.layers[layer][:, gene_idx]
        if scipy.sparse.issparse(data):
            return data.toarray().flatten()
        return data.flatten()

    # Loop through genes and create plots
    for i, gene in enumerate(genes_to_plot):
        try:
            # Get data for this gene
            before_data = get_data(original_layer, gene)
            after_data = get_data(normalized_layer, gene)

            # Original data plot
            sns.histplot(before_data, bins=50, kde=True, ax=axes[i, 0])
            axes[i, 0].set_title(f"{gene} - Before")
            axes[i, 0].set_ylabel("Frequency")
            axes[i, 0].text(
                0.05,
                0.95,
                f"Mean: {before_data.mean():.3f}\nStd: {before_data.std():.3f}\n"
                f"% zeros: {(before_data == 0).mean():.1%}",
                transform=axes[i, 0].transAxes,
                va="top",
            )

            # Normalized data plot
            sns.histplot(after_data, bins=50, kde=True, ax=axes[i, 1])
            suffix = " (Log-transformed)" if log_transformed else ""
            axes[i, 1].set_title(f"{gene} - After{suffix}")
            axes[i, 1].text(
                0.05,
                0.95,
                f"Mean: {after_data.mean():.3f}\nStd: {after_data.std():.3f}\n"
                f"% zeros: {(after_data == 0).mean():.1%}",
                transform=axes[i, 1].transAxes,
                va="top",
            )
        except Exception as e:
            log.warning(f"Failed to plot gene {gene}: {str(e)}")
            axes[i, 0].text(0.5, 0.5, f"Error plotting {gene}", ha="center")
            axes[i, 1].text(0.5, 0.5, f"Error plotting {gene}", ha="center")

    plt.tight_layout()

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(
            os.path.join(save_dir, "normalization_gene_comparison.png"), dpi=300
        )

    return fig


# ------ Main Functions ------
def normalize_data(
    adata: AnnData,
    method: Literal[
        "standard", "scran", "pearson_residuals", "sctransform", "clr"
    ] = "standard",
    layer: Optional[str] = "counts",
    output_layer: str = "log1p_norm",
    # Params for 'standard' method
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    # Common params
    log_transform: bool = True,
    plot: bool = True,
    save_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Normalize and log-transform single-cell RNA sequencing data.

    This function implements multiple normalization strategies to account for
    differences in sequencing depth and composition biases between cells.

    Args:
        adata: AnnData object containing single-cell data.
        method: Normalization method to use:
            - "standard": Standard library size normalization (Scanpy default)
            - "scran": Size factors using pooling (requires R, rpy2, and scran)
            - "pearson_residuals": Variance stabilizing using Pearson residuals
            - "sctransform": SCTransform method (requires R, rpy2, and Seurat)
            - "clr": Centered log-ratio transformation for compositional data
        layer: Name of the layer containing raw count data. If None, uses adata.X.
        target_sum: Total count each cell will have after normalization.
            If 1e6, this gives CPM normalization. If None, uses median total count.
        exclude_highly_expressed: Whether to exclude highly expressed genes from
            normalization factor calculation to avoid biases.
        max_fraction: When exclude_highly_expressed=True, genes with more counts than
            this fraction of the original total counts are excluded.
        log_transform: Whether to log-transform the data after normalization.
            Not applied for methods that already include log transformation.
        output_layer: Name of the layer to store normalized (and log-transformed) data.
        plot: Whether to generate distribution plots before and after normalization.
        save_dir: Directory to save plots. If None, plots are not saved.
        force: Whether to proceed even if the output layer already exists.

    Returns:
        AnnData object with normalized data in the specified output layer.

    Raises:
        ValueError: If parameters are invalid or required dependencies are missing.
        RuntimeError: If normalization fails due to computational issues.

    Examples:
        >>> # Standard normalization with log transformation
        >>> adata = normalize_data(adata, method="standard")
        >>>
        >>> # CLR normalization (useful for compositional data)
        >>> adata = normalize_data(adata, method="clr")
    """
    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata
    
    # --- Input Data Validation ---
    if layer and layer not in adata.layers:
        raise ValueError(f"Input layer '{layer}' not found in adata.layers.")
    if "counts" not in adata.layers and layer != "counts":
        adata.layers["counts"] = adata.X.copy()
        log.info("Saved current adata.X to adata.layers['counts'] for reference.")
    source_adata = sc.AnnData(adata.layers[layer].copy()) if layer else adata.copy()

    min_val = source_adata.X.min()
    if min_val < 0: # Check for zero or negative values that might cause problems
        log.warning("Input data contains negative values. Normalization may produce unexpected results.")
    elif np.all(source_adata.X < 1) and np.any(source_adata.X > 0):
        log.warning("Input data seems to be already normalized (all values < 1).")
    
    # Determine input data source
    if layer is not None and layer not in adata.layers:
        log.warning(
            f"Layer '{layer}' not found in adata.layers. Using adata.X instead."
        )
        input_data_source = adata.X
    else:
        input_data_source = adata.layers.get(layer, adata.X)

    # Calculate basic statistics for logging and validation
    total_counts = np.array(input_data_source.sum(axis=1)).flatten()

    log.info(
        f"Input data summary: median counts per cell = {np.median(total_counts):.1f}, "
        f"range = [{np.min(total_counts):.1f}, {np.max(total_counts):.1f}]"
    )
    log.info(f"Normalizing data from '{layer or 'adata.X'}' using method '{method}'.")

    # --- Method-Specific Normalization ---
    try:
        if method == "standard":
            log.info("Applying standard library size normalization...")
            if target_sum is not None and target_sum <= 0:
                raise ValueError("target_sum must be a positive number.")
            if not 0 <= max_fraction <= 1:
                raise ValueError("max_fraction must be between 0 and 1 (inclusive).")
            log.info(f"Parameters: target_sum={target_sum}, log_transform={log_transform}")
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
            log.warning(
                "The 'scran' method requires R, rpy2, and Bioconductor's scran package."
            )
            try:
                import rpy2
                from rpy2.robjects import numpy2ri, pandas2ri, r
                from rpy2.robjects.packages import importr

                pandas2ri.activate()
                numpy2ri.activate()

                # Check if the required R packages are installed
                try:
                    scran = importr("scran")
                    scuttle = importr("scuttle")
                    singlecellexperiment = importr("SingleCellExperiment")

                    log.info("Successfully loaded R packages for scran normalization.")

                    # Convert to sparse matrix if it's not already
                    X_mat = input_data_source.copy()
                    if not scipy.sparse.issparse(X_mat):
                        X_mat = scipy.sparse.csr_matrix(X_mat)

                    # Convert to R sparse matrix format
                    counts_r = pandas2ri.py2rpy(
                        pd.DataFrame(
                            X_mat.toarray(),
                            columns=adata.var_names,
                            index=adata.obs_names,
                        )
                    )

                    # Create SingleCellExperiment object
                    sce = r("""
                    function(counts) {
                        sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))
                        return(sce)
                    }
                    """)(counts_r)

                    # Run scran normalization
                    sce_norm = r("""
                    function(sce) {
                        set.seed(42)
                        clusters <- scran::quickCluster(sce)
                        sce <- scran::computeSumFactors(sce, clusters=clusters)
                        sce <- scuttle::logNormCounts(sce)
                        return(sce)
                    }
                    """)(sce)

                    # Extract normalized values
                    logcounts = r("function(sce) { return(logcounts(sce)) }")(sce_norm)
                    X_norm = numpy2ri.rpy2py(logcounts)

                    # scran already log-transforms the data
                    log_transform = False

                except Exception as r_err:
                    log.error(f"Failed to run scran normalization: {r_err}")
                    log.warning("Falling back to standard normalization.")
                    temp_adata = sc.AnnData(input_data_source.copy())
                    sc.pp.normalize_total(
                        temp_adata, target_sum=target_sum, inplace=True
                    )
                    X_norm = temp_adata.X.copy()

            except ImportError:
                log.warning(
                    "rpy2 is not installed. Falling back to standard normalization."
                )
                temp_adata = sc.AnnData(input_data_source.copy())
                sc.pp.normalize_total(temp_adata, target_sum=target_sum, inplace=True)
                X_norm = temp_adata.X.copy()

        elif method == "pearson_residuals":
            log.info(
                "Applying Pearson residuals normalization (variance stabilizing transformation)..."
            )
            X_norm = sc.experimental.pp.normalize_pearson_residuals(
                sc.AnnData(input_data_source.copy()), inplace=False
            )["X"]
            # Log transform is not needed for Pearson residuals
            log_transform = False

        elif method == "clr":
            log.info("Applying centered log-ratio (CLR) transformation...")
            # Apply pseudocount to handle zeros
            X_data = input_data_source.copy()

            # Handle sparse matrices
            is_sparse = scipy.sparse.issparse(X_data)
            if is_sparse:
                X_data = X_data.toarray()

            # Add pseudocount and calculate geometric mean
            X_data = X_data + 1  # Pseudocount
            log_X = np.log(X_data)
            geo_means = np.exp(np.mean(log_X, axis=1, keepdims=True))

            # Apply CLR transformation
            X_norm = log_X - np.log(geo_means)

            # Convert back to sparse if original was sparse
            if is_sparse:
                X_norm = scipy.sparse.csr_matrix(X_norm)

            # CLR already applies log transformation
            log_transform = False

        elif method == "sctransform":
            log.warning("sctransform method requires R, rpy2, and Seurat package.")
            try:
                import rpy2
                from rpy2.robjects import numpy2ri, pandas2ri, r
                from rpy2.robjects.packages import importr

                pandas2ri.activate()
                numpy2ri.activate()

                try:
                    seurat = importr("Seurat")
                    seuratobject = importr("SeuratObject")

                    log.info("Successfully loaded R packages for sctransform.")

                    # Convert to dense matrix if it's sparse
                    X_mat = input_data_source.copy()
                    if scipy.sparse.issparse(X_mat):
                        X_mat = X_mat.toarray()

                    # Convert to R matrix format
                    counts_r = numpy2ri.py2rpy(X_mat)
                    gene_names = numpy2ri.py2rpy(np.array(adata.var_names))
                    cell_names = numpy2ri.py2rpy(np.array(adata.obs_names))

                    # Run SCTransform
                    sct_result = r("""
                    function(counts, gene_names, cell_names) {
                        set.seed(42)
                        # Set row and column names
                        rownames(counts) <- gene_names
                        colnames(counts) <- cell_names
                        
                        # Create Seurat object
                        sobj <- SeuratObject::CreateSeuratObject(counts=t(counts))
                        
                        # Run SCTransform
                        sobj <- Seurat::SCTransform(sobj, verbose=FALSE)
                        
                        # Return normalized data
                        return(t(SeuratObject::GetAssayData(sobj, slot="scale.data")))
                    }
                    """)(counts_r, gene_names, cell_names)

                    X_norm = numpy2ri.rpy2py(sct_result)

                    # SCTransform applies its own normalization and scaling
                    log_transform = False

                except Exception as r_err:
                    log.error(f"Failed to run sctransform: {r_err}")
                    log.warning("Falling back to standard normalization.")
                    temp_adata = sc.AnnData(input_data_source.copy())
                    sc.pp.normalize_total(
                        temp_adata, target_sum=target_sum, inplace=True
                    )
                    X_norm = temp_adata.X.copy()

            except ImportError:
                log.warning(
                    "rpy2 is not installed. Falling back to standard normalization."
                )
                temp_adata = sc.AnnData(input_data_source.copy())
                sc.pp.normalize_total(temp_adata, target_sum=target_sum, inplace=True)
                X_norm = temp_adata.X.copy()

        else:
            raise ValueError(
                f"Unknown normalization method: {method}. "
                f"Choose from 'standard', 'scran', 'pearson_residuals', 'sctransform', 'clr'."
            )

    except Exception as e:
        log.error(f"Normalization failed: {str(e)}")
        raise RuntimeError(f"Failed to normalize data using {method} method: {str(e)}")

    # --- Store Results and Log Transform ---
    # Store normalized data before log transformation
    adata.layers["normalized"] = X_norm.copy()
    log.info("Stored pre-log normalized data in adata.layers['normalized']")

    # Apply log transform if needed
    if log_transform:
        log.info("Applying log1p transformation.")
        try:
            # sc.pp.log1p is robust for both sparse and dense matrices
            adata.layers[output_layer] = sc.pp.log1p(X_norm, copy=True)
        except Exception as e:
            log.error(f"Log transformation failed: {str(e)}")
            raise RuntimeError(f"Failed to apply log1p transformation: {str(e)}")
    else:
        adata.layers[output_layer] = X_norm.copy()

    # Calculate and log post-normalization statistics
    final_data = adata.layers[output_layer]
    final_sum = np.array(final_data.sum(axis=1)).flatten()

    log.info(f"Normalization complete. Final data in adata.layers['{output_layer}']")
    log.info(
        f"Final data summary: median sum per cell = {np.median(final_sum):.2f}, "
        f"range = [{np.min(final_sum):.2f}, {np.max(final_sum):.2f}]"
    )

    # --- Store parameters in the unified namespace ---
    if "scrnatk" not in adata.uns:
        adata.uns["scrnatk"] = {}
    if "preprocess" not in adata.uns["scrnatk"]:
        adata.uns["scrnatk"]["preprocess"] = {}
    adata.uns["scrnatk"]["preprocess"]["normalization"] = {
        "method": method,
        "target_sum": target_sum,
        "log_transform": log_transform,
        "output_layer": output_layer,
    }

    # --- Plotting ---
    if plot:
        try:
            log.info("Generating normalization comparison plots...")
            # Plot general cell-level distribution
            global_fig = _plot_normalization_global(
                adata,
                input_data_source,
                adata.layers[output_layer],
                method=method,
                log_transformed=log_transform,
                save_dir=save_dir,
            )

            # Plot gene-level distributions
            gene_fig = _plot_normalization_comparison(
                adata,
                "counts" if layer is None else layer,
                output_layer,
                log_transformed=log_transform,
                save_dir=save_dir,
            )

        except Exception as e:
            log.warning(f"Failed to generate normalization plots: {str(e)}")

    return adata




def regress_out(
    adata: AnnData,
    keys: List[str],
    layer: Optional[str] = "log1p_norm",
    n_jobs: Optional[int] = None,
    output_layer: str = "regressed_out",
    force: bool = False,
    calculate_variance_explained: bool = False,
) -> AnnData:
    """
    Regress out unwanted sources of variation from a specified layer.

    This function performs linear regression to remove the effect of specific
    variables (like cell cycle scores, mitochondrial percentage, etc.) from
    gene expression data.

    Args:
        adata: AnnData object
        keys: Variables to regress out (must be in adata.obs)
        layer: Layer to use as input for regression.
        n_jobs: Number of parallel jobs to use
        output_layer: Layer to store the regressed-out data.
        force: Whether to overwrite existing output_layer if it exists

    Returns:
        AnnData with regressed out variables

    Raises:
        ValueError: If keys are not found in adata.obs or layer not in adata.layers

    Examples:
        >>> # Regress out cell cycle effects
        >>> adata = regress_out(adata, keys=['S_score', 'G2M_score'])
        >>>
        >>> # Regress out multiple technical factors
        >>> adata = regress_out(adata, keys=['pct_counts_mt', 'total_counts'])
    """
    # Parameter validation
    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    # Check if keys exist in adata.obs
    missing_keys = [key for key in keys if key not in adata.obs]
    if missing_keys:
        raise ValueError(f"Keys not found in adata.obs: {', '.join(missing_keys)}")

    if layer not in adata.layers and layer is not None:
        raise ValueError(f"Layer '{layer}' not found in adata.layers.")# First, determine the name of the data source
    
    source_name = f"layer '{layer}'" if layer else "adata.X"
    
    log.info(f"Regressing out: {', '.join(keys)} from {source_name}")
    log.info(f"Using {n_jobs if n_jobs else 'all available'} processor cores.")

    # Calculate variance explained by each factor before regression
    var_explained = {}

    # Use the context manager to safely handle layers
    try:
        with use_layer_as_X(adata, layer):
            # Optionally calculate variance explained by each factor
            # (This is computationally expensive, so consider adding a flag to control it)
            if calculate_variance_explained:  # Set to True to enable variance analysis
                log.info("Calculating variance explained by each factor (this may be slow)...")
                from sklearn.linear_model import LinearRegression

                # Convert sparse matrix to dense if needed
                X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X

                # Calculate total variance for each gene
                gene_total_var = np.var(X, axis=0)

                for key in keys:
                    # Extract the covariate
                    covar = adata.obs[key].values.reshape(-1, 1)

                    # For each gene, calculate variance explained
                    gene_var_explained = []
                    for j in range(X.shape[1]):
                        if gene_total_var[j] > 0:  # Skip genes with no variance
                            model = LinearRegression().fit(covar, X[:, j])
                            y_pred = model.predict(covar)
                            explained_var = np.var(y_pred) / gene_total_var[j]
                            gene_var_explained.append(explained_var)
                        else:
                            gene_var_explained.append(0)

                    # Average across genes
                    var_explained[key] = np.mean(gene_var_explained)

                log.info("Variance explained by each factor:")
                for key, val in var_explained.items():
                    log.info(f"  - {key}: {val:.2%}")

            # The regress_out function modifies adata.X in place
            sc.pp.regress_out(adata, keys=keys, n_jobs=n_jobs)

            # Store the result from the modified adata.X into the output layer
            adata.layers[output_layer] = adata.X.copy()
    except Exception as e:
        log.error(f"Regression failed: {str(e)}")
        raise RuntimeError(f"Failed to regress out variables: {str(e)}")

    # Store regression information in uns
    if "regression" not in adata.uns:
        adata.uns["regression"] = {}

    adata.uns["regression"][output_layer] = {
        "regressed_variables": keys,
        "input_layer": layer,
        "variance_explained": var_explained,
    }

    log.info(
        f"Regression complete. Regressed data stored in adata.layers['{output_layer}']"
    )

    return adata
