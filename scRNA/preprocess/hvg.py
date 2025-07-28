"""
Functions for identifying highly variable genes in single-cell RNA-seq data.

This module provides methods for detecting genes with high expression variance
across cells, which are typically the most informative for downstream analysis.
It includes standard methods from Scanpy as well as custom approaches that
consider batch effects and technical characteristics.
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from ..utils.utils import use_layer_as_X

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = ["annotate_hvg", "select_hvg", "evaluate_hvg_stability", "plot_hvg_metrics"]

# --- Helper Functions ---
def _get_sample_specific_genes(
    adata: AnnData,
    sample_key: str,
    n_specific_genes: int,
    layer: Optional[str] = None,
    method: str = "t-test",
) -> List[str]:
    """
    Identifies top genes specific to each sample using differential expression analysis.

    This is a helper for the 'custom' HVG method to find genes that might
    drive batch effects rather than shared biological variation.

    Args:
        adata: AnnData object
        sample_key: Key in adata.obs that identifies sample groups
        n_specific_genes: Number of specific genes to identify per sample
        layer: Layer to use for analysis. If None, uses adata.X
        method: Method for differential expression analysis

    Returns:
        List of sample-specific gene names

    Raises:
        ValueError: If method is invalid or analysis fails
    """
    # Ensure there are multiple samples to compare
    n_samples = adata.obs[sample_key].nunique()
    if n_samples <= 1:
        log.info(
            f"Only one sample group found for '{sample_key}'. "
            "Skipping sample-specific gene identification."
        )
        return []

    log.info(f"Identifying sample-specific genes across {n_samples} groups...")

    try:
        # Use a copy to avoid modifying the original object's .uns field
        temp_adata = adata.copy()

        # Use the specified layer if provided
        if layer is not None:
            if layer not in temp_adata.layers:
                raise KeyError(f"Layer '{layer}' not found in AnnData object")
            temp_adata.X = temp_adata.layers[layer].copy()

        # rank_genes_groups requires log-transformed data
        # Check if data is already log-transformed
        if layer is None and "log1p" not in str(temp_adata.X[:5]):  # Heuristic check
            log.info("Log-transforming data for differential expression analysis")
            sc.pp.log1p(temp_adata)

        # Validate method
        valid_methods = ["t-test", "wilcoxon", "logreg"]
        if method not in valid_methods:
            log.warning(
                f"Method '{method}' not in {valid_methods}. Defaulting to 't-test'."
            )
            method = "t-test"

        # Run differential expression analysis
        sc.tl.rank_genes_groups(
            temp_adata,
            groupby=sample_key,
            method=method,
            n_genes=n_specific_genes,
            pts=True,  # Include fraction of cells expressing gene
        )

        # Extract gene names from the result
        specific_genes_df = pd.DataFrame(temp_adata.uns["rank_genes_groups"]["names"])
        specific_genes = list(np.unique(specific_genes_df.values.flatten()))
        log.info(f"Identified {len(specific_genes)} sample-specific genes")

        return specific_genes

    except Exception as e:
        log.error(f"Could not identify sample-specific genes: {str(e)}")
        log.exception("Detailed error:")
        return []


def _detect_gene_types(var_names: pd.Index) -> Dict[str, np.ndarray]:
    """
    Detects common gene types based on naming patterns.

    Args:
        var_names: Index of gene names to analyze

    Returns:
        Dictionary with gene type masks
    """
    # Convert to strings in case of other index types
    gene_names = var_names.astype(str)

    # Create masks for different gene types
    gene_types = {
        "mitochondrial": gene_names.str.match(r"^(MT-|mt-|MT\.|mt\.)")
        | gene_names.str.match(r"^(MTRNR|MTATP|MTND|MTCO|MTCYB)"),
        "ribosomal": gene_names.str.match(r"^(RP[SL]|Rp[sl])"),
        "hemoglobin": gene_names.str.match(r"^HB[^P]")
        | gene_names.str.contains(r"^hemoglobin", case=False),
        "heat_shock": gene_names.str.match(r"^HSP")
        | gene_names.str.contains(r"^heat shock", case=False),
        "immediate_early": np.isin(
            gene_names, ["FOS", "JUN", "JUNB", "EGR1", "NR4A1", "ZFP36"]
        ),
    }

    return gene_types


def _exclude_genes(
    adata: AnnData,
    hvg_mask: np.ndarray,
    exclude_types: List[str] = ["mitochondrial", "ribosomal"],
    gene_types: Optional[Dict[str, np.ndarray]] = None,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Exclude specific gene types from HVG selection.

    Args:
        adata: AnnData object
        hvg_mask: Boolean mask of current HVGs
        exclude_types: List of gene types to exclude
        gene_types: Pre-computed gene type masks (optional)

    Returns:
        Updated HVG mask and counts of excluded genes by type
    """
    # Detect gene types if not provided
    if gene_types is None:
        gene_types = _detect_gene_types(adata.var_names)

    # Initialize counters
    excluded_counts = {}

    # Create a copy of the mask to modify
    updated_mask = hvg_mask.copy()

    # Process each gene type to exclude
    for gene_type in exclude_types:
        if gene_type in gene_types:
            type_mask = gene_types[gene_type]
            # Count how many HVGs are of this type
            excluded_counts[gene_type] = (hvg_mask & type_mask).sum()
            # Remove these genes from HVGs
            updated_mask = updated_mask & ~type_mask
            log.info(f"Excluded {excluded_counts[gene_type]} {gene_type} genes")
        else:
            log.warning(f"Gene type '{gene_type}' not recognized")

    return updated_mask, excluded_counts


# --- Main Functions ---


def annotate_hvg(
    adata: AnnData,
    method: Literal["scanpy", "custom", "triku", "cell_ranger"] = "scanpy",
    # General Parameters
    layer: Optional[str] = "log1p_norm",
    n_top_genes: Optional[int] = 2000,
    output_key: Optional[str] = None,
    # Scanpy-specific args
    flavor: Literal[
        "seurat", "seurat_v3", "cell_ranger", "pearson_residuals"
    ] = "seurat_v3",
    batch_key: Optional[str] = None,
    min_disp: Optional[float] = 0.5,
    max_disp: Optional[float] = None,
    min_mean: Optional[float] = 0.0125,
    max_mean: Optional[float] = 3,
    span: Optional[float] = 0.3,
    # Custom-specific args
    sample_key: str = "sampleID",
    min_n_samples: int = 2,
    n_highly_expressed_genes: int = 50,
    #find_specific_genes: bool = True,
    n_specific_genes: int = 20,
    # Common Parameters
    exclude_gene_types: Optional[List[str]] = ["mitochondrial", "ribosomal"],
    plot: bool = False,
    save_path: Optional[str] = None,
    copy: bool = False,
    force: bool = False,
) -> AnnData:
    """
    Annotates highly variable genes (HVGs) using different methods.

    This function identifies genes with high expression variance across cells,
    which are typically the most informative for downstream analysis. It adds
    a boolean column to `adata.var` indicating HVGs.

    Args:
        adata: AnnData object with gene expression data.
        method: Method for HVG selection:
            - "scanpy": Use Scanpy's built-in methods (flavor parameter controls specifics)
            - "custom": Multi-criteria approach considering sample specificity
            - "triku": Uses k-nearest neighbors graph to identify informative genes
            - "cell_ranger": Approach used by 10x Genomics Cell Ranger
        
    --- General Parameters ---
        layer: Layer to use for HVG detection. If None, uses adata.X.
        n_top_genes: Number of top HVGs to select. Behavior depends on the method.
        output_key: Key to store results in adata.var. If None, auto-generated.

    --- Parameters for method='scanpy' or 'cell_ranger' ---
        flavor: The specific algorithm for 'scanpy': 'seurat', 'seurat_v3', or 'cell_ranger'.
            - "seurat": Original Seurat approach
            - "seurat_v3": Improved approach from Seurat v3
            - "cell_ranger": Similar to 10x Genomics approach
            - "pearson_residuals": Uses Pearson residuals variance
        batch_key: For 'scanpy' method, batch key for batch-aware HVG selection.
        min_disp, max_disp, min_mean, max_mean, span: Filtering and smoothing parameters.
            - min_disp, max_disp: Min/max dispersion thresholds for filtering.
            - min_mean, max_mean: Min/max mean expression thresholds for filtering.
            - span: Span parameter for LOWESS smoothing in dispersion calculation.

    --- Parameters for method='custom' ---
        sample_key: Key in adata.obs for per-sample HVG selection.
        min_n_samples: A gene must be an HVG in at least this many samples to be kept.
        n_highly_expressed_genes: Number of top expressed genes to exclude as noise.
        find_specific_genes: If True, find and exclude sample-specific (batch) genes.
                             Set to False to speed up the process if batch effects are minimal.
        n_specific_genes: Number of sample-specific genes to exclude per sample.

    --- Common Parameters ---
        exclude_gene_types: Types of genes to exclude from the final HVG list.
        plot: Whether to plot dispersion vs. mean expression.
        save_path: Path to save the plot.
        copy: Whether to return a copy of the AnnData object.
        force: Whether to overwrite existing annotations.

    Returns:
        AnnData object with HVG annotations in `adata.var`.

    Raises:
        ValueError: If method or parameters are invalid.
        KeyError: If required keys are not in adata.obs or adata.layers.

    Examples:
        >>> # Standard approach using Scanpy's method
        >>> adata = annotate_hvg(adata, method="scanpy", flavor="seurat_v3")
        >>>
        >>> # Custom approach considering batch effects
        >>> adata = annotate_hvg(
        ...     adata,
        ...     method="custom",
        ...     sample_key="batch",
        ...     min_n_samples=2
        ... )
    """
    # Handle copy
    if copy:
        adata = adata.copy()

    # Determine output key
    if output_key is None:
        output_key = f"highly_variable_{method}"
        if method == "scanpy":
            output_key = f"{output_key}_{flavor}"

    # Check if output already exists
    if output_key in adata.var and not force:
        n_existing = adata.var[output_key].sum()
        log.info(
            f"HVG annotations already exist under '{output_key}' with {n_existing} genes. "
            f"Use force=True to overwrite."
        )
        return adata

    # Check if layer exists
    if layer is not None and layer not in adata.layers:
        available = list(adata.layers.keys())
        raise KeyError(
            f"Layer '{layer}' not found in adata.layers. Available: {available}"
        )

    # Beginning log message
    log.info(f"Identifying highly variable genes using '{method}' method")
    if layer:
        log.info(f"Using expression data from layer: '{layer}'")
    else:
        log.info("Using expression data from adata.X")

    # Main processing block
    with use_layer_as_X(adata, layer):
        if method == "scanpy":
            log.info(f"Running HVG selection with 'scanpy' method (flavor: {flavor})")

            # Validate parameters
            if flavor in ["seurat_v3", "pearson_residuals"] and n_top_genes is None:
                raise ValueError(f"flavor='{flavor}' requires `n_top_genes` to be set")

            # Batch-aware processing
            if batch_key:
                if batch_key not in adata.obs:
                    raise KeyError(f"Batch key '{batch_key}' not found in adata.obs")
                log.info(
                    f"Using batch-aware HVG selection with key '{batch_key}' "
                    f"({adata.obs[batch_key].nunique()} batches)"
                )

            # Process parameters
            hvg_kwargs = {
                "flavor": flavor,
                "n_top_genes": n_top_genes,
                "min_mean": min_mean,
                "max_mean": max_mean,
                "min_disp": min_disp,
                "max_disp": max_disp,
                "span": span,
                "batch_key": batch_key,
                "inplace": True,
            }

            # Run Scanpy's HVG detection
            try:
                sc.pp.highly_variable_genes(adata, **hvg_kwargs)

                # Rename the default output for clarity
                adata.var[output_key] = adata.var["highly_variable"].copy()

                # Keep the metrics
                for metric in ["means", "dispersions", "dispersions_norm"]:
                    if metric in adata.var:
                        adata.var[f"{output_key}_{metric}"] = adata.var[metric].copy()

                # Clean up temporary columns if needed
                if output_key != "highly_variable":
                    del adata.var["highly_variable"]

            except Exception as e:
                log.error(f"HVG selection failed: {str(e)}")
                log.exception("Detailed error:")
                raise RuntimeError(f"Failed to compute highly variable genes: {str(e)}")

        elif method == "custom":
            log.info("Running HVG selection with 'custom' multi-criteria method")

            # Validate parameters
            if sample_key not in adata.obs.columns:
                raise KeyError(f"Sample key '{sample_key}' not found in adata.obs")

            samples = adata.obs[sample_key].unique()
            log.info(f"Found {len(samples)} samples using key '{sample_key}'")

            # 1. Find HVGs per sample
            hvg_masks = []
            sample_counts = []

            for sample in samples:
                sample_mask = adata.obs[sample_key] == sample
                sample_size = np.sum(sample_mask)
                sample_counts.append(sample_size)

                if sample_size > 10:  # Min cells required
                    log.info(f"Processing sample '{sample}' with {sample_size} cells")
                    sample_adata = adata[sample_mask].copy()

                    # Use the same parameters as the scanpy method
                    sc.pp.highly_variable_genes(
                        sample_adata,
                        n_top_genes=n_top_genes,
                        min_mean=min_mean,
                        max_mean=max_mean,
                        min_disp=min_disp,
                        max_disp=max_disp,
                        span=span,
                        inplace=True,
                    )

                    hvg_masks.append(sample_adata.var["highly_variable"])
                else:
                    log.warning(
                        f"Sample '{sample}' has only {sample_size} cells, "
                        f"which is fewer than the minimum 10 required. Skipping."
                    )

            if not hvg_masks:
                raise ValueError("No samples had enough cells to compute HVGs")

            # 2. Combine HVGs across samples
            combined_df = pd.concat(hvg_masks, axis=1)
            combined_df.columns = [f"sample_{i}" for i in range(len(hvg_masks))]

            # Count in how many samples each gene is highly variable
            sample_counts = combined_df.sum(axis=1)

            # Apply minimum sample threshold
            combined_hvgs = sample_counts >= min_n_samples
            log.info(
                f"Found {combined_hvgs.sum()} genes considered HVG in at least "
                f"{min_n_samples} samples"
            )

            # Store the per-sample information for reference
            adata.var[f"{output_key}_sample_count"] = sample_counts

            # 3. Identify genes to exclude

            # 3a. Highly expressed genes
            gene_expr = np.array(adata.X.sum(axis=0)).flatten()
            top_expr_indices = np.argsort(-gene_expr)[:n_highly_expressed_genes]
            top_expr_genes = adata.var_names[top_expr_indices]

            log.info(
                f"Identified {len(top_expr_genes)} highly expressed genes to exclude"
            )

            # 3b. Sample-specific genes
            specific_genes = _get_sample_specific_genes(
                adata, sample_key, n_specific_genes, layer=layer
            )

            log.info(
                f"Identified {len(specific_genes)} sample-specific genes to exclude"
            )

            # 3c. Create exclude set
            exclude_genes = set(top_expr_genes) | set(specific_genes)
            exclude_mask = adata.var_names.isin(exclude_genes)

            log.info(
                f"Total genes to exclude: {len(exclude_genes)} "
                f"({exclude_mask.sum()} found in dataset)"
            )

            # 4. Final HVG list
            final_hvg_mask = combined_hvgs & ~exclude_mask
            adata.var[output_key] = final_hvg_mask

            # Store additional metrics
            adata.var[f"{output_key}_highly_expressed"] = adata.var_names.isin(
                top_expr_genes
            )
            adata.var[f"{output_key}_sample_specific"] = adata.var_names.isin(
                specific_genes
            )

        elif method == "triku":
            try:
                import triku
            except ImportError:
                log.error(
                    "Triku method requires the 'triku' package. "
                    "Install with 'pip install triku'"
                )
                raise ImportError("Please install triku: pip install triku")

            log.info("Running HVG selection with 'triku' method")

            try:
                # Use Triku to identify HVGs
                result = triku.tl.triku(
                    adata, n_pcs=50, n_neighbors=15, return_all=True
                )

                # Store results
                adata.var[output_key] = result["highly_variable"]
                adata.var[f"{output_key}_score"] = result["score"]

                log.info(
                    f"Triku identified {result['highly_variable'].sum()} highly variable genes"
                )

            except Exception as e:
                log.error(f"Triku HVG selection failed: {str(e)}")
                log.exception("Detailed error:")
                raise RuntimeError(f"Failed to compute HVGs with triku: {str(e)}")

        elif method == "cell_ranger":
            log.info("Running HVG selection with 'cell_ranger' method")

            # This is a wrapper around Scanpy's implementation with cell_ranger flavor
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    flavor="cell_ranger",
                    n_top_genes=n_top_genes,
                    batch_key=batch_key,
                    inplace=True,
                )

                # Rename the default output for clarity
                adata.var[output_key] = adata.var["highly_variable"].copy()

                # Keep the metrics
                for metric in ["means", "dispersions", "dispersions_norm"]:
                    if metric in adata.var:
                        adata.var[f"{output_key}_{metric}"] = adata.var[metric].copy()

                # Clean up temporary columns if needed
                if output_key != "highly_variable":
                    del adata.var["highly_variable"]

            except Exception as e:
                log.error(f"Cell Ranger HVG selection failed: {str(e)}")
                raise RuntimeError(
                    f"Failed to compute HVGs with Cell Ranger method: {str(e)}"
                )

        else:
            available_methods = ["scanpy", "custom", "triku", "cell_ranger"]
            raise ValueError(
                f"Unknown method '{method}'. Choose from: {available_methods}"
            )

    # Common exclusion step for all methods
    if exclude_gene_types:
        log.info(f"Excluding gene types: {exclude_gene_types}")

        # Detect gene types
        gene_types = _detect_gene_types(adata.var_names)

        # Get current mask and update it
        current_mask = adata.var[output_key]
        updated_mask, excluded_counts = _exclude_genes(
            adata, current_mask, exclude_gene_types, gene_types
        )

        # Update the mask
        adata.var[output_key] = updated_mask

        # Log exclusion results
        for gene_type, count in excluded_counts.items():
            log.info(f"Excluded {count} {gene_type} genes")

    # Final count
    final_count = adata.var[output_key].sum()
    log.info(f"Final number of highly variable genes ({method} method): {final_count}")

    # Store metadata about the HVG selection
    if "hvg" not in adata.uns:
        adata.uns["hvg"] = {}

    adata.uns["hvg"][output_key] = {
        "method": method,
        "params": {
            "n_top_genes": n_top_genes,
            "layer": layer,
            "flavor": flavor if method == "scanpy" else None,
            "batch_key": batch_key,
            "min_n_samples": min_n_samples if method == "custom" else None,
            "exclude_gene_types": exclude_gene_types,
        },
        "n_hvgs": final_count,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    # Create plots if requested
    if plot:
        try:
            fig = plot_hvg_metrics(adata, output_key, save_path=save_path)
        except Exception as e:
            log.warning(f"Failed to generate HVG plots: {str(e)}")

    return adata


def select_hvg(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_top_genes: Optional[int] = None,
    subset: bool = True,
    sort_by: Optional[str] = "dispersions_norm",
    keep_raw: bool = True,
    layer_to_subset: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Selects a final set of HVGs and optionally subsets the AnnData object.

    This function can be used to refine HVG selection or to create a smaller
    AnnData object containing only HVGs for downstream analysis.

    Args:
        adata: AnnData object with HVG annotations.
        hvg_key: The column in `adata.var` to use (e.g., 'highly_variable_scanpy').
        n_top_genes: If provided, select the top N genes based on the sort_by column.
        subset: If True, subsets `adata` to the selected HVGs and returns it.
                If False, adds a new boolean mask to `adata.var` and returns None.
        sort_by: Column to use for sorting genes when selecting top N.
                 Common options: 'dispersions_norm', 'means', 'dispersions'.
        keep_raw: Whether to store the full dataset in `.raw` when subsetting.
        layer_to_subset: If provided, only subset this layer (for memory efficiency)
                        rather than all layers.
        copy: Whether to return a copy even when subset=False.

    Returns:
        A subsetted AnnData object if `subset=True`, otherwise the original
        (or a copy if copy=True).

    Raises:
        KeyError: If hvg_key is not found or sort_by column doesn't exist.
        ValueError: If invalid parameters are provided.

    Examples:
        >>> # Select top 1000 genes from previously annotated HVGs
        >>> adata_hvg = select_hvg(adata, hvg_key="highly_variable_scanpy", n_top_genes=1000)
        >>>
        >>> # Just update the standard 'highly_variable' flag without subsetting
        >>> adata = select_hvg(adata, hvg_key="highly_variable_custom", subset=False)
    """
    # Parameter validation
    if hvg_key not in adata.var:
        raise KeyError(
            f"HVG key '{hvg_key}' not found in `adata.var`. Please run `annotate_hvg` first."
        )

    if n_top_genes is not None and n_top_genes <= 0:
        raise ValueError("n_top_genes must be a positive integer")

    if n_top_genes is not None and sort_by is not None and sort_by not in adata.var:
        # Look for method-specific column names
        method_specific = f"{hvg_key}_{sort_by}"
        if method_specific in adata.var:
            sort_by = method_specific
            log.info(f"Using method-specific sorting column: '{sort_by}'")
        else:
            available_metrics = [
                col
                for col in adata.var.columns
                if any(x in col for x in ["mean", "disp", "score"])
            ]
            raise KeyError(
                f"Sorting column '{sort_by}' not found in `adata.var`. "
                f"Available metrics: {available_metrics}"
            )

    if layer_to_subset is not None and layer_to_subset not in adata.layers:
        raise KeyError(f"Layer '{layer_to_subset}' not found in adata.layers")

    # Create a copy if needed
    if copy:
        adata = adata.copy()

    # Store raw data if subsetting and it doesn't already exist
    if subset and keep_raw and adata.raw is None:
        log.info("Storing original data in adata.raw before subsetting")
        adata.raw = adata.copy()

    log.info(f"Processing HVGs from '{hvg_key}' column")

    # Get the initial mask
    final_hvg_mask = adata.var[hvg_key].copy()
    initial_hvg_count = final_hvg_mask.sum()

    log.info(f"Initial HVG count: {initial_hvg_count}")

    # Refine selection if n_top_genes is provided
    if n_top_genes is not None:
        if initial_hvg_count == 0:
            log.warning("No HVGs found in the initial mask. Cannot select top genes.")
        elif sort_by is None:
            log.warning("No sort_by column provided. Selecting HVGs randomly.")
            hvg_indices = np.where(final_hvg_mask)[0]
            if len(hvg_indices) > n_top_genes:
                selected_indices = np.random.choice(
                    hvg_indices, n_top_genes, replace=False
                )
                new_mask = np.zeros_like(final_hvg_mask)
                new_mask[selected_indices] = True
                final_hvg_mask = new_mask
        else:
            # Sort genes by the specified metric, but only consider those already marked as HVG
            try:
                metrics = adata.var.loc[final_hvg_mask, sort_by]

                # Determine if we should sort ascending or descending
                ascending = True  # Default for most metrics
                if any(x in sort_by for x in ["disp", "var", "score"]):
                    ascending = False  # Higher dispersion/variance is better

                # Get top genes
                if ascending:
                    top_genes = metrics.nsmallest(n_top_genes).index
                else:
                    top_genes = metrics.nlargest(n_top_genes).index

                # Create a new mask with only the top N genes
                final_hvg_mask = adata.var_names.isin(top_genes)
                log.info(
                    f"Selected top {final_hvg_mask.sum()} HVGs based on '{sort_by}'"
                )

            except Exception as e:
                log.error(f"Failed to select top genes: {str(e)}")
                log.exception("Detailed error:")

    # Update adata.var['highly_variable'] which is the default for downstream tools
    adata.var["highly_variable"] = final_hvg_mask
    log.info(f"Updated adata.var['highly_variable'] with {final_hvg_mask.sum()} genes")

    # Subset the AnnData object if requested
    if subset:
        log.info("Subsetting AnnData object to selected HVGs")

        if layer_to_subset is not None:
            # Only subset a specific layer to save memory
            log.info(f"Subsetting only layer '{layer_to_subset}' and X for efficiency")

            # Create a view with only the selected genes
            subset_view = adata[:, adata.var["highly_variable"]]

            # Create a new AnnData with only the selected layer
            import scipy.sparse

            if scipy.sparse.issparse(subset_view.layers[layer_to_subset]):
                X_new = subset_view.layers[layer_to_subset].copy()
            else:
                X_new = subset_view.layers[layer_to_subset].copy()

            # Create a new AnnData with the subsetted data
            from anndata import AnnData as AnnData_class

            adata_subset = AnnData_class(
                X=X_new,
                obs=subset_view.obs.copy(),
                var=subset_view.var.copy(),
                uns=subset_view.uns.copy(),
            )

            # Add the raw if it exists
            if subset_view.raw is not None:
                adata_subset.raw = subset_view.raw.copy()

            return adata_subset
        else:
            # Standard subsetting of all data
            return adata[:, adata.var["highly_variable"]].copy()

    return adata  # Return the updated AnnData (original or copy)


def evaluate_hvg_stability(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_bootstrap: int = 20,
    sample_fraction: float = 0.8,
    method: str = "scanpy",
    flavor: str = "seurat",
    n_top_genes: Optional[int] = 2000,
    layer: Optional[str] = None,
    random_state: Optional[int] = 42,
    plot: bool = True,
    save_path: Optional[str] = None,
) -> AnnData:
    """
    Evaluates the stability of HVG selection through bootstrap resampling.

    This function repeatedly samples cells and recalculates HVGs to determine
    how consistently each gene is selected, providing a measure of confidence
    in the gene selection.

    Args:
        adata: AnnData object with gene expression data
        hvg_key: Key in adata.var containing the current HVG selection
        n_bootstrap: Number of bootstrap samples to generate
        sample_fraction: Fraction of cells to include in each bootstrap sample
        method: HVG selection method to use in each bootstrap
        flavor: Flavor for the scanpy method
        n_top_genes: Number of top genes to select in each bootstrap
        layer: Layer to use for expression data
        random_state: Seed for random number generator
        plot: Whether to generate visualization of stability results
        save_path: Path to save the plot

    Returns:
        AnnData with added stability metrics in .var

    Examples:
        >>> # Evaluate stability of previously selected HVGs
        >>> adata = evaluate_hvg_stability(adata, hvg_key="highly_variable")
        >>>
        >>> # Plot genes by stability
        >>> sc.pl.scatter(
        ...     adata,
        ...     x="hvg_selection_frequency",
        ...     y="dispersions_norm",
        ...     color="highly_variable"
        ... )
    """
    import random

    # Parameter validation
    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")

    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be a positive integer")

    if not 0 < sample_fraction < 1:
        raise ValueError("sample_fraction must be between 0 and 1")

    # Get current HVGs
    current_hvgs = set(adata.var_names[adata.var[hvg_key]])
    log.info(f"Evaluating stability of {len(current_hvgs)} HVGs")

    # Set random seed for reproducibility
    if random_state is not None:
        random.seed(random_state)
        np.random.seed(random_state)

    # Storage for results
    gene_selection_count = {gene: 0 for gene in adata.var_names}
    n_cells_per_bootstrap = int(adata.n_obs * sample_fraction)

    log.info(
        f"Performing {n_bootstrap} bootstrap iterations with "
        f"{n_cells_per_bootstrap} cells each ({sample_fraction:.1%} of data)"
    )

    # Progress reporting variables
    report_interval = max(1, n_bootstrap // 10)

    # Perform bootstrap iterations
    for i in range(n_bootstrap):
        if i % report_interval == 0:
            log.info(f"Bootstrap iteration {i + 1}/{n_bootstrap}")

        # Random subsample of cells
        cell_indices = np.random.choice(
            adata.n_obs, size=n_cells_per_bootstrap, replace=False
        )

        # Create subsampled data
        bootstrap_adata = adata[cell_indices].copy()

        # Run HVG selection on the bootstrap sample
        annotate_hvg(
            bootstrap_adata,
            method=method,
            flavor=flavor,
            n_top_genes=n_top_genes,
            layer=layer,
            output_key="bootstrap_hvg",
            force=True,
        )

        # Get HVGs from this bootstrap
        bootstrap_hvgs = set(
            bootstrap_adata.var_names[bootstrap_adata.var["bootstrap_hvg"]]
        )

        # Update counts
        for gene in bootstrap_hvgs:
            if gene in gene_selection_count:
                gene_selection_count[gene] += 1

    # Calculate selection frequency for each gene
    selection_frequency = {
        gene: count / n_bootstrap for gene, count in gene_selection_count.items()
    }

    # Add to adata.var
    adata.var["hvg_selection_frequency"] = pd.Series(
        [selection_frequency.get(gene, 0) for gene in adata.var_names],
        index=adata.var_names,
    )

    # Calculate stability metrics
    stability_score = np.mean(
        [selection_frequency.get(gene, 0) for gene in current_hvgs]
    )
    top_quartile = np.quantile(
        [selection_frequency.get(gene, 0) for gene in current_hvgs], 0.75
    )
    bottom_quartile = np.quantile(
        [selection_frequency.get(gene, 0) for gene in current_hvgs], 0.25
    )

    log.info("HVG stability metrics:")
    log.info(f"  - Overall stability score: {stability_score:.3f}")
    log.info(f"  - Top 25% of HVGs selected with frequency >= {top_quartile:.3f}")
    log.info(f"  - Bottom 25% of HVGs selected with frequency <= {bottom_quartile:.3f}")

    # Store stability metrics in uns
    adata.uns["hvg_stability"] = {
        "overall_score": stability_score,
        "top_quartile": top_quartile,
        "bottom_quartile": bottom_quartile,
        "n_bootstrap": n_bootstrap,
        "sample_fraction": sample_fraction,
        "method": method,
    }

    # Create visualization if requested
    if plot:
        try:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # Histogram of selection frequencies
            sns.histplot(
                adata.var["hvg_selection_frequency"], bins=30, kde=True, ax=axes[0]
            )
            axes[0].set_title("HVG Selection Frequency Distribution")
            axes[0].set_xlabel("Selection Frequency")
            axes[0].set_ylabel("Number of Genes")

            # Add vertical line for average stability
            axes[0].axvline(
                stability_score,
                color="red",
                linestyle="--",
                label=f"Avg HVG Stability: {stability_score:.3f}",
            )
            axes[0].legend()

            # Scatter plot of mean vs stability
            if "means" in adata.var:
                x = "means"
            elif f"{hvg_key}_means" in adata.var:
                x = f"{hvg_key}_means"
            else:
                # Calculate means if not available
                if layer is None:
                    adata.var["temp_means"] = np.array(adata.X.mean(axis=0)).flatten()
                else:
                    adata.var["temp_means"] = np.array(
                        adata.layers[layer].mean(axis=0)
                    ).flatten()
                x = "temp_means"

            scatter = axes[1].scatter(
                adata.var[x],
                adata.var["hvg_selection_frequency"],
                c=adata.var[hvg_key].astype(int),
                alpha=0.6,
                cmap="coolwarm",
                s=10,
            )

            axes[1].set_title("HVG Stability vs. Mean Expression")
            axes[1].set_xlabel("Mean Expression")
            axes[1].set_ylabel("Selection Frequency")

            # Add colorbar legend
            cbar = plt.colorbar(scatter, ax=axes[1])
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(["Not HVG", "HVG"])

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                log.info(f"Saved stability plot to {save_path}")

            # Clean up temporary column if created
            if "temp_means" in adata.var:
                del adata.var["temp_means"]

        except Exception as e:
            log.warning(f"Failed to create stability plot: {str(e)}")

    return adata


def plot_hvg_metrics(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_top_genes: int = 20,
    metrics: Optional[List[str]] = None,
    show_gene_labels: bool = True,
    size_by_expr: bool = True,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Creates visualizations of HVG metrics to evaluate selection quality.

    This function generates plots showing the relationship between various
    metrics used in HVG selection, helping to understand and validate the
    gene selection process.

    Args:
        adata: AnnData object with HVG annotations
        hvg_key: Key in adata.var containing HVG selection
        n_top_genes: Number of top genes to label in the plot
        metrics: Specific metrics to plot (defaults to auto-detection)
        show_gene_labels: Whether to show gene names on the plot
        size_by_expr: Whether to size points by mean expression
        save_path: Path to save the figure

    Returns:
        matplotlib Figure object

    Examples:
        >>> # Plot metrics for default HVG selection
        >>> fig = plot_hvg_metrics(adata)
        >>> plt.show()
        >>>
        >>> # Plot custom HVG selection with more labeled genes
        >>> fig = plot_hvg_metrics(
        ...     adata,
        ...     hvg_key="highly_variable_custom",
        ...     n_top_genes=50
        ... )
    """
    # Parameter validation
    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")

    # Auto-detect available metrics
    available_metrics = {}
    for column in adata.var.columns:
        if any(x in column for x in ["mean", "disp", "var", "score"]):
            # Determine what type of metric this is
            if "mean" in column:
                available_metrics["mean"] = column
            elif any(x in column for x in ["disp", "var"]):
                available_metrics["dispersion"] = column
            elif "norm" in column:
                available_metrics["norm_dispersion"] = column
            elif "score" in column:
                available_metrics["score"] = column

    # Use specific metrics if provided
    if metrics is not None:
        for metric in metrics:
            if metric not in adata.var.columns:
                log.warning(f"Metric '{metric}' not found in adata.var columns")

    # Determine which metrics to plot
    if "norm_dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["norm_dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "score" in available_metrics:
        x = available_metrics["mean"] if "mean" in available_metrics else None
        y = available_metrics["score"]
        plot_type = "score"
    else:
        # Fallback: check for method-specific columns
        method_specific_x = f"{hvg_key}_means"
        method_specific_y = f"{hvg_key}_dispersions_norm"

        if method_specific_x in adata.var and method_specific_y in adata.var:
            x = method_specific_x
            y = method_specific_y
            plot_type = "dispersion_vs_mean"
        else:
            log.warning("Could not find appropriate metrics for plotting")
            # Calculate basic metrics if needed
            adata.var["_temp_mean"] = np.array(adata.X.mean(axis=0)).flatten()
            x = "_temp_mean"

            # Use hvg_selection_frequency if available
            if "hvg_selection_frequency" in adata.var:
                y = "hvg_selection_frequency"
                plot_type = "stability"
            else:
                # Can't create a meaningful plot without metrics
                log.error("Insufficient metrics available for plotting")
                raise ValueError("Cannot create HVG plot: no appropriate metrics found")

    # Create the figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Determine point sizes if requested
    if size_by_expr and "mean" in available_metrics:
        sizes = np.clip(adata.var[available_metrics["mean"]] * 20, 5, 200)
    else:
        sizes = 30

    # Create the scatter plot
    scatter = ax.scatter(
        adata.var[x],
        adata.var[y],
        s=sizes,
        c=adata.var[hvg_key].astype(int),
        cmap="coolwarm",
        alpha=0.7,
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Not HVG", "HVG"])

    # Set appropriate labels based on plot type
    if plot_type == "dispersion_vs_mean":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Dispersion (normalized)")
        ax.set_title("Highly Variable Genes: Dispersion vs. Mean")
        ax.set_xscale("log")
    elif plot_type == "score":
        ax.set_xlabel("Mean Expression" if x else "Gene Index")
        ax.set_ylabel("HVG Score")
        ax.set_title("Highly Variable Genes: Score Distribution")
    elif plot_type == "stability":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Selection Frequency")
        ax.set_title("HVG Selection Stability")

    # Label top genes if requested
    if show_gene_labels and n_top_genes > 0:
        # Get indices of HVGs sorted by the y metric
        hvg_mask = adata.var[hvg_key]
        if hvg_mask.sum() > 0:
            top_indices = adata.var.loc[hvg_mask, y].nlargest(n_top_genes).index

            for idx in top_indices:
                gene_name = idx
                ax.annotate(
                    gene_name,
                    (adata.var.loc[idx, x], adata.var.loc[idx, y]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    bbox=dict(
                        boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8
                    ),
                )

    # Add count information to the title
    n_hvgs = adata.var[hvg_key].sum()
    total_genes = len(adata.var)
    ax.set_title(
        f"{ax.get_title()}\n{n_hvgs} HVGs selected ({n_hvgs / total_genes:.1%} of {total_genes} genes)"
    )

    # Save the figure if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved HVG metrics plot to {save_path}")

    # Clean up temporary columns if created
    if "_temp_mean" in adata.var:
        del adata.var["_temp_mean"]

    return fig
