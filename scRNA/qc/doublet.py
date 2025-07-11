"""
Doublet detection for single-cell RNA-seq data.

This module provides functions for identifying potential doublet cells
using the Scrublet algorithm and custom filtering approaches.
"""

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import scrublet as scr
from typing import Dict, List, Optional

__all__ = [
    "is_doublet", 
    "filter_doublets"
]

def is_doublet(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    rate: float = 0.1,
    n_pcs: int = 30,
    threshold: float = 0.2,
    over_genes: float = 0.99,
    plot_umap: bool = True,
    save_dir: Optional[str] = None,
) -> sc.AnnData:
    """
    Identify and plot potential doublet cells using the scrublet package.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        rate (float, optional): Expected doublet rate. Defaults to 0.1.
        n_pcs (int, optional): Number of principal components to use. Defaults to 30.
        threshold (float, optional): Threshold for calling doublets. Defaults to 0.2.
        over_genes (float, optional): Quantile threshold for overexpressed genes. Defaults to 0.99.
        plot_umap (bool, optional): Whether to plot UMAP embedding with doublet scores. Defaults to True.
        save_dir (str, optional): Directory to save plots. If None, plots are not saved to disk. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with doublet scores and predictions added.
    """
    import gc
    
    adata.obs["doublet_scores"] = 0.0
    adata.obs["predicted_doublets"] = False
    adata.obs["predicted_doublets_final"] = False
    adata.obs["overexpressed_doublets"] = False

    # Get total number of samples for progress reporting
    samples = adata.obs[sample_key].unique()
    total_samples = len(samples)
    
    for i, sample in enumerate(samples):
        print(f"Processing sample {i+1}/{total_samples}: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        # Get the current sample's cell count and feature count
        n_cells, n_features = data.shape
        # Dynamically set n_pcs parameter
        actual_n_pcs = min(n_pcs, n_cells-1, n_features-1)
        
        # Initialize arrays for doublet scores and predictions
        doublet_scores = np.full(data.shape[0], np.nan)
        predicted_doublets = np.full(data.shape[0], False)
        final_doublets = np.full(data.shape[0], False)
        
        # Check if the sample has enough cells for doublet detection
        if n_cells < 10:
            print(f"  Warning: Sample {sample} has fewer than 10 cells ({n_cells}). Skipping doublet detection.")
            adata.obs.loc[data.obs.index, "doublet_scores"] = doublet_scores
            adata.obs.loc[data.obs.index, "predicted_doublets"] = predicted_doublets
            adata.obs.loc[data.obs.index, "predicted_doublets_final"] = final_doublets
            continue

        try:
            print(f"  Running Scrublet with {actual_n_pcs} PCs...")
            scrub = scr.Scrublet(data.X, expected_doublet_rate=rate)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(
                verbose=False, n_prin_comps=actual_n_pcs
            )
            final_doublets = scrub.call_doublets(threshold=threshold)
            
            print(f"  Doublet detection complete. Found {sum(final_doublets)} potential doublets.")
            
            if plot_umap:
                try:
                    scrub.set_embedding(
                        "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
                    )
                    fig = scrub.plot_embedding("UMAP", order_points=True)
                    
                    if save_dir:
                        import os
                        os.makedirs(save_dir, exist_ok=True)
                        plt.savefig(os.path.join(save_dir, f"{sample}_doublets_umap.png"), dpi=300)
                    plt.show()
                    plt.close()
                except Exception as e:
                    print(f"  Warning: Could not generate UMAP for doublets: {e}")
                    
        except Exception as e:
            print(f"  Error: Scrublet failed for sample {sample}: {e}")
            predicted_doublets = np.full(data.shape[0], False)
            final_doublets = np.full(data.shape[0], False)
            doublet_scores = np.full(data.shape[0], np.nan)
            
        adata.obs.loc[data.obs.index, "doublet_scores"] = doublet_scores
        adata.obs.loc[data.obs.index, "predicted_doublets"] = predicted_doublets
        adata.obs.loc[data.obs.index, "predicted_doublets_final"] = final_doublets

        # Force garbage collection
        gc.collect()
        
        # Print progress
        print(f"  Completed {i+1}/{total_samples} samples ({(i+1)/total_samples*100:.1f}%)")
        
    # Identify cells with overexpressed genes as potential doublets
    top_genes = np.quantile(adata.obs.n_genes_by_counts, over_genes)
    adata.obs["overexpressed_doublets"] = adata.obs["n_genes_by_counts"] > top_genes

    # Print the overall statistics for the entire adata object
    total_cells = adata.n_obs
    print(f"\nOverall statistics for the entire adata object:")
    print(f"Total number of cells: {total_cells}")

    potential_doublets = adata.obs["predicted_doublets"].sum()
    print(
        f"Potential doublet cells (based on doublet scores): {potential_doublets} ({potential_doublets / total_cells * 100:.2f}%)"
    )

    final_doublets = adata.obs["predicted_doublets_final"].sum()
    print(
        f"Potential doublet cells (based on threshold {threshold}): {final_doublets} ({final_doublets / total_cells * 100:.2f}%)"
    )

    overexpressed = adata.obs["overexpressed_doublets"].sum()
    print(
        f"Potential doublet cells (based on detected genes > {top_genes:.1f}): {overexpressed} ({overexpressed / total_cells * 100:.2f}%)"
    )

    combined_doublets = adata.obs.filter(regex="doublets").sum(axis=1).value_counts()
    for n_doublets, count in sorted(combined_doublets.items()):
        print(
            f"{n_doublets} types of doublets: {count} ({count / total_cells * 100:.2f}%)"
        )

    return adata

def filter_doublets(
    adata: sc.AnnData,
    only_predicted_final: bool = True,
    include_overexpressed: bool = False,
    min_score: Optional[float] = None,
    copy: bool = False
) -> sc.AnnData:
    """
    Filter out doublet cells from the AnnData object.
    
    This function removes cells identified as doublets by the is_doublet function.
    
    Args:
        adata (AnnData): AnnData object with doublet information (use is_doublet first)
        only_predicted_final (bool, optional): If True, only filter cells marked as 
                                              final predicted doublets. Defaults to True.
        include_overexpressed (bool, optional): If True, also filter cells with 
                                               overexpressed genes. Defaults to False.
        min_score (float, optional): If provided, filter cells with doublet scores 
                                    above this threshold. Defaults to None.
        copy (bool, optional): If True, return a filtered copy of the AnnData object.
                              If False, filter the AnnData object in place. Defaults to False.
    
    Returns:
        AnnData: Filtered AnnData object with doublet cells removed
        
    Example:
        >>> adata = sc.read_h5ad("data.h5ad")
        >>> adata = scRNA.qc.is_doublet(adata, sample_key="batch")
        >>> adata_filtered = scRNA.qc.filter_doublets(
        ...     adata, only_predicted_final=True, include_overexpressed=True
        ... )
    """
    # Check if doublet detection has been run
    required_cols = ["doublet_scores", "predicted_doublets", "predicted_doublets_final"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(
            f"Missing required columns in adata.obs: {missing_cols}. "
            f"Run is_doublet() first to identify doublets."
        )
    
    # Create a mask for cells to keep (non-doublets)
    keep_cells = np.ones(adata.n_obs, dtype=bool)
    
    # Filter based on final predicted doublets
    if only_predicted_final:
        keep_cells &= ~adata.obs["predicted_doublets_final"].values
        filter_type = "final predicted doublets"
    else:
        keep_cells &= ~adata.obs["predicted_doublets"].values
        filter_type = "all predicted doublets"
    
    # Optionally filter based on overexpressed genes
    if include_overexpressed:
        if "overexpressed_doublets" not in adata.obs.columns:
            print("Warning: 'overexpressed_doublets' column not found. Skipping this filter.")
        else:
            keep_cells &= ~adata.obs["overexpressed_doublets"].values
            filter_type += " and overexpressed genes"
    
    # Optionally filter based on doublet score threshold
    if min_score is not None:
        if min_score < 0 or min_score > 1:
            raise ValueError(f"min_score must be between 0 and 1, got {min_score}")
        
        keep_cells &= (adata.obs["doublet_scores"] < min_score).values
        filter_type += f" and score threshold {min_score}"
    
    # Get statistics before filtering
    cells_before = adata.n_obs
    cells_to_remove = np.sum(~keep_cells)
    
    # Create output object
    if copy:
        adata_out = adata[keep_cells].copy()
    else:
        adata_out = adata[keep_cells]
    
    # Print filtering statistics
    cells_after = adata_out.n_obs
    
    print(f"Doublet filtering based on {filter_type}:")
    print(f"  Cells before filtering: {cells_before}")
    print(f"  Cells after filtering: {cells_after}")
    print(f"  Removed doublets: {cells_to_remove} ({cells_to_remove/cells_before:.2%})")
    
    # Create a summary breakdown of removed cells
    if cells_to_remove > 0:
        if only_predicted_final:
            n_final = np.sum(adata.obs["predicted_doublets_final"])
            print(f"  - Final predicted doublets: {n_final} ({n_final/cells_before:.2%})")
        else:
            n_predicted = np.sum(adata.obs["predicted_doublets"])
            print(f"  - All predicted doublets: {n_predicted} ({n_predicted/cells_before:.2%})")
        
        if include_overexpressed and "overexpressed_doublets" in adata.obs.columns:
            n_overexp = np.sum(adata.obs["overexpressed_doublets"])
            print(f"  - Overexpressed gene doublets: {n_overexp} ({n_overexp/cells_before:.2%})")
        
        if min_score is not None:
            n_score = np.sum(adata.obs["doublet_scores"] >= min_score)
            print(f"  - Score-based doublets (>={min_score}): {n_score} ({n_score/cells_before:.2%})")
    
    # Optional visualization of filtering results
    try:
        import matplotlib.pyplot as plt
        
        if cells_to_remove > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Plot doublet score distribution
            ax.hist(adata.obs["doublet_scores"], bins=50, alpha=0.7, 
                   label=f"All cells (n={cells_before})")
            
            # Highlight removed cells
            if min_score is not None:
                ax.axvline(x=min_score, color='r', linestyle='--', 
                          label=f"Score threshold ({min_score})")
            
            removed_scores = adata.obs.loc[~keep_cells, "doublet_scores"]
            if not removed_scores.empty:
                ax.hist(removed_scores, bins=30, alpha=0.7, color='red',
                       label=f"Removed doublets (n={cells_to_remove})")
            
            ax.set_xlabel("Doublet score")
            ax.set_ylabel("Number of cells")
            ax.set_title("Distribution of doublet scores")
            ax.legend()
            
            plt.tight_layout()
            plt.show()
    except Exception as e:
        print(f"Couldn't generate visualization: {e}")
    
    return adata_out