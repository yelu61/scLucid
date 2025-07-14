"""
Doublet detection for single-cell RNA-seq data.

This module provides functions for identifying potential doublet cells
using the Scrublet algorithm and custom filtering approaches.
"""

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr

__all__ = ["is_doublet"]


def is_doublet(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    rate: float = 0.1,
    n_pcs: int = 30,
    threshold: float = 0.2,
    over_genes: float = 0.99,
    plot_umap: bool = True,
    save_dir: Optional[str] = None,
    show: bool = True,
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
        print(f"Processing sample {i + 1}/{total_samples}: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        # Get the current sample's cell count and feature count
        n_cells, n_features = data.shape
        # Dynamically set n_pcs parameter
        actual_n_pcs = min(n_pcs, n_cells - 1, n_features - 1)

        # Initialize arrays for doublet scores and predictions
        doublet_scores = np.full(data.shape[0], np.nan)
        predicted_doublets = np.full(data.shape[0], False)
        final_doublets = np.full(data.shape[0], False)

        # Check if the sample has enough cells for doublet detection
        if n_cells < 10:
            print(
                f"  Warning: Sample {sample} has fewer than 10 cells ({n_cells}). Skipping doublet detection."
            )
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

            print(
                f"  Doublet detection complete. Found {sum(final_doublets)} potential doublets."
            )

            if plot_umap:
                try:
                    scrub.set_embedding(
                        "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
                    )
                    fig = scrub.plot_embedding("UMAP", order_points=True)

                    if save_dir:
                        import os

                        os.makedirs(save_dir, exist_ok=True)
                        plt.savefig(
                            os.path.join(save_dir, f"{sample}_doublets_umap.png"),
                            dpi=300,
                        )
                    if show:
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
        print(
            f"  Completed {i + 1}/{total_samples} samples ({(i + 1) / total_samples * 100:.1f}%)"
        )

    # Identify cells with overexpressed genes as potential doublets
    top_genes = np.quantile(adata.obs.n_genes_by_counts, over_genes)
    adata.obs["overexpressed_doublets"] = adata.obs["n_genes_by_counts"] > top_genes

    # Print the overall statistics for the entire adata object
    total_cells = adata.n_obs
    print("\nOverall statistics for the entire adata object:")
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
