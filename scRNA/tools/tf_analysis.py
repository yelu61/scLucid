import infercnvpy as cnv
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Union
import os

def find_tumor(
    adata: ad.AnnData,
    alpha: int = 2,
    key_added: str = "tumor",
    copy: bool = False,
) -> ad.AnnData:
    """
    Identify tumor cells based on the "cnv_score" column in adata.obs.

    Args:
        adata (AnnData): AnnData object with a "cnv_score" column in adata.obs.
        alpha (float, optional): Number of standard deviations to use for identifying tumor cells. Defaults to 2.0.
        key_added (str, optional): Name of the column to be added to adata.obs indicating tumor cells. Defaults to "tumor".
        copy (bool, optional): Whether to return a copy of adata or modify the original object. Defaults to False.

    Returns:
        adata (AnnData): AnnData object with an added column indicating tumor cells.

    Raises:
        ValueError: If "cnv_score" is not present in adata.obs_keys().
    """
    if "cnv_score" not in adata.obs_keys(): 
        raise ValueError("cnv_score not in adata.obs_names, please run infercnvpy first")
    
    cnv_score = np.sort(adata.obs["cnv_score"].unique())
    tumor_threshold_index = np.diff(cnv_score).argmax() + 1
    tumor_cnv_scores = cnv_score[tumor_threshold_index:]
    
    if len(tumor_cnv_scores) == 1:
        min_tumor_cnv_score = tumor_cnv_scores[0]
    else:
        lower_tumor_cnv_score_threshold = tumor_cnv_scores.mean() - alpha * tumor_cnv_scores.std()
        min_tumor_cnv_score = tumor_cnv_scores[tumor_cnv_scores >= lower_tumor_cnv_score_threshold].min()
    
    adata.obs[key_added] = adata.obs["cnv_score"].map(lambda x: 1 if x >= min_tumor_cnv_score else 0)
    
    return adata if copy else adata

def run_cnv_analysis(
    adata: ad.AnnData,
    sample_key: str = "sampleID",
    ref_obs: str = "Main_celltype",
    ref_keys: Union[str, list] = "Immune",
    wins: int = 250,
    step: int = 1,
    plot_heatmap: bool = True,
    heatmap_groupby: str = "celltype", 
    plot_umap: bool = True,
    plot_tumor: bool = True,
    figsize: tuple = (12, 3),
    #save_dir: Optional[str] = None,
) -> ad.AnnData:
    """
    Perform copy number variation (CNV) analysis on single-cell data using infercnvpy.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        ref_obs (str, optional): Reference column for CNV analysis. Defaults to "Main_celltype".
        ref_keys (Union[str, list], optional): Cell type(s) to use as reference. Defaults to "Immune".
        wins (int, optional): Window size for CNV analysis. Defaults to 250.
        step (int, optional): Step size for CNV analysis. Defaults to 1.
        plot_heatmap (bool, optional): Whether to plot chromosome heatmap. Defaults to True.
        heatmap_groupby (str, optional): Column name in adata.obs to use for grouping cells. Defaults to "celltype".
        plot_umap (bool, optional): Whether to plot UMAP embedding with CNV scores and clusters. Defaults to True.
        plot_tumor (bool, optional): Whether to plot UMAP embedding with predicted tumor cells. Defaults to True.
        figsize (tuple, optional): Figure size for UMAP plots. Defaults to (12, 3).
        save_dir (str, optional): Directory path to save the generated plots. If None, plots will be shown but not saved. Defaults to None.

    Returns:
        adata (AnnData): AnnData object with CNV analysis results added.
    """    
    adata.obs["cnv_score"] = 0.0
    adata.obs["is_tumor"] = 0

    if isinstance(ref_keys, str):
        ref_keys = [ref_keys]
        
    for sample in adata.obs[sample_key].unique():
        print(f"Begin of CNV analysis for sample: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        cnv.tl.infercnv(
            data,
            reference_key=ref_obs,
            reference_cat=ref_keys,
            window_size=wins,
            key_added="cnv",
            step=step,
        )

        if plot_heatmap:
            cnv.pl.chromosome_heatmap(data, groupby=heatmap_groupby, save=f"cnv_heatmap_{sample}.png")
            cnv.tl.pca(data)
            cnv.pp.neighbors(data)
            cnv.tl.leiden(data)
            cnv.pl.chromosome_heatmap(data, groupby="cnv_leiden", dendrogram=True,
                                      save=f"cnv_heatmap_leiden_{sample}.png")
            cnv.tl.umap(data)
            cnv.tl.cnv_score(data)

        if plot_umap:
            fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize)
            cnv.pl.umap(
                data,
                color="cnv_leiden",
                legend_loc="on data",
                legend_fontoutline=2,
                ax=axes[0],
                show=False,
            )
            axes[0].set_title("UMAP (Leiden Clusters)")
            cnv.pl.umap(data, color="cnv_score", ax=axes[1], show=False)
            axes[1].set_title("UMAP (CNV Score)")
            cnv.pl.umap(data, color=heatmap_groupby, ax=axes[2], show=False)
            axes[2].set_title(f"UMAP ({heatmap_groupby})")
            plt.tight_layout()
            plt.savefig(f"umap_cnv_{sample}.png", bbox_inches="tight")

        data = find_tumor(data, key_added="is_tumor")
        
        if plot_tumor:
            cnv.pl.umap(data, color="is_tumor", save=f"umap_tumor_{sample}.png")

        adata.obs.loc[data.obs.index, "cnv_score"] = data.obs["cnv_score"]
        adata.obs.loc[data.obs.index, "is_tumor"] = data.obs["is_tumor"]

    return adata

# Example usage:
# adata = run_cnv_analysis(adata, ref_keys=["Immune", "Stromal"])