"""
Quality control module for single-cell RNA-seq data.
"""

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr

__all__ = ["calculate_qc_metric", "identify_outliers", "is_low_quality_cell", "is_doublet"]

def calculate_qc_metric(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    keys: list = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ],
    plot_violin: bool = True,
    plot_scatter: bool = True,
):
    """
    Calculate and plot QC metrics for each sample in the AnnData object.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        keys (list, optional): List of QC metrics to plot. Defaults to config.QC_METRICS.
        plot_violin (bool, optional): Whether to plot violin plots for QC metrics. Defaults to True.
        plot_scatter (bool, optional): Whether to plot scatter plot for total_counts vs n_genes_by_counts. Defaults to True.

    Returns:
        adata (AnnData): AnnData object with QC metrics added.
    """
    adata.obs["total_counts"] = 0
    adata.obs["log1p_total_counts"] = 0
    adata.obs["n_genes_by_counts"] = 0
    adata.obs["log1p_n_genes_by_counts"] = 0
    adata.obs["pct_counts_in_top_20_genes"] = 0
    adata.obs["pct_counts_mt"] = 0.0
    adata.obs["pct_counts_ribo"] = 0.0
    adata.obs["pct_counts_hb"] = 0.0

    for sample in adata.obs[sample_key].unique():
        print(f"Begin of QC metric calculation and QC plot for sample: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        # Calculate the QC covariates or metric
        data.var["mt"] = data.var_names.str.contains(
            r"^(MT-|mt-)", regex=True, na=False
        )  # mitochondrial genes
        data.var["ribo"] = data.var_names.str.contains(
            r"^(RPS|RPL|Rps|Rpl|Gm\d+)", regex=True, na=False
        )  # ribosomal genes
        data.var["hb"] = data.var_names.str.contains(
            r"^(HB|Hb)[^(P|p)]", regex=True, na=False
        )  # hemoglobin genes
        sc.pp.calculate_qc_metrics(
            data,
            qc_vars=["mt", "ribo", "hb"],
            inplace=True,
            percent_top=[20],
            log1p=True,
        )

        # Plot the QC covariates per sample to assess how well the respective QC metric separates cell populations
        if plot_violin:
            fig, axes = plt.subplots(nrows=1, ncols=len(keys), figsize=(12, 3))
            for i, ax in enumerate(axes.flat):
                sc.pl.violin(data, keys[i], ax=ax, jitter=0.4, show=False)
                ax.set_title(f"Sample: {sample}", fontsize=10)
            plt.tight_layout()
            plt.show()

        if plot_scatter:
            sc.pl.scatter(
                data,
                x="total_counts",
                y="n_genes_by_counts",
                color="pct_counts_mt",
                title=f"Sample: {sample}",
                show=False,
                legend_loc="right margin",
            )
        
        metric_cols = [
            "total_counts",
            "log1p_total_counts",
            "n_genes_by_counts",
            "log1p_n_genes_by_counts",
            "pct_counts_in_top_20_genes",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
        ]
        adata.obs.loc[data.obs.index, metric_cols] = data.obs[metric_cols]
        print("Done.")

    return adata


def identify_outliers(adata: sc.AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Identify outliers based on the given metric and number of median absolute deviations.

    Args:
        adata (AnnData): AnnData object to check for outliers.
        metric (str): The metric to use for outlier detection. Must be a valid column in adata.obs.
        nmads (int): Number of median absolute deviations for outlier detection.

    Returns:
        outliers (pandas.Series): Boolean mask indicating if a cell is an outlier or not.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Invalid metric '{metric}'. Must be a column in adata.obs.")

    values = adata.obs[metric]
    median = np.median(values)
    mad = median_abs_deviation(values)
    outliers = [abs(value - median) > nmads * mad for value in values]
    return pd.Series(outliers, index=adata.obs_names)


def is_low_quality_cell(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    min_genes: int = 200,
    nmad: int = 5,
    pc_mt: int = 20,
    pc_hb: int = 20,
):
    """
    Identify and filter out low-quality cells based on gene counts, total counts, mitochondrial percentage, and hemoglobin percentage.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        sample_key (str, optional): The key in adata.obs to identify different samples. Defaults to "sampleID".
        min_genes (int, optional): Minimum number of genes expressed in a cell. Defaults to 200.
        nmad (int, optional): Number of median absolute deviations for outlier detection. Defaults to 5.
        pc_mt (int, optional): Maximum percentage of mitochondrial counts allowed. Defaults to 20.
        pc_hb (int, optional): Maximum percentage of hemoglobin counts allowed. Defaults to 20.

    Returns:
        adata (AnnData): AnnData object with low-quality cells filtered out.
    """
    
    # Filter out cells with low gene counts
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    adata.obs["outlier"] = False
    adata.obs["mt_outlier"] = False
    adata.obs["hb_outlier"] = False

    for sample in adata.obs[sample_key].unique():
        print(f"QC of low quality cells for sample: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]
        
        # Identify outlier cells
        data.obs["outlier"] = (
            identify_outliers(data, metric="log1p_total_counts", nmads=nmad)
            | identify_outliers(data, metric="log1p_n_genes_by_counts", nmads=nmad)
            | identify_outliers(data, metric="pct_counts_in_top_20_genes", nmads=nmad)
        )
        print(
            f"Outlier cells (based on total counts, gene counts, and top gene percentage): {data.obs.outlier.sum()}"
        )

        # Identify cells with high mitochondrial percentage
        data.obs["mt_outlier"] = data.obs["pct_counts_mt"] > pc_mt
        print(
            f"Cells with high mitochondrial percentage (> {pc_mt}%): {data.obs.mt_outlier.sum()}"
        )

        # Identify cells with high hemoglobin percentage
        data.obs["hb_outlier"] = data.obs["pct_counts_hb"] > pc_hb
        print(
            f"Cells with high hemoglobin percentage (> {pc_hb}%): {data.obs.hb_outlier.sum()}"
        )

        adata.obs.loc[data.obs.index, "outlier"] = data.obs["outlier"]
        adata.obs.loc[data.obs.index, "mt_outlier"] = data.obs["mt_outlier"]
        adata.obs.loc[data.obs.index, "hb_outlier"] = data.obs["hb_outlier"]

    # Print the overall statistics for the entire adata object
    print("\nOverall statistics for the entire adata object:")
    total_cells = adata.n_obs
    print(f"Total number of cells: {total_cells}")

    outlier_cells = adata.obs.outlier.sum()
    print(
        f"Outlier cells (based on total counts, gene counts, and top gene percentage): {outlier_cells} ({outlier_cells / total_cells * 100:.2f}%)"
    )

    mt_outlier_cells = adata.obs.mt_outlier.sum()
    print(
        f"Cells with high mitochondrial percentage (> {pc_mt}%): {mt_outlier_cells} ({mt_outlier_cells / total_cells * 100:.2f}%)"
    )

    hb_outlier_cells = adata.obs.hb_outlier.sum()
    print(
        f"Cells with high hemoglobin percentage (> {pc_hb}%): {hb_outlier_cells} ({hb_outlier_cells / total_cells * 100:.2f}%)"
    )

    combined_outliers = adata.obs.filter(regex="outlier").sum(axis=1).value_counts()
    for n_outliers, count in combined_outliers.items():
        print(
            f"{n_outliers} types of outliers: {count} ({count / total_cells * 100:.2f}%)"
        )

    return adata


def is_doublet(
    adata: sc.AnnData,
    sample_key: str = "sampleID",
    rate: float = 0.1,
    n_pcs: int = 30,
    threshold: float = 0.2,
    over_genes: float = 0.99,
    plot_umap: bool = True,
):
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

    Returns:
        adata (AnnData): AnnData object with doublet scores and predictions added.
    """

    adata.obs["doublet_scores"] = 0.0
    adata.obs["predicted_doublets"] = False
    adata.obs["predicted_doublets_final"] = False
    adata.obs["overexpressed_doublets"] = False

    for sample in adata.obs[sample_key].unique():
        print(f"Begin of post doublets removal and QC plot for sample: {sample}")
        data = adata[adata.obs[sample_key] == sample, :]

        # 获取当前样本的细胞数和特征数
        n_cells, n_features = data.shape
        # 动态设置n_pcs参数
        n_pcs = min(n_pcs, n_cells, n_features)

        try:
            scrub = scr.Scrublet(data.X, expected_doublet_rate=rate)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(
                verbose=False, n_prin_comps=n_pcs
            )
            final_doublets = scrub.call_doublets(threshold=threshold)
        except Exception as e:
            print(f"Scrublet failed for sample {sample}: {e}")
            predicted_doublets = None
            final_doublets = None
            doublet_scores = np.full(data.shape[0], np.nan)
            
        if predicted_doublets is None:
            predicted_doublets = np.full(data.shape[0], True)
        if final_doublets is None:
            final_doublets = np.full(data.shape[0], True)
            
        if plot_umap:
            scrub.set_embedding(
                "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
            )
            scrub.plot_embedding("UMAP", order_points=True)

        adata.obs.loc[data.obs.index, "doublet_scores"] = doublet_scores
        adata.obs.loc[data.obs.index, "predicted_doublets"] = predicted_doublets.astype(bool)
        adata.obs.loc[data.obs.index, "predicted_doublets_final"] = final_doublets.astype(bool)

    # Identify cells with overexpressed genes as potential doublets
    top_genes = np.quantile(adata.obs.n_genes_by_counts, over_genes)
    adata.obs["overexpressed_doublets"] = adata.obs["n_genes_by_counts"] > top_genes

    # Print the overall statistics for the entire adata object
    total_cells = adata.n_obs
    print(f"\nTotal number of cells: {total_cells}")

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
        f"Potential doublet cells (based on detected genes > {top_genes}): {overexpressed} ({overexpressed / total_cells * 100:.2f}%)"
    )

    combined_doublets = adata.obs.filter(regex="doublets").sum(axis=1).value_counts()
    for n_doublets, count in combined_doublets.items():
        print(
            f"{n_doublets} types of doublets: {count} ({count / total_cells * 100:.2f}%)"
        )

    return adata