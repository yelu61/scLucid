"""
Cell-cell communication analysis using CellPhoneDB.

This module provides a wrapper for CellPhoneDB, a popular tool for inferring
ligand-receptor interactions between cell populations.
"""

import logging
import os
import subprocess
from typing import Dict, List, Optional

import anndata
import pandas as pd

log = logging.getLogger(__name__)


def run_cellphonedb(
    adata: anndata.AnnData,
    groupby: str,
    out_dir: str,
    iterations: int = 1000,
    threads: int = 4,
    copy: bool = False,
) -> anndata.AnnData:
    """
    Run CellPhoneDB analysis for cell-cell communication.
    Results are saved to adata.uns['sclucid']['cellphonedb']
    """
    if copy:
        adata = adata.copy()

    # --- 1. Input Validation and Setup ---
    log.info("Starting CellPhoneDB analysis...")
    if groupby not in adata.obs:
        raise ValueError(f"Groupby key '{groupby}' not found in adata.obs.")
    if adata.raw is None:
        log.warning(
            "adata.raw is not set. Using adata.X. It is highly recommended to use raw counts for CellPhoneDB."
        )

    os.makedirs(out_dir, exist_ok=True)

    # --- 2. Prepare Input Files ---
    meta_file = os.path.join(out_dir, "meta.txt")
    counts_file = os.path.join(out_dir, "counts.h5ad")

    # Prepare metadata file
    meta_df = adata.obs[[groupby]].copy()
    meta_df.reset_index(inplace=True)
    meta_df.columns = ["Cell", "cell_type"]
    meta_df.to_csv(meta_file, sep="\t", index=False)
    log.info(f"Metadata file written to: {meta_file}")

    # Prepare counts file (using adata.raw if available)
    if adata.raw is not None:
        counts_adata = anndata.AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs)
    else:
        counts_adata = adata
    counts_adata.write_h5ad(counts_file)
    log.info(f"Counts file written to: {counts_file}")

    # --- 3. Run CellPhoneDB Command-Line Tool ---
    log.info("Executing CellPhoneDB statistical analysis...")
    cmd = [
        "cellphonedb",
        "method",
        "statistical_analysis",
        meta_file,
        counts_file,
        "--output-path",
        out_dir,
        "--iterations",
        str(iterations),
        "--threads",
        str(threads),
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        log.info("CellPhoneDB analysis completed successfully.")
        log.debug(f"CellPhoneDB stdout:\n{result.stdout}")
    except FileNotFoundError:
        raise FileNotFoundError(
            "`cellphonedb` command not found. Please ensure CellPhoneDB is installed and in your PATH."
        )
    except subprocess.CalledProcessError as e:
        log.error("CellPhoneDB analysis failed.")
        log.error(f"Stderr:\n{e.stderr}")
        raise RuntimeError(f"CellPhoneDB failed with error: {e.stderr}")

    # --- 4. Load Results back into AnnData ---
    log.info("Loading results back into AnnData object.")
    results = {}
    result_files = ["pvalues.txt", "significant_means.txt", "means.txt", "deconvoluted.txt"]
    for f in result_files:
        path = os.path.join(out_dir, f)
        if os.path.exists(path):
            key_name = f.replace(".txt", "")
            results[key_name] = pd.read_csv(path, sep="\t")

    adata.uns.setdefault("sclucid", {})["cellphonedb"] = {
        "results": results,
        "params": {
            "groupby": groupby,
            "iterations": iterations,
            "threads": threads,
            "out_dir": out_dir,
        },
    }
    if "significant_means" in results:
        log.info(f"Found {len(results['significant_means'])} significant interactions.")
    return adata


def run_cellphonedb_batch(
    adatas: List[anndata.AnnData],
    groupby: str,
    out_dir: str,
    iterations: int = 1000,
    threads: int = 4,
    sample_ids: Optional[List[str]] = None,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Batch run CellPhoneDB for a list of AnnData objects.
    Returns dict[sample_id] = AnnData with CellPhoneDB results or None if failed.
    """
    results = {}
    if sample_ids is None:
        sample_ids = [f"sample{i+1}" for i in range(len(adatas))]
    for adata, sid in zip(adatas, sample_ids):
        sample_dir = os.path.join(out_dir, sid)
        try:
            results[sid] = run_cellphonedb(
                adata=adata,
                groupby=groupby,
                out_dir=sample_dir,
                iterations=iterations,
                threads=threads,
            )
            log.info(f"CellPhoneDB completed for {sid}")
        except Exception as e:
            log.error(f"CellPhoneDB failed for {sid}: {e}")
            results[sid] = None
    return results


def run_cellphonedb_by_group(
    adata: anndata.AnnData,
    groupby: str,
    split_by: str,
    out_dir: str,
    iterations: int = 1000,
    threads: int = 4,
    min_cells: int = 100,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Run CellPhoneDB for each group (e.g., sample/condition) in AnnData.
    Returns dict[group] = AnnData with CellPhoneDB results.
    """
    results = {}
    groups = adata.obs[split_by].unique()
    for group in groups:
        adata_sub = adata[adata.obs[split_by] == group].copy()
        if adata_sub.n_obs < min_cells:
            log.warning(f"Group {group} skipped (too few cells: {adata_sub.n_obs})")
            continue
        group_dir = os.path.join(out_dir, f"{group}")
        try:
            results[group] = run_cellphonedb(
                adata=adata_sub,
                groupby=groupby,
                out_dir=group_dir,
                iterations=iterations,
                threads=threads,
            )
            log.info(f"CellPhoneDB completed for group {group}")
        except Exception as e:
            log.error(f"CellPhoneDB failed for group {group}: {e}")
            results[group] = None
    return results


def summarize_cellphonedb(
    adata: anndata.AnnData,
    save_dir: Optional[str] = None,
    top_n: int = 30,
):
    """
    Summarize and export main CellPhoneDB result tables and top interactions.
    """
    db = adata.uns.get("sclucid", {}).get("cellphonedb", {})
    results = db.get("results", {})
    params = db.get("params", {})

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    # Export main tables
    for k, df in results.items():
        if not isinstance(df, pd.DataFrame):
            continue
        out_path = os.path.join(save_dir, f"{k}.csv") if save_dir else None
        df.to_csv(out_path or f"{k}.csv", index=False)

    # Export top interactions (by mean or pvalue)
    if "significant_means" in results:
        top_df = results["significant_means"].sort_values("mean", ascending=False).head(top_n)
        top_path = os.path.join(save_dir, "top_interactions.csv") if save_dir else None
        top_df.to_csv(top_path or "top_interactions.csv", index=False)

    # Export parameter summary
    if params and save_dir:
        with open(os.path.join(save_dir, "cellphonedb_params.txt"), "w") as f:
            for k, v in params.items():
                f.write(f"{k}: {v}\n")
    log.info(f"CellPhoneDB summary exported to: {save_dir or os.getcwd()}")


# 可选：自动主图/markdown分析报告等可依赖上面输出扩展
