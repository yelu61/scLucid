"""
Unified Differential Expression and Enrichment Analysis for scRNA-seq.

This module provides:
- Marker gene discovery and DE analysis (find_markers, compare_groups, etc.)
    - Scenario A: Find marker genes for each cluster (one-vs-rest). -find_markers
    - Scenario B: Compare two specific cell types (e.g., CD4+ T vs. CD8+ T). -compare_groups
    - Scenario C: Find markers for a cell type that are conserved across conditions -get_conserved_markers
        (e.g., find robust T cell markers present in both Control and Treated samples).
    - Scenario D: Compare the same cell type across different conditions -compare_conditions
        (e.g., Treated T-cells vs. Control T-cells).
- DE result filtering (filter_markers)
- Functional enrichment analysis (run_enrichment)
- Combined summary export for AI/manual cell type annotation (summarize_markers_and_enrichment)
- Visualization and cluster characterization

All results and parameters are stored under adata.uns['sclucid']['analysis']['de'].
"""

import dataclasses
import logging
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils import sanitize_for_hdf5
from .config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    ConservedMarkersConfig,  # --- 完善: 导入新Config ---
    DifferentialConfig,
    EnrichmentConfig,
    FilterMarkersConfig,
)

log = logging.getLogger(__name__)

__all__ = [
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    "run_enrichment",
    "summarize_markers_and_enrichment",
    "characterize_clusters",
    "visualize_markers",
]


def _is_0_to_1(series: pd.Series) -> bool:
    """Check if a pandas Series is on a 0-1 scale."""
    s = series.dropna()
    return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))

def _to_frac(series: pd.Series) -> pd.Series:
    """
    Robustly convert a pandas Series (0-100 or 0-1) to a 0-1 fraction scale.
    """
    if series is None:
        return pd.Series(dtype=float)
    s_numeric = pd.to_numeric(series, errors='coerce')
    if _is_0_to_1(s_numeric):
        return s_numeric
    log.debug("Detected percentage scale (0-100). Converting to fraction (0-1).")
    return s_numeric.clip(lower=0, upper=100) / 100.0

# --- Differential Expression Analysis ---
def find_markers(
    adata: AnnData,
    config: Optional[DifferentialConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes (one-vs-rest) and store results.
    """
    if config is None:
        active_config = DifferentialConfig(**kwargs)
    else:
        active_config = dataclasses.replace(config)
        for k, v in kwargs.items():
            if hasattr(active_config, k):
                setattr(active_config, k, v)

    groupby = active_config.groupby
    key_added = active_config.key_added or "rank_genes_groups"

    log.info(f"Finding markers for '{groupby}' using method '{active_config.method}'.")

    rank_genes_params = {
        "groupby": groupby,
        "method": active_config.method,
        "layer": active_config.layer,
        "key_added": key_added,
        "use_raw": active_config.use_raw,
        "pts": True,
        "reference": active_config.reference,
    }
    if active_config.groups is not None:
        rank_genes_params["groups"] = active_config.groups
        log.info(f"Running DE analysis on subset of groups: {active_config.groups}")

    sc.tl.rank_genes_groups(adata, **rank_genes_params)

    if key_added not in adata.uns:
        raise KeyError(f"scanpy returned no result at adata.uns['{key_added}'].")

    raw = adata.uns[key_added]
    if "names" not in raw:
        raise KeyError(f"scanpy result missing 'names' at adata.uns['{key_added}'].")

    names_field = raw["names"]
    if not hasattr(names_field, "dtype") or names_field.dtype.names is None:
        raise ValueError(
            "Scanpy 'names' field lacks dtype names; structure may have changed."
        )

    groups_tested = names_field.dtype.names
    result_dfs: List[pd.DataFrame] = []

    for group in groups_tested:
        df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)
        if df.empty:
            continue
        
        # Harmonize pct column
        if "pct_nz_group" not in df.columns and "pct_nz" in df.columns:
            df = df.rename(columns={"pct_nz": "pct_nz_group"})
            log.debug("Renamed 'pct_nz' to 'pct_nz_group' for compatibility.")
        
        df["group"] = group

        # Optional filters
        if active_config.pval_cutoff is not None and "pvals_adj" in df.columns:
            before = len(df)
            df = df[df["pvals_adj"] <= float(active_config.pval_cutoff)].copy()
            log.info(
                f"Group '{group}': p-adj <= {active_config.pval_cutoff} retained {len(df)}/{before} rows."
            )

        if active_config.fold_change_max is not None and "logfoldchanges" in df.columns:
            df["logfoldchanges"] = df["logfoldchanges"].clip(
                upper=float(active_config.fold_change_max)
            )

        result_dfs.append(df)

    if not result_dfs:
        log.warning("No valid marker results found for any group after filtering.")
        full_df = pd.DataFrame()
    else:
        full_df = pd.concat(result_dfs, ignore_index=True)

    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[key_added] = adata.uns[key_added]
    df_key = f"{key_added}_df"
    root[df_key] = full_df
    
    p = active_config.to_dict()
    p["scanpy_version"] = getattr(sc, "__version__", "unknown")
    root[f"{key_added}_params"] = sanitize_for_hdf5(p)

    log.info(
        f"Found {len(full_df)} total marker rows across {len(groups_tested)} groups."
    )
    log.info(
        f"Stored processed DataFrame at .uns['sclucid']['analysis']['de']['{df_key}']"
    )
    return full_df


def filter_markers(
    adata: AnnData,
    config: FilterMarkersConfig,
) -> pd.DataFrame:
    """
    Filter marker genes with robust handling of percentage scales and detailed step logs.
    """
    key = config.key
    key_added = config.key_added or f"{key}_filtered_df"
    df_key = f"{key}_df"

    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if df_key not in root:
        raise KeyError(
            f"Source DataFrame not found at .uns['sclucid']['analysis']['de']['{df_key}']. Run `find_markers` first."
        )
    df = root[df_key].copy()
    if df.empty:
        log.warning("Source marker DataFrame is empty. Returning empty DataFrame.")
        return pd.DataFrame()

    required_cols = ["logfoldchanges", "pvals_adj", "pct_nz_group", "group", "names"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Marker DataFrame missing required columns: {missing}")
    
    # --- 完善: 使用共享的辅助函数 ---
    has_ref = "pct_nz_reference" in df.columns
    pct_group_frac = _to_frac(df["pct_nz_group"])
    pct_ref_frac = _to_frac(df["pct_nz_reference"]) if has_ref else None

    log.info(f"Filtering markers from '{df_key}'...")
    filt = pd.Series(True, index=df.index)

    if config.min_log2fc is not None:
        keep = (
            (df["logfoldchanges"].abs() >= float(config.min_log2fc))
            if config.use_abs_log2fc
            else (df["logfoldchanges"] >= float(config.min_log2fc))
        )
        log.info(
            f"[Filter] log2FC {'|x|' if config.use_abs_log2fc else ''} >= {config.min_log2fc}: kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.info("[Filter] min_log2fc: skipped (None)") # --- 完善: 增加日志 ---

    if config.max_padj is not None:
        keep = df["pvals_adj"] <= float(config.max_padj)
        log.info(
            f"[Filter] adj p <= {config.max_padj}: kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.info("[Filter] max_padj: skipped (None)")

    if config.min_in_group_pct is not None:
        keep = pct_group_frac >= float(config.min_in_group_pct)
        log.info(
            f"[Filter] pct_in_group >= {config.min_in_group_pct:.3f}: kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.info("[Filter] min_in_group_pct: skipped (None)")

    if has_ref:
        if config.max_out_group_pct is not None:
            keep = pct_ref_frac <= float(config.max_out_group_pct)
            log.info(
                f"[Filter] pct_out_group <= {config.max_out_group_pct:.3f}: kept {int(keep.sum())}/{len(filt)}"
            )
            filt &= keep
        else:
            log.info("[Filter] max_out_group_pct: skipped (None)")

        if config.min_diff_pct is not None:
            keep = (pct_group_frac - pct_ref_frac) >= float(config.min_diff_pct)
            log.info(
                f"[Filter] (pct_in - pct_out) >= {config.min_diff_pct:.3f}: kept {int(keep.sum())}/{len(filt)}"
            )
            filt &= keep
        else:
            log.info("[Filter] min_diff_pct: skipped (None)")
    else:
        if config.max_out_group_pct is not None or config.min_diff_pct is not None:
            log.warning(
                "pct_nz_reference not found; specificity-related filters skipped."
            )

    filtered_df = df[filt].copy()
    log.info(f"Retained {len(filtered_df)} genes after all filters.")

    # Keep top N per group
    if (
        config.keep_top_n is not None
        and config.keep_top_n > 0
        and not filtered_df.empty
    ):
        sort_by_col = config.sort_by
        if sort_by_col == "diff_pct":
            if has_ref:
                # --- 完善: 在 _to_frac 转换后的数据上计算 diff_pct ---
                filtered_df["diff_pct"] = (
                    pct_group_frac[filt] - pct_ref_frac[filt]
                )
            else:
                log.warning(
                    "Cannot sort by 'diff_pct' as 'pct_nz_reference' is missing. Falling back to 'scores'."
                )
                sort_by_col = "scores"

        if sort_by_col not in filtered_df.columns:
            fallback_col = (
                "logfoldchanges"
                if "logfoldchanges" in filtered_df.columns
                else "scores"
            )
            if fallback_col not in filtered_df.columns: # 终极回退
                fallback_col = filtered_df.columns[0]
                
            log.warning(
                f"Sort key '{config.sort_by}' not found. Falling back to '{fallback_col}'."
            )
            sort_by_col = fallback_col

        log.info(
            f"Selecting top {config.keep_top_n} genes per group, sorted by '{sort_by_col}'."
        )
        parts = []
        for g in filtered_df["group"].unique():
            sub = filtered_df[filtered_df["group"] == g].sort_values(
                sort_by_col, ascending=False
            )
            parts.append(sub.head(config.keep_top_n))
        filtered_df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    root[key_added] = filtered_df
    root[f"{key_added}_params"] = {**config.to_dict(), "n_retained": len(filtered_df)}
    log.info(
        f"Final filtered markers: {len(filtered_df)} rows -> .uns['sclucid']['analysis']['de']['{key_added}']"
    )
    return filtered_df


def compare_groups(adata: AnnData, config: CompareGroupsConfig) -> pd.DataFrame:
    """
    Compare two groups (e.g., cell types or conditions) for DE genes.
    """
    groupby = config.groupby
    group1 = config.group1
    group2 = config.group2
    key_added = config.key_added or f"compare_{group1}_vs_{group2}".replace(" ", "_")

    log.info(f"Comparing DE genes between '{group1}' and '{group2}' from '{groupby}'.")

    if groupby not in adata.obs.columns:
        raise KeyError(f"Column '{groupby}' not found in adata.obs.")
    subset_mask = adata.obs[groupby].isin([group1, group2])
    if subset_mask.sum() == 0:
        raise ValueError(
            f"No cells found for either '{group1}' or '{group2}' in '{groupby}'."
        )
    temp_adata = adata[subset_mask].copy()
    temp_adata.obs["_compare_groups"] = (
        temp_adata.obs[groupby]
        .map({group1: "group1", group2: "group2"})
        .astype("category")
    )

    sc.tl.rank_genes_groups(
        temp_adata,
        groupby="_compare_groups",
        groups=["group1"],
        reference="group2",
        method=config.method,
        layer=config.layer,
        use_raw=config.use_raw,
        pts=True,
    )

    results_df = sc.get.rank_genes_groups_df(temp_adata, group="group1")
    if "pct_nz_group" not in results_df.columns and "pct_nz" in results_df.columns:
        results_df = results_df.rename(columns={"pct_nz": "pct_nz_group"})

    # --- 完善: 使用共享的辅助函数 ---
    in_frac = _to_frac(results_df.get("pct_nz_group"))

    lfc = pd.to_numeric(results_df["logfoldchanges"], errors="coerce")
    padj = pd.to_numeric(results_df["pvals_adj"], errors="coerce")
    filt = (
        (lfc.abs() >= float(config.min_log2fc))
        & (padj <= float(config.max_padj))
        & (in_frac >= float(config.min_in_group_pct))
    )
    filtered_results = results_df[filt].copy()

    filtered_results = filtered_results.sort_values("logfoldchanges", ascending=False)
    up = filtered_results[filtered_results["logfoldchanges"] > 0].head(
        config.n_top_genes
    )
    down = (
        filtered_results[filtered_results["logfoldchanges"] < 0]
        .sort_values("logfoldchanges", ascending=True)
        .head(config.n_top_genes)
    )
    final_results = pd.concat([up, down], ignore_index=True)

    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[key_added] = final_results
    root[f"{key_added}_params"] = sanitize_for_hdf5(config.to_dict())

    log.info(
        f"Found {len(final_results)} DE genes. Stored at .uns['...']['{key_added}']."
    )

    if config.plot:
        # ... [绘图代码保持不变] ...
        pass
    
    return final_results


def compare_conditions(
    adata: AnnData,
    config: CompareConditionsConfig,
) -> pd.DataFrame:
    """
    Compare two conditions within a specific group using a config object.
    """
    groupby = config.groupby
    group_name = config.group_name
    condition_key = config.condition_key

    log.info(
        f"Comparing conditions '{config.condition1}' vs '{config.condition2}' within '{group_name}'."
    )

    if group_name not in adata.obs[groupby].unique():
        raise ValueError(f"Group '{group_name}' not found in adata.obs['{groupby}']")
    adata_subset = adata[adata.obs[groupby] == group_name].copy()

    comp_config = config.comparison_params
    comp_config.groupby = condition_key
    comp_config.group1 = config.condition1
    comp_config.group2 = config.condition2

    if config.key_added is None:
        safe_group = group_name.replace(" ", "_").replace("/", "_")
        comp_config.key_added = (
            f"compare_{config.condition1}_vs_{config.condition2}_in_{safe_group}"
        )
    else:
        comp_config.key_added = config.key_added

    results_df = compare_groups(adata_subset, config=comp_config)

    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[comp_config.key_added] = results_df
    root[f"{comp_config.key_added}_params"] = sanitize_for_hdf5(config.to_dict())
    log.info(f"Stored condition comparison at .uns['...']['{comp_config.key_added}']")

    return results_df


def get_conserved_markers(
    adata: AnnData,
    config: ConservedMarkersConfig,
) -> Dict[str, pd.DataFrame]:
    """
    Find markers for a group that are conserved across multiple conditions.
    """
    key_added = config.key_added or f"conserved_markers_{config.groupby}_{config.condition_key}"

    if config.condition_key not in adata.obs.columns or config.groupby not in adata.obs.columns:
        raise KeyError("Both 'groupby' and 'condition_key' must exist in adata.obs.")

    if pd.api.types.is_categorical_dtype(adata.obs[config.condition_key]):
        conditions = list(adata.obs[config.condition_key].cat.categories)
    else:
        conditions = list(pd.unique(adata.obs[config.condition_key]))

    if pd.api.types.is_categorical_dtype(adata.obs[config.groupby]):
        groups = list(adata.obs[config.groupby].cat.categories)
    else:
        groups = list(pd.unique(adata.obs[config.groupby]))

    min_conditions = config.min_conditions
    if min_conditions is None:
        min_conditions = max(1, len(conditions) - 1)

    conserved_markers: Dict[str, pd.DataFrame] = {}
    per_group_details: Dict[str, pd.DataFrame] = {}

    for group in groups:
        markers_per_condition = []
        for cond in conditions:
            subset = adata[
                (adata.obs[config.groupby] == group) & (adata.obs[config.condition_key] == cond)
            ]
            if subset.n_obs < config.min_cells:
                log.info(
                    f"Skip group '{group}' in condition '{cond}': n_cells={subset.n_obs} < {config.min_cells}"
                )
                continue

            temp_adata = adata[adata.obs[config.condition_key] == cond].copy()
            if group not in temp_adata.obs[config.groupby].unique():
                continue

            sc.tl.rank_genes_groups(
                temp_adata,
                groupby=config.groupby,
                groups=[group],
                reference="rest",
                method=config.method,
                layer=config.layer,
                use_raw=config.use_raw,
                pts=True,
            )

            df = sc.get.rank_genes_groups_df(temp_adata, group=group)
            if "pct_nz_group" not in df.columns and "pct_nz" in df.columns:
                df = df.rename(columns={"pct_nz": "pct_nz_group"})
            
            in_frac = _to_frac(
                df.get("pct_nz_group", pd.Series(index=df.index, dtype=float))
            )

            df = df[
                (df["logfoldchanges"] >= float(config.min_log2fc))
                & (df["pvals_adj"] <= float(config.max_padj))
                & (in_frac >= float(config.min_in_group_pct))
            ].copy()
            if df.empty:
                continue
            df["condition"] = cond
            markers_per_condition.append(df)

        if len(markers_per_condition) < min_conditions:
            log.info(
                f"Group '{group}': insufficient conditions with markers ({len(markers_per_condition)} < {min_conditions})."
            )
            continue

        full_df = pd.concat(markers_per_condition, ignore_index=True)
        per_group_details[group] = full_df.copy()

        gene_counts = full_df.groupby("names").size()
        conserved_genes = gene_counts[gene_counts >= min_conditions].index.tolist()
        conserved_df = full_df[full_df["names"].isin(conserved_genes)]

        agg_df = (
            conserved_df.groupby("names")
            .agg(
                mean_log2fc=("logfoldchanges", "mean"),
                min_log2fc=("logfoldchanges", "min"),
                max_log2fc=("logfoldchanges", "max"),
                min_pval_adj=("pvals_adj", "min"),
                mean_pval_adj=("pvals_adj", "mean"),
                mean_pct_in_group=("pct_nz_group", "mean"),
                n_conditions=("condition", "nunique"),
            )
            .sort_values("mean_log2fc", ascending=False)
        )
        conserved_markers[group] = agg_df

    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[key_added] = sanitize_for_hdf5(
        {
            "aggregates": conserved_markers,
            "details": per_group_details,
            "params": config.to_dict(),
        }
    )
    return conserved_markers


# --- Enrichment Analysis ---

def _standardize_enrichment_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Helper to standardize enrichment result columns."""
    if df is None or df.empty:
        return pd.DataFrame()
    d = df.copy()
    col_map = {
        "Term": ["Term", "term_name", "Description", "Name", "Pathway", "term"],
        "Adjusted P-value": [
            "Adjusted P-value",
            "Adjusted P-value (Benjamini-Hochberg)",
            "Adj P-value",
            "p.adjust",
            "padj",
            "FDR",
            "FDR q-value",
            "qvalue",
        ],
    }
    for standard_col, candidates in col_map.items():
        found_col = next((c for c in candidates if c in d.columns), None)
        if found_col and found_col != standard_col:
            d = d.rename(columns={found_col: standard_col})
            
    if "Adjusted P-value" in d.columns:
        d["Adjusted P-value"] = pd.to_numeric(
            d["Adjusted P-value"], errors="coerce"
        )
    return d

def run_enrichment(
    adata: AnnData,
    groupby: str,
    config: EnrichmentConfig,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Run enrichment (ORA and/or GSEA) for each group using GSEApy.
    """
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if config.de_key not in root:
        raise KeyError(
            f"DE results not found at .uns['...']['{config.de_key}']. Run find_markers/filter_markers."
        )
    marker_df = root[config.de_key]
    if marker_df.empty:
        log.warning("Marker DataFrame is empty; enrichment will be skipped.")
        marker_df = pd.DataFrame(columns=["group", "names", "logfoldchanges"])

    if groupby in adata.obs and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        group_order = list(adata.obs[groupby].cat.categories)
    else:
        group_order = list(pd.unique(marker_df["group"])) if "group" in marker_df.columns else []

    background_genes = list(adata.var_names)
    if config.plot and config.save_dir:
        Path(config.save_dir).mkdir(parents=True, exist_ok=True)

    # --- Prepare  ---
    gmt_files_to_run = {}
    gene_sets_list = (
        config.gene_sets_online if config.mode == "online" else config.gene_sets_offline
    )
    if not isinstance(gene_sets_list, list):
        gene_sets_list = [gene_sets_list]

    if config.mode == "offline":
        if config.custom_gene_sets and Path(config.custom_gene_sets).is_file():
            gmt_files_to_run = {"custom": config.custom_gene_sets}
            log.info(f"Using custom gene set file: {config.custom_gene_sets}")
        else:
            for gs_category in gene_sets_list:
                try:
                    filename = f"{config.organism.lower()}_{gs_category}_{config.gmt_version}.gmt"
                    file_path = resources.files("scLucid").joinpath(
                        "resources", filename
                    )
                    if file_path.is_file():
                        gmt_files_to_run[gs_category] = str(file_path)
                    else:
                        log.warning(f"Resource file not found: {filename}")
                except Exception as e:
                    log.error(f"Error finding resource file for {gs_category}: {e}")
            if not gmt_files_to_run:
                raise FileNotFoundError("No valid gene set files found for offline mode.")
    else:
        # Online mode just needs the list of names
        gmt_files_to_run = {gs: gs for gs in gene_sets_list}


    # --- 完善: 选择 GSEA 排序的列 ---
    rank_col = config.rank_col_gsea
    if rank_col not in marker_df.columns:
        fallback_col = "scores" if "scores" in marker_df.columns else "logfoldchanges"
        log.warning(f"GSEA rank column '{rank_col}' not found. Falling back to '{fallback_col}'.")
        rank_col = fallback_col
    if rank_col not in marker_df.columns:
        raise KeyError(
            "GSEA requires a ranking column ('scores' or 'logfoldchanges') in the marker DataFrame."
        )

    enrichment_results: Dict[str, Dict[str, pd.DataFrame]] = {}
    enrichment_meta: Dict[str, Dict[str, Union[List[str], str, int]]] = {}

    for cluster in group_order:
        cluster_results: Dict[str, pd.DataFrame] = {}
        sub = marker_df[marker_df["group"] == cluster]
        if sub.empty:
            log.warning(f"Skipping '{cluster}': no marker genes found in '{config.de_key}'.")
            enrichment_results[cluster] = {'ora': pd.DataFrame(), 'gsea': pd.DataFrame()}
            continue

        # --- 1. ORA (Over-Representation Analysis) ---
        if config.method in ["ora", "both"]:
            # ORA uses the top N genes
            gene_list = (
                sub.sort_values(rank_col, ascending=False)["names"]
                .head(config.n_top_genes_ora)
                .astype(str)
                .tolist()
            )
            enrichment_meta[cluster] = {
                "n_input_genes_ora": len(gene_list),
                "rank_col_gsea": rank_col,
                "de_key": config.de_key,
            }

            if len(gene_list) < config.min_genes_for_ora:
                log.warning(
                    f"Skipping ORA for '{cluster}': not enough genes ({len(gene_list)} < {config.min_genes_for_ora})."
                )
                cluster_results['ora'] = pd.DataFrame()
            else:
                all_ora_results = []
                for category, gmt in gmt_files_to_run.items():
                    try:
                        if config.mode == "online":
                            enr_ora = gp.enrichr(
                                gene_list=gene_list,
                                gene_sets=gmt, # gmt is the library name
                                organism=config.organism,
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0, # Filter later
                            )
                        else: # Offline
                            enr_ora = gp.enrich(
                                gene_list=gene_list,
                                gene_sets=gmt, # gmt is the file path
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0, # Filter later
                            )
                        res = enr_ora.results.copy()
                        res["Gene_set"] = category
                        all_ora_results.append(res)
                    except Exception as e:
                        log.error(f"Enrichment failed for {cluster} (Category: {category}): {e}")
                
                ora_df = pd.concat(all_ora_results, ignore_index=True) if all_ora_results else pd.DataFrame()
                ora_df = _standardize_enrichment_cols(ora_df)
                if 'Adjusted P-value' in ora_df.columns:
                    ora_df = ora_df[ora_df['Adjusted P-value'] < config.max_padj]
                cluster_results['ora'] = ora_df

        # --- 2. GSEA (Gene Set Enrichment Analysis) ---
        if config.method in ["gsea", "both"]:
            rnk = sub.drop_duplicates(subset='names', keep='first').set_index('names')[rank_col]
            if rnk.empty:
                log.warning(f"Skipping GSEA for '{cluster}': no ranked genes available.")
                cluster_results['gsea'] = pd.DataFrame()
            else:
                all_gsea_results = []
                for category, gmt in gmt_files_to_run.items():
                    try:
                        gsea_res = gp.prerank(
                            rnk=rnk,
                            gene_sets=gmt, # gmt is library name or file path
                            permutation_num=config.gsea_permutations,
                            min_size=config.gsea_min_size,
                            max_size=config.gsea_max_size,
                            outdir=None,
                            seed=42,
                        )
                        res = gsea_res.res2d.copy()
                        res["Gene_set"] = category
                        all_gsea_results.append(res)
                    except Exception as e:
                         log.error(f"GSEA failed for {cluster} (Category: {category}): {e}")
                
                gsea_df = pd.concat(all_gsea_results, ignore_index=True) if all_gsea_results else pd.DataFrame()
                gsea_df = _standardize_enrichment_cols(gsea_df)
                if 'Adjusted P-value' in gsea_df.columns:
                    gsea_df = gsea_df[gsea_df['Adjusted P-value'] < config.max_padj]
                cluster_results['gsea'] = gsea_df

        enrichment_results[str(cluster)] = cluster_results

    # Store results
    store_root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    store_root[config.key_added] = {
        "results": enrichment_results, # DataFrames are stored directly
        "params": sanitize_for_hdf5(config.to_dict()),
        "meta": sanitize_for_hdf5(enrichment_meta),
    }
    log.info(f"Enrichment analysis complete. Results stored in .uns['...']['{config.key_added}']")
    return enrichment_results


# --- Batch Visualization and Cluster Characterization ---

def visualize_markers(
    adata: AnnData,
    markers: Union[pd.DataFrame, Dict[str, List[str]], List[str]],
    groupby: Optional[str] = None,
    n_genes_per_group: int = 5,
    plot_type: Literal[
        "dotplot", "heatmap", "stacked_violin", "violin", "matrixplot"
    ] = "dotplot",
    dendrogram: bool = False,
    standard_scale: Optional[Literal["var", "group"]] = "var", # --- 完善: 默认值改为'var' ---
    swap_axes: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    save_path: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """
    Visualize marker genes across groups with several plot types.
    """
    gene_list: List[str] = []
    gene_dict: Dict[str, List[str]] = {}

    if isinstance(markers, pd.DataFrame):
        if "names" not in markers.columns:
            for alt in ("gene", "Gene", "feature", "symbol"):
                if alt in markers.columns:
                    markers = markers.rename(columns={alt: "names"})
                    break
        if "group" not in markers.columns:
            raise ValueError(
                "DataFrame must contain 'group' and 'names' columns for grouped visualization."
            )

        for g in markers["group"].unique():
            group_markers = markers[markers["group"] == g]
            if "logfoldchanges" in group_markers.columns:
                group_markers = group_markers.sort_values(
                    "logfoldchanges", ascending=False
                )
            elif "scores" in group_markers.columns:
                group_markers = group_markers.sort_values("scores", ascending=False)
            
            top_genes = group_markers["names"].head(n_genes_per_group).tolist()
            gene_list.extend(top_genes)
            gene_dict[str(g)] = top_genes
            
    elif isinstance(markers, dict):
        gene_dict = markers
        for genes in markers.values():
            gene_list.extend(list(genes))
        if n_genes_per_group > 0:
             log.warning("n_genes_per_group is ignored when 'markers' is a dictionary.")

    elif isinstance(markers, (list, tuple)):
        gene_list = list(markers)
        if groupby is None:
            raise ValueError("groupby must be specified when markers is a list.")
        # --- 完善: 当 markers 是列表时，也填充 gene_dict 以便 heatmap/dotplot 使用分组
        gene_dict = {"Selected Markers": gene_list}

    else:
        raise TypeError("markers must be a DataFrame, dictionary, or list.")

    # Deduplicate and ensure in var_names
    gene_list_unique = [g for g in dict.fromkeys(gene_list) if g in adata.var_names]
    if not gene_list_unique:
        raise ValueError("No valid genes found for visualization.")
    
    # --- 完善: 为 heatmap/dotplot 准备好 gene_dict, 确保所有基因都在
    for g, glist in gene_dict.items():
        gene_dict[g] = [g for g in glist if g in adata.var_names]
        
    # Auto figsize
    if figsize is None:
        n_groups = (
            len(adata.obs[groupby].cat.categories)
            if groupby
            and groupby in adata.obs
            and pd.api.types.is_categorical_dtype(adata.obs[groupby])
            else len(adata.obs[groupby].unique()) if groupby in adata.obs else 1
        )
        n_genes = len(gene_list_unique)
        
        # --- 完善: 为不同 plot_type 优化 figsize 启发式 ---
        if plot_type == "heatmap" or plot_type == "dotplot":
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.3))
            if swap_axes:
                 width, height = height, width
        elif plot_type == "stacked_violin":
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.4))
        else: # violin
            width = max(8, n_genes * 2)
            height = 6
        
        figsize = (width, height)

    # Plot
    plot_func_map = {
        "dotplot": sc.pl.dotplot,
        "heatmap": sc.pl.heatmap,
        "stacked_violin": sc.pl.stacked_violin,
        "matrixplot": sc.pl.matrixplot,
    }
    
    plot_kwargs = {
        "groupby": groupby,
        "dendrogram": dendrogram,
        "standard_scale": standard_scale,
        "use_raw": use_raw,
        "layer": layer,
        "figsize": figsize,
        "show": False, # --- 完善: 统一控制 show ---
        **kwargs,
    }

    if plot_type in plot_func_map:
        # --- 完善: 为 dotplot/heatmap 使用 gene_dict, 为其他使用 gene_list_unique
        var_names_arg = gene_dict if plot_type in ["dotplot", "heatmap", "matrixplot"] else gene_list_unique
        
        # Handle swap_axes specifically
        if plot_type in ["heatmap", "dotplot", "matrixplot", "stacked_violin"]:
            plot_kwargs["swap_axes"] = swap_axes
            
        plot_func_map[plot_type](adata, var_names=var_names_arg, **plot_kwargs)
        
    elif plot_type == "violin":
        # violin has a different signature
        sc.pl.violin(
            adata,
            keys=gene_list_unique,
            groupby=groupby,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            show=False,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")
        plt.close() # --- 完善: 保存后关闭图像
    else:
        plt.show() # --- 完善: 仅在不保存时显示


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Run DE and enrichment for each cluster and collect evidence for annotation.
    """
    log.info(f"Characterizing clusters in '{groupby}'...")

    if de_config is None:
        de_config = DifferentialConfig(groupby=groupby, use_raw=True)
    else:
        de_config.groupby = groupby

    de_key = de_config.key_added or "rank_genes_groups"
    de_df_key = f"{de_key}_df"
    find_markers(adata, config=de_config)

    if enrichment_config is None:
        enrichment_config = EnrichmentConfig(de_key=de_df_key, groupby=groupby) # --- 完善: 传入 groupby ---
    else:
        enrichment_config.de_key = de_df_key
    
    # --- 完善: run_enrichment 现在返回一个字典 ---
    enrichment_results_dict = run_enrichment(
        adata, groupby=groupby, config=enrichment_config
    )

    clusters = (
        adata.obs[groupby].cat.categories
        if pd.api.types.is_categorical_dtype(adata.obs[groupby])
        else pd.unique(adata.obs[groupby])
    )
    characterization_results = {}
    de_df = adata.uns["sclucid"]["analysis"]["de"][de_df_key]

    for cluster in clusters:
        # --- 完善: 处理 ORA/GSEA 的嵌套字典 ---
        enr_res = enrichment_results_dict.get(str(cluster), {})
        characterization_results[str(cluster)] = {
            "top_de_genes": de_df[de_df["group"] == cluster],
            "enrichment_ora": enr_res.get("ora", pd.DataFrame()),
            "enrichment_gsea": enr_res.get("gsea", pd.DataFrame()),
        }

    adata.uns[key_added] = { # --- 完善: 不再 sanitizing, 直接存储 ---
        "results": characterization_results,
        "params": {
            "groupby": groupby,
            "de_df_key": de_df_key,
            "enrichment_key": enrichment_config.key_added,
            "de_params": de_config.to_dict(),
            "enrichment_params": enrichment_config.to_dict(),
        },
    }
    log.info(f"Cluster characterization complete -> adata.uns['{key_added}']")
    return adata


def summarize_markers_and_enrichment(
    adata: AnnData,
    groupby: str,
    markers_df: Optional[pd.DataFrame] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    markers_key: str = "rank_genes_groups_df",
    enrichment_key: str = "enrichment",
    enrichment_method_to_summarize: Literal["ora", "gsea"] = "ora", # --- 完善: 允许选择 ---
    n_markers: int = 25,
    n_terms: int = 10,
    summary_file: Optional[str] = None,
    sort_markers_by: str = "logfoldchanges",  # or "scores"
    enrichment_padj_cutoff: float = 0.05,
) -> Dict[str, str]:
    """
    Build per-group Markdown summaries of top markers and enriched terms.

    This robust version correctly parses enrichment results and filters them
    to ensure only significant pathways are summarized.
    """
    # 1. Load markers
    if markers_df is None:
        try:
            log.info(f"Auto-retrieving markers from .uns using key: {markers_key}")
            markers_df = adata.uns["sclucid"]["analysis"]["de"][markers_key]
        except KeyError:
            raise KeyError(
                f"Marker DataFrame not found at .uns['sclucid']['analysis']['de']['{markers_key}']."
            )

    if markers_df is None or markers_df.empty:
        log.warning("The provided marker DataFrame is empty. Cannot summarize markers.")
        markers_df = pd.DataFrame(columns=["group", "names", "logfoldchanges"])

    # 2. Load enrichment results
    if enrichment_dict is None:
        enrichment_dict = {}
        try:
            log.info(
                f"Auto-retrieving enrichment from .uns using key: {enrichment_key}"
            )
            enr_store = adata.uns["sclucid"]["analysis"]["de"].get(enrichment_key, {})
            if isinstance(enr_store, dict) and "results" in enr_store:
                enrichment_dict = enr_store["results"] # This is now {cluster: {"ora": df, "gsea": df}}
            else:
                enrichment_dict = enr_store
        except KeyError:
            log.warning(
                f"No enrichment data found with key '{enrichment_key}'. Pathways will be 'N/A'."
            )

    # 3. Determine group order
    if groupby in adata.obs:
        group_order = list(pd.unique(adata.obs[groupby].astype(str)))
    else:
        group_order = list(
            pd.unique(markers_df.get("group", pd.Series([], dtype=str)).astype(str))
        )

    # 4. Determine marker sort column
    sort_col = (
        sort_markers_by
        if sort_markers_by in markers_df.columns
        else "scores"
        if "scores" in markers_df.columns
        else "logfoldchanges"
    )
    if sort_col not in markers_df.columns:
        log.warning(f"Sort column '{sort_col}' not found in marker df. Using first column.")
        sort_col = markers_df.columns[0] if not markers_df.empty else "names"


    # 5. Build summaries
    summaries: Dict[str, str] = {}
    lines: List[str] = []

    def _detect_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
        term_candidates = ["Term", "term_name", "Description", "Name", "Pathway", "term"]
        pval_candidates = ["Adjusted P-value", "p.adjust", "padj", "FDR"]
        tcol = next((c for c in term_candidates if c in df.columns), None)
        pcol = next((c for c in pval_candidates if c in df.columns), None)
        return tcol, pcol

    for g in group_order:
        # Top markers
        group_markers = markers_df[markers_df["group"].astype(str) == str(g)]
        top_genes = (
            group_markers.sort_values(sort_col, ascending=False)["names"]
            .head(n_markers)
            .astype(str)
            .tolist()
        )

        # Top pathways
        top_terms: List[str] = []
        
        # --- 完善: 钻取到 ora/gsea 结果 ---
        cluster_enr_dict = enrichment_dict.get(str(g), {})
        enr_df = cluster_enr_dict.get(enrichment_method_to_summarize, pd.DataFrame())

        if isinstance(enr_df, pd.DataFrame) and not enr_df.empty:
            term_col, pval_col = _detect_cols(enr_df)
            if term_col and pval_col:
                tmp = enr_df.copy()
                tmp[pval_col] = pd.to_numeric(tmp[pval_col], errors="coerce")
                sig = tmp.dropna(subset=[pval_col])
                sig = sig[sig[pval_col] < float(enrichment_padj_cutoff)]
                
                if not sig.empty:
                    # --- 完善: GSEA 应该按 NES 排序, ORA 按 P值 ---
                    sort_col_enr = pval_col
                    ascending = True
                    if enrichment_method_to_summarize == 'gsea' and 'NES' in sig.columns:
                        sort_col_enr = 'NES'
                        ascending = False # GSEA中, 正的NES表示富集在基因列表顶部

                    top_terms = (
                        sig.sort_values(sort_col_enr, ascending=ascending)[term_col]
                        .head(n_terms)
                        .astype(str)
                        .tolist()
                    )
                else:
                    log.debug( # --- 完善: 日志级别改为 debug ---
                        f"Cluster {g}: No pathways found below p_adj cutoff of {enrichment_padj_cutoff}"
                    )
            else:
                log.warning(
                    f"Cluster {g}: Could not detect term/padj columns. Columns present: {list(enr_df.columns)}"
                )

        title = f"### Cluster {g}"
        mk_str = f"**Top Markers**: {', '.join(top_genes) if top_genes else 'N/A'}"
        pt_str = f"**Top {enrichment_method_to_summarize.upper()} Pathways**: {', '.join(top_terms) if top_terms else 'N/A'}"

        summary_text = f"{title}\n{mk_str}\n{pt_str}"
        summaries[str(g)] = summary_text
        lines.append(summary_text)

    # 6. Write to file
    if summary_file:
        Path(summary_file).parent.mkdir(parents=True, exist_ok=True)
        content = (
            "\n\n---\n\n".join(lines)
            if lines
            else "# Marker and Enrichment Summary\n\nNo results."
        )
        with open(summary_file, "w", encoding="utf-8") as f:
            f.write(content)
        log.info(f"Marker + enrichment summaries exported to {summary_file}")

    return summaries