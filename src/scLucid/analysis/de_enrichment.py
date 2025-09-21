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


# --- Differential Expression Analysis ---
def find_markers(
    adata: AnnData,
    config: Optional[DifferentialConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes (one-vs-rest) and store both raw scanpy output and a concatenated DataFrame.

    Enhancements:
    - Robustness to Scanpy structure changes.
    - Optional post-hoc p-value cutoff and FC clipping.
    - Trace saved with scanpy version and parameters.
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

    # Run scanpy rank genes
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
            log.warning("Renamed 'pct_nz' to 'pct_nz_group' for compatibility.")
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
    root[key_added] = adata.uns[key_added]  # raw scanpy result
    df_key = f"{key_added}_df"
    root[df_key] = full_df
    # Params trace
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

    Enhancements:
    - Detect pct scale (0–1 or 0–100) and convert to fraction.
    - More robust sorting fallback and top-N selection logging.
    - Trace saved under .uns with filter counts.
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

    def _is_0_to_1(series: pd.Series) -> bool:
        s = series.dropna()
        return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))

    def _to_frac(series: pd.Series) -> pd.Series:
        return series if _is_0_to_1(series) else series.clip(lower=0, upper=100) / 100.0

    has_ref = "pct_nz_reference" in df.columns
    pct_group_frac = _to_frac(df["pct_nz_group"])
    pct_ref_frac = _to_frac(df["pct_nz_reference"]) if has_ref else None

    log.info("Using fraction scale (0–1) for percentage-based filters.")

    filt = pd.Series(True, index=df.index)

    if config.min_log2fc is not None:
        keep = (
            (df["logfoldchanges"].abs() >= float(config.min_log2fc))
            if config.use_abs_log2fc
            else (df["logfoldchanges"] >= float(config.min_log2fc))
        )
        log.info(
            f"[Filter] log2FC {'|x|' if config.use_abs_log2fc else ''} >= {config.min_log2fc}: kept {int(keep.sum())}/{len(keep)}"
        )
        filt &= keep

    if config.max_padj is not None:
        keep = df["pvals_adj"] <= float(config.max_padj)
        log.info(
            f"[Filter] adj p <= {config.max_padj}: kept {int(keep.sum())}/{len(keep)}"
        )
        filt &= keep

    if config.min_in_group_pct is not None:
        keep = pct_group_frac >= float(config.min_in_group_pct)
        log.info(
            f"[Filter] pct_in_group >= {config.min_in_group_pct:.3f}: kept {int(keep.sum())}/{len(keep)}"
        )
        filt &= keep

    if has_ref:
        if config.max_out_group_pct is not None:
            keep = pct_ref_frac <= float(config.max_out_group_pct)
            log.info(
                f"[Filter] pct_out_group <= {config.max_out_group_pct:.3f}: kept {int(keep.sum())}/{len(keep)}"
            )
            filt &= keep
        if config.min_diff_pct is not None:
            keep = (pct_group_frac - pct_ref_frac) >= float(config.min_diff_pct)
            log.info(
                f"[Filter] (pct_in - pct_out) >= {config.min_diff_pct:.3f}: kept {int(keep.sum())}/{len(keep)}"
            )
            filt &= keep
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
            if "pct_nz_reference" in filtered_df.columns:
                filtered_df["diff_pct"] = (
                    filtered_df["pct_nz_group"] - filtered_df["pct_nz_reference"]
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
                else filtered_df.columns[0]
            )
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

    # Store
    root[key_added] = filtered_df
    # Save params and counts
    root[f"{key_added}_params"] = {**config.to_dict(), "n_retained": len(filtered_df)}
    log.info(
        f"Final filtered markers: {len(filtered_df)} rows -> .uns['sclucid']['analysis']['de']['{key_added}']"
    )
    return filtered_df


def compare_groups(adata: AnnData, config: CompareGroupsConfig) -> pd.DataFrame:
    """
    Compare two groups (e.g., cell types or conditions) for DE genes using a config object.
    Improvements:
    - Preserve raw/layers by operating on a subset copy (no new empty AnnData).
    - Unify pct scale to 0–1 before applying thresholds.
    - More explicit up/down selection avoiding tail() pitfalls.
    - Store params for reproducibility.
    """
    groupby = config.groupby
    group1 = config.group1
    group2 = config.group2
    key_added = config.key_added or f"compare_{group1}_vs_{group2}".replace(" ", "_")

    log.info(f"Comparing DE genes between '{group1}' and '{group2}' from '{groupby}'.")

    # --- Subset while preserving raw & layers ---
    if groupby not in adata.obs.columns:
        raise KeyError(f"Column '{groupby}' not found in adata.obs.")
    subset_mask = adata.obs[groupby].isin([group1, group2])
    if subset_mask.sum() == 0:
        raise ValueError(
            f"No cells found for either '{group1}' or '{group2}' in '{groupby}'."
        )
    temp_adata = adata[subset_mask].copy()  # keeps .raw and .layers
    temp_adata.obs["_compare_groups"] = (
        temp_adata.obs[groupby]
        .map({group1: "group1", group2: "group2"})
        .astype("category")
    )

    # --- Run DE on the subset copy ---
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
    # Harmonize pct column
    if "pct_nz_group" not in results_df.columns and "pct_nz" in results_df.columns:
        results_df = results_df.rename(columns={"pct_nz": "pct_nz_group"})
        log.warning("Renamed 'pct_nz' to 'pct_nz_group' for compatibility.")

    # --- Scale to fraction for filtering ---
    def _is_0_to_1(s: pd.Series) -> bool:
        s = s.dropna()
        return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))

    def _to_frac(s: pd.Series) -> pd.Series:
        if s is None:
            return pd.Series(index=results_df.index, dtype=float)
        s = pd.to_numeric(s, errors="coerce")
        if _is_0_to_1(s):
            return s
        return s.clip(lower=0, upper=100) / 100.0

    in_frac = _to_frac(results_df.get("pct_nz_group"))

    # --- Apply thresholds ---
    lfc = pd.to_numeric(results_df["logfoldchanges"], errors="coerce")
    padj = pd.to_numeric(results_df["pvals_adj"], errors="coerce")
    filt = (
        (lfc.abs() >= float(config.min_log2fc))
        & (padj <= float(config.max_padj))
        & (in_frac >= float(config.min_in_group_pct))
    )
    filtered_results = results_df[filt].copy()

    # --- Split upregulated/downregulated and take top N by magnitude ---
    filtered_results = filtered_results.sort_values("logfoldchanges", ascending=False)
    up = filtered_results[filtered_results["logfoldchanges"] > 0].head(
        config.n_top_genes
    )
    down = (
        filtered_results[filtered_results["logfoldchanges"] < 0]
        .sort_values("logfoldchanges", ascending=True)  # most negative first
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

    # --- Volcano plot (unchanged aesthetics) ---
    if config.plot:
        plt.figure(figsize=(10, 7))
        plt.scatter(
            results_df["logfoldchanges"],
            -np.log10(pd.to_numeric(results_df["pvals_adj"], errors="coerce").clip(1e-300)),
            alpha=0.3,
            s=10,
            color="grey",
            label="Not significant",
            rasterized=True,
        )
        plt.scatter(
            up["logfoldchanges"],
            -np.log10(pd.to_numeric(up["pvals_adj"], errors="coerce").clip(1e-300)),
            alpha=0.7,
            s=30,
            color="red",
            label=f"Higher in {group1}",
        )
        plt.scatter(
            down["logfoldchanges"],
            -np.log10(pd.to_numeric(down["pvals_adj"], errors="coerce").clip(1e-300)),
            alpha=0.7,
            s=30,
            color="blue",
            label=f"Higher in {group2}",
        )
        plt.axvline(x=config.min_log2fc, color="grey", linestyle="--", alpha=0.5)
        plt.axvline(x=-config.min_log2fc, color="grey", linestyle="--", alpha=0.5)
        plt.axhline(
            y=-np.log10(config.max_padj), color="grey", linestyle="--", alpha=0.5
        )
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-log10(Adjusted p-value)")
        plt.title(f"Differential Expression: {group1} vs {group2}")
        plt.legend()
        plt.tight_layout()
        if config.save_dir:
            Path(config.save_dir).mkdir(parents=True, exist_ok=True)
            safe_key = str(key_added).replace("/", "_").replace(" ", "_")
            plt.savefig(Path(config.save_dir) / f"{safe_key}_volcano.png", dpi=300)
        plt.show()

    return final_results


def compare_conditions(
    adata: AnnData,
    config: CompareConditionsConfig,
) -> pd.DataFrame:
    """
    Compare two conditions within a specific group using a config object.
    Improvement:
    - Store params for reproducibility.
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
    groupby: str,
    condition_key: str,
    method: str = "wilcoxon",
    min_cells: int = 10,
    min_conditions: Optional[int] = None,
    min_log2fc: float = 0.5,
    max_padj: float = 0.05,
    min_in_group_pct: float = 0.25,
    layer: Optional[str] = None,
    use_raw: bool = False,
    key_added: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Find markers for a group that are conserved across multiple conditions.
    Improvements:
    - Unify pct scale (0–1).
    - Store per-condition detail alongside aggregates.
    - More informative logging.
    """
    if key_added is None:
        key_added = f"conserved_markers_{groupby}_{condition_key}"

    if condition_key not in adata.obs.columns or groupby not in adata.obs.columns:
        raise KeyError("Both 'groupby' and 'condition_key' must exist in adata.obs.")

    if pd.api.types.is_categorical_dtype(adata.obs[condition_key]):
        conditions = list(adata.obs[condition_key].cat.categories)
    else:
        conditions = list(pd.unique(adata.obs[condition_key]))

    if pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        groups = list(adata.obs[groupby].cat.categories)
    else:
        groups = list(pd.unique(adata.obs[groupby]))

    if min_conditions is None:
        min_conditions = max(1, len(conditions) - 1)

    def _is_0_to_1(s: pd.Series) -> bool:
        s = s.dropna()
        return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))

    def _to_frac(s: pd.Series) -> pd.Series:
        if _is_0_to_1(s):
            return s
        return s.clip(lower=0, upper=100) / 100.0

    conserved_markers: Dict[str, pd.DataFrame] = {}
    per_group_details: Dict[str, pd.DataFrame] = {}

    for group in groups:
        markers_per_condition = []
        for cond in conditions:
            subset = adata[
                (adata.obs[groupby] == group) & (adata.obs[condition_key] == cond)
            ]
            if subset.n_obs < min_cells:
                log.info(
                    f"Skip group '{group}' in condition '{cond}': n_cells={subset.n_obs} < {min_cells}"
                )
                continue

            temp_adata = adata[adata.obs[condition_key] == cond].copy()
            # Ensure group exists within this condition subset
            if group not in temp_adata.obs[groupby].unique():
                continue

            sc.tl.rank_genes_groups(
                temp_adata,
                groupby=groupby,
                groups=[group],
                reference="rest",
                method=method,
                layer=layer,
                use_raw=use_raw,
                pts=True,
            )

            df = sc.get.rank_genes_groups_df(temp_adata, group=group)
            if "pct_nz_group" not in df.columns and "pct_nz" in df.columns:
                df = df.rename(columns={"pct_nz": "pct_nz_group"})
            # Unify to fraction
            in_frac = _to_frac(
                df.get("pct_nz_group", pd.Series(index=df.index, dtype=float))
            )

            df = df[
                (df["logfoldchanges"] >= float(min_log2fc))
                & (df["pvals_adj"] <= float(max_padj))
                & (in_frac >= float(min_in_group_pct))
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
    params_dict = sanitize_for_hdf5(
        {
            "groupby": groupby,
            "condition_key": condition_key,
            "method": method,
            "min_cells": min_cells,
            "min_conditions": min_conditions,
            "min_log2fc": min_log2fc,
            "max_padj": max_padj,
            "min_in_group_pct": min_in_group_pct,
            "layer": layer,
            "use_raw": use_raw,
        }
    )

    root[key_added] = sanitize_for_hdf5(
        {
            "aggregates": conserved_markers,
            "details": per_group_details,
            "params": params_dict,
        }
    )
    return conserved_markers


# --- Enrichment Analysis ---


def run_enrichment(
    adata: AnnData,
    groupby: str,
    config: EnrichmentConfig,
) -> Dict[str, pd.DataFrame]:
    """
    Run enrichment for each group using GSEApy.
    Improvements:
    - Flexible marker ranking (scores preferred if available).
    - Standardize result column names (Term, Adjusted P-value) to ease downstream usage.
    - Record the gene list used per group for reproducibility.
    - Preserve group order from adata.obs[groupby] if categorical.
    """
    de_key = config.de_key
    mode = config.mode
    organism = config.organism
    gene_sets_to_use = (
        config.gene_sets_online if mode == "online" else config.gene_sets_offline
    )
    gmt_version = config.gmt_version
    custom_gene_sets = config.custom_gene_sets
    n_top_genes = config.n_top_genes
    key_added = config.key_added
    min_genes_for_enrichment = config.min_genes_for_enrichment
    max_padj = config.max_padj
    plot = config.plot
    save_dir = config.save_dir
    n_plot_terms = config.n_plot_terms

    log.info(f"Running enrichment for '{groupby}' in '{mode}' mode.")

    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if de_key not in root:
        raise KeyError(
            f"DE results not found at .uns['sclucid']['analysis']['de']['{de_key}']. Run find_markers/filter_markers first."
        )
    marker_df = root[de_key]
    if marker_df.empty:
        log.warning(
            "Marker DataFrame is empty; enrichment will be skipped for all groups."
        )
        marker_df = pd.DataFrame(columns=["group", "names", "logfoldchanges"])

    # Determine group order
    if groupby in adata.obs and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        group_order = list(adata.obs[groupby].cat.categories)
    else:
        group_order = (
            list(pd.unique(marker_df["group"])) if "group" in marker_df.columns else []
        )

    background_genes = list(adata.var_names)
    if plot and save_dir:
        Path(save_dir).mkdir(parents=True, exist_ok=True)

    # Prepare gene sets for offline
    gmt_files_to_run = {}
    if mode == "offline":
        if custom_gene_sets and Path(custom_gene_sets).is_file():
            gmt_files_to_run = {"custom": custom_gene_sets}
            log.info(f"Using custom gene set file: {custom_gene_sets}")
        else:
            for gs_category in (
                gene_sets_to_use
                if isinstance(gene_sets_to_use, list)
                else [gene_sets_to_use]
            ):
                try:
                    filename = f"{organism.lower()}_{gs_category}_{gmt_version}.gmt"
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
                raise FileNotFoundError(
                    "No valid gene set files found for offline mode."
                )

    # Helper to standardize enrichment result columns
    def _standardize_enrichment_cols(df: pd.DataFrame) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame()
        d = df.copy()
        # Standardize term column
        for c in ["Term", "term_name", "Description", "Name", "Pathway", "term"]:
            if c in d.columns:
                if c != "Term":
                    d = d.rename(columns={c: "Term"})
                break
        # Standardize adjusted p column
        for c in [
            "Adjusted P-value",
            "Adjusted P-value (Benjamini-Hochberg)",
            "Adj P-value",
            "p.adjust",
            "padj",
            "FDR",
            "FDR q-value",
            "qvalue",
        ]:
            if c in d.columns:
                if c != "Adjusted P-value":
                    d = d.rename(columns={c: "Adjusted P-value"})
                break
        # Ensure numeric type for Adjusted P-value if present
        if "Adjusted P-value" in d.columns:
            d["Adjusted P-value"] = pd.to_numeric(
                d["Adjusted P-value"], errors="coerce"
            )

        return d

    # Choose ranking column for top gene selection
    rank_col = (
        "scores"
        if (config.prefer_score_for_enrichment and "scores" in marker_df.columns)
        else "logfoldchanges"
        if "logfoldchanges" in marker_df.columns
        else None
    )
    if rank_col is None:
        raise KeyError(
            "marker_df must contain either 'scores' or 'logfoldchanges' to rank genes for enrichment."
        )

    enrichment_results: Dict[str, pd.DataFrame] = {}
    enrichment_meta: Dict[str, Dict[str, Union[List[str], str, int]]] = {}

    for cluster in group_order:
        try:
            sub = marker_df[marker_df["group"] == cluster]
            if sub.empty:
                log.info(f"Skipping group '{cluster}': no markers.")
                enrichment_results[cluster] = pd.DataFrame()
                continue

            gene_list = (
                sub.sort_values(rank_col, ascending=False)["names"]
                .head(n_top_genes)
                .astype(str)
                .tolist()
            )
            enrichment_meta[cluster] = {
                "n_input_genes": len(gene_list),
                "rank_col": rank_col,
                "de_key": de_key,
            }

            if len(gene_list) < min_genes_for_enrichment:
                log.warning(
                    f"Skipping '{cluster}': not enough genes ({len(gene_list)} < {min_genes_for_enrichment})."
                )
                enrichment_results[cluster] = pd.DataFrame()
                continue

            if mode == "online":
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=gene_sets_to_use,
                    organism=organism,
                    background=len(background_genes),
                    outdir=None,
                    cutoff=max_padj,
                )
                results = enr.results
            else:
                all_offline_results = []
                for category, gmt_file in gmt_files_to_run.items():
                    enr = gp.enrich(
                        gene_list=gene_list,
                        gene_sets=gmt_file,
                        background=len(background_genes),
                        outdir=None,
                        cutoff=max_padj,
                    )
                    res = enr.results.copy()
                    res["Gene_set"] = category
                    all_offline_results.append(res)
                results = (
                    pd.concat(all_offline_results, ignore_index=True)
                    if all_offline_results
                    else pd.DataFrame()
                )

            results = _standardize_enrichment_cols(results)
            sorted_results = (
                results.sort_values("Adjusted P-value")
                if not results.empty
                else pd.DataFrame()
            )
            enrichment_results[cluster] = sorted_results

            if plot and not sorted_results.empty:
                top_pathways = sorted_results.head(n_plot_terms)
                plt.figure(figsize=(10, max(4, len(top_pathways) * 0.4)))
                y_labels = [
                    t[:70] + "..." if isinstance(t, str) and len(t) > 70 else t
                    for t in top_pathways["Term"]
                ]
                plt.barh(
                    y_labels,
                    -np.log10(top_pathways["Adjusted P-value"]),
                    color="skyblue",
                )
                plt.xlabel("-log10(Adjusted P-value)")
                plt.title(f"Top {len(top_pathways)} Enriched Pathways for {cluster}")
                plt.gca().invert_yaxis()
                plt.tight_layout()
                if save_dir:
                    safe_cluster = str(cluster).replace("/", "_").replace(" ", "_")
                    plt.savefig(
                        Path(save_dir) / f"{safe_cluster}_enrichment.png",
                        dpi=300,
                        bbox_inches="tight",
                    )
                plt.close()

        except Exception as e:
            log.warning(f"Enrichment failed for group '{cluster}': {str(e)}")
            enrichment_results[cluster] = pd.DataFrame()

    # Store results without sanitizing DataFrames; only params/meta sanitized
    store_root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    store_root[key_added] = {
        "results": enrichment_results,  # keep DataFrames as-is
        "params": sanitize_for_hdf5(config.to_dict()),
        "meta": sanitize_for_hdf5(enrichment_meta),
    }
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
    standard_scale: Optional[Literal["var", "group"]] = None,
    swap_axes: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    save_path: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """
    Visualize marker genes across groups with several plot types.
    Improvements:
    - Robust gene extraction from df/dict/list.
    - Safer default figsize with min/max bounds.
    """
    gene_list: List[str] = []

    if isinstance(markers, pd.DataFrame):
        if "names" not in markers.columns:
            # Try to remap common alternatives
            for alt in ("gene", "Gene", "feature", "symbol"):
                if alt in markers.columns:
                    markers = markers.rename(columns={alt: "names"})
                    log.info(f"Renamed '{alt}' to 'names' for plotting.")
                    break
        if "group" not in markers.columns:
            raise ValueError(
                "DataFrame must contain 'group' and 'names' columns for grouped visualization."
            )

        if n_genes_per_group > 0:
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
        else:
            gene_list = markers["names"].tolist()

    elif isinstance(markers, dict):
        if n_genes_per_group > 0:
            for genes in markers.values():
                gene_list.extend(list(genes)[:n_genes_per_group])
        else:
            for genes in markers.values():
                gene_list.extend(list(genes))

    elif isinstance(markers, (list, tuple)):
        gene_list = list(markers)
        if groupby is None:
            raise ValueError("groupby must be specified when markers is a list.")
    else:
        raise TypeError("markers must be a DataFrame, dictionary, or list.")

    # Deduplicate and ensure in var_names
    gene_list = [g for g in dict.fromkeys(gene_list) if g in adata.var_names]
    if not gene_list:
        raise ValueError("No valid genes found for visualization.")

    # Auto figsize
    if figsize is None:
        n_groups = (
            len(adata.obs[groupby].cat.categories)
            if groupby
            and groupby in adata.obs
            and pd.api.types.is_categorical_dtype(adata.obs[groupby])
            else 1
        )
        n_genes = len(gene_list)
        width = max(6, min(16, n_genes * 0.5))
        height = max(4, min(12, n_groups * 0.4))
        # If swapping axes, flip heuristic a bit
        if swap_axes:
            width, height = (
                max(6, min(16, n_groups * 0.5)),
                max(4, min(12, n_genes * 0.4)),
            )
        figsize = (width, height)

    # Plot
    if plot_type == "dotplot":
        sc.pl.dotplot(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "heatmap":
        sc.pl.heatmap(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "stacked_violin":
        sc.pl.stacked_violin(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "violin":
        sc.pl.violin(
            adata,
            keys=gene_list,
            groupby=groupby,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    elif plot_type == "matrixplot":
        sc.pl.matrixplot(
            adata,
            var_names=gene_list,
            groupby=groupby,
            dendrogram=dendrogram,
            standard_scale=standard_scale,
            swap_axes=swap_axes,
            use_raw=use_raw,
            layer=layer,
            figsize=figsize,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Run DE and enrichment for each cluster and collect evidence for annotation.
    Improvements:
    - Store used keys and params for reproducibility.
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
        enrichment_config = EnrichmentConfig(de_key=de_df_key)
    else:
        enrichment_config.de_key = de_df_key

    enrichment_results = run_enrichment(
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
        characterization_results[cluster] = {
            "top_de_genes": de_df[de_df["group"] == cluster],
            "enrichment": enrichment_results.get(cluster, pd.DataFrame()),
        }

    adata.uns[key_added] = sanitize_for_hdf5(
        {
            "results": characterization_results,
            "params": {
                "groupby": groupby,
                "de_df_key": de_df_key,
                "enrichment_key": enrichment_config.key_added,
                "de_params": de_config.to_dict(),
                "enrichment_params": enrichment_config.to_dict(),
            },
        }
    )
    log.info(f"Cluster characterization complete -> adata.uns['{key_added}']")
    return adata


# --- Marker + Enrichment Summary for AI/manual annotation ---


# In de_enrichment.py


def summarize_markers_and_enrichment(
    adata: AnnData,
    groupby: str,
    markers_df: Optional[pd.DataFrame] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    markers_key: str = "rank_genes_groups_df",
    enrichment_key: str = "enrichment",
    n_markers: int = 25,
    n_terms: int = 10,
    summary_file: Optional[str] = None,
    sort_markers_by: str = "logfoldchanges",  # or "scores"
    enrichment_padj_cutoff: float = 0.05,  #
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
                enrichment_dict = enr_store["results"]
            else:
                enrichment_dict = enr_store  # Fallback for older format
        except KeyError:
            log.warning(
                f"No enrichment data found with key '{enrichment_key}'. Pathways will be 'N/A'."
            )

    # 2.1 Attempt to auto-deserialize any non-DataFrame entries (compat with older sanitized storage)
    def _to_df(obj) -> pd.DataFrame:
        if isinstance(obj, pd.DataFrame):
            return obj
        if obj is None:
            return pd.DataFrame()
        # Common sanitized patterns: JSON string, dict-of-lists, list-of-records
        if isinstance(obj, str):
            # try JSON
            try:
                return pd.read_json(obj)
            except Exception:
                try:
                    # sometimes stored as repr(list[dict])
                    tmp = eval(obj, {"__builtins__": {}})  # only if trusted
                    return pd.DataFrame(tmp)
                except Exception:
                    return pd.DataFrame()
        if isinstance(obj, dict):
            # dict-of-lists or dict like DataFrame.to_dict()
            try:
                return pd.DataFrame(obj)
            except Exception:
                # dict with 'data' or 'records'
                for k in ("data", "records"):
                    if k in obj:
                        try:
                            return pd.DataFrame(obj[k])
                        except Exception:
                            pass
                return pd.DataFrame()
        if isinstance(obj, (list, tuple)):
            try:
                return pd.DataFrame(obj)
            except Exception:
                return pd.DataFrame()
        return pd.DataFrame()

    enrichment_dict = {str(k): _to_df(v) for k, v in (enrichment_dict or {}).items()}

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

    # 5. Build summaries
    summaries: Dict[str, str] = {}
    lines: List[str] = []

    # Helper to detect term and p-adj columns robustly
    def _detect_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
        term_candidates = [
            "Term",
            "term_name",
            "Description",
            "Name",
            "Pathway",
            "term",
        ]
        pval_candidates = [
            "Adjusted P-value",
            "Adjusted P-value (Benjamini-Hochberg)",
            "Adj P-value",
            "p.adjust",
            "padj",
            "FDR",
            "FDR q-value",
            "qvalue",
        ]
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
        enr_df = enrichment_dict.get(str(g), pd.DataFrame())

        if isinstance(enr_df, pd.DataFrame) and not enr_df.empty:
            term_col, pval_col = _detect_cols(enr_df)
            if term_col and pval_col:
                tmp = enr_df.copy()
                tmp[pval_col] = pd.to_numeric(tmp[pval_col], errors="coerce")
                sig = tmp[tmp[pval_col] < float(enrichment_padj_cutoff)]
                if not sig.empty:
                    top_terms = (
                        sig.sort_values(pval_col, ascending=True)[term_col]
                        .head(n_terms)
                        .astype(str)
                        .tolist()
                    )
                else:
                    log.info(
                        f"Cluster {g}: No pathways found below p_adj cutoff of {enrichment_padj_cutoff}"
                    )
            else:
                log.warning(
                    f"Cluster {g}: Could not detect term/padj columns. Columns present: {list(enr_df.columns)}"
                )

        title = f"### Cluster {g}"
        mk_str = f"**Top Markers**: {', '.join(top_genes) if top_genes else 'N/A'}"
        pt_str = f"**Top Pathways**: {', '.join(top_terms) if top_terms else 'N/A'}"

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
