"""
Differential Expression and Enrichment Analysis Module for scLucid.

This module provides comprehensive differential expression analysis and functional
enrichment capabilities for single-cell RNA-seq data, combining the robustness of
enterprise-grade engineering with user-friendly interfaces.

Key Features:
- Marker gene discovery (one-vs-rest, pairwise comparisons)
- Conserved marker identification across conditions
- Flexible filtering with detailed parameter tracking
- ORA and GSEA enrichment analysis (online/offline modes)
- Publication-quality visualizations
- Batch processing workflows
- Complete result traceability and HDF5 compatibility

All results are stored in adata.uns['sclucid']['analysis']['de'] with full
parameter provenance for reproducibility.
"""

import dataclasses
import logging
import pickle
import re
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text
from anndata import AnnData
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from ..utils import sanitize_for_hdf5
from .config import (
    BaseConfig,
    CompareConditionsConfig,
    CompareGroupsConfig,
    ConservedMarkersConfig,
    DifferentialConfig,
    EnrichmentConfig,
    FilterMarkersConfig,
)

log = logging.getLogger(__name__)

__all__ = [
    # Core DE functions
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    # Enrichment
    "run_enrichment",
    # Batch analysis
    "batch_celltype_deg_enrichment",
    # Advanced analysis
    "characterize_clusters",
    "summarize_markers_and_enrichment",
    # Visualization
    "visualize_markers",
    "plot_volcano",
    "plot_multi_cluster_deg",
    # Utilities
    "export_enrichment_results",
    "ResultManager",
    "save_results",
    "load_results",
]


# ==================== Helper Functions ====================


def _safe_filename(s: str) -> str:
    """Convert string to filesystem-safe filename."""
    return re.sub(r"[^0-9a-zA-Z\u4e00-\u9fa5\-_\.]", "_", s)


def _is_0_to_1(series: pd.Series) -> bool:
    """
    Check if a pandas Series is on a 0-1 scale.

    Args:
        series: Input Series

    Returns:
        True if all non-null values are in [0, 1]
    """
    s = series.dropna()
    return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))


def _to_frac(series: pd.Series) -> pd.Series:
    """
    Convert Series from 0-100 or 0-1 range to 0-1 fraction scale.

    Handles:
    - None/empty Series
    - Already scaled 0-1 data
    - Percentage scale (0-100)
    - Invalid values (converted to NaN)

    Args:
        series: Input percentage/fraction Series

    Returns:
        Series scaled to [0, 1] range
    """
    if series is None:
        return pd.Series(dtype=float)

    s_numeric = pd.to_numeric(series, errors="coerce")

    if _is_0_to_1(s_numeric):
        return s_numeric

    log.debug("Detected percentage scale (0-100). Converting to fraction (0-1).")
    return s_numeric.clip(lower=0, upper=100) / 100.0


def _standardize_pct_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize percentage column names across different Scanpy versions.

    Converts:
    - pct.1 / pct_nz_group -> pct_nz_group (0-1 scale)
    - pct.2 / pct_nz_reference -> pct_nz_reference (0-1 scale)

    Args:
        df: DataFrame with percentage columns

    Returns:
        DataFrame with standardized column names and values
    """
    df = df.copy()

    # Handle group column
    if "pct_nz_group" in df.columns:
        df["pct_nz_group"] = _to_frac(df["pct_nz_group"])
    elif "pct.1" in df.columns:
        df["pct_nz_group"] = _to_frac(df["pct.1"])
        df.drop(columns=["pct.1"], inplace=True)

    # Handle reference column
    if "pct_nz_reference" in df.columns:
        df["pct_nz_reference"] = _to_frac(df["pct_nz_reference"])
    elif "pct.2" in df.columns:
        df["pct_nz_reference"] = _to_frac(df["pct.2"])
        df.drop(columns=["pct.2"], inplace=True)

    return df


def _standardize_enrichment_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize enrichment result column names across GSEApy/Enrichr outputs.

    Args:
        df: Enrichment results DataFrame

    Returns:
        DataFrame with standardized column names
    """
    if df is None or df.empty:
        return pd.DataFrame()

    df = df.copy()

    rename_map = {
        "Adjusted P-value": "pval_adj",
        "P-value": "pval",
        "Overlap": "overlap",
        "Genes": "genes",
        "Gene_set": "gene_set",
        "Term": "term",
        "NES": "nes",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer_pval",
    }

    df.rename(columns=rename_map, inplace=True)

    # Ensure pval_adj exists
    if "pval_adj" not in df.columns and "fdr" in df.columns:
        df["pval_adj"] = df["fdr"]

    return df


def _store_results(
    adata: AnnData,
    key: str,
    results: Union[pd.DataFrame, Dict],
    config: BaseConfig,
    result_type: str = "de",
) -> None:
    """
    Store analysis results to adata.uns with full provenance.

    Creates structure:
    adata.uns['sclucid']['analysis'][result_type][key] = results
    adata.uns['sclucid']['analysis'][result_type][f'{key}_config'] = config_dict

    Args:
        adata: AnnData object
        key: Storage key for results
        results: Results DataFrame or dictionary
        config: Configuration object used
        result_type: Type of analysis ('de', 'enrichment', etc.)
    """
    # Initialize nested structure
    if "sclucid" not in adata.uns:
        adata.uns["sclucid"] = {}
    if "analysis" not in adata.uns["sclucid"]:
        adata.uns["sclucid"]["analysis"] = {}
    if result_type not in adata.uns["sclucid"]["analysis"]:
        adata.uns["sclucid"]["analysis"][result_type] = {}

    # Store results
    adata.uns["sclucid"]["analysis"][result_type][key] = results

    # Store configuration with HDF5 compatibility
    config_dict = {
        k: v
        for k, v in config.__dict__.items()
        if not k.startswith("_") and v is not None
    }
    adata.uns["sclucid"]["analysis"][result_type][f"{key}_config"] = sanitize_for_hdf5(
        config_dict
    )

    if config.verbose:
        log.info(
            f"Results stored: adata.uns['sclucid']['analysis']['{result_type}']['{key}']"
        )


# ==================== Core Differential Expression Functions ====================


def find_markers(
    adata: AnnData,
    config: Optional[DifferentialConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes using one-vs-rest differential expression analysis.

    This function:
    1. Runs Scanpy's rank_genes_groups
    2. Extracts and standardizes results for all groups
    3. Stores both raw and processed results
    4. Returns complete unfiltered DataFrame

    Subsequent filtering should be done via filter_markers() for flexibility.

    Args:
        adata: AnnData object
        config: DifferentialConfig object with analysis parameters
        **kwargs: Additional parameters to override config

    Returns:
        Complete marker gene DataFrame with columns:
        - names: Gene names
        - scores: Statistical scores
        - logfoldchanges: Log2 fold changes
        - pvals: P-values
        - pvals_adj: Adjusted p-values
        - pct_nz_group: Expression % in group (0-1 scale)
        - pct_nz_reference: Expression % in reference (0-1 scale)
        - group: Cluster/group identifier

    Example:
        >>> config = DifferentialConfig(
        ...     groupby="leiden",
        ...     method="wilcoxon",
        ...     use_raw=True
        ... )
        >>> markers = find_markers(adata, config)
        >>> # Then filter:
        >>> filter_config = FilterMarkersConfig(min_log2fc=1.0, max_padj=0.01)
        >>> filtered = filter_markers(adata, filter_config)

    Notes:
        - Results stored at: adata.uns['sclucid']['analysis']['de']['{key}_df']
        - Raw scanpy output at: adata.uns['sclucid']['analysis']['de']['{key}']
        - Parameters at: adata.uns['sclucid']['analysis']['de']['{key}_params']
    """
    # Configuration handling with immutability
    if config is None:
        active_config = DifferentialConfig(**kwargs)
    else:
        active_config = dataclasses.replace(config)
        for k, v in kwargs.items():
            if hasattr(active_config, k):
                setattr(active_config, k, v)

    groupby = active_config.groupby
    key_added = active_config.key_added or "rank_genes_groups"

    if active_config.verbose:
        log.info(
            f"Finding markers: groupby='{groupby}', method='{active_config.method}'"
        )

    # Build Scanpy parameters
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
        log.info(f"Analyzing subset of groups: {active_config.groups}")

    # Run differential expression
    sc.tl.rank_genes_groups(adata, **rank_genes_params)

    # Robust result extraction with validation
    if key_added not in adata.uns:
        raise KeyError(
            f"Scanpy returned no result at adata.uns['{key_added}']. "
            "This may indicate an issue with the input data or parameters."
        )

    raw = adata.uns[key_added]
    if "names" not in raw:
        raise KeyError(
            f"Scanpy result missing 'names' field at adata.uns['{key_added}']. "
            "The structure may have changed in newer Scanpy versions."
        )

    names_field = raw["names"]
    if not hasattr(names_field, "dtype") or names_field.dtype.names is None:
        raise ValueError(
            "Scanpy 'names' field lacks structured dtype. "
            "Cannot extract group-wise results."
        )

    groups_tested = names_field.dtype.names
    result_dfs: List[pd.DataFrame] = []

    # Extract results for each group
    for group in groups_tested:
        df = sc.get.rank_genes_groups_df(adata, key=key_added, group=group)
        if df.empty:
            log.warning(f"No results for group '{group}'")
            continue

        # Harmonize column names (Scanpy version compatibility)
        if "pct_nz_group" not in df.columns and "pct_nz" in df.columns:
            df = df.rename(columns={"pct_nz": "pct_nz_group"})
            log.debug("Renamed 'pct_nz' to 'pct_nz_group' for compatibility")

        df["group"] = group

        # Optional in-function filtering (light touch)
        if active_config.pval_cutoff is not None and "pvals_adj" in df.columns:
            before = len(df)
            df = df[df["pvals_adj"] <= float(active_config.pval_cutoff)].copy()
            if active_config.verbose:
                log.info(
                    f"Group '{group}': p_adj <= {active_config.pval_cutoff} "
                    f"retained {len(df)}/{before} genes"
                )

        if active_config.fold_change_max is not None and "logfoldchanges" in df.columns:
            df["logfoldchanges"] = df["logfoldchanges"].clip(
                upper=float(active_config.fold_change_max)
            )

        result_dfs.append(df)

    # Combine and standardize
    if not result_dfs:
        log.warning("No valid marker results found for any group after filtering")
        full_df = pd.DataFrame()
    else:
        full_df = pd.concat(result_dfs, ignore_index=True)
        full_df = _standardize_pct_columns(full_df)

    # Store with provenance
    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )

    root[key_added] = adata.uns[key_added]  # Raw Scanpy output
    df_key = f"{key_added}_df"
    root[df_key] = full_df  # Processed DataFrame

    # Parameter tracking
    params = active_config.to_dict()
    params["scanpy_version"] = getattr(sc, "__version__", "unknown")
    root[f"{key_added}_params"] = sanitize_for_hdf5(params)

    if active_config.verbose:
        log.info(
            f"Found {len(full_df)} total markers across {len(groups_tested)} groups"
        )
        log.info(f"Results stored at .uns['sclucid']['analysis']['de']['{df_key}']")
        log.info("Use filter_markers() for advanced filtering")

    return full_df


def filter_markers(
    adata: AnnData,
    config: FilterMarkersConfig,
) -> pd.DataFrame:
    """
    Filter marker genes with comprehensive criteria and detailed logging.

    Supports filtering by:
    - Statistical significance (p-value)
    - Effect size (log fold change)
    - Expression prevalence (% cells expressing)
    - Specificity (difference in % between groups)
    - Top N selection per group

    Args:
        adata: AnnData object containing DE results
        config: FilterMarkersConfig with filtering parameters

    Returns:
        Filtered marker DataFrame

    Example:
        >>> config = FilterMarkersConfig(
        ...     key="rank_genes_groups",
        ...     min_log2fc=1.0,
        ...     max_padj=0.01,
        ...     min_in_group_pct=0.25,
        ...     max_out_group_pct=0.50,
        ...     min_diff_pct=0.15,
        ...     keep_top_n=50,
        ...     sort_by="scores"
        ... )
        >>> filtered_markers = filter_markers(adata, config)

    Notes:
        - Input: adata.uns['sclucid']['analysis']['de']['{key}_df']
        - Output: adata.uns['sclucid']['analysis']['de']['{key_added}']
    """
    key = config.key
    key_added = config.key_added or f"{key}_filtered_df"
    df_key = f"{key}_df"

    # Load source data
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if df_key not in root:
        raise KeyError(
            f"Source DataFrame not found at "
            f".uns['sclucid']['analysis']['de']['{df_key}']. "
            "Run find_markers() first."
        )

    df = root[df_key].copy()
    if df.empty:
        log.warning("Source marker DataFrame is empty. Returning empty DataFrame.")
        return pd.DataFrame()

    # Validate required columns
    required_cols = ["logfoldchanges", "pvals_adj", "pct_nz_group", "group", "names"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Marker DataFrame missing required columns: {missing}")

    # Standardize percentage columns
    df = _standardize_pct_columns(df)

    has_ref = "pct_nz_reference" in df.columns
    pct_group_frac = df["pct_nz_group"]
    pct_ref_frac = df["pct_nz_reference"] if has_ref else None

    log.info(f"Filtering markers from '{df_key}'...")
    filt = pd.Series(True, index=df.index)

    # Filter 1: Log2 Fold Change
    if config.min_log2fc is not None:
        lfc = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        if config.use_abs_log2fc:
            keep = lfc.abs() >= float(config.min_log2fc)
        else:
            keep = lfc >= float(config.min_log2fc)
        keep = keep.fillna(False)

        log.info(
            f"[Filter] log2FC {'|x|' if config.use_abs_log2fc else ''} "
            f">= {config.min_log2fc}: kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.debug("[Filter] min_log2fc: skipped (None)")

    # Filter 2: Adjusted P-value
    if config.max_padj is not None:
        padj = pd.to_numeric(df["pvals_adj"], errors="coerce")
        keep = (padj <= float(config.max_padj)).fillna(False)

        log.info(
            f"[Filter] adj_p <= {config.max_padj}: kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.debug("[Filter] max_padj: skipped (None)")

    # Filter 3: In-group expression prevalence
    if config.min_in_group_pct is not None:
        keep = pct_group_frac >= float(config.min_in_group_pct)

        log.info(
            f"[Filter] pct_in_group >= {config.min_in_group_pct:.3f}: "
            f"kept {int(keep.sum())}/{len(filt)}"
        )
        filt &= keep
    else:
        log.debug("[Filter] min_in_group_pct: skipped (None)")

    # Filter 4 & 5: Specificity filters (require reference group)
    if has_ref:
        if config.max_out_group_pct is not None:
            keep = pct_ref_frac <= float(config.max_out_group_pct)

            log.info(
                f"[Filter] pct_out_group <= {config.max_out_group_pct:.3f}: "
                f"kept {int(keep.sum())}/{len(filt)}"
            )
            filt &= keep
        else:
            log.debug("[Filter] max_out_group_pct: skipped (None)")

        if config.min_diff_pct is not None:
            diff_pct = pct_group_frac - pct_ref_frac
            keep = diff_pct >= float(config.min_diff_pct)

            log.info(
                f"[Filter] (pct_in - pct_out) >= {config.min_diff_pct:.3f}: "
                f"kept {int(keep.sum())}/{len(filt)}"
            )
            filt &= keep
        else:
            log.debug("[Filter] min_diff_pct: skipped (None)")
    else:
        if config.max_out_group_pct is not None or config.min_diff_pct is not None:
            log.warning(
                "'pct_nz_reference' not found; specificity-related filters skipped"
            )

    filtered_df = df[filt].copy()
    log.info(f"Retained {len(filtered_df)} genes after all filters")

    # Post-filter: Keep top N per group
    if (
        config.keep_top_n is not None
        and config.keep_top_n > 0
        and not filtered_df.empty
    ):
        sort_by_col = config.sort_by

        # Handle special case: diff_pct
        if sort_by_col == "diff_pct":
            if has_ref:
                filtered_df["diff_pct"] = pct_group_frac[filt] - pct_ref_frac[filt]
            else:
                log.warning(
                    "Cannot sort by 'diff_pct' without 'pct_nz_reference'. "
                    "Falling back to 'scores'"
                )
                sort_by_col = "scores"

        # Fallback if sort column missing
        if sort_by_col not in filtered_df.columns:
            fallback_col = (
                "logfoldchanges"
                if "logfoldchanges" in filtered_df.columns
                else "scores"
                if "scores" in filtered_df.columns
                else filtered_df.columns[0]
            )
            log.warning(
                f"Sort key '{config.sort_by}' not found. "
                f"Falling back to '{fallback_col}'"
            )
            sort_by_col = fallback_col

        log.info(
            f"Selecting top {config.keep_top_n} genes per group, "
            f"sorted by '{sort_by_col}'"
        )

        parts = []
        for g in filtered_df["group"].unique():
            sub = filtered_df[filtered_df["group"] == g].sort_values(
                sort_by_col, ascending=False
            )
            parts.append(sub.head(config.keep_top_n))

        filtered_df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    # Store results
    root[key_added] = filtered_df
    root[f"{key_added}_params"] = {**config.to_dict(), "n_retained": len(filtered_df)}

    log.info(
        f"Final filtered markers: {len(filtered_df)} rows -> "
        f".uns['sclucid']['analysis']['de']['{key_added}']"
    )

    return filtered_df


def compare_groups(
    adata: AnnData, config: Optional[CompareGroupsConfig] = None, **kwargs
) -> pd.DataFrame:
    """
    Compare two specific groups (e.g., cell types, conditions) for DE genes.

    Combines the robustness of careful validation with the convenience of
    integrated visualization.

    Args:
        adata: AnnData object
        config: CompareGroupsConfig with comparison parameters
        **kwargs: Parameters to override config

    Returns:
        Filtered DE genes DataFrame with top N up/down-regulated genes

    Example:
        >>> config = CompareGroupsConfig(
        ...     groupby="celltype",
        ...     group1="CD4_T",
        ...     group2="CD8_T",
        ...     min_log2fc=1.0,
        ...     max_padj=0.01,
        ...     n_top_genes=50,
        ...     plot=True
        ... )
        >>> degs = compare_groups(adata, config)

    Notes:
        - Automatically selects top N up and down-regulated genes
        - Optionally generates volcano plot
        - Results stored with full parameter provenance
    """
    # Configuration handling
    if config is None:
        config = CompareGroupsConfig(**kwargs)
    else:
        config = dataclasses.replace(config, **kwargs)

    groupby = config.groupby
    group1 = config.group1
    group2 = config.group2
    key_added = config.key_added or f"compare_{group1}_vs_{group2}".replace(" ", "_")

    if config.verbose:
        log.info(f"Comparing DE genes: '{group1}' vs '{group2}' in '{groupby}'")

    # Input validation
    if groupby not in adata.obs.columns:
        raise KeyError(f"Column '{groupby}' not found in adata.obs")

    subset_mask = adata.obs[groupby].isin([group1, group2])
    if subset_mask.sum() == 0:
        raise ValueError(
            f"No cells found for either '{group1}' or '{group2}' in '{groupby}'"
        )

    # Create temporary subset with standardized group labels
    temp_adata = adata[subset_mask].copy()
    temp_adata.obs["_cmp_grp"] = (
        temp_adata.obs[groupby].map({group1: "grp1", group2: "grp2"}).astype("category")
    )

    # Run differential expression
    sc.tl.rank_genes_groups(
        temp_adata,
        groupby="_cmp_grp",
        groups=["grp1"],
        reference="grp2",
        method=config.method,
        layer=config.layer,
        use_raw=config.use_raw,
        pts=True,
        tie_correct=True,
    )

    # Extract and standardize results
    df = sc.get.rank_genes_groups_df(temp_adata, group="grp1")
    df = _standardize_pct_columns(df)

    # Apply filters
    lfc = pd.to_numeric(df["logfoldchanges"], errors="coerce")
    padj = pd.to_numeric(df["pvals_adj"], errors="coerce")
    pct_in = _to_frac(df.get("pct_nz_group", pd.Series(1, index=df.index)))

    filt = (
        (lfc.abs() >= float(config.min_log2fc))
        & (padj <= float(config.max_padj))
        & (pct_in >= float(config.min_pct))
    )

    filtered = df[filt].copy()

    # Select top N up and down-regulated genes
    up = filtered[filtered["logfoldchanges"] > 0].head(config.n_top_genes)
    down = (
        filtered[filtered["logfoldchanges"] < 0]
        .sort_values("logfoldchanges", ascending=True)
        .head(config.n_top_genes)
    )
    final = pd.concat([up, down], ignore_index=True)

    if config.verbose:
        log.info(f"Found {len(final)} DE genes ({len(up)} up, {len(down)} down)")

    # Store results
    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[key_added] = final
    root[f"{key_added}_params"] = sanitize_for_hdf5(config.to_dict())

    log.info(f"Results stored at .uns['sclucid']['analysis']['de']['{key_added}']")

    # Optional visualization
    if config.plot and not filtered.empty:
        save_path = None
        if config.save_dir:
            Path(config.save_dir).mkdir(parents=True, exist_ok=True)
            save_path = str(Path(config.save_dir) / f"{key_added}_volcano.pdf")

        plot_volcano(
            filtered,
            title=f"{group1} vs {group2}",
            subtitle=f"Differential Expression Analysis (n={temp_adata.n_obs} cells)",
            top_n_up=config.n_top_genes,
            top_n_down=config.n_top_genes,
            lfc_threshold=config.min_log2fc,
            pval_threshold=config.max_padj,
            savepath=save_path,
        )

    return final


def compare_conditions(
    adata: AnnData, config: Optional[CompareConditionsConfig] = None, **kwargs
) -> pd.DataFrame:
    """
    Compare two conditions within a specific cell type/group.

    This is a specialized wrapper around compare_groups() that first
    subsets to a specific cell type, then compares conditions.

    Args:
        adata: AnnData object
        config: CompareConditionsConfig
        **kwargs: Override parameters

    Returns:
        Filtered DE genes DataFrame

    Example:
        >>> config = CompareConditionsConfig(
        ...     groupby="celltype",
        ...     group_name="T_cells",
        ...     condition_key="treatment",
        ...     condition1="Treated",
        ...     condition2="Control",
        ...     min_log2fc=0.5,
        ...     max_padj=0.05
        ... )
        >>> degs = compare_conditions(adata, config)
    """
    if config is None:
        config = CompareConditionsConfig(**kwargs)
    else:
        config = dataclasses.replace(config, **kwargs)

    groupby = config.groupby
    group_name = config.group_name
    condition_key = config.condition_key

    log.info(
        f"Comparing conditions '{config.condition1}' vs '{config.condition2}' "
        f"within '{group_name}'"
    )

    # Validate group exists
    if group_name not in adata.obs[groupby].unique():
        raise ValueError(f"Group '{group_name}' not found in adata.obs['{groupby}']")

    # Subset to specific cell type
    adata_subset = adata[adata.obs[groupby] == group_name].copy()

    # Create comparison config
    comp_config = CompareGroupsConfig(
        groupby=condition_key,
        group1=config.condition1,
        group2=config.condition2,
        method=config.method,
        min_log2fc=config.min_log2fc,
        max_padj=config.max_padj,
        min_pct=config.min_pct,
        n_top_genes=config.n_top_genes,
        layer=config.layer,
        use_raw=config.use_raw,
        key_added=config.key_added
        or (
            f"compare_{config.condition1}_vs_{config.condition2}_"
            f"in_{_safe_filename(group_name)}"
        ),
        plot=config.plot,
        save_dir=config.save_dir,
        verbose=config.verbose,
    )

    # Run comparison on subset
    results_df = compare_groups(adata_subset, config=comp_config)

    # Add metadata
    results_df["celltype"] = group_name

    # Store in parent AnnData
    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[comp_config.key_added] = results_df
    root[f"{comp_config.key_added}_params"] = sanitize_for_hdf5(config.to_dict())

    log.info(
        f"Condition comparison complete: {len(results_df)} DE genes. "
        f"Stored at .uns['...']['{comp_config.key_added}']"
    )

    return results_df


def get_conserved_markers(
    adata: AnnData, config: Optional[ConservedMarkersConfig] = None, **kwargs
) -> Dict[str, pd.DataFrame]:
    """
    Find markers for each group that are conserved across multiple conditions.

    Strategy:
    1. For each group, run DE in each condition separately
    2. Identify genes that are significant in >= min_conditions
    3. Aggregate statistics (mean/min/max fold change, etc.)
    4. Return conserved markers per group

    Args:
        adata: AnnData object
        config: ConservedMarkersConfig
        **kwargs: Override parameters

    Returns:
        Dictionary mapping group names to conserved marker DataFrames

    Example:
        >>> config = ConservedMarkersConfig(
        ...     groupby="celltype",
        ...     condition_key="batch",
        ...     min_log2fc=1.0,
        ...     max_padj=0.01,
        ...     min_cells=20,
        ...     min_conditions=2
        ... )
        >>> conserved = get_conserved_markers(adata, config)
        >>> # Access T cell conserved markers:
        >>> t_cell_markers = conserved["T_cells"]

    Notes:
        - Useful for finding robust markers across batches/samples
        - Requires sufficient cells per group per condition
    """
    if config is None:
        config = ConservedMarkersConfig(**kwargs)
    else:
        config = dataclasses.replace(config, **kwargs)

    key_added = config.key_added or (
        f"conserved_markers_{config.groupby}_{config.condition_key}"
    )

    # Validate columns exist
    if config.condition_key not in adata.obs.columns:
        raise KeyError(f"Condition key '{config.condition_key}' not in adata.obs")
    if config.groupby not in adata.obs.columns:
        raise KeyError(f"Groupby key '{config.groupby}' not in adata.obs")

    # Get unique conditions and groups
    conditions = (
        list(adata.obs[config.condition_key].cat.categories)
        if pd.api.types.is_categorical_dtype(adata.obs[config.condition_key])
        else list(pd.unique(adata.obs[config.condition_key]))
    )

    groups = (
        list(adata.obs[config.groupby].cat.categories)
        if pd.api.types.is_categorical_dtype(adata.obs[config.groupby])
        else list(pd.unique(adata.obs[config.groupby]))
    )

    min_conditions = config.min_conditions or max(1, len(conditions) - 1)

    if config.verbose:
        log.info(
            f"Finding conserved markers across {len(conditions)} conditions, "
            f"requiring significance in >={min_conditions} conditions"
        )

    conserved_markers: Dict[str, pd.DataFrame] = {}

    for group in groups:
        markers_per_condition = []

        # Run DE in each condition
        for cond in conditions:
            subset = adata[
                (adata.obs[config.groupby] == group)
                & (adata.obs[config.condition_key] == cond)
            ]

            if subset.n_obs < config.min_cells:
                log.debug(
                    f"Skip group '{group}' in condition '{cond}': "
                    f"n_cells={subset.n_obs} < {config.min_cells}"
                )
                continue

            # Create temporary adata for this condition
            temp_adata = adata[adata.obs[config.condition_key] == cond].copy()

            if group not in temp_adata.obs[config.groupby].unique():
                continue

            # Run DE
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
            df = _standardize_pct_columns(df)

            # Filter for significance
            pct_in = _to_frac(df.get("pct_nz_group", pd.Series(1, index=df.index)))

            df = df[
                (df["logfoldchanges"] >= float(config.min_log2fc))
                & (df["pvals_adj"] <= float(config.max_padj))
                & (pct_in >= float(config.min_pct))
            ].copy()

            if df.empty:
                continue

            df["condition"] = cond
            markers_per_condition.append(df)

        # Check if we have enough conditions
        if len(markers_per_condition) < min_conditions:
            log.info(
                f"Group '{group}': insufficient conditions with markers "
                f"({len(markers_per_condition)} < {min_conditions}). Skipping."
            )
            continue

        # Combine results across conditions
        full_df = pd.concat(markers_per_condition, ignore_index=True)

        # Count how many conditions each gene appears in
        gene_counts = full_df.groupby("names").size()
        conserved_genes = gene_counts[gene_counts >= min_conditions].index.tolist()

        conserved_df = full_df[full_df["names"].isin(conserved_genes)]

        # Aggregate statistics
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

        if config.verbose:
            log.info(
                f"Group '{group}': found {len(agg_df)} conserved markers "
                f"across {agg_df['n_conditions'].min()}-{agg_df['n_conditions'].max()} conditions"
            )

    # Store results
    root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    root[key_added] = sanitize_for_hdf5(
        {
            "aggregates": conserved_markers,
            "params": config.to_dict(),
        }
    )

    log.info(
        f"Conserved marker analysis complete: {len(conserved_markers)} groups. "
        f"Results stored at .uns['...']['{key_added}']"
    )

    return conserved_markers


# ==================== Enrichment Analysis ====================


def run_enrichment(
    adata: AnnData,
    groupby: Optional[str] = None,
    config: Optional[EnrichmentConfig] = None,
    **kwargs,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Run functional enrichment analysis (ORA and/or GSEA) using GSEApy.

    Supports:
    - Online mode: Direct API calls to Enrichr/GSEA servers
    - Offline mode: Local GMT files (faster, more stable)
    - Multiple gene set categories per analysis
    - Both ORA (Over-Representation) and GSEA (Prerank) methods

    Args:
        adata: AnnData object with DE results
        groupby: Optional grouping column (usually inferred from DE results)
        config: EnrichmentConfig object
        **kwargs: Override config parameters

    Returns:
        Nested dictionary:
        {
            "cluster1": {
                "ora": DataFrame with ORA results,
                "gsea": DataFrame with GSEA results
            },
            "cluster2": {...},
            ...
        }

    Example:
        >>> config = EnrichmentConfig(
        ...     de_key="rank_genes_groups_filtered_df",
        ...     mode="offline",
        ...     method="both",
        ...     organism="human",
        ...     gene_sets_offline=["hallmark", "go_bp", "reactome"]
        ... )
        >>> enr_results = run_enrichment(adata, groupby="leiden", config=config)
        >>>
        >>> # Access cluster 0 ORA results:
        >>> cluster0_ora = enr_results["0"]["ora"]

    Notes:
        - Offline mode requires GMT files in scLucid/resources/
        - Results stored at: adata.uns['sclucid']['analysis']['de']['{key_added}']
        - Background gene list automatically set to all genes in adata
    """
    # Configuration
    if config is None:
        config = EnrichmentConfig(**kwargs)
    else:
        config = dataclasses.replace(config, **kwargs)

    # Load DE results
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if config.de_key not in root:
        raise KeyError(
            f"DE results not found at .uns['sclucid']['analysis']['de']['{config.de_key}']. "
            "Run find_markers() or filter_markers() first."
        )

    marker_df = root[config.de_key]
    if marker_df.empty:
        log.warning("Marker DataFrame is empty. Skipping enrichment analysis.")
        return {}

    # Determine groups
    if groupby and "group" in marker_df.columns:
        group_order = list(pd.unique(marker_df["group"]))
    else:
        group_order = ["all"]
        marker_df = marker_df.copy()
        marker_df["group"] = "all"

    if config.verbose:
        log.info(
            f"Running {config.method.upper()} enrichment for {len(group_order)} groups "
            f"in {config.mode} mode"
        )

    # Background genes
    background_genes = list(adata.var_names)

    # Prepare gene sets
    gmt_files_to_run = {}
    gene_sets_list = (
        config.gene_sets_online if config.mode == "online" else config.gene_sets_offline
    )

    if not isinstance(gene_sets_list, list):
        gene_sets_list = [gene_sets_list]

    if config.mode == "offline":
        # Custom GMT file takes priority
        if config.custom_gene_sets and Path(config.custom_gene_sets).is_file():
            gmt_files_to_run = {"custom": config.custom_gene_sets}
            log.info(f"Using custom gene set file: {config.custom_gene_sets}")
        else:
            # Load built-in GMT files
            for gs_category in gene_sets_list:
                try:
                    filename = f"{config.organism.lower()}_{gs_category}_{config.gmt_version}.gmt"
                    file_path = resources.files("scLucid").joinpath(
                        "resources", filename
                    )

                    if file_path.is_file():
                        gmt_files_to_run[gs_category] = str(file_path)
                        log.debug(f"Loaded GMT file: {filename}")
                    else:
                        log.warning(f"GMT file not found: {filename}")

                except Exception as e:
                    log.error(f"Error loading GMT file for '{gs_category}': {e}")

            if not gmt_files_to_run:
                raise FileNotFoundError(
                    f"No valid GMT files found for offline mode. "
                    f"Searched for: {gene_sets_list}"
                )
    else:
        # Online mode: gene set names are used directly
        gmt_files_to_run = {gs: gs for gs in gene_sets_list}

    # Select GSEA ranking column
    rank_col = config.rank_col_gsea
    if config.prefer_score_for_enrichment and "scores" in marker_df.columns:
        rank_col = "scores"
        log.debug("Using 'scores' for GSEA ranking (prefer_score_for_enrichment=True)")
    elif rank_col not in marker_df.columns:
        fallback = "scores" if "scores" in marker_df.columns else "logfoldchanges"
        log.warning(
            f"GSEA rank column '{rank_col}' not found. Falling back to '{fallback}'"
        )
        rank_col = fallback

    if rank_col not in marker_df.columns:
        raise KeyError(
            f"GSEA requires a ranking column ('{rank_col}') in the marker DataFrame"
        )

    # Run enrichment for each group
    enrichment_results: Dict[str, Dict[str, pd.DataFrame]] = {}

    for cluster in group_order:
        cluster_results: Dict[str, pd.DataFrame] = {}
        sub = marker_df[marker_df["group"] == cluster]

        if sub.empty:
            log.warning(
                f"Skipping '{cluster}': no marker genes found in '{config.de_key}'"
            )
            continue

        # === ORA (Over-Representation Analysis) ===
        if config.method in ["ora", "both"]:
            gene_list = (
                sub.sort_values(rank_col, ascending=False)["names"]
                .head(config.n_top_genes_ora)
                .astype(str)
                .tolist()
            )

            if len(gene_list) < config.min_genes_for_ora:
                log.warning(
                    f"Skipping ORA for '{cluster}': only {len(gene_list)} genes "
                    f"(< {config.min_genes_for_ora})"
                )
            else:
                all_ora_results = []

                for category, gmt in gmt_files_to_run.items():
                    try:
                        if config.mode == "online":
                            enr_ora = gp.enrichr(
                                gene_list=gene_list,
                                gene_sets=gmt,
                                organism=config.organism,
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0,
                            )
                        else:
                            enr_ora = gp.enrich(
                                gene_list=gene_list,
                                gene_sets=gmt,
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0,
                            )

                        res = enr_ora.results.copy()
                        res["gene_set_category"] = category
                        all_ora_results.append(res)

                    except Exception as e:
                        log.error(
                            f"ORA failed for cluster '{cluster}', "
                            f"category '{category}': {e}"
                        )

                if all_ora_results:
                    ora_df = pd.concat(all_ora_results, ignore_index=True)
                    ora_df = _standardize_enrichment_cols(ora_df)

                    if "pval_adj" in ora_df.columns:
                        ora_df = ora_df[ora_df["pval_adj"] < config.cutoff_pval]

                    cluster_results["ora"] = ora_df
                else:
                    cluster_results["ora"] = pd.DataFrame()

        # === GSEA (Gene Set Enrichment Analysis) ===
        if config.method in ["gsea", "both"]:
            # Build ranked gene list
            rnk = (
                sub.drop_duplicates(subset="names", keep="first")
                .set_index("names")[rank_col]
                .sort_values(ascending=False)
            )

            if rnk.empty:
                log.warning(f"Skipping GSEA for '{cluster}': no ranked genes available")
            else:
                all_gsea_results = []

                for category, gmt in gmt_files_to_run.items():
                    try:
                        gsea_res = gp.prerank(
                            rnk=rnk,
                            gene_sets=gmt,
                            permutation_num=config.gsea_permutations,
                            min_size=config.gsea_min_size,
                            max_size=config.gsea_max_size,
                            outdir=None,
                            seed=42,
                            processes=4,
                        )

                        res = gsea_res.res2d.copy()
                        res["gene_set_category"] = category
                        all_gsea_results.append(res)

                    except Exception as e:
                        log.error(
                            f"GSEA failed for cluster '{cluster}', "
                            f"category '{category}': {e}"
                        )

                if all_gsea_results:
                    gsea_df = pd.concat(all_gsea_results, ignore_index=True)
                    gsea_df = _standardize_enrichment_cols(gsea_df)

                    if "pval_adj" in gsea_df.columns:
                        gsea_df = gsea_df[gsea_df["pval_adj"] < config.cutoff_pval]

                    cluster_results["gsea"] = gsea_df
                else:
                    cluster_results["gsea"] = pd.DataFrame()

        enrichment_results[str(cluster)] = cluster_results

        # Save individual results
        if config.save_dir:
            Path(config.save_dir).mkdir(parents=True, exist_ok=True)
            for method, res_df in cluster_results.items():
                if not res_df.empty:
                    output_path = (
                        Path(config.save_dir)
                        / f"{_safe_filename(str(cluster))}_{method}_enrichment.csv"
                    )
                    res_df.to_csv(output_path, index=False)

    # Store results
    store_root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("de", {})
    )
    store_root[config.key_added] = {
        "results": enrichment_results,
        "params": sanitize_for_hdf5(config.to_dict()),
    }

    if config.verbose:
        log.info(
            f"Enrichment analysis complete: {len(enrichment_results)} groups analyzed. "
            f"Results stored at .uns['...']['{config.key_added}']"
        )

    return enrichment_results


def export_enrichment_results(
    adata: AnnData,
    enrichment_key: str = "enrichment",
    output_path: str = "enrichment_results.xlsx",
) -> None:
    """
    Export enrichment results to Excel with separate sheets per cluster/method.

    Args:
        adata: AnnData with enrichment results
        enrichment_key: Key for enrichment results in adata.uns
        output_path: Output Excel file path

    Example:
        >>> run_enrichment(adata, groupby="leiden", config=enr_config)
        >>> export_enrichment_results(adata, output_path="enrichment.xlsx")
    """
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})

    if enrichment_key not in root:
        raise KeyError(
            f"Enrichment results '{enrichment_key}' not found. "
            "Run run_enrichment() first."
        )

    enr_store = root[enrichment_key]
    if "results" not in enr_store:
        raise ValueError("Enrichment store missing 'results' key")

    enrichment_results = enr_store["results"]

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for cluster, methods_dict in enrichment_results.items():
            for method, df in methods_dict.items():
                if not df.empty:
                    sheet_name = f"{_safe_filename(cluster)}_{method}"[
                        :31
                    ]  # Excel limit
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

    log.info(f"Enrichment results exported to {output_path}")


# ==================== Batch Analysis ====================


def batch_celltype_deg_enrichment(
    adata: AnnData,
    celltype_col: str,
    condition_col: str,
    condition1: str,
    condition2: str,
    outdir: str,
    celltypes: Optional[List[str]] = None,
    min_cells: int = 20,
    # DE parameters
    de_method: str = "wilcoxon",
    min_log2fc: float = 0.5,
    max_padj: float = 0.05,
    min_pct: float = 0.1,
    # Enrichment parameters
    run_enrichment_analysis: bool = True,
    gene_sets: Optional[List[str]] = None,
    organism: str = "human",
    enrichment_mode: str = "online",
    # Visualization
    plot_volcano_charts: bool = True,
    plot_enrichment_charts: bool = True,
    # Other
    save_pickle: bool = True,
    verbose: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Batch DEG and enrichment analysis across multiple cell types.

    For each cell type:
    1. Subset cells
    2. Run condition comparison (condition1 vs condition2)
    3. Optional: Run enrichment on up/down-regulated genes
    4. Optional: Generate volcano plots
    5. Save results

    Args:
        adata: AnnData object
        celltype_col: Column in obs with cell type labels
        condition_col: Column in obs with condition labels
        condition1: First condition (numerator in fold change)
        condition2: Second condition (reference/denominator)
        outdir: Output directory for results
        celltypes: Specific cell types to analyze (None = all)
        min_cells: Minimum cells required per cell type

        de_method: DE test method
        min_log2fc: Minimum log2 fold change threshold
        max_padj: Maximum adjusted p-value threshold
        min_pct: Minimum expression percentage

        run_enrichment_analysis: Whether to run enrichment
        gene_sets: Gene sets for enrichment (default: GO_BP)
        organism: Species for enrichment
        enrichment_mode: 'online' or 'offline'

        plot_volcano_charts: Generate volcano plots
        plot_enrichment_charts: Generate enrichment bar plots

        save_pickle: Save complete results as pickle
        verbose: Verbose logging

    Returns:
        Dictionary mapping cell types to results:
        {
            "T_cells": {
                "degs": DataFrame,
                "sig_degs": DataFrame (filtered),
                "enr_up": Enrichr object,
                "enr_down": Enrichr object,
                "n_cells": int
            },
            ...
        }

    Example:
        >>> results = batch_celltype_deg_enrichment(
        ...     adata,
        ...     celltype_col="celltype",
        ...     condition_col="treatment",
        ...     condition1="Treated",
        ...     condition2="Control",
        ...     outdir="./deg_results",
        ...     min_log2fc=1.0,
        ...     max_padj=0.01
        ... )
    """
    # Create output directory
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Determine cell types
    if celltypes is None:
        celltypes = adata.obs[celltype_col].unique().tolist()

    if verbose:
        log.info(
            f"Batch DEG analysis: {len(celltypes)} cell types, "
            f"{condition1} vs {condition2}"
        )

    # Default gene sets
    if gene_sets is None:
        gene_sets = ["GO_Biological_Process_2023"]

    results = {}

    for celltype in celltypes:
        safe_celltype = _safe_filename(celltype)

        if verbose:
            log.info(f"\n{'=' * 60}")
            log.info(f"Processing: {celltype}")
            log.info(f"{'=' * 60}")

        # Subset data
        adata_sub = adata[adata.obs[celltype_col] == celltype].copy()
        adata_sub = adata_sub[
            adata_sub.obs[condition_col].isin([condition1, condition2])
        ]

        n_cells = adata_sub.n_obs

        if n_cells < min_cells:
            log.warning(f"Skipping {celltype}: only {n_cells} cells (< {min_cells})")
            continue

        # Differential expression
        de_config = CompareConditionsConfig(
            groupby=celltype_col,
            group_name=celltype,
            condition_key=condition_col,
            condition1=condition1,
            condition2=condition2,
            method=de_method,
            min_log2fc=min_log2fc,
            max_padj=max_padj,
            min_pct=min_pct,
            save_dir=outdir,
            plot=False,
            verbose=verbose,
        )

        try:
            degs = compare_conditions(adata_sub, config=de_config)
        except Exception as e:
            log.error(f"DE analysis failed for {celltype}: {e}")
            continue

        # Save DEG table
        deg_table_path = (
            Path(outdir) / f"{safe_celltype}_DEG_{condition1}_vs_{condition2}.csv"
        )
        degs.to_csv(deg_table_path, index=False)

        # Volcano plot
        if plot_volcano_charts and len(degs) > 0:
            volcano_path = (
                Path(outdir)
                / f"{safe_celltype}_volcano_{condition1}_vs_{condition2}.pdf"
            )

            try:
                plot_volcano(
                    degs_df=degs,
                    title=f"DEG: {celltype}",
                    subtitle=f"{condition1} vs {condition2} (n={n_cells} cells)",
                    lfc_threshold=min_log2fc,
                    pval_threshold=max_padj,
                    savepath=str(volcano_path),
                )
            except Exception as e:
                log.warning(f"Volcano plot failed for {celltype}: {e}")

        # Enrichment analysis
        enr_up, enr_down = None, None
        sig_degs = None

        if run_enrichment_analysis and len(degs) > 0:
            sig_degs = degs[
                (degs["pvals_adj"] < max_padj)
                & (degs["logfoldchanges"].abs() > min_log2fc)
            ]

            up_genes = sig_degs[sig_degs["logfoldchanges"] > min_log2fc][
                "names"
            ].tolist()

            down_genes = sig_degs[sig_degs["logfoldchanges"] < -min_log2fc][
                "names"
            ].tolist()

            # Up-regulated enrichment
            if len(up_genes) > 5:
                try:
                    enr_up = gp.enrichr(
                        gene_list=up_genes,
                        gene_sets=gene_sets,
                        organism=organism,
                        outdir=None,
                    )

                    enr_up_df = enr_up.results
                    up_path = (
                        Path(outdir)
                        / f"{safe_celltype}_enrichment_up_{condition1}_vs_{condition2}.csv"
                    )
                    enr_up_df.to_csv(up_path, index=False)

                    # Visualization
                    if plot_enrichment_charts and not enr_up_df.empty:
                        plt.figure(figsize=(10, 8))
                        gp.barplot(
                            enr_up_df,
                            title=f"{celltype} UP ({condition1})",
                            top_term=20,
                            cutoff=1,
                        )
                        plt.tight_layout()
                        plt.savefig(
                            Path(outdir)
                            / f"{safe_celltype}_enrichment_up_{condition1}_vs_{condition2}.pdf",
                            dpi=300,
                        )
                        plt.close()

                except Exception as e:
                    log.warning(f"Enrichment (up) failed for {celltype}: {e}")

            # Down-regulated enrichment
            if len(down_genes) > 5:
                try:
                    enr_down = gp.enrichr(
                        gene_list=down_genes,
                        gene_sets=gene_sets,
                        organism=organism,
                        outdir=None,
                    )

                    enr_down_df = enr_down.results
                    down_path = (
                        Path(outdir)
                        / f"{safe_celltype}_enrichment_down_{condition1}_vs_{condition2}.csv"
                    )
                    enr_down_df.to_csv(down_path, index=False)

                    # Visualization
                    if plot_enrichment_charts and not enr_down_df.empty:
                        plt.figure(figsize=(10, 8))
                        gp.barplot(
                            enr_down_df,
                            title=f"{celltype} DOWN ({condition2})",
                            top_term=20,
                            cutoff=1,
                        )
                        plt.tight_layout()
                        plt.savefig(
                            Path(outdir)
                            / f"{safe_celltype}_enrichment_down_{condition1}_vs_{condition2}.pdf",
                            dpi=300,
                        )
                        plt.close()

                except Exception as e:
                    log.warning(f"Enrichment (down) failed for {celltype}: {e}")

        # Collect results
        results[celltype] = {
            "degs": degs,
            "sig_degs": sig_degs,
            "enr_up": enr_up,
            "enr_down": enr_down,
            "n_cells": n_cells,
        }

    # Save summary
    if save_pickle:
        pickle_path = (
            Path(outdir) / f"all_DEG_enrichment_{condition1}_vs_{condition2}.pkl"
        )
        with open(pickle_path, "wb") as f:
            pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
        log.info(f"Complete results saved: {pickle_path}")

    # Summary statistics
    if verbose:
        log.info(f"\n{'=' * 60}")
        log.info("Summary:")
        log.info(f"{'=' * 60}")
        for ct, res in results.items():
            n_degs = len(res["degs"]) if res["degs"] is not None else 0
            log.info(f"{ct}: {n_degs} DEGs, {res['n_cells']} cells")

    return results


# ==================== Advanced Analysis Functions ====================


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Comprehensive cluster characterization: DE + enrichment in one step.

    This is a convenience function that:
    1. Runs find_markers() for each cluster
    2. Runs run_enrichment() on the markers
    3. Stores combined results for easy access

    Args:
        adata: AnnData object
        groupby: Column for clustering (e.g., 'leiden')
        de_config: Optional DE configuration
        enrichment_config: Optional enrichment configuration
        key_added: Key for storing combined results

    Returns:
        Modified AnnData with results in adata.uns[key_added]

    Example:
        >>> adata = characterize_clusters(
        ...     adata,
        ...     groupby="leiden",
        ...     de_config=DifferentialConfig(use_raw=True),
        ...     enrichment_config=EnrichmentConfig(mode="offline")
        ... )
        >>> # Access cluster 0 characterization:
        >>> cluster0 = adata.uns["cluster_characterization"]["results"]["0"]
        >>> cluster0_markers = cluster0["top_de_genes"]
        >>> cluster0_pathways = cluster0["enrichment_ora"]
    """
    log.info(f"Characterizing clusters in '{groupby}'...")

    # Set up DE config
    if de_config is None:
        de_config = DifferentialConfig(groupby=groupby, use_raw=True)
    else:
        de_config.groupby = groupby

    de_key = de_config.key_added or "rank_genes_groups"
    de_df_key = f"{de_key}_df"

    # Run DE
    find_markers(adata, config=de_config)

    # Set up enrichment config
    if enrichment_config is None:
        enrichment_config = EnrichmentConfig(
            de_key=de_df_key, mode="offline", method="ora"
        )
    else:
        enrichment_config.de_key = de_df_key

    # Run enrichment
    enrichment_results_dict = run_enrichment(
        adata, groupby=groupby, config=enrichment_config
    )

    # Combine results
    clusters = (
        adata.obs[groupby].cat.categories
        if pd.api.types.is_categorical_dtype(adata.obs[groupby])
        else pd.unique(adata.obs[groupby])
    )

    de_df = adata.uns["sclucid"]["analysis"]["de"][de_df_key]
    characterization_results = {}

    for cluster in clusters:
        enr_res = enrichment_results_dict.get(str(cluster), {})

        characterization_results[str(cluster)] = {
            "top_de_genes": de_df[de_df["group"] == cluster],
            "enrichment_ora": enr_res.get("ora", pd.DataFrame()),
            "enrichment_gsea": enr_res.get("gsea", pd.DataFrame()),
        }

    # Store combined results
    adata.uns[key_added] = {
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
    enrichment_method_to_summarize: Literal["ora", "gsea"] = "ora",
    n_markers: int = 25,
    n_terms: int = 10,
    summary_file: Optional[str] = None,
    sort_markers_by: str = "logfoldchanges",
    enrichment_padj_cutoff: float = 0.05,
) -> Dict[str, str]:
    """
    Generate human-readable Markdown summaries for AI/manual annotation.

    Creates per-cluster summaries with:
    - Top marker genes
    - Top enriched pathways

    Ideal for:
    - Feeding to LLMs for automated annotation
    - Quick manual review
    - Documentation

    Args:
        adata: AnnData with DE and enrichment results
        groupby: Grouping column
        markers_df: Optional marker DataFrame (auto-loaded if None)
        enrichment_dict: Optional enrichment dict (auto-loaded if None)
        markers_key: Key for markers in adata.uns
        enrichment_key: Key for enrichment in adata.uns
        enrichment_method_to_summarize: 'ora' or 'gsea'
        n_markers: Number of top markers to include
        n_terms: Number of top terms to include
        summary_file: Optional output file path
        sort_markers_by: Column for sorting markers
        enrichment_padj_cutoff: P-value cutoff for pathways

    Returns:
        Dictionary mapping cluster names to markdown summaries

    Example:
        >>> summaries = summarize_markers_and_enrichment(
        ...     adata,
        ...     groupby="leiden",
        ...     summary_file="cluster_summaries.md"
        ... )
        >>> print(summaries["0"])
        ### Cluster 0
        **Top Markers**: CD3D, CD3E, CD3G, IL7R, TRAC
        **Top ORA Pathways**: T cell activation, T cell differentiation, ...
    """
    # Load markers
    if markers_df is None:
        try:
            log.info(f"Auto-retrieving markers from '{markers_key}'")
            markers_df = adata.uns["sclucid"]["analysis"]["de"][markers_key]
        except KeyError:
            raise KeyError(
                f"Marker DataFrame not found at "
                f".uns['sclucid']['analysis']['de']['{markers_key}']"
            )

    if markers_df is None or markers_df.empty:
        log.warning("Marker DataFrame is empty. Cannot summarize.")
        markers_df = pd.DataFrame(columns=["group", "names", "logfoldchanges"])

    # Load enrichment
    if enrichment_dict is None:
        enrichment_dict = {}
        try:
            log.info(f"Auto-retrieving enrichment from '{enrichment_key}'")
            enr_store = adata.uns["sclucid"]["analysis"]["de"].get(enrichment_key, {})

            if isinstance(enr_store, dict) and "results" in enr_store:
                enrichment_dict = enr_store["results"]
            else:
                enrichment_dict = enr_store

        except KeyError:
            log.warning("No enrichment data found. Pathways will be 'N/A'")

    # Determine group order
    if groupby in adata.obs:
        group_order = list(pd.unique(adata.obs[groupby].astype(str)))
    else:
        group_order = list(
            pd.unique(markers_df.get("group", pd.Series([], dtype=str)).astype(str))
        )

    # Determine sort column
    sort_col = (
        sort_markers_by
        if sort_markers_by in markers_df.columns
        else "scores"
        if "scores" in markers_df.columns
        else "logfoldchanges"
    )

    if sort_col not in markers_df.columns:
        log.warning(f"Sort column '{sort_col}' not found. Using first column.")
        sort_col = markers_df.columns[0] if not markers_df.empty else "names"

    # Helper: Detect term/pval columns
    def _detect_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
        term_candidates = [
            "Term",
            "term_name",
            "term",
            "Description",
            "Name",
            "Pathway",
        ]
        pval_candidates = ["pval_adj", "Adjusted P-value", "p.adjust", "padj", "FDR"]

        tcol = next((c for c in term_candidates if c in df.columns), None)
        pcol = next((c for c in pval_candidates if c in df.columns), None)

        return tcol, pcol

    # Build summaries
    summaries: Dict[str, str] = {}
    lines: List[str] = []

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
                    # GSEA: sort by NES (descending)
                    # ORA: sort by p-value (ascending)
                    sort_col_enr = pval_col
                    ascending = True

                    if (
                        enrichment_method_to_summarize == "gsea"
                        and "nes" in sig.columns
                    ):
                        sort_col_enr = "nes"
                        ascending = False

                    top_terms = (
                        sig.sort_values(sort_col_enr, ascending=ascending)[term_col]
                        .head(n_terms)
                        .astype(str)
                        .tolist()
                    )
                else:
                    log.debug(
                        f"Cluster {g}: No pathways below p_adj={enrichment_padj_cutoff}"
                    )
            else:
                log.warning(
                    f"Cluster {g}: Could not detect term/pval columns. "
                    f"Columns: {list(enr_df.columns)}"
                )

        # Format summary
        title = f"### Cluster {g}"
        mk_str = f"**Top Markers**: {', '.join(top_genes) if top_genes else 'N/A'}"
        pt_str = (
            f"**Top {enrichment_method_to_summarize.upper()} Pathways**: "
            f"{', '.join(top_terms) if top_terms else 'N/A'}"
        )

        summary_text = f"{title}\n{mk_str}\n{pt_str}"
        summaries[str(g)] = summary_text
        lines.append(summary_text)

    # Write to file
    if summary_file:
        Path(summary_file).parent.mkdir(parents=True, exist_ok=True)
        content = (
            "\n\n---\n\n".join(lines)
            if lines
            else "# Marker and Enrichment Summary\n\nNo results."
        )

        with open(summary_file, "w", encoding="utf-8") as f:
            f.write(content)

        log.info(f"Summaries exported to {summary_file}")

    return summaries


# ==================== Visualization Functions ====================


def visualize_markers(
    adata: AnnData,
    markers: Union[pd.DataFrame, Dict[str, List[str]], List[str]],
    groupby: Optional[str] = None,
    n_genes_per_group: int = 5,
    plot_type: Literal[
        "dotplot", "heatmap", "stacked_violin", "violin", "matrixplot"
    ] = "dotplot",
    dendrogram: bool = False,
    standard_scale: Optional[Literal["var", "group"]] = "var",
    swap_axes: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    save_path: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """
    Visualize marker genes with automatic formatting and error handling.

    Supports multiple input formats:
    - DataFrame with 'group' and 'names' columns
    - Dictionary mapping groups to gene lists
    - Simple list of genes

    Args:
        adata: AnnData object
        markers: Marker genes (DataFrame/dict/list)
        groupby: Grouping column (required for list input)
        n_genes_per_group: Top N genes per group (for DataFrame input)
        plot_type: Visualization type
        dendrogram: Add dendrogram
        standard_scale: Standardization ('var' or 'group')
        swap_axes: Swap x/y axes
        layer: Data layer to use
        use_raw: Use .raw
        save_path: Output file path
        figsize: Figure size (auto-calculated if None)
        **kwargs: Additional arguments to Scanpy plotting functions

    Example:
        >>> # From DataFrame
        >>> visualize_markers(
        ...     adata,
        ...     markers=filtered_markers,
        ...     groupby="leiden",
        ...     plot_type="dotplot"
        ... )
        >>>
        >>> # From dictionary
        >>> marker_dict = {
        ...     "T_cells": ["CD3D", "CD3E", "CD3G"],
        ...     "B_cells": ["CD19", "MS4A1", "CD79A"]
        ... }
        >>> visualize_markers(adata, markers=marker_dict)
    """
    gene_list: List[str] = []
    gene_dict: Dict[str, List[str]] = {}

    # Parse input markers
    if isinstance(markers, pd.DataFrame):
        # Standardize column names
        if "names" not in markers.columns:
            for alt in ("gene", "Gene", "feature", "symbol"):
                if alt in markers.columns:
                    markers = markers.rename(columns={alt: "names"})
                    break

        if "group" not in markers.columns:
            raise ValueError(
                "DataFrame must contain 'group' and 'names' columns "
                "for grouped visualization"
            )

        # Extract top genes per group
        for g in markers["group"].unique():
            group_markers = markers[markers["group"] == g]

            # Sort by logfoldchanges or scores
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

    elif isinstance(markers, (list, tuple)):
        gene_list = list(markers)
        if groupby is None:
            raise ValueError("groupby must be specified when markers is a list")
        gene_dict = {"Selected Markers": gene_list}

    else:
        raise TypeError("markers must be a DataFrame, dictionary, or list")

    # Deduplicate and validate
    gene_list_unique = [g for g in dict.fromkeys(gene_list) if g in adata.var_names]
    if not gene_list_unique:
        raise ValueError("No valid genes found in adata.var_names")

    # Prepare gene_dict for grouped plots
    for g, glist in gene_dict.items():
        gene_dict[g] = [gene for gene in glist if gene in adata.var_names]

    # Auto-calculate figsize
    if figsize is None:
        if groupby and groupby in adata.obs:
            n_groups = (
                len(adata.obs[groupby].cat.categories)
                if pd.api.types.is_categorical_dtype(adata.obs[groupby])
                else len(adata.obs[groupby].unique())
            )
        else:
            n_groups = 1

        n_genes = len(gene_list_unique)

        if plot_type in ["heatmap", "dotplot", "matrixplot"]:
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.3))
            if swap_axes:
                width, height = height, width
        elif plot_type == "stacked_violin":
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.4))
        else:  # violin
            width = max(8, n_genes * 2)
            height = 6

        figsize = (width, height)

    # Plot
    plot_kwargs = {
        "groupby": groupby,
        "dendrogram": dendrogram,
        "standard_scale": standard_scale,
        "use_raw": use_raw,
        "layer": layer,
        "figsize": figsize,
        "show": False,
        **kwargs,
    }

    try:
        if plot_type == "dotplot":
            sc.pl.dotplot(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "heatmap":
            sc.pl.heatmap(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "stacked_violin":
            sc.pl.stacked_violin(
                adata, var_names=gene_list_unique, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "matrixplot":
            sc.pl.matrixplot(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "violin":
            sc.pl.violin(adata, keys=gene_list_unique, **plot_kwargs)
        else:
            raise ValueError(f"Unknown plot type: {plot_type}")

    except Exception as e:
        log.error(f"Plotting failed: {e}")
        raise

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")
        plt.close()
    else:
        plt.show()


def plot_volcano(
    degs_df: pd.DataFrame,
    title: str,
    subtitle: Optional[str] = None,
    top_n_up: int = 15,
    top_n_down: int = 15,
    genes_to_highlight: Optional[List[str]] = None,
    lfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    palette: Optional[Dict[str, str]] = None,
    savepath: Optional[str] = None,
    figsize: tuple = (12, 12),
    dpi: int = 300,
) -> None:
    """
    Publication-quality volcano plot with intelligent label placement.

    Features:
    - Smart label selection based on ranking score (|LFC| * -log10(p))
    - adjustText for anti-overlap
    - Statistical summary box
    - Custom gene highlighting

    Args:
        degs_df: DE results DataFrame
        title: Main title
        subtitle: Subtitle (e.g., sample info)
        top_n_up: Number of top up-regulated genes to label
        top_n_down: Number of top down-regulated genes to label
        genes_to_highlight: Additional genes to highlight
        lfc_threshold: Log2 fold change threshold
        pval_threshold: Adjusted p-value threshold
        palette: Color scheme
        savepath: Output file path
        figsize: Figure size
        dpi: Resolution

    Example:
        >>> plot_volcano(
        ...     degs_df,
        ...     title="T cells: Treated vs Control",
        ...     subtitle="n=1234 cells",
        ...     top_n_up=20,
        ...     genes_to_highlight=["CD3D", "CD4"]
        ... )
    """
    df = degs_df.copy()

    # Calculate -log10(p-adj)
    df["-log10_pvals_adj"] = -np.log10(df["pvals_adj"].astype(float) + 1e-300)

    # Categorize genes
    df["status"] = "Not significant"
    df.loc[
        (df["logfoldchanges"] > lfc_threshold) & (df["pvals_adj"] < pval_threshold),
        "status",
    ] = "Up-regulated"
    df.loc[
        (df["logfoldchanges"] < -lfc_threshold) & (df["pvals_adj"] < pval_threshold),
        "status",
    ] = "Down-regulated"

    # Default palette
    if palette is None:
        palette = {
            "Up-regulated": "#d62728",
            "Down-regulated": "#1f77b4",
            "Not significant": "#cccccc",
        }

    # Create figure
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=figsize)

    # Plot points (layered for z-order)
    for status, color, alpha, size in [
        ("Not significant", palette["Not significant"], 0.4, 15),
        ("Up-regulated", palette["Up-regulated"], 0.8, 30),
        ("Down-regulated", palette["Down-regulated"], 0.8, 30),
    ]:
        mask = df["status"] == status
        ax.scatter(
            df.loc[mask, "logfoldchanges"],
            df.loc[mask, "-log10_pvals_adj"],
            s=size,
            alpha=alpha,
            c=color,
            label=status,
            ec="none",
            zorder=2 if status != "Not significant" else 1,
        )

    # Smart label selection
    df["ranking_score"] = np.abs(df["logfoldchanges"]) * df["-log10_pvals_adj"]

    up_genes = df[df["status"] == "Up-regulated"].nlargest(top_n_up, "ranking_score")
    down_genes = df[df["status"] == "Down-regulated"].nlargest(
        top_n_down, "ranking_score"
    )

    genes_to_label_df = pd.concat([up_genes, down_genes])

    # Add custom highlights
    if genes_to_highlight:
        highlight_df = df[df["names"].isin(genes_to_highlight)]
        genes_to_label_df = pd.concat(
            [genes_to_label_df, highlight_df]
        ).drop_duplicates(subset=["names"])

    # Add labels
    texts = []
    for _, row in genes_to_label_df.iterrows():
        txt = ax.text(
            row["logfoldchanges"],
            row["-log10_pvals_adj"],
            row["names"],
            fontsize=10,
            zorder=3,
        )
        texts.append(txt)

    # adjustText for anti-overlap
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.7),
            expand_points=(2.0, 2.0),
            expand_text=(1.3, 1.3),
            force_points=(0.3, 0.6),
            force_text=(0.5, 1.0),
            lim=1000,
            precision=0.01,
        )

    # Threshold lines
    ax.axhline(
        y=-np.log10(pval_threshold),
        color="grey",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax.axvline(x=lfc_threshold, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(x=-lfc_threshold, color="grey", linestyle="--", linewidth=1, alpha=0.7)

    # Statistical summary
    num_up = (df["status"] == "Up-regulated").sum()
    num_down = (df["status"] == "Down-regulated").sum()

    ax.text(
        0.02,
        0.98,
        f"Up: {num_up}\nDown: {num_down}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.8, ec="none"),
    )

    # Titles and labels
    fig.suptitle(title, fontsize=20, weight="bold", y=0.98)
    if subtitle:
        ax.set_title(subtitle, fontsize=14, pad=10)

    ax.set_xlabel("Log2 Fold Change", fontsize=14, weight="bold")
    ax.set_ylabel("-log10(Adjusted P-value)", fontsize=14, weight="bold")

    # Legend
    ax.legend(loc="upper right", frameon=True, fontsize=12, shadow=True)

    # Clean up
    sns.despine(ax=ax)
    plt.tight_layout()

    # Save or show
    if savepath:
        plt.savefig(savepath, dpi=dpi, bbox_inches="tight")
        log.info(f"Volcano plot saved: {savepath}")
        plt.close()
    else:
        plt.show()


def plot_multi_cluster_deg(
    df: pd.DataFrame,
    highlight_genes: Optional[List[str]] = None,
    pval_cutoff: float = 0.01,
    logfc_threshold: float = 1.0,
    top_n: int = 5,
    point_size_by_pval: bool = False,
    add_colored_bottom: bool = True,
    cluster_color_dict: Optional[Dict] = None,
    out_path: Optional[str] = None,
    figsize: Optional[tuple] = None,
    dpi: int = 300,
) -> None:
    """
    Multi-cluster differential expression overview plot.

    Shows:
    - All genes across all clusters in a single view
    - Significance-based coloring
    - Smart label placement for top genes
    - Optional colored bottom strip for cluster identification

    Args:
        df: DE DataFrame with columns: group, names, logfoldchanges, pvals_adj
        highlight_genes: Genes to highlight in green
        pval_cutoff: P-value cutoff for significance
        logfc_threshold: Log fold change threshold
        top_n: Top N genes to label per cluster (up and down)
        point_size_by_pval: Scale point size by -log10(p)
        add_colored_bottom: Add colored cluster strip at bottom
        cluster_color_dict: Custom cluster colors
        out_path: Output file path
        figsize: Figure size (auto-calculated if None)
        dpi: Resolution

    Example:
        >>> plot_multi_cluster_deg(
        ...     markers_df,
        ...     top_n=10,
        ...     highlight_genes=["CD3D", "CD19"],
        ...     out_path="cluster_overview.pdf"
        ... )
    """
    # Standardize column names
    if "names" in df.columns and "Gene" not in df.columns:
        df = df.rename(columns={"names": "Gene"})
    if "logfoldchanges" in df.columns and "avg_logFC" not in df.columns:
        df = df.rename(columns={"logfoldchanges": "avg_logFC"})
    if "group" in df.columns and "Cluster" not in df.columns:
        df = df.rename(columns={"group": "Cluster"})

    # Sort clusters
    try:
        clusters = sorted(df["Cluster"].unique(), key=int)
    except (ValueError, TypeError):
        clusters = sorted(df["Cluster"].unique())

    x_pos = np.arange(len(clusters))
    cluster_map = dict(zip(clusters, x_pos))

    # Colors
    if cluster_color_dict:
        color_map = cluster_color_dict
    else:
        cluster_colors = plt.cm.Spectral(np.linspace(0, 1, len(clusters)))
        color_map = dict(zip(clusters, cluster_colors))

    # Auto figsize
    if figsize is None:
        fig_width = max(16, len(clusters) * 1.8)
        fig_height = max(8, 8 + top_n * 0.2)
        figsize = (fig_width, fig_height)

    fig, ax = plt.subplots(figsize=figsize)

    texts = []
    points_coords = []

    for c in clusters:
        sub = df[df["Cluster"] == c].copy()
        idx = cluster_map[c]

        y = sub["avg_logFC"].values
        sig = sub["pvals_adj"].values < pval_cutoff

        up = (y > logfc_threshold) & sig
        down = (y < -logfc_threshold) & sig
        ns = ~sig

        # -log10(p) for sizing
        sub["neg_log_p"] = -np.log10(np.clip(sub["pvals_adj"], 1e-10, 1))

        # X jitter
        x = np.full(len(sub), idx)
        x_jitter = x + np.random.uniform(-0.45, 0.45, len(sub))

        # Point sizes
        base_size = 5
        if point_size_by_pval:
            sizes_ns = base_size * np.ones(sum(ns))
            sizes_up = base_size + 5 * sub.loc[up, "neg_log_p"]
            sizes_down = base_size + 5 * sub.loc[down, "neg_log_p"]
        else:
            sizes_ns = base_size
            sizes_up = base_size * 1.6
            sizes_down = base_size * 1.6

        # Plot points
        ax.scatter(x_jitter[ns], y[ns], c="#cccccc", s=sizes_ns, alpha=0.4, zorder=1)
        ax.scatter(x_jitter[up], y[up], c="#d62728", s=sizes_up, alpha=0.8, zorder=2)
        ax.scatter(
            x_jitter[down], y[down], c="#1f77b4", s=sizes_down, alpha=0.8, zorder=2
        )

        # Smart labeling
        sub["ranking_score"] = np.abs(sub["avg_logFC"]) * sub["neg_log_p"]

        # Top up
        top_up = (
            sub[up]
            .nlargest(top_n, "ranking_score")
            .sort_values("avg_logFC", ascending=False)
        )

        for j, (_, row) in enumerate(top_up.iterrows()):
            x_offset = [-0.25, 0, 0.25][j % 3]
            y_offset = 0.15

            txt = ax.text(
                idx + x_offset,
                row["avg_logFC"] + y_offset,
                row["Gene"],
                fontsize=6.5,
                ha="center",
                va="bottom",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
                zorder=3,
            )
            texts.append(txt)
            points_coords.append((idx, row["avg_logFC"]))

        # Top down
        top_down = (
            sub[down]
            .nlargest(top_n, "ranking_score")
            .sort_values("avg_logFC", ascending=True)
        )

        for j, (_, row) in enumerate(top_down.iterrows()):
            x_offset = [-0.25, 0, 0.25][j % 3]
            y_offset = -0.15

            txt = ax.text(
                idx + x_offset,
                row["avg_logFC"] + y_offset,
                row["Gene"],
                fontsize=6.5,
                ha="center",
                va="top",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
                zorder=3,
            )
            texts.append(txt)
            points_coords.append((idx, row["avg_logFC"]))

        # Highlight genes
        if highlight_genes:
            high_sub = sub[sub["Gene"].isin(highlight_genes)]
            for _, row in high_sub.iterrows():
                va = "bottom" if row["avg_logFC"] > 0 else "top"
                y_off = 0.2 if row["avg_logFC"] > 0 else -0.2

                txt = ax.text(
                    idx,
                    row["avg_logFC"] + y_off,
                    row["Gene"],
                    fontsize=7.5,
                    fontweight="bold",
                    color="green",
                    ha="center",
                    va=va,
                    bbox=dict(
                        boxstyle="round,pad=0.25",
                        facecolor="yellow",
                        edgecolor="green",
                        alpha=0.6,
                        linewidth=1.5,
                    ),
                    zorder=4,
                )
                texts.append(txt)
                points_coords.append((idx, row["avg_logFC"]))

    # adjustText
    if texts:
        adjust_text(
            texts,
            x=[p[0] for p in points_coords],
            y=[p[1] for p in points_coords],
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.4),
            expand_points=(2.5, 2.5),
            expand_text=(1.5, 1.8),
            expand_objects=(1.5, 1.5),
            force_points=(0.4, 0.8),
            force_text=(0.6, 1.2),
            force_objects=(0.4, 0.6),
            lim=2000,
            precision=0.001,
            only_move={"points": "xy", "text": "xy"},
            avoid_self=True,
            avoid_points=True,
            avoid_text=True,
            autoalign="xy",
        )

    # Threshold lines
    ax.axhline(logfc_threshold, ls="--", c="black", alpha=0.5, linewidth=1)
    ax.axhline(-logfc_threshold, ls="--", c="black", alpha=0.5, linewidth=1)
    ax.axhline(0, ls="--", c="gray", linewidth=0.8)

    # Colored bottom strip
    if add_colored_bottom:
        ylim = ax.get_ylim()
        dy = (ylim[1] - ylim[0]) * 0.035
        y_margin = (ylim[1] - ylim[0]) * 0.15
        ax.set_ylim(ylim[0] - dy, ylim[1] + y_margin)

        for i, c in enumerate(clusters):
            color = color_map.get(c, "gray")
            ax.add_patch(
                Rectangle(
                    (i - 0.5, ylim[0] - dy),
                    1,
                    dy,
                    color=color,
                    edgecolor="white",
                    linewidth=0.5,
                    clip_on=False,
                    zorder=0,
                )
            )

            # Auto text color
            if isinstance(color, str) and color.startswith("#"):
                rgb = [int(color.lstrip("#")[k : k + 2], 16) / 255 for k in (0, 2, 4)]
                text_color = "white" if np.mean(rgb) < 0.5 else "black"
            else:
                text_color = "black"

            ax.text(
                i,
                ylim[0] - dy / 2,
                str(c),
                ha="center",
                va="center",
                fontsize=9,
                color=text_color,
                weight="bold",
            )

        ax.set_xticks([])
        ax.set_xlabel("")
    else:
        ax.set_xticks(x_pos)
        ax.set_xticklabels(clusters, rotation=45, ha="right")
        ax.set_xlabel("Cluster")

    # Labels and title
    ax.set_ylabel("Average Log2 Fold Change", fontsize=12, weight="bold")
    ax.set_title(
        "Differential Expression per Cluster", fontsize=14, weight="bold", pad=20
    )

    # Grid
    ax.grid(True, ls="--", alpha=0.2, linewidth=0.5)

    # Spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    # Legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"Sig Up (P < {pval_cutoff})",
            markerfacecolor="#d62728",
            markersize=8,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"Sig Down (P < {pval_cutoff})",
            markerfacecolor="#1f77b4",
            markersize=8,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Non-Sig",
            markerfacecolor="#cccccc",
            markersize=8,
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=10,
    )

    plt.tight_layout()

    if out_path:
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
        log.info(f"Multi-cluster DEG plot saved: {out_path}")
        plt.close()
    else:
        plt.show()


# ==================== Result Management ====================


class ResultManager:
    """Unified result saving and loading manager."""

    SUPPORTED_FORMATS = {
        "csv": ".csv",
        "tsv": ".tsv",
        "excel": ".xlsx",
        "pickle": ".pkl",
        "parquet": ".parquet",
    }

    def __init__(self, base_dir: Union[str, Path]):
        """
        Initialize result manager.

        Args:
            base_dir: Base directory for storing results
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.de_dir = self.base_dir / "differential_expression"
        self.enrichment_dir = self.base_dir / "enrichment"
        self.plots_dir = self.base_dir / "plots"

        for dir_path in [self.de_dir, self.enrichment_dir, self.plots_dir]:
            dir_path.mkdir(exist_ok=True)

    def save_deg_results(
        self,
        results: pd.DataFrame,
        name: str,
        formats: Union[str, List[str]] = "csv",
        metadata: Optional[Dict] = None,
    ) -> Dict[str, Path]:
        """
        Save DE results in multiple formats.

        Args:
            results: Results DataFrame
            name: Base filename (without extension)
            formats: Format(s) to save ('csv', 'tsv', 'excel', 'pickle', 'parquet')
            metadata: Optional metadata dictionary

        Returns:
            Dictionary mapping formats to saved file paths
        """
        if isinstance(formats, str):
            formats = [formats]

        saved_paths = {}

        for fmt in formats:
            if fmt not in self.SUPPORTED_FORMATS:
                log.warning(f"Unsupported format '{fmt}', skipping")
                continue

            file_path = self.de_dir / f"{name}{self.SUPPORTED_FORMATS[fmt]}"

            try:
                if fmt == "csv":
                    results.to_csv(file_path, index=False)

                elif fmt == "tsv":
                    results.to_csv(file_path, index=False, sep="\t")

                elif fmt == "excel":
                    with pd.ExcelWriter(file_path, engine="openpyxl") as writer:
                        results.to_excel(writer, sheet_name="DEGs", index=False)
                        if metadata:
                            meta_df = pd.DataFrame([metadata]).T
                            meta_df.columns = ["Value"]
                            meta_df.to_excel(writer, sheet_name="Metadata")

                elif fmt == "pickle":
                    data_to_save = {"results": results, "metadata": metadata or {}}
                    with open(file_path, "wb") as f:
                        pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)

                elif fmt == "parquet":
                    results.to_parquet(file_path, index=False)

                saved_paths[fmt] = file_path
                log.info(f"Saved {fmt.upper()}: {file_path}")

            except Exception as e:
                log.error(f"Failed to save {fmt.upper()}: {e}")

        return saved_paths

    def load_deg_results(self, name: str, format: str = "csv") -> pd.DataFrame:
        """
        Load DE results.

        Args:
            name: Base filename (without extension)
            format: Format to load

        Returns:
            Results DataFrame
        """
        if format not in self.SUPPORTED_FORMATS:
            raise ValueError(f"Unsupported format '{format}'")

        file_path = self.de_dir / f"{name}{self.SUPPORTED_FORMATS[format]}"

        if not file_path.exists():
            raise FileNotFoundError(f"Result file not found: {file_path}")

        try:
            if format == "csv":
                return pd.read_csv(file_path)
            elif format == "tsv":
                return pd.read_csv(file_path, sep="\t")
            elif format == "excel":
                return pd.read_excel(file_path, sheet_name="DEGs")
            elif format == "pickle":
                with open(file_path, "rb") as f:
                    data = pickle.load(f)
                    return data["results"] if isinstance(data, dict) else data
            elif format == "parquet":
                return pd.read_parquet(file_path)

        except Exception as e:
            log.error(f"Failed to load {format.upper()}: {e}")
            raise


def save_results(
    results: Union[pd.DataFrame, Dict],
    name: str,
    outdir: str,
    formats: Union[str, List[str]] = "csv",
    **kwargs,
) -> Dict:
    """
    Convenience function for saving results.

    Args:
        results: Results to save
        name: Base filename
        outdir: Output directory
        formats: Format(s) to save
        **kwargs: Additional metadata

    Returns:
        Dictionary of saved file paths
    """
    manager = ResultManager(outdir)
    return manager.save_deg_results(results, name, formats, **kwargs)


def load_results(name: str, outdir: str, format: str = "csv") -> pd.DataFrame:
    """
    Convenience function for loading results.

    Args:
        name: Base filename
        outdir: Directory containing results
        format: Format to load

    Returns:
        Results DataFrame
    """
    manager = ResultManager(outdir)
    return manager.load_deg_results(name, format)
