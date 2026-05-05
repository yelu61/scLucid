"""
Core differential expression analysis functions.

This module provides the main DE analysis functions:
- find_markers: One-vs-rest marker gene discovery
- filter_markers: Filter DE results by criteria
- compare_groups: Pairwise group comparisons
- compare_conditions: Compare conditions within cell types
- get_conserved_markers: Find conserved markers across conditions
"""

import logging
from typing import Dict, List, Optional

import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...base_config import apply_config_overrides
from ...utils.helpers import sanitize_for_hdf5
from ..config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    ConservedMarkersConfig,
    DifferentialConfig,
    FilterMarkersConfig,
)
from .de_plots import plot_volcano
from .de_utils import _safe_filename
from .scanpy_compat import _to_frac, standardize_pct_columns as _standardize_pct_columns
from importlib.metadata import PackageNotFoundError, version

log = logging.getLogger(__name__)


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
    if config is None:
        active_config = DifferentialConfig()
        active_config = apply_config_overrides(active_config, **kwargs)
    else:
        active_config = apply_config_overrides(config, **kwargs)

    groupby = active_config.groupby
    key_added = active_config.key_added or "rank_genes_groups"

    if active_config.verbose:
        log.info(f"Finding markers: groupby='{groupby}', method='{active_config.method}'")

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
            "Scanpy 'names' field lacks structured dtype. " "Cannot extract group-wise results."
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
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})

    root[key_added] = adata.uns[key_added]  # Raw Scanpy output
    df_key = f"{key_added}_df"
    root[df_key] = full_df  # Processed DataFrame

    # Parameter tracking
    params = active_config.to_dict()
    params["scanpy_version"] = version("scanpy")
    root[f"{key_added}_params"] = sanitize_for_hdf5(params)

    if active_config.verbose:
        log.info(f"Found {len(full_df)} total markers across {len(groups_tested)} groups")
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

        log.info(f"[Filter] adj_p <= {config.max_padj}: kept {int(keep.sum())}/{len(filt)}")
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
            log.warning("'pct_nz_reference' not found; specificity-related filters skipped")

    filtered_df = df[filt].copy()
    log.info(f"Retained {len(filtered_df)} genes after all filters")

    # Post-filter: Keep top N per group
    if config.keep_top_n is not None and config.keep_top_n > 0 and not filtered_df.empty:
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
                else "scores" if "scores" in filtered_df.columns else filtered_df.columns[0]
            )
            log.warning(
                f"Sort key '{config.sort_by}' not found. " f"Falling back to '{fallback_col}'"
            )
            sort_by_col = fallback_col

        log.info(
            f"Selecting top {config.keep_top_n} genes per group, " f"sorted by '{sort_by_col}'"
        )

        parts = []
        for g in filtered_df["group"].unique():
            sub = filtered_df[filtered_df["group"] == g].sort_values(sort_by_col, ascending=False)
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
        config = config.model_copy(update=kwargs)

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
        raise ValueError(f"No cells found for either '{group1}' or '{group2}' in '{groupby}'")

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
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
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
        config = config.model_copy(update=kwargs)

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
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
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
        config = config.model_copy(update=kwargs)

    key_added = config.key_added or (f"conserved_markers_{config.groupby}_{config.condition_key}")

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
                (adata.obs[config.groupby] == group) & (adata.obs[config.condition_key] == cond)
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
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
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
