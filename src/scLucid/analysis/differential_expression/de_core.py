"""
Core differential expression analysis functions.

This module provides the main DE analysis functions:
- find_markers: One-vs-rest cell-level marker gene discovery
- filter_markers: Filter DE results by criteria
- compare_groups: Exploratory pairwise cell-level group comparisons
- compare_conditions: Exploratory cell-level condition comparisons within cell types
- get_conserved_markers: Find conserved markers across conditions
"""

import logging
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse, stats

from ...base_config import apply_config_overrides
from ...utils.helpers import sanitize_for_hdf5
from ..config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    ConservedMarkersConfig,
    DifferentialConfig,
    FilterMarkersConfig,
    PseudobulkDEConfig,
)
from .de_plots import plot_volcano
from .de_utils import _hierarchical_fdr_correction, _safe_filename
from .scanpy_compat import _to_frac
from .scanpy_compat import standardize_pct_columns as _standardize_pct_columns

log = logging.getLogger(__name__)


# ==================== Core Differential Expression Functions ====================


def _tag_cell_level_de_result(
    df: pd.DataFrame,
    *,
    inference_level: str,
    analysis_intent: str,
) -> pd.DataFrame:
    """Add explicit inference semantics to cell-level DE result tables."""
    if df.empty:
        return df
    tagged = df.copy()
    tagged["inference_level"] = inference_level
    tagged["claim_level"] = (
        "exploratory_marker_screen"
        if inference_level == "cell_level_marker_discovery"
        else "exploratory_hypothesis_generation"
    )
    tagged["analysis_intent"] = analysis_intent
    tagged["valid_for_publication_inference"] = False
    tagged["pseudoreplication_warning"] = True
    tagged["recommended_formal_inference_api"] = "run_pseudobulk_de"
    tagged["de_warning"] = (
        "Cell-level rank_genes_groups treats cells as independent observations. "
        "Use sample-level pseudobulk DE for publication-grade condition inference."
    )
    return tagged


def _tag_pseudobulk_de_result(df: pd.DataFrame) -> pd.DataFrame:
    """Add explicit claim semantics to pseudobulk DE result tables."""
    if df.empty:
        return df
    tagged = df.copy()
    if "inference_level" not in tagged:
        tagged["inference_level"] = "sample_level"
    if "valid_for_publication_inference" not in tagged:
        tagged["valid_for_publication_inference"] = False

    def _claim(row: pd.Series) -> str:
        inference = str(row.get("inference_level", ""))
        valid = bool(row.get("valid_for_publication_inference", False))
        if inference == "sample_level" and valid:
            return "replicate_aware_sample_level_condition_inference"
        if inference == "descriptive_single_sample":
            return "descriptive_effect_size_only"
        if inference == "exploratory_cell_level":
            return "exploratory_hypothesis_generation"
        return "review_required_condition_inference"

    tagged["claim_level"] = tagged.apply(_claim, axis=1)
    tagged["recommended_formal_inference_api"] = "run_pseudobulk_de"
    tagged["de_review_note"] = tagged["claim_level"].map(
        {
            "replicate_aware_sample_level_condition_inference": (
                "Sample-level pseudobulk result with biological replicates; still review design, covariates, and model assumptions."
            ),
            "descriptive_effect_size_only": (
                "Single-sample or unreplicated pseudobulk result; use effect sizes descriptively, not as formal DE."
            ),
            "exploratory_hypothesis_generation": (
                "Cell-level fallback treats cells as observations; use only for screening/hypothesis generation."
            ),
            "review_required_condition_inference": (
                "Condition inference status is ambiguous; review sample replication and method metadata."
            ),
        }
    )
    return tagged


def find_markers(
    adata: AnnData,
    config: Optional[DifferentialConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Find marker genes using one-vs-rest cell-level differential expression.

    This function is intended for marker discovery and cluster/cell-type
    characterization. It is not a publication-grade condition-level inference
    test because cells from the same biological sample are not independent
    replicates. Use ``run_pseudobulk_de`` for condition DE when sample
    replicates are available.

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
            "Scanpy 'names' field lacks structured dtype. Cannot extract group-wise results."
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
        full_df = _tag_cell_level_de_result(
            full_df,
            inference_level="cell_level_marker_discovery",
            analysis_intent="marker_discovery",
        )

    # Store with provenance
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})

    root[key_added] = adata.uns[key_added]  # Raw Scanpy output
    df_key = f"{key_added}_df"
    root[df_key] = full_df  # Processed DataFrame

    # Parameter tracking
    params = active_config.to_dict()
    params["scanpy_version"] = version("scanpy")
    params["inference_level"] = "cell_level_marker_discovery"
    params["claim_level"] = "exploratory_marker_screen"
    params["valid_for_publication_inference"] = False
    params["de_warning"] = (
        "find_markers is for marker discovery; use run_pseudobulk_de for formal condition DE."
    )
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
    if not filtered_df.empty and "claim_level" not in filtered_df.columns:
        filtered_df = _tag_cell_level_de_result(
            filtered_df,
            inference_level="cell_level_marker_discovery",
            analysis_intent="filtered_marker_discovery",
        )
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
                    "Cannot sort by 'diff_pct' without 'pct_nz_reference'. Falling back to 'scores'"
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
            log.warning(f"Sort key '{config.sort_by}' not found. Falling back to '{fallback_col}'")
            sort_by_col = fallback_col

        log.info(f"Selecting top {config.keep_top_n} genes per group, sorted by '{sort_by_col}'")

        parts = []
        for g in filtered_df["group"].unique():
            sub = filtered_df[filtered_df["group"] == g].sort_values(sort_by_col, ascending=False)
            parts.append(sub.head(config.keep_top_n))

        filtered_df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    # Store results
    root[key_added] = filtered_df
    root[f"{key_added}_params"] = {
        **config.to_dict(),
        "n_retained": len(filtered_df),
        "claim_level": "exploratory_marker_screen",
        "valid_for_publication_inference": False,
    }

    log.info(
        f"Final filtered markers: {len(filtered_df)} rows -> "
        f".uns['sclucid']['analysis']['de']['{key_added}']"
    )

    return filtered_df


def compare_groups(
    adata: AnnData, config: Optional[CompareGroupsConfig] = None, **kwargs
) -> pd.DataFrame:
    """
    Compare two specific groups (e.g., cell types, conditions) for DE genes at cell level.

    Combines the robustness of careful validation with the convenience of
    integrated visualization. This is an exploratory cell-level comparison,
    suitable for marker discovery or cluster characterization. For
    condition-level biological inference, use ``run_pseudobulk_de`` so
    biological samples, not cells, define replicates.

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
        ...     plot=True,
        ...     acknowledge_exploratory=True,
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

    if not config.acknowledge_exploratory:
        raise ValueError(
            "compare_groups performs exploratory cell-level DE that treats cells as independent "
            "observations and is not valid for formal condition inference. "
            "Set acknowledge_exploratory=True in CompareGroupsConfig to proceed after reviewing "
            "this warning, or use run_pseudobulk_de for sample-level inference."
        )

    if config.verbose:
        log.info(f"Comparing DE genes: '{group1}' vs '{group2}' in '{groupby}'")
    log.warning(
        "compare_groups runs exploratory cell-level DE. For condition-level "
        "inference, use run_pseudobulk_de with sample_col and biological replicates."
    )

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
    final = _tag_cell_level_de_result(
        final,
        inference_level="exploratory_cell_level",
        analysis_intent="exploratory_pairwise_cell_level_de",
    )

    if config.verbose:
        log.info(f"Found {len(final)} DE genes ({len(up)} up, {len(down)} down)")

    # Store results
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    root[key_added] = final
    params = config.to_dict()
    params["inference_level"] = "exploratory_cell_level"
    params["claim_level"] = "exploratory_hypothesis_generation"
    params["valid_for_publication_inference"] = False
    params["de_warning"] = (
        "compare_groups is cell-level exploratory DE; use run_pseudobulk_de for formal condition DE."
    )
    params["recommended_formal_inference_api"] = "run_pseudobulk_de"
    root[f"{key_added}_params"] = sanitize_for_hdf5(params)

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
    Compare two conditions within a specific cell type/group at cell level.

    This is a specialized wrapper around compare_groups() that first
    subsets to a specific cell type, then compares conditions. It is
    exploratory and not formal condition inference because cells are not
    biological replicates. Prefer ``run_pseudobulk_de`` for publication-grade
    condition DE.

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
        ...     max_padj=0.05,
        ...     acknowledge_exploratory=True,
        ... )
        >>> degs = compare_conditions(adata, config)
    """
    if config is None:
        config = CompareConditionsConfig(**kwargs)
    else:
        config = config.model_copy(update=kwargs)

    if not config.acknowledge_exploratory:
        raise ValueError(
            "compare_conditions performs exploratory cell-level DE that treats cells as independent "
            "observations and is not valid for formal condition inference. "
            "Set acknowledge_exploratory=True in CompareConditionsConfig to proceed after reviewing "
            "this warning, or use run_pseudobulk_de for sample-level inference."
        )

    groupby = config.groupby
    group_name = config.group_name
    condition_key = config.condition_key

    log.info(
        f"Comparing conditions '{config.condition1}' vs '{config.condition2}' within '{group_name}'"
    )
    log.warning(
        "compare_conditions runs exploratory cell-level DE. For condition-level "
        "inference, use run_pseudobulk_de with sample_col and biological replicates."
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
        or (f"compare_{config.condition1}_vs_{config.condition2}_in_{_safe_filename(group_name)}"),
        plot=config.plot,
        save_dir=config.save_dir,
        verbose=config.verbose,
        acknowledge_exploratory=True,
    )

    # Run comparison on subset
    results_df = compare_groups(adata_subset, config=comp_config)

    # Add metadata
    results_df["celltype"] = group_name

    # Store in parent AnnData
    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    root[comp_config.key_added] = results_df
    params = config.to_dict()
    params["inference_level"] = "exploratory_cell_level"
    params["claim_level"] = "exploratory_hypothesis_generation"
    params["valid_for_publication_inference"] = False
    params["de_warning"] = (
        "compare_conditions is cell-level exploratory DE; use run_pseudobulk_de for formal condition DE."
    )
    params["recommended_formal_inference_api"] = "run_pseudobulk_de"
    root[f"{comp_config.key_added}_params"] = sanitize_for_hdf5(params)

    log.info(
        f"Condition comparison complete: {len(results_df)} DE genes. "
        f"Stored at .uns['...']['{comp_config.key_added}']"
    )

    return results_df


def _get_expression_matrix(adata: AnnData, layer: Optional[str], use_raw: bool):
    """Return expression matrix and gene names for pseudobulk aggregation."""
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        return adata.layers[layer], pd.Index(adata.var_names)
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is not available")
        return adata.raw.X, pd.Index(adata.raw.var_names)
    return adata.X, pd.Index(adata.var_names)


def _aggregate_counts_by_sample(
    adata: AnnData,
    sample_col: str,
    condition_key: str,
    layer: Optional[str],
    use_raw: bool,
    min_cells_per_sample: int,
    covariate_cols: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate cell-level counts to sample-level pseudobulk counts."""
    if sample_col not in adata.obs:
        raise KeyError(f"Column '{sample_col}' not found in adata.obs")
    if condition_key not in adata.obs:
        raise KeyError(f"Column '{condition_key}' not found in adata.obs")
    covariate_cols = list(dict.fromkeys(covariate_cols or []))
    missing_covariates = [col for col in covariate_cols if col not in adata.obs]
    if missing_covariates:
        raise KeyError(f"Covariate column(s) not found in adata.obs: {missing_covariates}")

    X, var_names = _get_expression_matrix(adata, layer=layer, use_raw=use_raw)
    sample_values = adata.obs[sample_col].astype(str)
    samples = list(pd.unique(sample_values))
    rows = []
    meta_rows = []

    for sample in samples:
        mask = (sample_values == sample).to_numpy()
        n_cells = int(mask.sum())
        if n_cells < min_cells_per_sample:
            continue
        block = X[mask]
        summed = (
            np.asarray(block.sum(axis=0)).ravel() if sparse.issparse(block) else block.sum(axis=0)
        )
        rows.append(summed)
        cond_values = adata.obs.loc[mask, condition_key].astype(str).unique()
        if len(cond_values) != 1:
            raise ValueError(
                f"Sample '{sample}' has multiple conditions in '{condition_key}': {cond_values}"
            )
        meta_row = {
            "sample": sample,
            condition_key: cond_values[0],
            "n_cells": n_cells,
            "library_size": float(np.sum(summed)),
        }
        for covariate in covariate_cols:
            cov_values = adata.obs.loc[mask, covariate].astype(str).unique()
            if len(cov_values) != 1:
                raise ValueError(
                    f"Sample '{sample}' has multiple values in covariate '{covariate}': {cov_values}"
                )
            meta_row[covariate] = cov_values[0]
        meta_rows.append(meta_row)

    if not rows:
        return pd.DataFrame(columns=var_names), pd.DataFrame()

    counts = pd.DataFrame(rows, index=[row["sample"] for row in meta_rows], columns=var_names)
    meta = pd.DataFrame(meta_rows).set_index("sample")
    return counts, meta


def _benjamini_hochberg(pvals: pd.Series, method: str = "fdr_bh") -> pd.Series:
    """Adjust p-values with statsmodels when available, with a local BH fallback."""
    pvals = pd.to_numeric(pvals, errors="coerce")
    valid = pvals.notna()
    adjusted = pd.Series(np.nan, index=pvals.index, dtype=float)
    if valid.sum() == 0:
        return adjusted

    try:
        from statsmodels.stats.multitest import multipletests

        _, padj, _, _ = multipletests(pvals[valid].to_numpy(), method=method)
    except Exception:
        order = np.argsort(pvals[valid].to_numpy())
        ranked = pvals[valid].to_numpy()[order]
        n = len(ranked)
        padj_ordered = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
        padj = np.empty_like(padj_ordered)
        padj[order] = np.clip(padj_ordered, 0, 1)
    adjusted.loc[valid] = padj
    return adjusted


def _run_welch_logcpm_de(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    condition_key: str,
    condition1: str,
    condition2: str,
    min_counts: int,
    pseudocount: float,
    p_adjust_method: str,
) -> pd.DataFrame:
    """Run a conservative Python pseudobulk fallback on logCPM values."""
    selected = meta[meta[condition_key].astype(str).isin([condition1, condition2])].copy()
    if selected.empty:
        return pd.DataFrame()

    selected_counts = counts.loc[selected.index]
    keep = selected_counts.sum(axis=0) >= min_counts
    selected_counts = selected_counts.loc[:, keep]
    if selected_counts.empty:
        return pd.DataFrame()

    lib_sizes = selected_counts.sum(axis=1)
    zero_lib = lib_sizes == 0
    if zero_lib.any():
        log.warning(
            f"{zero_lib.sum()} sample(s) have zero library size after filtering; "
            "excluding them from CPM calculation."
        )
        selected_counts = selected_counts.loc[~zero_lib]
        lib_sizes = lib_sizes.loc[~zero_lib]
        selected = selected.loc[selected_counts.index]
        if selected.empty:
            return pd.DataFrame()
    cpm = selected_counts.div(lib_sizes, axis=0) * 1e6
    logcpm = np.log2(cpm + pseudocount)

    group1 = logcpm.loc[selected[condition_key].astype(str) == condition1]
    group2 = logcpm.loc[selected[condition_key].astype(str) == condition2]
    n1, n2 = len(group1), len(group2)
    if n1 == 0 or n2 == 0:
        return pd.DataFrame()

    records = []
    for gene in logcpm.columns:
        x1 = group1[gene].astype(float)
        x2 = group2[gene].astype(float)
        mean1 = float(x1.mean())
        mean2 = float(x2.mean())
        log2fc = mean2 - mean1
        if n1 >= 2 and n2 >= 2:
            stat, pval = stats.ttest_ind(x2, x1, equal_var=False, nan_policy="omit")
            method = "welch_logcpm_n2" if min(n1, n2) == 2 else "welch_logcpm_n3plus"
        else:
            stat, pval = np.nan, np.nan
            method = "insufficient_replicates"
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": log2fc,
                "log2fc": log2fc,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "mean_logcpm_condition1": mean1,
                "mean_logcpm_condition2": mean2,
                "base_mean": float(selected_counts[gene].mean()),
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": n1,
                "n_samples_condition2": n2,
                "method": method,
            }
        )

    result = pd.DataFrame(records)
    result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
    result["padj"] = result["pvals_adj"]
    return result.sort_values(["pvals_adj", "pvals"], na_position="last").reset_index(drop=True)


def _run_linear_model_logcpm_de(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    condition_key: str,
    condition1: str,
    condition2: str,
    min_counts: int,
    pseudocount: float,
    p_adjust_method: str,
    covariates: Optional[List[str]] = None,
    robust_cov_type: str = "HC3",
) -> pd.DataFrame:
    """Run sample-level OLS on logCPM values with optional covariates."""
    try:
        import statsmodels.formula.api as smf
    except ImportError as exc:
        raise ImportError("statsmodels is required for linear_model_logcpm") from exc

    covariates = list(dict.fromkeys(covariates or []))
    selected = meta[meta[condition_key].astype(str).isin([condition1, condition2])].copy()
    if selected.empty:
        return pd.DataFrame()
    missing_covariates = [cov for cov in covariates if cov not in selected.columns]
    if missing_covariates:
        raise KeyError(
            f"Covariate column(s) missing from pseudobulk metadata: {missing_covariates}"
        )

    selected_counts = counts.loc[selected.index]
    keep = selected_counts.sum(axis=0) >= min_counts
    selected_counts = selected_counts.loc[:, keep]
    if selected_counts.empty:
        return pd.DataFrame()

    lib_sizes = selected_counts.sum(axis=1)
    nonzero_lib = lib_sizes > 0
    if not nonzero_lib.all():
        selected_counts = selected_counts.loc[nonzero_lib]
        lib_sizes = lib_sizes.loc[nonzero_lib]
        selected = selected.loc[selected_counts.index]
        if selected.empty:
            return pd.DataFrame()

    selected["__condition"] = pd.Categorical(
        selected[condition_key].astype(str),
        categories=[condition1, condition2],
        ordered=False,
    )
    design_terms = ["C(__condition)"]
    model_covariates: List[str] = []
    for idx, covariate in enumerate(covariates):
        if selected[covariate].nunique(dropna=True) < 2:
            log.info(
                "Skipping covariate '%s' in pseudobulk linear model because it has <2 levels.",
                covariate,
            )
            continue
        safe_col = f"__cov_{idx}"
        selected[safe_col] = selected[covariate].astype(str)
        design_terms.append(f"C({safe_col})")
        model_covariates.append(covariate)

    cpm = selected_counts.div(lib_sizes, axis=0) * 1e6
    logcpm = np.log2(cpm + pseudocount)
    n1 = int((selected["__condition"].astype(str) == condition1).sum())
    n2 = int((selected["__condition"].astype(str) == condition2).sum())
    if n1 < 2 or n2 < 2:
        return pd.DataFrame()

    formula = "value ~ " + " + ".join(design_terms)
    term = f"C(__condition)[T.{condition2}]"
    records = []
    for gene in logcpm.columns:
        model_df = selected.copy()
        model_df["value"] = logcpm[gene].astype(float)
        try:
            fit = smf.ols(formula, data=model_df).fit()
            result_fit = (
                fit
                if robust_cov_type == "nonrobust"
                else fit.get_robustcov_results(cov_type=robust_cov_type)
            )
        except Exception as exc:
            log.warning("Linear model failed for gene '%s': %s", gene, exc)
            continue
        exog_names = list(getattr(result_fit.model, "exog_names", []))
        if term not in exog_names:
            continue
        term_idx = exog_names.index(term)
        coef = float(np.asarray(result_fit.params)[term_idx])
        pval = float(np.asarray(result_fit.pvalues)[term_idx])
        stat = float(np.asarray(result_fit.tvalues)[term_idx])
        ci = np.asarray(result_fit.conf_int())[term_idx]
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": coef,
                "log2fc": coef,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "ci_lower": float(ci[0]),
                "ci_upper": float(ci[1]),
                "mean_logcpm_condition1": float(
                    logcpm.loc[selected["__condition"].astype(str) == condition1, gene].mean()
                ),
                "mean_logcpm_condition2": float(
                    logcpm.loc[selected["__condition"].astype(str) == condition2, gene].mean()
                ),
                "base_mean": float(selected_counts[gene].mean()),
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": n1,
                "n_samples_condition2": n2,
                "method": "linear_model_logcpm",
                "covariance_type": robust_cov_type,
                "design_formula": formula,
                "design_covariates": ",".join(model_covariates),
            }
        )

    result = pd.DataFrame(records)
    if result.empty:
        return result
    result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
    result["padj"] = result["pvals_adj"]
    result["model_warning"] = (
        "Python sample-level logCPM linear model with "
        f"{robust_cov_type} covariance; use as replicate-aware inference "
        "with explicit covariates, not as a replacement for full edgeR/limma/DESeq2 validation."
    )
    return result.sort_values(["pvals_adj", "pvals"], na_position="last").reset_index(drop=True)


def _make_pseudobulk_adata(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    condition_key: str,
    covariates: Optional[List[str]] = None,
    min_counts: int = 0,
) -> AnnData:
    """Convert sample-level counts and metadata into an AnnData for statsmodels backends."""
    covariates = list(dict.fromkeys(covariates or []))
    obs = meta.copy()
    if sample_col not in obs.columns:
        obs[sample_col] = obs.index.astype(str)
    if condition_key not in obs.columns:
        raise KeyError(f"Condition key '{condition_key}' not found in pseudobulk metadata")
    for cov in covariates:
        if cov not in obs.columns:
            raise KeyError(f"Covariate '{cov}' not found in pseudobulk metadata")
    filtered_counts = counts
    if min_counts > 0:
        keep = counts.sum(axis=0) >= min_counts
        filtered_counts = counts.loc[:, keep]
    return AnnData(
        X=np.asarray(filtered_counts.values, dtype=float),
        obs=obs,
        var=pd.DataFrame(index=filtered_counts.columns),
    )


def _run_statsmodels_glm_de(
    pseudobulk_adata: AnnData,
    design_col: str,
    sample_col: str,
    covariates: Optional[List[str]] = None,
    family: str = "NegativeBinomial",
    condition1: Optional[str] = None,
    condition2: Optional[str] = None,
    p_adjust_method: str = "fdr_bh",
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """
    Run per-gene GLM on pseudobulk counts using statsmodels.

    Supports NegativeBinomial (raw counts, log link), Gamma (counts + pseudocount,
    log link), and Gaussian (logCPM, identity link with HC3 sandwich covariance).
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    covariates = list(dict.fromkeys(covariates or []))
    df = pseudobulk_adata.obs.copy()
    if sample_col not in df.columns:
        df[sample_col] = df.index.astype(str)

    cond_values = pd.unique(df[design_col].astype(str))
    if condition1 is None and condition2 is None:
        if len(cond_values) != 2:
            raise ValueError(
                "Exactly two conditions required when condition1/condition2 are not specified"
            )
        condition1, condition2 = sorted(cond_values)
    elif condition1 is None or condition2 is None:
        raise ValueError("Provide both condition1 and condition2 or neither")

    selected_mask = df[design_col].astype(str).isin([condition1, condition2])
    df = df[selected_mask].copy()
    X = pseudobulk_adata.X[selected_mask.to_numpy()]
    if hasattr(X, "toarray"):
        X = X.toarray()

    df["__condition"] = pd.Categorical(
        df[design_col].astype(str),
        categories=[condition1, condition2],
        ordered=False,
    )
    terms = ["C(__condition)"]
    valid_covariates: List[str] = []
    for idx, cov in enumerate(covariates):
        if cov in df.columns and df[cov].nunique(dropna=True) >= 2:
            safe = f"__cov_{idx}"
            df[safe] = df[cov].astype(str)
            terms.append(f"C({safe})")
            valid_covariates.append(cov)

    formula = "value ~ " + " + ".join(terms)
    term = f"C(__condition)[T.{condition2}]"

    family_norm = str(family).strip().lower()
    if family_norm == "negativebinomial":
        sm_family = sm.families.NegativeBinomial(alpha=1.0)
        use_logcpm = False
        log_link = True
    elif family_norm == "gamma":
        sm_family = sm.families.Gamma(link=sm.families.links.Log())
        use_logcpm = False
        log_link = True
    elif family_norm == "gaussian":
        sm_family = sm.families.Gaussian()
        use_logcpm = True
        log_link = False
    else:
        raise ValueError(f"Unsupported GLM family: {family}")

    counts_df = pd.DataFrame(X, index=df.index, columns=pseudobulk_adata.var_names)
    if use_logcpm:
        lib_sizes = counts_df.sum(axis=1)
        nonzero = lib_sizes > 0
        if not nonzero.all():
            counts_df = counts_df.loc[nonzero]
            df = df.loc[nonzero]
        cpm = counts_df.div(lib_sizes.loc[nonzero], axis=0) * 1e6
        response = np.log2(cpm + pseudocount)
    else:
        response = counts_df + pseudocount

    n1 = int((df["__condition"].astype(str) == condition1).sum())
    n2 = int((df["__condition"].astype(str) == condition2).sum())

    records = []
    for gene in pseudobulk_adata.var_names:
        if gene not in response.columns:
            continue
        model_df = df.copy()
        model_df["value"] = response[gene].astype(float)
        if model_df["value"].notna().sum() < 3:
            continue
        try:
            if family_norm == "gaussian":
                fit = smf.glm(formula, data=model_df, family=sm_family).fit(cov_type="HC3")
            else:
                fit = smf.glm(formula, data=model_df, family=sm_family).fit()
        except Exception as exc:
            log.warning("Statsmodels GLM failed for gene '%s': %s", gene, exc)
            continue
        exog_names = list(getattr(fit.model, "exog_names", []))
        if term not in exog_names:
            continue
        term_idx = exog_names.index(term)
        coef = float(np.asarray(fit.params)[term_idx])
        pval = float(np.asarray(fit.pvalues)[term_idx])
        stat = float(np.asarray(fit.tvalues)[term_idx])
        if log_link:
            log2fc = float(coef / np.log(2))
        else:
            log2fc = float(coef)
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logFC": log2fc,
                "log2fc": log2fc,
                "logfoldchanges": log2fc,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": n1,
                "n_samples_condition2": n2,
                "method": f"statsmodels_glm_{family_norm}",
                "glm_family": family_norm,
                "design_formula": formula,
                "design_covariates": ",".join(valid_covariates),
            }
        )

    result = pd.DataFrame(records)
    if result.empty:
        return result
    result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
    result["padj"] = result["pvals_adj"]
    return result.sort_values(["pvals_adj", "pvals"], na_position="last").reset_index(drop=True)


def _run_statsmodels_gee_de(
    pseudobulk_adata: AnnData,
    design_col: str,
    sample_col: str,
    covariates: Optional[List[str]] = None,
    family: str = "NegativeBinomial",
    condition1: Optional[str] = None,
    condition2: Optional[str] = None,
    p_adjust_method: str = "fdr_bh",
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """
    Run per-gene GEE on pseudobulk counts using statsmodels.

    Uses an exchangeable working correlation within ``sample_col`` to account
    for within-sample correlation. This is a Python-native sample-aware
    alternative to mixed models for pseudobulk data.
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.genmod.cov_struct import Exchangeable, Independence

    covariates = list(dict.fromkeys(covariates or []))
    df = pseudobulk_adata.obs.copy()
    if sample_col not in df.columns:
        df[sample_col] = df.index.astype(str)

    cond_values = pd.unique(df[design_col].astype(str))
    if condition1 is None and condition2 is None:
        if len(cond_values) != 2:
            raise ValueError(
                "Exactly two conditions required when condition1/condition2 are not specified"
            )
        condition1, condition2 = sorted(cond_values)
    elif condition1 is None or condition2 is None:
        raise ValueError("Provide both condition1 and condition2 or neither")

    selected_mask = df[design_col].astype(str).isin([condition1, condition2])
    df = df[selected_mask].copy()
    X = pseudobulk_adata.X[selected_mask.to_numpy()]
    if hasattr(X, "toarray"):
        X = X.toarray()

    df["__condition"] = pd.Categorical(
        df[design_col].astype(str),
        categories=[condition1, condition2],
        ordered=False,
    )
    terms = ["C(__condition)"]
    valid_covariates: List[str] = []
    for idx, cov in enumerate(covariates):
        if cov in df.columns and df[cov].nunique(dropna=True) >= 2:
            safe = f"__cov_{idx}"
            df[safe] = df[cov].astype(str)
            terms.append(f"C({safe})")
            valid_covariates.append(cov)

    formula = "value ~ " + " + ".join(terms)
    term = f"C(__condition)[T.{condition2}]"

    family_norm = str(family).strip().lower()
    if family_norm == "negativebinomial":
        sm_family = sm.families.NegativeBinomial(alpha=1.0)
        log_link = True
    elif family_norm == "gamma":
        sm_family = sm.families.Gamma(link=sm.families.links.Log())
        log_link = True
    elif family_norm == "gaussian":
        sm_family = sm.families.Gaussian()
        log_link = False
    else:
        raise ValueError(f"Unsupported GEE family: {family}")

    counts_df = pd.DataFrame(X, index=df.index, columns=pseudobulk_adata.var_names)
    response = counts_df + pseudocount

    n1 = int((df["__condition"].astype(str) == condition1).sum())
    n2 = int((df["__condition"].astype(str) == condition2).sum())

    records = []
    for gene in pseudobulk_adata.var_names:
        if gene not in response.columns:
            continue
        model_df = df.copy()
        model_df["value"] = response[gene].astype(float)
        if model_df["value"].notna().sum() < 3:
            continue
        # Try exchangeable first; fall back to independence if it fails
        # (e.g. one observation per sample).
        fit = None
        for cov_struct in [Exchangeable(), Independence()]:
            try:
                fit = smf.gee(
                    formula,
                    data=model_df,
                    groups=model_df[sample_col],
                    family=sm_family,
                    cov_struct=cov_struct,
                ).fit()
                break
            except Exception as exc:
                log.debug("GEE %s failed for gene '%s': %s", type(cov_struct).__name__, gene, exc)
                continue
        if fit is None:
            log.warning("Statsmodels GEE failed for gene '%s'", gene)
            continue
        exog_names = list(getattr(fit.model, "exog_names", []))
        if term not in exog_names:
            continue
        term_idx = exog_names.index(term)
        coef = float(np.asarray(fit.params)[term_idx])
        pval = float(np.asarray(fit.pvalues)[term_idx])
        stat = float(np.asarray(fit.tvalues)[term_idx])
        if log_link:
            log2fc = float(coef / np.log(2))
        else:
            log2fc = float(coef)
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logFC": log2fc,
                "log2fc": log2fc,
                "logfoldchanges": log2fc,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": n1,
                "n_samples_condition2": n2,
                "method": f"statsmodels_gee_{family_norm}",
                "gee_family": family_norm,
                "design_formula": formula,
                "design_covariates": ",".join(valid_covariates),
            }
        )

    result = pd.DataFrame(records)
    if result.empty:
        return result
    result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
    result["padj"] = result["pvals_adj"]
    return result.sort_values(["pvals_adj", "pvals"], na_position="last").reset_index(drop=True)


def _run_descriptive_single_sample_pseudobulk(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    condition_key: str,
    condition1: str,
    condition2: str,
    min_counts: int,
    pseudocount: float,
) -> pd.DataFrame:
    """Return effect-size-only pseudobulk summaries when formal replication is absent."""
    selected = meta[meta[condition_key].astype(str).isin([condition1, condition2])].copy()
    if selected.empty:
        return pd.DataFrame()

    selected_counts = counts.loc[selected.index]
    keep = selected_counts.sum(axis=0) >= min_counts
    selected_counts = selected_counts.loc[:, keep]
    if selected_counts.empty:
        return pd.DataFrame()

    lib_sizes = selected_counts.sum(axis=1)
    nonzero_lib = lib_sizes > 0
    selected_counts = selected_counts.loc[nonzero_lib]
    lib_sizes = lib_sizes.loc[nonzero_lib]
    selected = selected.loc[selected_counts.index]
    if selected.empty:
        return pd.DataFrame()

    cpm = selected_counts.div(lib_sizes, axis=0) * 1e6
    logcpm = np.log2(cpm + pseudocount)
    group1_mask = selected[condition_key].astype(str) == condition1
    group2_mask = selected[condition_key].astype(str) == condition2
    group1 = logcpm.loc[group1_mask]
    group2 = logcpm.loc[group2_mask]
    counts1 = selected_counts.loc[group1_mask]
    counts2 = selected_counts.loc[group2_mask]
    n1, n2 = len(group1), len(group2)
    if n1 == 0 or n2 == 0:
        return pd.DataFrame()

    records = []
    for gene in logcpm.columns:
        mean1 = float(group1[gene].mean())
        mean2 = float(group2[gene].mean())
        log2fc = mean2 - mean1
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": log2fc,
                "log2fc": log2fc,
                "scores": np.nan,
                "statistic": np.nan,
                "pvals": np.nan,
                "pval": np.nan,
                "pvals_adj": np.nan,
                "padj": np.nan,
                "mean_logcpm_condition1": mean1,
                "mean_logcpm_condition2": mean2,
                "mean_counts_condition1": float(counts1[gene].mean()),
                "mean_counts_condition2": float(counts2[gene].mean()),
                "base_mean": float(selected_counts[gene].mean()),
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": n1,
                "n_samples_condition2": n2,
                "method": "descriptive_pseudobulk",
            }
        )

    result = pd.DataFrame(records)
    return result.sort_values("logfoldchanges", key=lambda s: s.abs(), ascending=False).reset_index(
        drop=True
    )


def _run_pydeseq2_de(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    condition_key: str,
    condition1: str,
    condition2: str,
    min_counts: int,
    p_adjust_method: str,
) -> pd.DataFrame:
    """Run DESeq2 on pseudobulk counts via pydeseq2."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ImportError("pydeseq2 is not installed") from exc

    selected = meta[meta[condition_key].astype(str).isin([condition1, condition2])].copy()
    if selected.empty:
        return pd.DataFrame()

    selected_counts = counts.loc[selected.index].round().astype(int)
    keep = selected_counts.sum(axis=0) >= min_counts
    selected_counts = selected_counts.loc[:, keep]
    if selected_counts.empty:
        return pd.DataFrame()

    design_col = "__condition"
    selected = pd.DataFrame(index=selected.index)
    selected[design_col] = pd.Categorical(
        meta.loc[selected.index, condition_key].astype(str),
        categories=[condition1, condition2],
        ordered=False,
    )

    dds = DeseqDataSet(
        counts=selected_counts,
        metadata=selected,
        design=f"~{design_col}",
        quiet=True,
        n_cpus=1,
    )
    dds.deseq2()
    deseq_stats = DeseqStats(
        dds,
        contrast=[design_col, condition2, condition1],
        quiet=True,
        n_cpus=1,
    )
    deseq_stats.summary()
    result = deseq_stats.results_df.copy()

    rename_map = {
        "log2FoldChange": "logfoldchanges",
        "pvalue": "pvals",
        "padj": "pvals_adj",
        "stat": "scores",
        "baseMean": "base_mean",
    }
    result = result.rename(columns=rename_map)
    result["names"] = result.index.astype(str)
    result["gene"] = result["names"]
    if "logfoldchanges" in result:
        result["log2fc"] = result["logfoldchanges"]
    if "scores" in result:
        result["statistic"] = result["scores"]
    if "pvals" in result:
        result["pval"] = result["pvals"]
    if "pvals_adj" not in result and "pvals" in result:
        result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
    if "pvals_adj" in result:
        result["padj"] = result["pvals_adj"]

    n1 = int((selected[design_col].astype(str) == condition1).sum())
    n2 = int((selected[design_col].astype(str) == condition2).sum())
    result["condition1"] = condition1
    result["condition2"] = condition2
    result["contrast"] = f"{condition2}_vs_{condition1}"
    result["direction"] = f"{condition2} - {condition1}"
    result["n_samples_condition1"] = n1
    result["n_samples_condition2"] = n2
    result["method"] = "deseq2"

    sort_cols = [col for col in ["pvals_adj", "pvals"] if col in result.columns]
    if sort_cols:
        result = result.sort_values(sort_cols, na_position="last")
    return result.reset_index(drop=True)


def run_pseudobulk_de(
    adata: AnnData,
    config: Optional[PseudobulkDEConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Run sample-level pseudobulk DE for one or more condition contrasts.

    This is the preferred path for condition-level DEG when biological sample
    replicates are available. If a contrast has only one sample in either
    condition and ``fallback_to_cell_level=True``, it falls back to the legacy
    cell-level ``compare_conditions``/``compare_groups`` path and marks the
    result as exploratory.
    """
    if config is None:
        active_config = PseudobulkDEConfig(**kwargs)
    else:
        active_config = config.model_copy(update=kwargs)

    groupby = active_config.groupby
    if groupby is not None and groupby not in adata.obs:
        raise KeyError(f"Column '{groupby}' not found in adata.obs")

    if groupby is None:
        group_names: List[Optional[str]] = [None]
    elif active_config.group_names is not None:
        group_names = list(active_config.group_names)
    else:
        group_names = list(pd.unique(adata.obs[groupby].astype(str)))

    tasks: List[Tuple[Optional[str], str, str]] = [
        (group_name, condition1, condition2)
        for group_name in group_names
        for condition1, condition2 in active_config.contrasts
    ]
    design_covariates = list(dict.fromkeys(active_config.design_covariates or []))
    if active_config.block_col:
        design_covariates = list(dict.fromkeys([*design_covariates, active_config.block_col]))

    def _run_one(task: Tuple[Optional[str], str, str]) -> pd.DataFrame:
        group_name, condition1, condition2 = task
        if group_name is None:
            adata_sub = adata
            group_label = "all"
        else:
            adata_sub = adata[adata.obs[groupby].astype(str) == str(group_name)].copy()
            group_label = str(group_name)

        counts, meta = _aggregate_counts_by_sample(
            adata_sub,
            sample_col=active_config.sample_col,
            condition_key=active_config.condition_key,
            layer=active_config.layer,
            use_raw=active_config.use_raw,
            min_cells_per_sample=active_config.min_cells_per_sample,
            covariate_cols=design_covariates,
        )
        if meta.empty:
            return pd.DataFrame()

        selected = meta[
            meta[active_config.condition_key].astype(str).isin([condition1, condition2])
        ]
        n1 = int((selected[active_config.condition_key].astype(str) == condition1).sum())
        n2 = int((selected[active_config.condition_key].astype(str) == condition2).sum())

        if min(n1, n2) < active_config.min_samples_per_condition:
            log.warning(
                f"Skipping {group_label} {condition2} vs {condition1}: "
                f"replicates are {n1}/{n2}, below min_samples_per_condition"
            )
            return pd.DataFrame()

        has_formal_replicates = min(n1, n2) >= 2
        should_use_cell_fallback = active_config.method == "cell_level_fallback" or (
            active_config.method == "auto"
            and min(n1, n2) < 2
            and active_config.fallback_to_cell_level
        )
        should_use_descriptive_single_sample = (
            active_config.method == "auto"
            and min(n1, n2) < 2
            and not active_config.fallback_to_cell_level
            and active_config.single_sample_mode == "descriptive"
        )
        should_use_linear_model = (
            active_config.method == "linear_model_logcpm"
            or (active_config.method == "auto" and bool(design_covariates))
        ) and min(n1, n2) >= 2
        should_use_statsmodels_glm = active_config.method == "statsmodels_glm" and min(n1, n2) >= 2
        should_use_statsmodels_gee = active_config.method == "statsmodels_gee" and min(n1, n2) >= 2
        should_try_deseq2 = (
            active_config.method in {"auto", "deseq2"}
            and not should_use_linear_model
            and not should_use_statsmodels_glm
            and not should_use_statsmodels_gee
            and min(n1, n2) >= 2
        )
        should_use_welch = (
            active_config.method in {"auto", "welch_logcpm"}
            and not should_use_linear_model
            and not should_use_statsmodels_glm
            and not should_use_statsmodels_gee
            and min(n1, n2) >= 2
        )

        if should_use_descriptive_single_sample:
            result = _run_descriptive_single_sample_pseudobulk(
                counts,
                meta,
                condition_key=active_config.condition_key,
                condition1=condition1,
                condition2=condition2,
                min_counts=active_config.min_counts,
                pseudocount=active_config.pseudocount,
            )
            if not result.empty:
                result["pseudobulk_warning"] = (
                    "Only one biological sample in at least one condition; "
                    "returned descriptive pseudobulk effect sizes without formal p-values."
                )
        elif should_use_cell_fallback:
            if (
                active_config.method != "cell_level_fallback"
                and not active_config.fallback_to_cell_level
            ):
                log.warning(
                    f"Skipping {group_label} {condition2} vs {condition1}: "
                    "cell-level fallback is disabled"
                )
                return pd.DataFrame()

            if groupby is None:
                fallback_config = CompareGroupsConfig(
                    groupby=active_config.condition_key,
                    group1=condition2,
                    group2=condition1,
                    min_log2fc=0.0,
                    max_padj=1.0,
                    n_top_genes=active_config.n_genes
                    if hasattr(active_config, "n_genes")
                    else 5000,
                    layer=active_config.layer,
                    use_raw=active_config.use_raw,
                    plot=False,
                    acknowledge_exploratory=True,
                )
                result = compare_groups(adata_sub, fallback_config)
            else:
                fallback_config = CompareConditionsConfig(
                    groupby=groupby,
                    group_name=str(group_name),
                    condition_key=active_config.condition_key,
                    condition1=condition2,
                    condition2=condition1,
                    min_log2fc=0.0,
                    max_padj=1.0,
                    n_top_genes=5000,
                    layer=active_config.layer,
                    use_raw=active_config.use_raw,
                    plot=False,
                    acknowledge_exploratory=True,
                )
                result = compare_conditions(adata, fallback_config)
            result = result.copy()
            result["method"] = "cell_level_fallback"
            if min(n1, n2) < 2:
                result["pseudobulk_warning"] = "Only one sample in at least one condition"
            else:
                result["pseudobulk_warning"] = (
                    "Cell-level fallback forced by method='cell_level_fallback'"
                )
        elif should_use_linear_model:
            result = _run_linear_model_logcpm_de(
                counts,
                meta,
                condition_key=active_config.condition_key,
                condition1=condition1,
                condition2=condition2,
                min_counts=active_config.min_counts,
                pseudocount=active_config.pseudocount,
                p_adjust_method=active_config.p_adjust_method,
                covariates=design_covariates,
                robust_cov_type=active_config.robust_cov_type,
            )
        elif should_use_statsmodels_glm:
            pb_adata = _make_pseudobulk_adata(
                counts,
                meta,
                sample_col=active_config.sample_col,
                condition_key=active_config.condition_key,
                covariates=design_covariates,
                min_counts=active_config.min_counts,
            )
            result = _run_statsmodels_glm_de(
                pb_adata,
                design_col=active_config.condition_key,
                sample_col=active_config.sample_col,
                covariates=design_covariates,
                family=active_config.statsmodels_glm_family,
                condition1=condition1,
                condition2=condition2,
                p_adjust_method=active_config.p_adjust_method,
                pseudocount=active_config.pseudocount,
            )
        elif should_use_statsmodels_gee:
            pb_adata = _make_pseudobulk_adata(
                counts,
                meta,
                sample_col=active_config.sample_col,
                condition_key=active_config.condition_key,
                covariates=design_covariates,
                min_counts=active_config.min_counts,
            )
            result = _run_statsmodels_gee_de(
                pb_adata,
                design_col=active_config.condition_key,
                sample_col=active_config.sample_col,
                covariates=design_covariates,
                family=active_config.statsmodels_glm_family,
                condition1=condition1,
                condition2=condition2,
                p_adjust_method=active_config.p_adjust_method,
                pseudocount=active_config.pseudocount,
            )
        elif should_try_deseq2:
            try:
                result = _run_pydeseq2_de(
                    counts,
                    meta,
                    condition_key=active_config.condition_key,
                    condition1=condition1,
                    condition2=condition2,
                    min_counts=active_config.min_counts,
                    p_adjust_method=active_config.p_adjust_method,
                )
            except Exception as exc:
                if active_config.method == "deseq2":
                    log.warning(
                        f"Skipping {group_label} {condition2} vs {condition1}: "
                        f"DESeq2 failed ({exc})"
                    )
                    return pd.DataFrame()
                log.warning(
                    f"DESeq2 failed for {group_label} {condition2} vs {condition1}; "
                    f"falling back to Welch logCPM ({exc})"
                )
                result = _run_welch_logcpm_de(
                    counts,
                    meta,
                    condition_key=active_config.condition_key,
                    condition1=condition1,
                    condition2=condition2,
                    min_counts=active_config.min_counts,
                    pseudocount=active_config.pseudocount,
                    p_adjust_method=active_config.p_adjust_method,
                )
                if not result.empty:
                    result["pseudobulk_warning"] = "DESeq2 failed; used Welch logCPM fallback"
        elif should_use_welch:
            result = _run_welch_logcpm_de(
                counts,
                meta,
                condition_key=active_config.condition_key,
                condition1=condition1,
                condition2=condition2,
                min_counts=active_config.min_counts,
                pseudocount=active_config.pseudocount,
                p_adjust_method=active_config.p_adjust_method,
            )
        else:
            log.warning(
                f"Skipping {group_label} {condition2} vs {condition1}: "
                f"method='{active_config.method}' requires >=2 samples per condition"
            )
            return pd.DataFrame()

        if result.empty:
            return result
        result["group"] = group_label
        result["groupby"] = groupby or "all"
        result["condition1"] = condition1
        result["condition2"] = condition2
        result["contrast"] = f"{condition2}_vs_{condition1}"
        result["direction"] = f"{condition2} - {condition1}"
        result["design_covariates"] = ",".join(design_covariates)
        result["block_col"] = active_config.block_col or ""
        if "n_samples_condition1" not in result:
            result["n_samples_condition1"] = n1
        if "n_samples_condition2" not in result:
            result["n_samples_condition2"] = n2
        result["replicate_status"] = (
            "replicated" if has_formal_replicates else "single_sample_or_unreplicated"
        )
        result["inference_level"] = "sample_level"
        result["valid_for_publication_inference"] = bool(has_formal_replicates)
        result["pseudoreplication_warning"] = False
        if result["method"].astype(str).eq("descriptive_pseudobulk").all():
            result["inference_level"] = "descriptive_single_sample"
            result["valid_for_publication_inference"] = False
        elif result["method"].astype(str).eq("cell_level_fallback").all():
            result["inference_level"] = "exploratory_cell_level"
            result["valid_for_publication_inference"] = False
            result["pseudoreplication_warning"] = True
        if not has_formal_replicates and "pseudobulk_warning" not in result:
            result["pseudobulk_warning"] = "Insufficient biological replicates for formal inference"
        return _tag_pseudobulk_de_result(result)

    if active_config.n_jobs > 1 and len(tasks) > 1:
        with ThreadPoolExecutor(max_workers=active_config.n_jobs) as pool:
            result_parts = list(pool.map(_run_one, tasks))
    else:
        result_parts = [_run_one(task) for task in tasks]

    result_parts = [part for part in result_parts if part is not None and not part.empty]
    results = pd.concat(result_parts, ignore_index=True) if result_parts else pd.DataFrame()

    if active_config.hierarchical_correction and not results.empty and "pvals" in results.columns:
        group_cols = []
        if "group" in results.columns:
            group_cols.append("group")
        if "contrast" in results.columns:
            group_cols.append("contrast")
        results = _hierarchical_fdr_correction(
            results,
            pval_col="pvals",
            group_cols=group_cols or ["contrast"],
            global_method="benjamini_yekutieli",
        )

    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    root[active_config.key_added] = results
    params = active_config.to_dict()
    params["recommended_primary_method"] = "sample_level_pseudobulk"
    params["cell_level_fallback_policy"] = "exploratory_only"
    params["claim_level"] = (
        "replicate_aware_sample_level_condition_inference"
        if not results.empty and results["valid_for_publication_inference"].all()
        else "review_required_condition_inference"
    )
    root[f"{active_config.key_added}_params"] = sanitize_for_hdf5(params)
    return results


def run_mixedlm_de(
    adata: AnnData,
    sample_col: str,
    condition_col: str,
    condition1: Optional[str] = None,
    condition2: Optional[str] = None,
    covariates: Optional[List[str]] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    p_adjust_method: str = "fdr_bh",
    max_genes: Optional[int] = 5000,
    key_added: str = "mixedlm_de",
) -> pd.DataFrame:
    """
    Cell-level sample-aware differential expression using statsmodels MixedLM.

    For each gene, fits a linear mixed model with a random intercept per sample
    and tests the condition coefficient. This accounts for within-sample
    correlation but is still cell-level and therefore tagged as not valid for
    publication inference on its own.

    Args:
        adata: AnnData object.
        sample_col: Column in ``adata.obs`` containing sample identifiers.
        condition_col: Column in ``adata.obs`` containing condition labels.
        condition1: Reference condition. If None, inferred from the data.
        condition2: Test condition. If None, inferred from the data.
        covariates: Additional fixed-effect covariates.
        layer: Layer to use for expression values.
        use_raw: Use ``adata.raw``.
        p_adjust_method: Multiple-testing method.
        max_genes: Maximum number of genes to test. Warns if exceeded.
        key_added: Key for storing results in ``adata.uns``.

    Returns:
        DataFrame with logFC, pvalue, padj and sample-aware tags.
    """
    from statsmodels.regression.mixed_linear_model import MixedLM

    if sample_col not in adata.obs.columns:
        raise KeyError(f"Column '{sample_col}' not found in adata.obs")
    if condition_col not in adata.obs.columns:
        raise KeyError(f"Column '{condition_col}' not found in adata.obs")

    covariates = list(dict.fromkeys(covariates or []))
    for cov in covariates:
        if cov not in adata.obs.columns:
            raise KeyError(f"Covariate '{cov}' not found in adata.obs")

    X, var_names = _get_expression_matrix(adata, layer=layer, use_raw=use_raw)
    if sparse.issparse(X):
        X = np.asarray(X.toarray())
    else:
        X = np.asarray(X)

    # Log1p-transform counts to put coefficients on approximately log scale.
    if X.max() > 50 and layer is None and not use_raw:
        endog_matrix = np.log1p(X)
    else:
        endog_matrix = X

    n_genes = len(var_names)
    if max_genes is not None and n_genes > max_genes:
        log.warning(
            "MixedLM DE requested for %d genes, exceeding max_genes=%d. "
            "Testing only the first %d genes; increase max_genes to test more.",
            n_genes,
            max_genes,
            max_genes,
        )
        endog_matrix = endog_matrix[:, :max_genes]
        var_names = var_names[:max_genes]
        n_genes = max_genes

    obs = adata.obs.copy()
    cond_values = pd.unique(obs[condition_col].astype(str))
    if condition1 is None and condition2 is None:
        if len(cond_values) != 2:
            raise ValueError(
                "Exactly two conditions required when condition1/condition2 are not specified"
            )
        condition1, condition2 = sorted(cond_values)
    elif condition1 is None or condition2 is None:
        raise ValueError("Provide both condition1 and condition2 or neither")

    selected_mask = obs[condition_col].astype(str).isin([condition1, condition2]).to_numpy()
    if selected_mask.sum() < 3:
        raise ValueError("Need at least 3 cells for MixedLM DE")

    obs = obs.loc[selected_mask].copy()
    obs["__condition"] = pd.Categorical(
        obs[condition_col].astype(str),
        categories=[condition1, condition2],
        ordered=False,
    )
    groups = obs[sample_col].astype(str)
    n_samples = groups.nunique()
    if n_samples < 2:
        raise ValueError("Need at least 2 samples for MixedLM DE")

    exog = pd.DataFrame(
        {
            "intercept": 1.0,
            "condition": (obs["__condition"].astype(str) == condition2).astype(float),
        }
    )
    for idx, cov in enumerate(covariates):
        if obs[cov].nunique(dropna=True) >= 2:
            safe = f"__cov_{idx}"
            obs[safe] = obs[cov].astype(str)
            dummies = pd.get_dummies(obs[safe], prefix=safe, drop_first=True)
            exog = pd.concat([exog, dummies.astype(float)], axis=1)

    term = "condition"
    records = []
    for gene_idx, gene in enumerate(var_names):
        endog = endog_matrix[selected_mask, gene_idx].astype(float)
        try:
            model = MixedLM(endog, exog, groups=groups)
            fit = model.fit(reml=False)
        except Exception as exc:
            log.warning("MixedLM failed for gene '%s': %s", gene, exc)
            continue
        exog_names = list(fit.model.exog_names)
        if term not in exog_names:
            continue
        term_idx = exog_names.index(term)
        coef = float(np.asarray(fit.params)[term_idx])
        pval = float(np.asarray(fit.pvalues)[term_idx])
        stat = float(np.asarray(fit.tvalues)[term_idx])
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logFC": coef,
                "log2fc": coef,
                "logfoldchanges": coef,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "condition1": condition1,
                "condition2": condition2,
                "contrast": f"{condition2}_vs_{condition1}",
                "direction": f"{condition2} - {condition1}",
                "n_samples": n_samples,
                "n_cells": int(selected_mask.sum()),
                "method": "statsmodels_mixedlm",
                "model": "statsmodels_mixedlm",
                "inference_level": "cell_level_sample_aware",
                "valid_for_publication_inference": False,
            }
        )

    result = pd.DataFrame(records)
    if not result.empty:
        result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=p_adjust_method)
        result["padj"] = result["pvals_adj"]
        result = result.sort_values(["pvals_adj", "pvals"], na_position="last").reset_index(
            drop=True
        )

    root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    root[key_added] = result
    root[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "sample_col": sample_col,
            "condition_col": condition_col,
            "condition1": condition1,
            "condition2": condition2,
            "covariates": covariates,
            "layer": layer,
            "use_raw": use_raw,
            "p_adjust_method": p_adjust_method,
            "max_genes": max_genes,
            "inference_level": "cell_level_sample_aware",
            "valid_for_publication_inference": False,
            "model": "statsmodels_mixedlm",
        }
    )
    return result


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
