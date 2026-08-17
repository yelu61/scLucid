"""
Cell type proportion analysis workflow (Pseudo-bulk method).

This module orchestrates Pseudo-bulk proportion analysis by combining
statistical testing and visualization functions from submodules.

Main workflow:
- Compute proportions
- Run statistical tests
- Generate visualizations
- Export results

For detailed statistical and plotting functions, see:
- stats.py: Statistical tests and effect sizes
- plots.py: Visualization functions
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from .config import ProportionConfig
from .plots import (
    plot_batch_effect,
    plot_box_summary,
    plot_cell_counts,
    plot_celltype_correlation,
    plot_celltype_variability,
    plot_composition,
    plot_composition_pca,
    plot_composition_transform_heatmap,
    plot_diff_stats,
    plot_effect_size_volcano,
    plot_individual_boxplots,
    plot_paired_proportion_shifts,
    plot_proportion_bar,
    plot_proportion_heatmap,
    plot_proportion_shifts,
    plot_proportion_timeseries,
    plot_proportion_with_ci,
)
from .stats import (
    compute_celltype_proportion,
    export_analysis_data,
    run_statistical_test,
)

log = logging.getLogger(__name__)


# ================= Helper Functions =================


def _set_publication_style():
    """Set publication-quality plotting style."""
    plt.style.use("seaborn-v0_8-paper")

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.linewidth": 1.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _natural_sort_key(text):
    """
    Generate key for natural sorting (handles numbers in strings).

    Parameters
    ----------
    text : str
        Text string to sort

    Returns:
    -------
    list
        Key for sorting
    """
    import re

    return [int(c) if c.isdigit() else c.lower() for c in re.split("([0-9]+)", text)]


# ================= Main Workflow =================


_REQUIRES_PAIRED_METHODS = {
    "clr-paired-t-test",
    "clr-paired-wilcoxon",
    "paired-t-test",
    "paired-wilcoxon",
}
_PAIRING_AWARE_METHODS = {
    *_REQUIRES_PAIRED_METHODS,
    "ancom-like-clr",
}


def _validate_proportion_metadata(
    adata: AnnData,
    config: ProportionConfig,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate sample-level metadata without silently collapsing conflicts."""
    experimental_unit_col = (
        config.experimental_unit_col or config.pairing_col or config.sample_col
    )
    metadata_cols = list(
        dict.fromkeys(
            col
            for col in (
                config.sample_col,
                config.condition_col,
                experimental_unit_col,
                config.pairing_col,
                config.batch_col,
                config.timepoint_col,
            )
            if col
        )
    )
    required_cols = list(
        dict.fromkeys([config.sample_col, config.condition_col, config.celltype_col, *metadata_cols])
    )
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Proportion metadata column(s) not found in adata.obs: {missing}")

    complete_cols = list(dict.fromkeys([*metadata_cols, config.celltype_col]))
    for col in complete_cols:
        values = adata.obs[col]
        missing_mask = values.isna() | values.astype(str).str.strip().eq("")
        if bool(missing_mask.any()):
            raise ValueError(
                f"Proportion analysis requires complete non-empty metadata in '{col}'; "
                f"found {int(missing_mask.sum())} invalid cell(s)."
            )

    cell_meta = adata.obs[metadata_cols].copy()
    conflicts: dict[str, int] = {}
    for col in metadata_cols:
        if col == config.sample_col:
            continue
        by_sample = cell_meta.groupby(config.sample_col, observed=True)[col].nunique(dropna=False)
        n_conflicts = int((by_sample > 1).sum())
        if n_conflicts:
            conflicts[col] = n_conflicts
    if conflicts:
        details = ", ".join(f"{col}={count}" for col, count in conflicts.items())
        raise ValueError(
            f"Each '{config.sample_col}' must map to exactly one condition/unit/covariate; "
            f"conflicting sample mappings: {details}."
        )

    sample_meta = cell_meta.drop_duplicates(subset=[config.sample_col]).set_index(
        config.sample_col
    )
    condition_values = sample_meta[config.condition_col].astype(str)
    conditions = list(pd.unique(condition_values))
    if len(conditions) < 2:
        raise ValueError(
            "Proportion comparison requires at least two conditions; "
            f"found {conditions}. Use descriptive proportions without a statistical test."
        )

    duplicate_unit_condition = (
        sample_meta.reset_index()
        .groupby([experimental_unit_col, config.condition_col], observed=True)[config.sample_col]
        .nunique()
    )
    duplicate_unit_condition = duplicate_unit_condition[duplicate_unit_condition > 1]
    if not duplicate_unit_condition.empty:
        raise ValueError(
            "Multiple sample rows map to the same experimental-unit/condition combination. "
            "Consolidate technical replicates before proportion inference."
        )

    unit_condition_counts = (
        sample_meta.reset_index()
        .groupby(experimental_unit_col, observed=True)[config.condition_col]
        .nunique()
    )
    repeated_units = int((unit_condition_counts > 1).sum())
    if repeated_units and not config.pairing_col:
        raise ValueError(
            f"Experimental units in '{experimental_unit_col}' occur in multiple conditions, "
            "but pairing_col is not configured. A paired design must be explicit."
        )
    if (
        repeated_units
        and config.pairing_col
        and config.pairing_col != experimental_unit_col
    ):
        unit_pair_rows = sample_meta.reset_index()[
            [experimental_unit_col, config.pairing_col]
        ].drop_duplicates()
        pairs_per_unit = unit_pair_rows.groupby(experimental_unit_col, observed=True)[
            config.pairing_col
        ].nunique(dropna=False)
        if bool((pairs_per_unit > 1).any()):
            raise ValueError(
                f"pairing_col='{config.pairing_col}' is not stable within each "
                "experimental unit."
            )
        repeated_unit_values = unit_condition_counts[unit_condition_counts > 1].index
        repeated_unit_pairs = unit_pair_rows[
            unit_pair_rows[experimental_unit_col].isin(repeated_unit_values)
        ]
        if repeated_unit_pairs[config.pairing_col].nunique(dropna=False) != repeated_units:
            raise ValueError(
                f"pairing_col='{config.pairing_col}' must uniquely identify each repeated "
                "experimental unit; a shared/coarser value does not model pairing."
            )

    unit_counts = {
        str(condition): int(frame[experimental_unit_col].nunique())
        for condition, frame in sample_meta.reset_index().groupby(
            config.condition_col, observed=True
        )
    }
    complete_pairs = 0
    if config.pairing_col:
        pair_condition_counts = (
            sample_meta.reset_index()
            .groupby(config.pairing_col, observed=True)[config.condition_col]
            .nunique()
        )
        complete_pairs = int((pair_condition_counts == len(conditions)).sum())

    design = {
        "status": "READY",
        "sample_col": config.sample_col,
        "condition_col": config.condition_col,
        "experimental_unit_col": experimental_unit_col,
        "pairing_col": config.pairing_col,
        "batch_col": config.batch_col,
        "n_samples": int(sample_meta.shape[0]),
        "conditions": [str(value) for value in conditions],
        "samples_per_condition": {
            str(key): int(value) for key, value in condition_values.value_counts().items()
        },
        "experimental_units_per_condition": unit_counts,
        "n_complete_pairs": complete_pairs,
        "replicate_basis": "complete_pairs" if config.pairing_col else experimental_unit_col,
        "mapping_conflicts": conflicts,
    }
    return sample_meta, design


def _validate_proportion_method_design(
    config: ProportionConfig,
    sample_meta: pd.DataFrame,
    design: dict[str, Any],
) -> None:
    """Fail closed when the selected test cannot identify the configured design."""
    conditions = design["conditions"]
    if config.test_method in _REQUIRES_PAIRED_METHODS:
        if not config.pairing_col:
            raise ValueError(f"{config.test_method} requires pairing_col")
        if len(conditions) != 2:
            raise ValueError(f"{config.test_method} requires exactly two conditions")
    elif config.test_method == "ancom-like-clr" and len(conditions) != 2:
        raise ValueError("ancom-like-clr requires exactly two conditions")
    elif (
        config.pairing_col
        and design["n_complete_pairs"]
        and config.test_method not in _PAIRING_AWARE_METHODS
    ):
        if config.legacy_exploratory:
            design["status"] = "REVIEW"
            design.setdefault("warnings", []).append(
                f"Repeated observations in '{config.pairing_col}' are ignored by "
                f"exploratory test_method='{config.test_method}'."
            )
        else:
            raise ValueError(
                f"Repeated/paired observations were detected in '{config.pairing_col}', but "
                f"test_method='{config.test_method}' does not model pairing."
            )

    if config.batch_col:
        n_batches = int(sample_meta[config.batch_col].nunique(dropna=True))
        if n_batches > 1 and config.test_method != "clr-ols":
            raise ValueError(
                f"batch_col='{config.batch_col}' has multiple levels but "
                f"test_method='{config.test_method}' would ignore it. Use clr-ols or remove "
                "the covariate after documenting why adjustment is unnecessary."
            )
        if config.test_method == "clr-ols" and n_batches < 2:
            raise ValueError("clr-ols requires batch_col with at least two observed levels")
        if config.test_method == "clr-ols":
            condition = pd.get_dummies(
                sample_meta[config.condition_col].astype(str), drop_first=True, dtype=float
            )
            batch = pd.get_dummies(
                sample_meta[config.batch_col].astype(str), drop_first=True, dtype=float
            )
            matrix = np.column_stack(
                [np.ones(sample_meta.shape[0]), condition.to_numpy(), batch.to_numpy()]
            )
            if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
                raise ValueError(
                    "The proportion design matrix is rank deficient: condition and batch "
                    "are not separately identifiable."
                )


def _auto_configure_analysis(adata: AnnData, config: ProportionConfig) -> ProportionConfig:
    """
    Automatically configure test method and plot types based on data characteristics.

    Logic:
    - N=1 per group: Force chi-square; disable boxplots
    - N=2 per group: Prefer DESeq2/t-test; enable basic plots
    - N>=3 per group: Prefer Wilcoxon/DESeq2; enable boxplots, volcano
    - Paired data: Prefer paired tests
    - Multi-group (>2): Prefer ANOVA

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : ProportionConfig
        Configuration object

    Returns:
    -------
    ProportionConfig
        Auto-configured settings
    """
    from copy import deepcopy

    config = deepcopy(config)
    if not config.auto_configure:
        log.info("Proportion auto-configuration is disabled; preserving method and plots.")
        return config

    # Extract metadata
    if config.condition_col not in adata.obs:
        log.warning(f"Condition column '{config.condition_col}' not found. Skipping auto-config.")
        return config

    sample_meta = adata.obs[[config.sample_col, config.condition_col]].drop_duplicates()
    condition_counts = sample_meta[config.condition_col].value_counts()
    n_groups = len(condition_counts)
    min_reps = condition_counts.min()
    max_reps = condition_counts.max()

    is_paired = False
    if config.pairing_col and config.pairing_col in adata.obs.columns:
        pair_counts = (
            adata.obs[[config.sample_col, config.pairing_col]]
            .drop_duplicates()[config.pairing_col]
            .value_counts()
        )
        if (pair_counts > 1).all():
            is_paired = True

    log.info(
        f"Auto-config detected: {n_groups} groups, min reps={min_reps}, max reps={max_reps}, paired={is_paired}"
    )

    # Auto-select test method. Prefer compositional sample-level tests by default.
    suggested_method = config.test_method

    if n_groups > 2:
        if config.test_method not in {"anova", "kruskal"}:
            suggested_method = "anova"
    elif is_paired:
        suggested_method = "clr-paired-wilcoxon" if min_reps >= 5 else "clr-paired-t-test"
    elif config.batch_col and config.batch_col in adata.obs.columns:
        suggested_method = "clr-ols"
    elif min_reps == 1:
        log.warning(
            "Detected N=1 in at least one group. Formal proportion inference will be "
            "descriptive/underpowered; p-values may be NaN."
        )
        suggested_method = "clr-t-test"
    elif min_reps == 2:
        if config.test_method in {"wilcoxon", "clr-wilcoxon"}:
            log.info("N=2 is too small for Wilcoxon power. Suggesting 'clr-t-test'.")
            suggested_method = "clr-t-test"

    # Update method
    if suggested_method != config.test_method:
        log.warning(
            f"Auto-config suggests '{suggested_method}' instead of '{config.test_method}' "
            f"based on data characteristics (n_groups={n_groups}, min_reps={min_reps}). "
            f"Set config.auto_configure=False to disable."
        )
        config.test_method = suggested_method

    # Auto-select plot types
    current_plots = set(config.plot_types)

    # N=1 specific adjustments
    if min_reps == 1:
        if "box" in current_plots:
            log.info("Removing 'box' plot (N=1 per group makes boxplots trivial).")
            current_plots.remove("box")
        current_plots.add("bar")
        current_plots.add("diff")

    # High N adjustments
    if min_reps >= 5:
        current_plots.add("box")
        current_plots.add("volcano")

    # Multi-group adjustments
    if n_groups > 2:
        current_plots.add("heatmap")

    config.plot_types = list(current_plots)
    return config


def celltype_proportion_analysis(
    adata: AnnData, config: ProportionConfig
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for cell type proportion analysis.

    This function orchestrates the complete analysis workflow:
    1. Compute cell type proportions
    2. Run statistical tests
    3. Generate visualizations
    4. Export results

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell data object
    config : ProportionConfig
        Configuration object containing analysis parameters

    Returns:
    -------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    stat_df : pd.DataFrame
        Statistical test results

    Examples:
    --------
    >>> from scLucid.analysis import ProportionConfig, celltype_proportion_analysis
    >>> config = ProportionConfig(
    ...     celltype_col='cell_type',
    ...     sample_col='sample_id',
    ...     condition_col='condition',
    ...     test_method='wilcoxon',
    ...     plot_types=['bar', 'box', 'heatmap', 'volcano'],
    ...     out_dir='./proportion_analysis'
    ... )
    >>> prop_df, stat_df = celltype_proportion_analysis(adata, config)
    """
    # Set publication style
    _set_publication_style()

    # Freeze sample/condition/experimental-unit mappings before any automatic
    # method selection can collapse metadata with drop_duplicates().
    sample_meta, design = _validate_proportion_metadata(adata, config)

    # Auto-configure
    config = _auto_configure_analysis(adata, config)
    _validate_proportion_method_design(config, sample_meta, design)
    design["test_method"] = config.test_method
    design["min_samples_per_condition"] = int(config.min_samples_per_condition)

    out_dir = Path(config.out_dir) if config.out_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Compute proportions
    log.info("Computing cell type proportions...")
    prop_df = compute_celltype_proportion(
        adata, celltype_col=config.celltype_col, sample_col=config.sample_col, normalize=True
    )

    # 2. Run statistical tests
    stat_df = pd.DataFrame()
    if config.condition_col:
        log.info(f"Running statistical tests ({config.test_method})...")

        sample_to_cond = sample_meta.loc[prop_df.index, config.condition_col]
        sample_to_pair = (
            sample_meta.loc[prop_df.index, config.pairing_col] if config.pairing_col else None
        )
        sample_to_batch = sample_meta.loc[prop_df.index, config.batch_col] if config.batch_col else None

        stat_df = run_statistical_test(
            prop_df,
            condition_col=config.condition_col,
            test_method=config.test_method,
            sample_to_cond=sample_to_cond,
            sample_to_pair=sample_to_pair,
            sample_to_batch=sample_to_batch,
            composition_pseudocount=config.composition_pseudocount,
            legacy_exploratory=config.legacy_exploratory,
        )

        if not stat_df.empty:
            unit_counts = design["experimental_units_per_condition"]
            if config.pairing_col:
                min_reps = int(design["n_complete_pairs"])
            else:
                min_reps = min(unit_counts.values()) if unit_counts else 0
            formal = min_reps >= int(config.min_samples_per_condition)
            stat_df["replicate_status"] = (
                "replicated" if formal else "single_sample_or_unreplicated"
            )
            semantically_valid = stat_df.get(
                "valid_for_publication_inference",
                pd.Series(False, index=stat_df.index),
            ).astype(bool)
            if "pval" in stat_df.columns:
                estimable = np.isfinite(pd.to_numeric(stat_df["pval"], errors="coerce"))
            else:
                estimable = pd.Series(False, index=stat_df.index)
            stat_df["valid_for_publication_inference"] = (
                semantically_valid & bool(formal) & estimable
            )
            stat_df["experimental_unit_col"] = design["experimental_unit_col"]
            stat_df["replicate_basis"] = design["replicate_basis"]
            stat_df["n_complete_pairs"] = int(design["n_complete_pairs"])
            stat_df["min_independent_units"] = int(min_reps)
            if not formal:
                stat_df["inference_level"] = "descriptive_sample_level"
                stat_df["claim_level"] = "descriptive_effect_size_only"
                stat_df["proportion_warning"] = (
                    "Insufficient biological replicates for publication-level "
                    "cell-type proportion inference."
                )
            design["replicate_requirement_met"] = bool(formal)
            design["valid_for_publication_inference"] = bool(
                stat_df["valid_for_publication_inference"].all()
            )
        else:
            design["replicate_requirement_met"] = False
            design["valid_for_publication_inference"] = False

        # Add effect sizes
        if config.test_method in ["t-test", "wilcoxon", "paired-t-test", "paired-wilcoxon"]:
            from .stats import _add_effect_sizes

            stat_df = _add_effect_sizes(stat_df, prop_df, sample_to_cond, method="cohens_d")

    # 3. Generate plots
    if config.plot_types:
        log.info(f"Generating {len(config.plot_types)} plots...")

        # Prepare data
        sample_meta = sample_meta.reindex(prop_df.index)
        condition = sample_meta[config.condition_col] if config.condition_col else None
        pair = sample_meta[config.pairing_col] if config.pairing_col else None
        batch = sample_meta[config.batch_col] if config.batch_col else None
        timepoints = sample_meta[config.timepoint_col] if config.timepoint_col else None
        ct_palette = config.ct_palette
        condition_palette = config.condition_palette

        for plot_type in config.plot_types:
            try:
                if plot_type == "counts":
                    plot_cell_counts(
                        adata,
                        celltype_col=config.celltype_col,
                        sample_col=config.sample_col,
                        group_col=condition.name if condition is not None else None,
                        palette=ct_palette,
                        out_dir=out_dir,
                    )

                elif plot_type == "bar":
                    sample_order = sorted(prop_df.index, key=_natural_sort_key)
                    plot_proportion_bar(
                        prop_df, sample_order=sample_order, palette=ct_palette, out_dir=out_dir
                    )

                elif plot_type in {"bar_composition", "composition"}:
                    if condition is not None:
                        plot_composition(prop_df, condition=condition, palette=ct_palette, out_dir=out_dir)

                elif plot_type == "box":
                    if condition is not None:
                        plot_box_summary(
                            prop_df,
                            condition=condition,
                            palette=condition_palette,
                            out_dir=out_dir,
                        )

                elif plot_type in {"individual_box", "individual_boxplots"}:
                    if condition is not None and not stat_df.empty:
                        plot_individual_boxplots(
                            prop_df,
                            condition=condition,
                            stat_df=stat_df,
                            palette=condition_palette,
                            out_dir=out_dir,
                        )

                elif plot_type in {"ci", "proportion_ci"}:
                    if condition is not None:
                        plot_proportion_with_ci(
                            prop_df,
                            condition=condition,
                            palette=condition_palette,
                            out_dir=out_dir,
                        )

                elif plot_type == "heatmap":
                    celltype_order = stat_df["cell_type"].values if not stat_df.empty else None
                    plot_proportion_heatmap(
                        prop_df,
                        celltype_order=celltype_order,
                        cluster_samples=True,
                        out_dir=out_dir,
                    )

                elif plot_type in {"clr_heatmap", "composition_heatmap"}:
                    plot_composition_transform_heatmap(prop_df, transform="clr", out_dir=out_dir)

                elif plot_type == "correlation":
                    plot_celltype_correlation(prop_df, out_dir=out_dir)

                elif plot_type == "volcano":
                    if not stat_df.empty:
                        plot_effect_size_volcano(stat_df, out_dir=out_dir)

                elif plot_type == "diff":
                    if condition is not None and not stat_df.empty:
                        plot_diff_stats(prop_df, stat_df, condition=condition, out_dir=out_dir)

                elif plot_type == "shift":
                    if condition is not None:
                        conditions = condition.dropna().astype(str).unique()
                        if len(conditions) == 2:
                            shift_df = prop_df.copy()
                            shift_df[condition.name or config.condition_col] = condition.astype(str)
                            plot_proportion_shifts(
                                shift_df,
                                condition_col=condition.name or config.condition_col,
                                condition1=conditions[0],
                                condition2=conditions[1],
                                palette=ct_palette,
                                out_dir=out_dir,
                            )

                elif plot_type in {"paired_shift", "paired_shifts"}:
                    if condition is not None and pair is not None:
                        conditions = condition.dropna().astype(str).unique()
                        if len(conditions) == 2:
                            plot_paired_proportion_shifts(
                                prop_df,
                                condition=condition.astype(str),
                                pair=pair.astype(str),
                                condition1=conditions[0],
                                condition2=conditions[1],
                                palette=condition_palette,
                                out_dir=out_dir,
                            )

                elif plot_type in {"composition_pca", "pca"}:
                    if condition is not None:
                        plot_composition_pca(
                            prop_df,
                            condition=condition,
                            palette=condition_palette,
                            out_dir=out_dir,
                        )

                elif plot_type == "variability":
                    plot_celltype_variability(prop_df, out_dir=out_dir)

                elif plot_type == "timeseries":
                    if timepoints is not None:
                        # Plot top varying cell types
                        celltype_var = prop_df.var(axis=0)
                        top_celltypes = celltype_var.nlargest(3).index.tolist()

                        for celltype in top_celltypes:
                            plot_proportion_timeseries(
                                prop_df,
                                timepoints=timepoints,
                                celltype=celltype,
                                group_col=condition,
                                palette=condition_palette,
                                out_dir=out_dir,
                            )

                elif plot_type == "batch_pca":
                    if batch is not None:
                        plot_batch_effect(
                            prop_df, batch=batch, method="pca", palette=condition_palette, out_dir=out_dir
                        )

            except Exception as e:
                log.error(f"Failed to generate {plot_type} plot: {e}")

    # 4. Export data
    if out_dir:
        export_analysis_data(prop_df, stat_df, out_dir)
        log.info(f"Analysis complete. Results saved to {out_dir}")

    proportion_ns = adata.uns.setdefault("sclucid", {}).setdefault("proportion", {})
    proportion_ns["design"] = design
    return prop_df, stat_df
