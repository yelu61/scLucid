"""
Cell type proportion statistical analysis.

This module provides statistical methods for analyzing cell type proportions,
including:
- Proportion computation from count matrices
- Multiple statistical tests (DESeq2, t-test, Wilcoxon, ANOVA, paired tests)
- Effect size calculation (Cohen's d, Cliff's Delta)
- Data export utilities
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

log = logging.getLogger(__name__)


def compute_celltype_proportion(
    adata: AnnData,
    celltype_col: str = "cell_type",
    sample_col: str = "sample_id",
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Compute cell type proportions and counts per sample.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    celltype_col : str
        Column in adata.obs containing cell type labels
    sample_col : str
        Column in adata.obs containing sample identifiers
    normalize : bool
        If True, return proportions; if False, return counts

    Returns:
    -------
    pd.DataFrame
        DataFrame with samples as rows and cell types as columns
    """
    # Extract relevant columns
    df = adata.obs[[sample_col, celltype_col]].copy()

    # Count cells per sample per cell type
    count_df = df.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)

    if normalize:
        # Calculate proportions; guard against zero-row sums
        row_sums = count_df.sum(axis=1)
        zero_sum_mask = row_sums == 0
        if zero_sum_mask.any():
            log.warning(
                f"{zero_sum_mask.sum()} sample(s) have zero total cells; "
                "proportions for these rows will be set to 0."
            )
            row_sums = row_sums.replace(0, np.nan)
        prop_df = count_df.div(row_sums, axis=0).fillna(0.0)
        return prop_df

    return count_df


# ================= Statistical Tests =================


def _run_deseq2(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    condition_col: str,
) -> pd.DataFrame:
    """
    Run DESeq2 differential abundance analysis.

    Parameters
    ----------
    count_df : pd.DataFrame
        Raw count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition
    condition_col : str
        Name of condition column

    Returns:
    -------
    pd.DataFrame
        DESeq2 results with p-values and log2 fold changes
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        log.warning("pydeseq2 not installed. Install with: pip install pydeseq2")
        return pd.DataFrame()

    sample_to_cond = sample_to_cond.reindex(count_df.index).dropna().astype(str)
    conditions = list(pd.unique(sample_to_cond))
    if len(conditions) != 2:
        log.warning("DESeq2 requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    condition1, condition2 = conditions
    counts = count_df.loc[sample_to_cond.index].round().astype(int)
    counts = counts.loc[:, counts.sum(axis=0) >= 10]
    if counts.empty:
        log.warning("DESeq2 skipped because no cell types passed min count filtering.")
        return pd.DataFrame()

    design_col = "__condition"
    metadata = pd.DataFrame(
        {
            design_col: pd.Categorical(
                sample_to_cond,
                categories=[condition1, condition2],
                ordered=False,
            )
        },
        index=sample_to_cond.index,
    )

    try:
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
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
    except Exception as exc:
        log.warning(f"DESeq2 failed: {exc}")
        return pd.DataFrame()

    res = deseq_stats.results_df.copy()
    res = res.rename(
        columns={
            "baseMean": "base_mean",
            "log2FoldChange": "log2fc",
            "stat": "statistic",
            "pvalue": "pval",
        }
    )
    res["cell_type"] = res.index.astype(str)
    res["condition1"] = condition1
    res["condition2"] = condition2
    res["mean_diff"] = (
        counts.loc[sample_to_cond == condition2].mean(axis=0)
        - counts.loc[sample_to_cond == condition1].mean(axis=0)
    ).reindex(res.index).to_numpy()
    res["direction"] = f"{condition2} - {condition1}"
    res["method"] = "deseq2"
    return res.reset_index(drop=True)


def _run_ttest(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run independent samples t-test for each cell type.

    Parameters
    ----------
    count_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition

    Returns:
    -------
    pd.DataFrame
        Test results with p-values and statistics
    """
    conditions = list(sample_to_cond.dropna().unique())
    if len(conditions) != 2:
        log.warning("t-test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    condition1, condition2 = conditions
    results = []
    for celltype in count_df.columns:
        group1 = count_df.loc[sample_to_cond == condition1, celltype]
        group2 = count_df.loc[sample_to_cond == condition2, celltype]

        # Statistic and effect direction are condition2 - condition1.
        stat, pval = stats.ttest_ind(group2, group1, equal_var=False)

        mean_diff = group2.mean() - group1.mean()

        results.append(
            {
                "cell_type": celltype,
                "condition1": condition1,
                "condition2": condition2,
                "statistic": stat,
                "pval": pval,
                "mean_diff": mean_diff,
                "direction": f"{condition2} - {condition1}",
            }
        )

    return pd.DataFrame(results)


def _run_wilcoxon(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run Mann-Whitney U test for each cell type.

    Parameters
    ----------
    count_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition

    Returns:
    -------
    pd.DataFrame
        Test results with p-values and statistics
    """
    conditions = list(sample_to_cond.dropna().unique())
    if len(conditions) != 2:
        log.warning("Wilcoxon test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    condition1, condition2 = conditions
    results = []
    for celltype in count_df.columns:
        group1 = count_df.loc[sample_to_cond == condition1, celltype]
        group2 = count_df.loc[sample_to_cond == condition2, celltype]

        # Perform Wilcoxon rank-sum test
        stat, pval = stats.mannwhitneyu(group2, group1, alternative="two-sided")

        mean_diff = group2.mean() - group1.mean()

        results.append(
            {
                "cell_type": celltype,
                "condition1": condition1,
                "condition2": condition2,
                "statistic": stat,
                "pval": pval,
                "mean_diff": mean_diff,
                "direction": f"{condition2} - {condition1}",
            }
        )

    return pd.DataFrame(results)


def _run_anova(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run one-way ANOVA for each cell type.

    Parameters
    ----------
    count_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition

    Returns:
    -------
    pd.DataFrame
        Test results with p-values and F-statistics
    """
    conditions = sample_to_cond.unique()
    if len(conditions) < 3:
        log.warning("ANOVA requires 3+ conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for celltype in count_df.columns:
        groups = [count_df.loc[sample_to_cond == cond, celltype] for cond in conditions]

        # Perform one-way ANOVA
        stat, pval = stats.f_oneway(*groups)

        results.append({"cell_type": celltype, "statistic": stat, "pval": pval})

    return pd.DataFrame(results)


def _run_kruskal(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run Kruskal-Wallis H-test for each cell type.

    This is the non-parametric multi-condition counterpart to one-way ANOVA.
    """
    conditions = sample_to_cond.unique()
    if len(conditions) < 2:
        log.warning("Kruskal-Wallis test requires 2+ conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for celltype in count_df.columns:
        groups = [count_df.loc[sample_to_cond == cond, celltype].dropna() for cond in conditions]
        groups = [group for group in groups if len(group) > 0]
        if len(groups) < 2:
            continue

        stat, pval = stats.kruskal(*groups)

        results.append({"cell_type": celltype, "statistic": stat, "pval": pval})

    return pd.DataFrame(results)


def _run_contingency_test(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run chi-square contingency table test.

    Parameters
    ----------
    count_df : pd.DataFrame
        Raw count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition

    Returns:
    -------
    pd.DataFrame
        Test results with chi-square statistics and p-values
    """
    conditions = list(sample_to_cond.dropna().unique())
    if len(conditions) < 2:
        log.warning("Chi-square test requires 2+ conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    observed = pd.DataFrame(index=conditions, columns=count_df.columns, dtype=float)
    for cond in conditions:
        cond_samples = sample_to_cond[sample_to_cond == cond].index
        observed.loc[cond] = count_df.loc[cond_samples].sum(axis=0)

    if observed.to_numpy().sum() == 0:
        log.warning("Chi-square test skipped because the contingency table is empty.")
        return pd.DataFrame()

    stat, pval, dof, expected = stats.chi2_contingency(observed)
    expected_df = pd.DataFrame(expected, index=observed.index, columns=observed.columns)
    # Guard against zero expected values in residuals
    with np.errstate(divide="ignore", invalid="ignore"):
        residuals = (observed - expected_df) / np.sqrt(expected_df)
        contributions = ((observed - expected_df) ** 2 / expected_df).sum(axis=0)
    residuals = residuals.fillna(0.0).replace([np.inf, -np.inf], 0.0)
    contributions = contributions.fillna(0.0).replace([np.inf, -np.inf], 0.0)

    results = []
    for celltype in observed.columns:
        row = {
            "cell_type": celltype,
            "statistic": float(contributions[celltype]),
            "pval": float(pval),
            "overall_statistic": float(stat),
            "overall_pval": float(pval),
            "dof": int(dof),
            "method_note": "Per-cell statistics are chi-square contributions; p-values are global.",
        }
        for cond in conditions:
            row[f"observed_{cond}"] = float(observed.loc[cond, celltype])
            row[f"expected_{cond}"] = float(expected_df.loc[cond, celltype])
            row[f"std_residual_{cond}"] = float(residuals.loc[cond, celltype])
        results.append(row)

    return pd.DataFrame(results)


def _run_paired_test(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    sample_to_pair: pd.Series,
    test_type: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Run paired statistical test for each cell type.

    Parameters
    ----------
    count_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Mapping from sample to condition (must have exactly 2)
    sample_to_pair : pd.Series
        Mapping from sample to pairing identifier (e.g., patient)
    test_type : str
        Type of test: 'wilcoxon' or 't-test'

    Returns:
    -------
    pd.DataFrame
        Test results with p-values and statistics
    """
    conditions = list(sample_to_cond.dropna().unique())
    if len(conditions) != 2:
        log.warning("Paired test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()
    condition1, condition2 = conditions

    results = []

    for celltype in count_df.columns:
        # Get paired samples
        pairs = []
        for pair_id in sample_to_pair.unique():
            pair_samples = sample_to_pair[sample_to_pair == pair_id].index

            # Check if we have both conditions for this pair
            if len(pair_samples) == 2:
                cond1_val = count_df.loc[
                    pair_samples[sample_to_cond.loc[pair_samples] == condition1], celltype
                ].values

                cond2_val = count_df.loc[
                    pair_samples[sample_to_cond.loc[pair_samples] == condition2], celltype
                ].values

                if len(cond1_val) > 0 and len(cond2_val) > 0:
                    pairs.append((cond1_val[0], cond2_val[0]))

        if len(pairs) < 3:
            log.warning(f"Insufficient pairs for {celltype}: {len(pairs)}. " "Skipping.")
            continue

        # Extract paired values
        group1 = np.array([p[0] for p in pairs])
        group2 = np.array([p[1] for p in pairs])

        # Perform test
        if test_type == "wilcoxon":
            stat, pval = stats.wilcoxon(group2, group1)
        else:  # paired t-test
            stat, pval = stats.ttest_rel(group2, group1)

        mean_diff = np.mean(group2 - group1)

        results.append(
            {
                "cell_type": celltype,
                "condition1": condition1,
                "condition2": condition2,
                "statistic": stat,
                "pval": pval,
                "mean_diff": mean_diff,
                "direction": f"{condition2} - {condition1}",
                "n_pairs": len(pairs),
            }
        )

    return pd.DataFrame(results)


def run_statistical_test(
    count_df: pd.DataFrame,
    condition_col: str,
    test_method: str = "wilcoxon",
    sample_to_cond: Optional[pd.Series] = None,
    sample_to_pair: Optional[pd.Series] = None,
    multiple_testing_correction: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Run statistical tests for differential cell type abundance.

    Parameters
    ----------
    count_df : pd.DataFrame
        Proportion or count matrix (samples × cell types)
        Must have index matching sample_to_cond
    condition_col : str
        Name of condition column
    test_method : str
        Statistical method: 'deseq2', 't-test', 'wilcoxon', 'anova',
        'paired-t-test', 'paired-wilcoxon'
    sample_to_cond : pd.Series, optional
        Mapping from sample to condition
    sample_to_pair : pd.Series, optional
        Mapping from sample to pairing identifier (for paired tests)
    multiple_testing_correction : str
        Method for multiple testing correction (see statsmodels)

    Returns:
    -------
    pd.DataFrame
        Test results with p-values, adjusted p-values, and statistics
    """
    if sample_to_cond is None:
        # Assume index is sample_id and need to map from count_df
        log.warning("sample_to_cond not provided. Using count_df index.")
        sample_to_cond = pd.Series(index=count_df.index, data=range(len(count_df)))

    # Dispatch to appropriate test function
    if test_method == "deseq2":
        res_df = _run_deseq2(count_df, sample_to_cond, condition_col)
    elif test_method == "t-test":
        res_df = _run_ttest(count_df, sample_to_cond)
    elif test_method == "wilcoxon":
        res_df = _run_wilcoxon(count_df, sample_to_cond)
    elif test_method == "anova":
        res_df = _run_anova(count_df, sample_to_cond)
    elif test_method == "kruskal":
        res_df = _run_kruskal(count_df, sample_to_cond)
    elif test_method == "chi-square":
        res_df = _run_contingency_test(count_df, sample_to_cond)
    elif test_method == "paired-t-test":
        if sample_to_pair is None:
            raise ValueError("sample_to_pair required for paired tests")
        res_df = _run_paired_test(count_df, sample_to_cond, sample_to_pair, "t-test")
    elif test_method == "paired-wilcoxon":
        if sample_to_pair is None:
            raise ValueError("sample_to_pair required for paired tests")
        res_df = _run_paired_test(count_df, sample_to_cond, sample_to_pair, "wilcoxon")
    else:
        raise ValueError(f"Unknown test method: {test_method}")

    if res_df.empty:
        return res_df

    # Multiple testing correction
    if "pval" in res_df.columns and multiple_testing_correction:
        try:
            from statsmodels.stats.multitest import multipletests

            _, res_df["padj"], _, _ = multipletests(
                res_df["pval"], method=multiple_testing_correction
            )
        except ImportError:
            log.warning(
                "statsmodels not installed. Skipping correction. "
                "Install with: pip install statsmodels"
            )
            res_df["padj"] = res_df["pval"]

    # Sort by adjusted p-value
    if "padj" in res_df.columns:
        res_df = res_df.sort_values("padj")
    elif "pval" in res_df.columns:
        res_df = res_df.sort_values("pval")

    return res_df


# ================= Effect Size =================


def _calculate_effect_size(group1: pd.Series, group2: pd.Series, method: str = "cohens_d") -> float:
    """
    Calculate effect size between two groups.

    Parameters
    ----------
    group1, group2 : pd.Series
        Data values for two groups
    method : str
        Effect size method: 'cohens_d' or 'cliffs_delta'

    Returns:
    -------
    float
        Effect size value
    """
    if method == "cohens_d":
        # Cohen's d
        n1, n2 = len(group1), len(group2)
        var1, var2 = group1.var(), group2.var()

        # Pooled standard deviation
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

        # Cohen's d
        d = (group1.mean() - group2.mean()) / pooled_std
        return d

    elif method == "cliffs_delta":
        # Cliff's Delta (non-parametric)
        n1, n2 = len(group1), len(group2)

        # Count pairwise comparisons
        greater = 0
        less = 0

        for x in group1:
            for y in group2:
                if x > y:
                    greater += 1
                elif x < y:
                    less += 1

        # Cliff's Delta
        delta = (greater - less) / (n1 * n2)
        return delta

    else:
        raise ValueError(f"Unknown effect size method: {method}")


def _add_effect_sizes(
    res_df: pd.DataFrame,
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    method: str = "cohens_d",
) -> pd.DataFrame:
    """
    Add effect sizes to results DataFrame.

    Parameters
    ----------
    res_df : pd.DataFrame
        Results from run_statistical_test
    count_df : pd.DataFrame
        Proportion matrix
    sample_to_cond : pd.Series
        Mapping from sample to condition
    method : str
        Effect size method

    Returns:
    -------
    pd.DataFrame
        Results DataFrame with effect size column added
    """
    if {"condition1", "condition2"}.issubset(res_df.columns):
        conditions = [res_df["condition1"].iloc[0], res_df["condition2"].iloc[0]]
    else:
        conditions = list(sample_to_cond.dropna().unique())

    if len(conditions) != 2:
        log.warning(
            "Effect size calculation requires exactly 2 conditions. " f"Got {len(conditions)}."
        )
        return res_df

    effect_sizes = []

    for _, row in res_df.iterrows():
        celltype = row["cell_type"]

        group1 = count_df.loc[sample_to_cond == conditions[0], celltype]
        group2 = count_df.loc[sample_to_cond == conditions[1], celltype]

        es = _calculate_effect_size(group2, group1, method)
        effect_sizes.append(es)

    res_df[f"effect_size_{method}"] = effect_sizes

    return res_df


# ================= Data Export =================


def export_analysis_data(
    prop_df: pd.DataFrame,
    stat_df: pd.DataFrame,
    out_dir: Union[str, Path],
    prefix: str = "proportion",
):
    """
    Export analysis results to CSV files.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix
    stat_df : pd.DataFrame
        Statistical test results
    out_dir : str or Path
        Output directory
    prefix : str
        Prefix for output files
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Export proportion matrix
    prop_df.to_csv(out_dir / f"{prefix}_matrix.csv")

    # Export statistical results
    if not stat_df.empty:
        stat_df.to_csv(out_dir / f"{prefix}_stats.csv", index=False)

    log.info(f"Exported analysis data to {out_dir}")
