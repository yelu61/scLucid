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
        # Calculate proportions
        prop_df = count_df.div(count_df.sum(axis=1), axis=0)
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
        from pydeseq2 import DESeq2, Preprocessing
    except ImportError:
        log.warning("pydeseq2 not installed. Install with: pip install pydeseq2")
        return pd.DataFrame()

    # Prepare data for DESeq2
    # Transpose to have cell types as rows, samples as columns
    counts = count_df.T
    counts = counts.astype(int)

    # Create metadata DataFrame
    metadata = pd.DataFrame({condition_col: sample_to_cond})

    # Filter low counts
    counts_filtered = Preprocessing.filter_genes(counts, min_cells=1, min_counts=10)

    # Run DESeq2
    dds = DESeq2(counts=counts_filtered, metadata=metadata, design_factors=[condition_col])
    dds.run_deseq()

    # Get results
    res = dds.results_df

    return res


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
    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning("t-test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for celltype in count_df.columns:
        group1 = count_df.loc[sample_to_cond == conditions[0], celltype]
        group2 = count_df.loc[sample_to_cond == conditions[1], celltype]

        # Perform t-test
        stat, pval = stats.ttest_ind(group1, group2, equal_var=False)

        # Calculate mean difference
        mean_diff = group1.mean() - group2.mean()

        results.append(
            {"cell_type": celltype, "statistic": stat, "pval": pval, "mean_diff": mean_diff}
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
    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning("Wilcoxon test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for celltype in count_df.columns:
        group1 = count_df.loc[sample_to_cond == conditions[0], celltype]
        group2 = count_df.loc[sample_to_cond == conditions[1], celltype]

        # Perform Wilcoxon rank-sum test
        stat, pval = stats.mannwhitneyu(group1, group2, alternative="two-sided")

        # Calculate mean difference
        mean_diff = group1.mean() - group2.mean()

        results.append(
            {"cell_type": celltype, "statistic": stat, "pval": pval, "mean_diff": mean_diff}
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
    conditions = sample_to_cond.unique()

    results = []
    for celltype in count_df.columns:
        # Create contingency table
        contingency = pd.DataFrame()

        for cond in conditions:
            cond_samples = sample_to_cond[sample_to_cond == cond].index
            counts = count_df.loc[cond_samples, celltype].sum()
            contingency.loc[celltype, cond] = counts

        # Perform chi-square test
        stat, pval, dof, expected = stats.chi2_contingency(contingency)

        results.append({"cell_type": celltype, "statistic": stat, "pval": pval, "dof": dof})

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
    conditions = sorted(sample_to_cond.unique())
    if len(conditions) != 2:
        log.warning("Paired test requires exactly 2 conditions. " f"Got {len(conditions)}.")
        return pd.DataFrame()

    results = []

    for celltype in count_df.columns:
        # Get paired samples
        pairs = []
        for pair_id in sample_to_pair.unique():
            pair_samples = sample_to_pair[sample_to_pair == pair_id].index

            # Check if we have both conditions for this pair
            if len(pair_samples) == 2:
                cond1_val = count_df.loc[
                    pair_samples[sample_to_cond.loc[pair_samples] == conditions[0]], celltype
                ].values

                cond2_val = count_df.loc[
                    pair_samples[sample_to_cond.loc[pair_samples] == conditions[1]], celltype
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
            stat, pval = stats.wilcoxon(group1, group2)
        else:  # paired t-test
            stat, pval = stats.ttest_rel(group1, group2)

        # Calculate mean difference
        mean_diff = np.mean(group1 - group2)

        results.append(
            {
                "cell_type": celltype,
                "statistic": stat,
                "pval": pval,
                "mean_diff": mean_diff,
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
    conditions = sorted(sample_to_cond.unique())

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

        es = _calculate_effect_size(group1, group2, method)
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
