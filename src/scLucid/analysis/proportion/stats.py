"""
Cell type proportion statistical analysis.

This module provides statistical methods for analyzing cell type proportions,
including:
- Proportion computation from count matrices
- CLR-transformed sample-level tests for compositional proportions
- Optional legacy raw-proportion tests for exploratory summaries
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


def _tag_proportion_result(res_df: pd.DataFrame, *, test_method: str) -> pd.DataFrame:
    """Add explicit inference semantics to cell-type proportion result tables."""
    if res_df.empty:
        return res_df
    tagged = res_df.copy()
    legacy_raw_methods = {
        "t-test",
        "wilcoxon",
        "anova",
        "kruskal",
        "paired-t-test",
        "paired-wilcoxon",
    }
    if test_method in legacy_raw_methods:
        tagged["inference_level"] = "exploratory_legacy_proportion"
        tagged["claim_level"] = "exploratory_hypothesis_generation"
        tagged["valid_for_publication_inference"] = False
        tagged["model_type"] = "raw_proportion_legacy_test"
        tagged["compositional_data_warning"] = (
            "Raw cell-type proportions are compositional; prefer CLR-transformed "
            "sample-level tests or a compositional model for formal inference."
        )
    elif test_method == "chi-square":
        tagged["inference_level"] = "global_contingency_screen"
        tagged["claim_level"] = "exploratory_global_composition_screen"
        tagged["valid_for_publication_inference"] = False
        tagged["model_type"] = "global_chi_square_contingency_contribution"
        tagged["compositional_data_warning"] = (
            "Chi-square reports a global contingency test; per-cell-type rows are contributions, "
            "not independent cell-type-specific p-values."
        )
    elif test_method == "deseq2":
        tagged["inference_level"] = tagged.get("inference_level", "sample_level")
        tagged["claim_level"] = "sample_level_compositional_count_model"
        tagged["valid_for_publication_inference"] = True
        tagged["model_type"] = "sample_level_deseq2_differential_abundance"
    elif test_method.startswith("clr-"):
        tagged["inference_level"] = tagged.get("inference_level", "sample_level")
        tagged["claim_level"] = "sample_level_clr_compositional_inference"
        tagged["valid_for_publication_inference"] = tagged.get(
            "valid_for_publication_inference",
            True,
        )
        tagged["model_type"] = "sample_level_clr_test"
        if "compositional_data_warning" not in tagged:
            tagged["compositional_data_warning"] = (
                "Cell type proportions are compositional; inference was run on sample-level "
                "CLR-transformed values."
            )
    elif test_method == "ancom-like-clr":
        tagged["inference_level"] = "sample_level"
        tagged["claim_level"] = "ancom_like_clr_heuristic"
        tagged["valid_for_publication_inference"] = False
        tagged["model_type"] = "sample_level_ancom_like_clr_heuristic"
        tagged["compositional_data_warning"] = (
            "ANCOM-like CLR heuristic is a simplified Python approximation of "
            "ANCOM-style testing; it is not the full ANCOM-BC model."
        )
    else:
        tagged["claim_level"] = "review_required_proportion_inference"
        tagged["valid_for_publication_inference"] = False
        tagged["model_type"] = "unknown_proportion_test"

    tagged["recommended_formal_inference"] = "sample_level_clr_or_compositional_model"
    tagged["proportion_review_note"] = tagged["claim_level"].map(
        {
            "sample_level_clr_compositional_inference": (
                "Sample-level CLR result; review replicate depth, pairing/batch design, and compositional assumptions."
            ),
            "sample_level_compositional_count_model": (
                "Sample-level count model result; review design and biological replicate balance."
            ),
            "exploratory_hypothesis_generation": (
                "Legacy raw-proportion test; use for screening only, not formal compositional inference."
            ),
            "exploratory_global_composition_screen": (
                "Global contingency screen; per-cell-type rows summarize contributions to one global test."
            ),
            "ancom_like_clr_heuristic": (
                "Simplified Python approximation of ANCOM-style CLR testing; "
                "not the full ANCOM-BC model. Use for exploratory screening only."
            ),
            "review_required_proportion_inference": (
                "Proportion inference semantics are ambiguous; review method and sample design."
            ),
        }
    )
    return tagged


def composition_transform(
    prop_df: pd.DataFrame,
    *,
    method: str = "clr",
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    """Transform sample-level compositions for valid compositional testing."""
    if method == "none":
        return prop_df.astype(float).copy()
    if method != "clr":
        raise ValueError(f"Unknown composition transform: {method}")
    values = prop_df.astype(float).copy()
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError("CLR composition transform requires finite values.")
    if (values < 0).to_numpy(dtype=bool).any():
        raise ValueError("CLR composition transform requires non-negative values.")

    row_sums = values.sum(axis=1)
    nonzero_rows = row_sums > 0
    if not nonzero_rows.all():
        log.warning(
            "%d sample(s) have zero total composition; CLR values for those rows will be set to 0.",
            int((~nonzero_rows).sum()),
        )

    positive_sums = row_sums[nonzero_rows]
    if positive_sums.empty:
        values.loc[:, :] = 0.0
    elif np.allclose(positive_sums, 1.0, atol=0.01, rtol=0.0) and values.max().max() <= 1.01:
        values.loc[nonzero_rows] = values.loc[nonzero_rows].div(positive_sums, axis=0)
    elif np.allclose(positive_sums, 100.0, atol=1.0, rtol=0.0) and values.max().max() <= 100.0:
        values.loc[nonzero_rows] = values.loc[nonzero_rows].div(100.0)
    else:
        values.loc[nonzero_rows] = values.loc[nonzero_rows].div(positive_sums, axis=0)
    values.loc[~nonzero_rows] = 0.0
    transformed = np.log(values + float(pseudocount))
    transformed = transformed.sub(transformed.mean(axis=1), axis=0)
    return transformed


def _mean_diff_ci(
    group1: pd.Series,
    group2: pd.Series,
    *,
    paired: bool = False,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Return condition2-condition1 mean difference and t-based confidence interval."""
    group1 = pd.to_numeric(group1, errors="coerce").dropna()
    group2 = pd.to_numeric(group2, errors="coerce").dropna()
    if paired:
        diffs = (group2.to_numpy() - group1.to_numpy()).astype(float)
    else:
        diffs = np.r_[
            group2.to_numpy(dtype=float) - float(group1.mean()),
        ]
        mean_diff = float(group2.mean() - group1.mean())
        se = np.sqrt(group1.var(ddof=1) / len(group1) + group2.var(ddof=1) / len(group2))
        df = min(len(group1), len(group2)) - 1
        if len(group1) < 2 or len(group2) < 2 or not np.isfinite(se) or se == 0:
            return mean_diff, np.nan, np.nan
        tcrit = stats.t.ppf(1 - alpha / 2, df=max(1, df))
        return mean_diff, float(mean_diff - tcrit * se), float(mean_diff + tcrit * se)

    if diffs.size < 2:
        return float(np.mean(diffs)) if diffs.size else np.nan, np.nan, np.nan
    mean_diff = float(np.mean(diffs))
    se = float(stats.sem(diffs, nan_policy="omit"))
    if not np.isfinite(se) or se == 0:
        return mean_diff, np.nan, np.nan
    tcrit = stats.t.ppf(1 - alpha / 2, df=max(1, diffs.size - 1))
    return mean_diff, float(mean_diff - tcrit * se), float(mean_diff + tcrit * se)


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


def _run_clr_sample_level_test(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    *,
    test_type: str = "t-test",
    sample_to_pair: Optional[pd.Series] = None,
    sample_to_batch: Optional[pd.Series] = None,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    """Run sample-level tests on CLR-transformed compositions."""
    conditions = list(sample_to_cond.dropna().astype(str).unique())
    if len(conditions) != 2:
        log.warning("CLR sample-level tests require exactly 2 conditions. Got %s.", len(conditions))
        return pd.DataFrame()
    condition1, condition2 = conditions
    aligned = count_df.loc[sample_to_cond.dropna().index]
    sample_to_cond = sample_to_cond.loc[aligned.index].astype(str)
    clr_df = composition_transform(aligned, method="clr", pseudocount=pseudocount)

    results = []
    for celltype in clr_df.columns:
        group1 = clr_df.loc[sample_to_cond == condition1, celltype]
        group2 = clr_df.loc[sample_to_cond == condition2, celltype]
        if len(group1) < 1 or len(group2) < 1:
            continue

        n1, n2 = len(group1), len(group2)
        stat = np.nan
        pval = np.nan
        n_pairs = np.nan
        ci_lower = np.nan
        ci_upper = np.nan
        mean_diff = float(group2.mean() - group1.mean())

        if test_type == "paired-t-test" or test_type == "paired-wilcoxon":
            if sample_to_pair is None:
                raise ValueError("sample_to_pair required for paired CLR tests")
            pairs = []
            pair_series = sample_to_pair.loc[clr_df.index].astype(str)
            for pair_id in pair_series.dropna().unique():
                pair_samples = pair_series[pair_series == pair_id].index
                conds = sample_to_cond.loc[pair_samples]
                vals1 = clr_df.loc[pair_samples[conds == condition1], celltype].to_numpy()
                vals2 = clr_df.loc[pair_samples[conds == condition2], celltype].to_numpy()
                if vals1.size and vals2.size:
                    pairs.append((float(vals1[0]), float(vals2[0])))
            n_pairs = len(pairs)
            if n_pairs >= 2:
                paired1 = pd.Series([p[0] for p in pairs], dtype=float)
                paired2 = pd.Series([p[1] for p in pairs], dtype=float)
                if test_type == "paired-wilcoxon":
                    stat, pval = stats.wilcoxon(paired2, paired1)
                else:
                    stat, pval = stats.ttest_rel(paired2, paired1)
                mean_diff, ci_lower, ci_upper = _mean_diff_ci(paired1, paired2, paired=True)
        elif test_type == "wilcoxon":
            if n1 >= 2 and n2 >= 2:
                stat, pval = stats.mannwhitneyu(group2, group1, alternative="two-sided")
                mean_diff, ci_lower, ci_upper = _mean_diff_ci(group1, group2)
        elif test_type == "ols":
            if sample_to_batch is None:
                log.warning("CLR OLS requested without sample_to_batch; falling back to t-test.")
                if n1 >= 2 and n2 >= 2:
                    stat, pval = stats.ttest_ind(group2, group1, equal_var=False)
                    mean_diff, ci_lower, ci_upper = _mean_diff_ci(group1, group2)
            else:
                try:
                    import statsmodels.formula.api as smf
                except ImportError:
                    log.warning("statsmodels not installed; CLR OLS cannot be run.")
                else:
                    model_df = pd.DataFrame(
                        {
                            "value": clr_df[celltype],
                            "condition": sample_to_cond,
                            "batch": sample_to_batch.loc[clr_df.index].astype(str),
                        }
                    ).dropna()
                    term = f"C(condition)[T.{condition2}]"
                    if model_df["condition"].nunique() == 2 and model_df["batch"].nunique() > 1:
                        fit = smf.ols("value ~ C(condition) + C(batch)", data=model_df).fit()
                        if term in fit.params:
                            stat = float(fit.tvalues[term])
                            pval = float(fit.pvalues[term])
                            mean_diff = float(fit.params[term])
                            ci = fit.conf_int().loc[term]
                            ci_lower, ci_upper = float(ci.iloc[0]), float(ci.iloc[1])
        else:
            if n1 >= 2 and n2 >= 2:
                stat, pval = stats.ttest_ind(group2, group1, equal_var=False)
                mean_diff, ci_lower, ci_upper = _mean_diff_ci(group1, group2)

        results.append(
            {
                "cell_type": celltype,
                "condition1": condition1,
                "condition2": condition2,
                "statistic": float(stat) if np.isfinite(stat) else np.nan,
                "pval": float(pval) if np.isfinite(pval) else np.nan,
                "mean_diff": mean_diff,
                "effect_size": mean_diff,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": int(n1),
                "n_samples_condition2": int(n2),
                "n_pairs": n_pairs,
                "method": f"clr-{test_type}",
                "transform": "clr",
                "inference_level": "sample_level",
                "compositional_data_warning": (
                    "Cell type proportions are compositional; inference was run on "
                    "sample-level CLR-transformed values."
                ),
            }
        )

    return pd.DataFrame(results)


def _run_ancom_like_clr_test(
    count_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_col: str,
    sample_col: str,
    fdr: float = 0.05,
    *,
    test_type: str = "t-test",
    pair_col: Optional[str] = None,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    """
    Simplified ANCOM-style CLR differential abundance test.

    This is a Python-native heuristic inspired by ANCOM V2. It runs a
    two-sample test on CLR-transformed compositions and computes an ANCOM-like
    W-statistic from pairwise associations across cell types. It is **not** the
    full ANCOM-BC model and should be treated as exploratory.

    Parameters
    ----------
    count_df : pd.DataFrame
        Raw or proportional count matrix (samples x cell types).
    metadata_df : pd.DataFrame
        Sample-level metadata containing at least ``sample_col`` and
        ``condition_col``. Optionally ``pair_col`` for paired designs.
    condition_col : str
        Column in ``metadata_df`` defining the two conditions to compare.
    sample_col : str
        Column in ``metadata_df`` identifying samples; must match ``count_df`` index.
    fdr : float
        False-discovery threshold for significance flags.
    test_type : str
        "t-test" (default) or "wilcoxon" for the per-cell-type test.
    pair_col : Optional[str]
        Column in ``metadata_df`` defining pairs; if provided, a paired test is run.
    pseudocount : float
        Small pseudo-count added before the CLR transform.

    Returns:
    -------
    pd.DataFrame
        Result table with per-cell-type statistics, W-statistic, and inference tags.
    """
    if sample_col not in metadata_df.columns:
        raise ValueError(f"sample_col '{sample_col}' not found in metadata_df")
    if condition_col not in metadata_df.columns:
        raise ValueError(f"condition_col '{condition_col}' not found in metadata_df")
    if pair_col is not None and pair_col not in metadata_df.columns:
        raise ValueError(f"pair_col '{pair_col}' not found in metadata_df")

    sample_meta = metadata_df.set_index(sample_col)
    common_idx = count_df.index.intersection(sample_meta.index)
    if common_idx.empty:
        log.warning("No matching samples between count_df and metadata_df")
        return pd.DataFrame()

    aligned = count_df.loc[common_idx]
    sample_meta = sample_meta.loc[common_idx]

    # CLR transform after adding pseudo-count
    prop_df = aligned.div(aligned.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    clr_df = np.log(prop_df + float(pseudocount))
    clr_df = clr_df.sub(clr_df.mean(axis=1), axis=0)

    conditions = sample_meta[condition_col].dropna().astype(str).unique()
    if len(conditions) != 2:
        log.warning("ANCOM-like CLR test requires exactly 2 conditions. Got %s.", len(conditions))
        return pd.DataFrame()
    condition1, condition2 = conditions

    # Pairwise correlation W-statistic (simplified ANCOM V2 heuristic)
    w_stats = pd.Series(0, index=clr_df.columns, dtype=int)
    n_ct = len(clr_df.columns)
    if n_ct > 1:
        corr_pvals = pd.DataFrame(np.nan, index=clr_df.columns, columns=clr_df.columns)
        for i, ct_i in enumerate(clr_df.columns):
            for j, ct_j in enumerate(clr_df.columns):
                if i >= j:
                    continue
                x = clr_df[ct_i].to_numpy(dtype=float)
                y = clr_df[ct_j].to_numpy(dtype=float)
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() < 3:
                    continue
                with np.errstate(invalid="ignore"):
                    _, p = stats.pearsonr(x[mask], y[mask])
                corr_pvals.loc[ct_i, ct_j] = p
                corr_pvals.loc[ct_j, ct_i] = p

        unique_mask = np.triu(np.ones((n_ct, n_ct), dtype=bool), k=1)
        unique_pvals = corr_pvals.where(unique_mask).stack().dropna()
        if not unique_pvals.empty:
            try:
                from statsmodels.stats.multitest import multipletests

                _, adj, _, _ = multipletests(unique_pvals.values, alpha=fdr, method="fdr_bh")
                sig_pairs = unique_pvals.index[adj < fdr]
            except ImportError:
                log.warning("statsmodels not installed; using uncorrected p-values for W-statistic.")
                sig_pairs = unique_pvals.index[unique_pvals.values < fdr]
            for ct_i, ct_j in sig_pairs:
                w_stats.loc[ct_i] += 1
                w_stats.loc[ct_j] += 1

    # Per-cell-type differential test on CLR values
    results = []
    for celltype in clr_df.columns:
        group1 = clr_df.loc[sample_meta[condition_col] == condition1, celltype]
        group2 = clr_df.loc[sample_meta[condition_col] == condition2, celltype]
        n1, n2 = len(group1), len(group2)
        stat = np.nan
        pval = np.nan
        n_pairs = np.nan

        if pair_col is not None:
            pairs = []
            for pair_id in sample_meta[pair_col].dropna().astype(str).unique():
                pair_samples = sample_meta[sample_meta[pair_col] == pair_id].index
                conds = sample_meta.loc[pair_samples, condition_col]
                vals1 = clr_df.loc[pair_samples[conds == condition1], celltype].to_numpy()
                vals2 = clr_df.loc[pair_samples[conds == condition2], celltype].to_numpy()
                if vals1.size and vals2.size:
                    pairs.append((float(vals1[0]), float(vals2[0])))
            n_pairs = len(pairs)
            if n_pairs >= 2:
                paired1 = np.array([p[0] for p in pairs], dtype=float)
                paired2 = np.array([p[1] for p in pairs], dtype=float)
                if test_type == "wilcoxon":
                    stat, pval = stats.wilcoxon(paired2, paired1)
                else:
                    stat, pval = stats.ttest_rel(paired2, paired1)
        elif n1 >= 2 and n2 >= 2:
            if test_type == "wilcoxon":
                stat, pval = stats.mannwhitneyu(group2, group1, alternative="two-sided")
            else:
                stat, pval = stats.ttest_ind(group2, group1, equal_var=False)

        mean_diff = float(group2.mean() - group1.mean())
        results.append(
            {
                "cell_type": celltype,
                "condition1": condition1,
                "condition2": condition2,
                "statistic": float(stat) if np.isfinite(stat) else np.nan,
                "pval": float(pval) if np.isfinite(pval) else np.nan,
                "mean_diff": mean_diff,
                "effect_size": mean_diff,
                "direction": f"{condition2} - {condition1}",
                "n_samples_condition1": int(n1),
                "n_samples_condition2": int(n2),
                "n_pairs": n_pairs,
                "w_statistic": int(w_stats.loc[celltype]),
                "method": "ancom-like-clr",
                "transform": "clr",
            }
        )

    res_df = pd.DataFrame(results)
    if res_df.empty:
        return res_df

    # Multiple-testing correction and significance flag
    try:
        from statsmodels.stats.multitest import multipletests

        valid = res_df["pval"].notna() & np.isfinite(res_df["pval"])
        res_df["padj"] = np.nan
        if valid.any():
            _, adj, _, _ = multipletests(res_df.loc[valid, "pval"], alpha=fdr, method="fdr_bh")
            res_df.loc[valid, "padj"] = adj
    except ImportError:
        log.warning("statsmodels not installed; skipping multiple-testing correction.")
        res_df["padj"] = res_df["pval"]

    res_df["significant"] = res_df["padj"].lt(fdr).fillna(False)

    return _tag_proportion_result(res_df, test_method="ancom-like-clr")


def run_statistical_test(
    count_df: pd.DataFrame,
    condition_col: str,
    test_method: str = "wilcoxon",
    sample_to_cond: Optional[pd.Series] = None,
    sample_to_pair: Optional[pd.Series] = None,
    sample_to_batch: Optional[pd.Series] = None,
    multiple_testing_correction: str = "fdr_bh",
    composition_pseudocount: float = 1e-6,
    legacy_exploratory: bool = False,
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
        Statistical method. Prefer 'clr-t-test', 'clr-wilcoxon',
        'clr-paired-t-test', 'clr-paired-wilcoxon', or 'clr-ols' for
        sample-level compositional inference. Raw-proportion 't-test',
        'wilcoxon', 'anova', 'kruskal', 'paired-t-test', and
        'paired-wilcoxon' are retained for legacy exploratory summaries.
    sample_to_cond : pd.Series, optional
        Mapping from sample to condition
    sample_to_pair : pd.Series, optional
        Mapping from sample to pairing identifier (for paired tests)
    multiple_testing_correction : str
        Method for multiple testing correction (see statsmodels)
    legacy_exploratory : bool
        If True, allow raw-proportion tests that violate compositional
        assumptions. Required for legacy 't-test', 'wilcoxon', 'anova',
        'kruskal', 'paired-t-test', and 'paired-wilcoxon' methods.

    Returns:
    -------
    pd.DataFrame
        Test results with p-values, adjusted p-values, and statistics
    """
    legacy_raw_methods = {
        "t-test",
        "wilcoxon",
        "anova",
        "kruskal",
        "paired-t-test",
        "paired-wilcoxon",
    }
    if test_method in legacy_raw_methods and not legacy_exploratory:
        raise ValueError(
            "Raw proportion tests are exploratory and violate compositional assumptions. "
            "Set legacy_exploratory=True to acknowledge."
        )

    if sample_to_cond is None:
        # Assume index is sample_id and need to map from count_df
        log.warning("sample_to_cond not provided. Using count_df index.")
        sample_to_cond = pd.Series(index=count_df.index, data=range(len(count_df)))

    # Dispatch to appropriate test function
    if test_method == "deseq2":
        res_df = _run_deseq2(count_df, sample_to_cond, condition_col)
    elif test_method == "clr-t-test":
        res_df = _run_clr_sample_level_test(
            count_df, sample_to_cond, test_type="t-test", pseudocount=composition_pseudocount
        )
    elif test_method == "clr-wilcoxon":
        res_df = _run_clr_sample_level_test(
            count_df, sample_to_cond, test_type="wilcoxon", pseudocount=composition_pseudocount
        )
    elif test_method == "clr-paired-t-test":
        if sample_to_pair is None:
            raise ValueError("sample_to_pair required for paired tests")
        res_df = _run_clr_sample_level_test(
            count_df,
            sample_to_cond,
            test_type="paired-t-test",
            sample_to_pair=sample_to_pair,
            pseudocount=composition_pseudocount,
        )
    elif test_method == "clr-paired-wilcoxon":
        if sample_to_pair is None:
            raise ValueError("sample_to_pair required for paired tests")
        res_df = _run_clr_sample_level_test(
            count_df,
            sample_to_cond,
            test_type="paired-wilcoxon",
            sample_to_pair=sample_to_pair,
            pseudocount=composition_pseudocount,
        )
    elif test_method == "clr-ols":
        res_df = _run_clr_sample_level_test(
            count_df,
            sample_to_cond,
            test_type="ols",
            sample_to_batch=sample_to_batch,
            pseudocount=composition_pseudocount,
        )
    elif test_method == "ancom-like-clr":
        metadata_df = pd.DataFrame(
            {
                "sample_id": count_df.index,
                condition_col: sample_to_cond.reindex(count_df.index).values,
            }
        )
        pair_col = None
        if sample_to_pair is not None:
            metadata_df["pair_id"] = sample_to_pair.reindex(count_df.index).values
            pair_col = "pair_id"
        res_df = _run_ancom_like_clr_test(
            count_df,
            metadata_df,
            condition_col=condition_col,
            sample_col="sample_id",
            pair_col=pair_col,
            pseudocount=composition_pseudocount,
        )
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

    res_df = _tag_proportion_result(res_df, test_method=test_method)

    # Multiple testing correction
    if "pval" in res_df.columns and multiple_testing_correction:
        try:
            from statsmodels.stats.multitest import multipletests

            valid = res_df["pval"].notna() & np.isfinite(res_df["pval"])
            res_df["padj"] = np.nan
            if valid.any():
                _, adjusted, _, _ = multipletests(
                    res_df.loc[valid, "pval"], method=multiple_testing_correction
                )
                res_df.loc[valid, "padj"] = adjusted
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
