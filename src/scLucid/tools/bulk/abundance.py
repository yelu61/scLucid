"""Differential abundance and clinical association for deconvolved bulk proportions."""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests

from .config import BulkAbundanceConfig, BulkClinicalAssociationConfig


def _adjust_pvalues(pvals: pd.Series, method: str = "fdr_bh") -> pd.Series:
    valid = pvals.notna() & np.isfinite(pvals)
    adjusted = pd.Series(np.nan, index=pvals.index)
    if valid.any():
        _, qvals, _, _ = multipletests(pvals[valid].values, method=method)
        adjusted.loc[valid] = qvals
    return adjusted


def run_bulk_abundance_test(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    config: Optional[BulkAbundanceConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """Test for differential abundance of deconvolved cell-type proportions.

    Parameters
    ----------
    proportions_df
        Samples x cell types DataFrame.
    metadata_df
        Sample metadata indexed by sample ID.
    config
        Abundance test configuration.
    **kwargs
        Config overrides.

    Returns
    -------
    pd.DataFrame
        Per-cell-type statistics with inference tags.
    """
    if config is None:
        config = BulkAbundanceConfig(**kwargs)
    else:
        config = config.model_copy(update=kwargs)

    data = proportions_df.join(metadata_df, how="inner")
    group1_samples = data[data[config.group_col] == config.group1].index
    group2_samples = data[data[config.group_col] == config.group2].index

    n1 = len(group1_samples)
    n2 = len(group2_samples)
    replicate_requirement_met = n1 >= 2 and n2 >= 2

    results = []
    for cell_type in proportions_df.columns:
        scores1 = data.loc[group1_samples, cell_type].dropna()
        scores2 = data.loc[group2_samples, cell_type].dropna()

        if len(scores1) < 2 or len(scores2) < 2:
            continue

        if config.method == "wilcoxon":
            stat, pval = mannwhitneyu(scores1, scores2, alternative="two-sided")
        else:
            stat, pval = ttest_ind(scores1, scores2)

        mean1 = float(scores1.mean())
        mean2 = float(scores2.mean())
        results.append(
            {
                "cell_type": cell_type,
                "statistic": stat,
                "pval": pval,
                "mean_abundance_group1": mean1,
                "mean_abundance_group2": mean2,
                "log2fc_abundance": np.log2(mean2 / (mean1 + 1e-10)),
                "n_group1": len(scores1),
                "n_group2": len(scores2),
                "method": config.method,
            }
        )

    result = pd.DataFrame(results)
    if result.empty:
        return result

    result["pvals_adj"] = _adjust_pvalues(result["pval"], method=config.p_adjust_method)
    result["inference_level"] = "sample_level" if replicate_requirement_met else "descriptive_sample_level"
    result["valid_for_publication_inference"] = replicate_requirement_met
    result["replicate_requirement_met"] = replicate_requirement_met
    result["result_warning"] = (
        None if replicate_requirement_met else "Insufficient biological replicates for formal abundance inference."
    )
    return result.sort_values("pval")


def correlate_abundance_with_clinical(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    config: Optional[BulkClinicalAssociationConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """Correlate inferred cell-type proportions with a continuous clinical variable.

    Parameters
    ----------
    proportions_df
        Samples x cell types DataFrame.
    metadata_df
        Sample metadata indexed by sample ID.
    config
        Clinical association configuration.
    **kwargs
        Config overrides.

    Returns
    -------
    pd.DataFrame
        Per-cell-type correlation results with inference tags.
    """
    if config is None:
        config = BulkClinicalAssociationConfig(**kwargs)
    else:
        config = config.model_copy(update=kwargs)

    data = proportions_df.join(metadata_df, how="inner")
    results = []

    for cell_type in proportions_df.columns:
        subset = data[[cell_type, config.clinical_variable]].dropna()
        if len(subset) < config.min_samples:
            continue

        if config.method == "pearson":
            corr, pval = pearsonr(subset[cell_type], subset[config.clinical_variable])
        else:
            corr, pval = spearmanr(subset[cell_type], subset[config.clinical_variable])

        results.append(
            {
                "cell_type": cell_type,
                "clinical_variable": config.clinical_variable,
                "correlation_coefficient": corr,
                "pval": pval,
                "n_samples": len(subset),
                "method": config.method,
            }
        )

    result = pd.DataFrame(results)
    if result.empty:
        return result

    result["pvals_adj"] = _adjust_pvalues(result["pval"])
    result["inference_level"] = "exploratory_trait_association"
    result["valid_for_publication_inference"] = False
    result["result_warning"] = "Observational correlation; not causal inference."
    return result.sort_values("pval")
