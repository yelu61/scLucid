"""Tests for differential abundance of deconvolved proportions."""

import numpy as np
import pandas as pd

from scLucid.tools.bulk import run_bulk_abundance_test


def test_run_bulk_abundance_test_wilcoxon():
    np.random.seed(0)
    proportions = pd.DataFrame(
        {
            "T": np.concatenate([np.random.uniform(0.3, 0.5, 3), np.random.uniform(0.5, 0.7, 3)]),
            "B": np.concatenate([np.random.uniform(0.3, 0.5, 3), np.random.uniform(0.1, 0.3, 3)]),
        },
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    metadata = pd.DataFrame(
        {"group": ["A", "A", "A", "B", "B", "B"]},
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    result = run_bulk_abundance_test(
        proportions, metadata, group_col="group", group1="A", group2="B", method="wilcoxon"
    )
    assert not result.empty
    assert "pval" in result.columns
    assert "pvals_adj" in result.columns
    assert result["valid_for_publication_inference"].all()
    assert result["replicate_requirement_met"].all()


def test_run_bulk_abundance_test_insufficient_replicates():
    # Group A has only one sample -> not enough replicates for formal inference
    proportions = pd.DataFrame(
        {"T": [0.4, 0.5, 0.55, 0.45], "B": [0.5, 0.4, 0.45, 0.5]},
        index=["S1", "S2", "S3", "S4"],
    )
    metadata = pd.DataFrame(
        {"group": ["A", "B", "B", "B"]},
        index=["S1", "S2", "S3", "S4"],
    )
    result = run_bulk_abundance_test(
        proportions, metadata, group_col="group", group1="A", group2="B", method="ttest"
    )
    # With a single sample in group A the function skips all cell types.
    assert result.empty


def test_run_bulk_abundance_test_low_replicate_returns_descriptive():
    # Group A has two samples and group B has three; scipy can compute a t-test,
    # but the function marks replicate_requirement_met as False if the minimum
    # of two per group is not satisfied for both.
    proportions = pd.DataFrame(
        {"T": [0.4, 0.42, 0.5, 0.55, 0.45], "B": [0.5, 0.48, 0.4, 0.45, 0.5]},
        index=["S1", "S2", "S3", "S4", "S5"],
    )
    metadata = pd.DataFrame(
        {"group": ["A", "A", "B", "B", "B"]},
        index=["S1", "S2", "S3", "S4", "S5"],
    )
    result = run_bulk_abundance_test(
        proportions, metadata, group_col="group", group1="A", group2="B", method="ttest"
    )
    assert not result.empty
    assert result["replicate_requirement_met"].all()
