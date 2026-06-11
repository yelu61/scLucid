"""Tests for clinical association with deconvolved proportions."""

import numpy as np
import pandas as pd

from scLucid.tools.bulk import correlate_abundance_with_clinical


def test_correlate_abundance_with_clinical():
    np.random.seed(0)
    n = 10
    proportions = pd.DataFrame(
        {
            "T": np.random.uniform(0.1, 0.6, n),
            "B": np.random.uniform(0.1, 0.6, n),
        },
        index=[f"S{i}" for i in range(n)],
    )
    # Create a clinical variable that correlates with T cell proportion
    metadata = pd.DataFrame(
        {
            "survival_months": proportions["T"] * 10 + np.random.normal(0, 0.5, n),
        },
        index=[f"S{i}" for i in range(n)],
    )
    result = correlate_abundance_with_clinical(
        proportions, metadata, clinical_variable="survival_months", method="pearson"
    )
    assert not result.empty
    assert "correlation_coefficient" in result.columns
    assert "pval" in result.columns
    assert not result["valid_for_publication_inference"].all()
    assert result["inference_level"].iloc[0] == "exploratory_trait_association"
