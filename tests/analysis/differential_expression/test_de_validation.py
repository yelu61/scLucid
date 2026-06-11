"""Tests for bulk-pseudobulk concordance."""

import numpy as np
import pandas as pd

from scLucid.analysis.differential_expression.de_validation import (
    compare_bulk_vs_pseudobulk_de,
)


def test_compare_bulk_vs_pseudobulk_de():
    genes = [f"G{i}" for i in range(20)]
    bulk = pd.DataFrame(
        {
            "gene": genes,
            "log2fc": np.concatenate([np.ones(10), -np.ones(10)]),
            "pvals_adj": np.concatenate([np.zeros(10), np.ones(10)]),
        }
    )
    pseudo = pd.DataFrame(
        {
            "gene": genes,
            "log2fc": np.concatenate([np.ones(10), -np.ones(10)]),
            "pvals_adj": np.concatenate([np.zeros(10), np.ones(10)]),
        }
    )
    result = compare_bulk_vs_pseudobulk_de(bulk, pseudo)
    assert result["log2fc_spearman"] == pytest.approx(1.0, abs=1e-6)
    assert result["directional_concordance"] == pytest.approx(1.0, abs=1e-6)
    assert result["jaccard_index"] == pytest.approx(1.0, abs=1e-6)


import pytest
