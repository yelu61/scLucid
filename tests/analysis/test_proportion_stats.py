"""Tests for cell-type proportion statistics."""

import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.analysis.proportion.config import ProportionConfig
from scLucid.analysis.proportion.pseudobulk import celltype_proportion_analysis
from scLucid.analysis.proportion.stats import run_statistical_test


def test_proportion_config_accepts_kruskal():
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="kruskal",
    )

    assert cfg.test_method == "kruskal"


def test_run_statistical_test_kruskal():
    count_df = pd.DataFrame(
        {
            "T": [0.1, 0.2, 0.8, 0.9, 0.4, 0.5],
            "B": [0.6, 0.5, 0.2, 0.1, 0.3, 0.2],
        },
        index=["s1", "s2", "s3", "s4", "s5", "s6"],
    )
    sample_to_cond = pd.Series(
        ["A", "A", "B", "B", "C", "C"],
        index=count_df.index,
    )

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="kruskal",
        sample_to_cond=sample_to_cond,
    )

    assert set(result["cell_type"]) == {"T", "B"}
    assert {"statistic", "pval", "padj"}.issubset(result.columns)


def test_pseudobulk_kruskal_uses_sample_level_metadata():
    adata = AnnData(X=np.ones((6, 1)))
    adata.obs["sample"] = ["s1", "s1", "s2", "s2", "s3", "s3"]
    adata.obs["condition"] = ["A", "A", "B", "B", "C", "C"]
    adata.obs["cell_type"] = ["T", "B", "T", "T", "B", "B"]

    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="kruskal",
        auto_configure=True,
        plot_types=[],
    )

    prop_df, stat_df = celltype_proportion_analysis(adata, cfg)

    assert list(prop_df.index) == ["s1", "s2", "s3"]
    assert not stat_df.empty
    assert set(stat_df["cell_type"]) == {"T", "B"}
