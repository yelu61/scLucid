"""Tests for the analysis validation scaffold runners."""

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from validation.analysis.run_annotation_accuracy_benchmark import (
    _accuracy_rows,
    _confusion_rows,
    _major_lineage,
    _prepare,
)
from validation.analysis.run_proportion_consistency_benchmark import _extract_direction
from validation.analysis.run_pseudobulk_de_type1_error_benchmark import (
    _fdr_at_alpha,
    _generate_null_adata,
)


class TestAnalysisValidationScaffold:
    def test_major_lineage_maps_t_cells(self):
        assert _major_lineage("CD8 T cell") == "lymphoid"
        assert _major_lineage("macrophage") == "myeloid"
        assert _major_lineage("fibroblast") == "stromal"

    def test_accuracy_rows_compute_agreement(self):
        ref = pd.Series(["A", "A", "B", "B"])
        pred = pd.Series(["A", "A", "B", "C"])
        rows = {r["metric"]: r["value"] for r in _accuracy_rows(ref, pred)}
        assert rows["exact_label_accuracy"] == 0.75
        assert 0.0 <= rows["major_lineage_accuracy"] <= 1.0

    def test_confusion_rows_best_match(self):
        ref = pd.Series(["A", "A", "B", "B"])
        pred = pd.Series(["A", "A", "B", "C"])
        rows = _confusion_rows(ref, pred)
        by_ref = {r["reference_label"]: r for r in rows}
        assert by_ref["A"]["best_match_fraction"] == 1.0
        assert by_ref["B"]["best_match_fraction"] == 0.5

    def test_generate_null_adata_has_required_obs(self):
        adata = _generate_null_adata(n_genes=100, n_samples=4, seed=0)
        assert set(adata.obs.columns) >= {"sample", "condition", "cell_type"}
        assert "counts" in adata.layers

    def test_fdr_at_alpha_on_uniform_pvals(self):
        pvals = pd.Series(np.linspace(0.001, 0.999, 100))
        fdr = _fdr_at_alpha(pvals, alpha=0.05)
        assert abs(fdr - 0.05) < 0.02

    def test_extract_direction_handles_tuple_result(self):
        prop_df = pd.DataFrame({"sample": ["S1", "S2"], "A": [0.4, 0.6]})
        stat_df = pd.DataFrame(
            [{"cell_type": "A", "log2_fold_change": 0.5}, {"cell_type": "B", "log2_fold_change": -0.2}]
        )
        assert _extract_direction((prop_df, stat_df), "A") == 1
        assert _extract_direction((prop_df, stat_df), "B") == -1
