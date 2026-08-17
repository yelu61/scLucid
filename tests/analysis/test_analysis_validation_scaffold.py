"""Tests for the analysis validation scaffold runners."""

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from validation.analysis.run_annotation_accuracy_benchmark import (
    _accuracy_rows,
    _confusion_rows,
    _major_lineage,
)
from validation.analysis.run_inference_contract_benchmark import run as run_inference_contract
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

    def test_real_inference_contract_runner_records_ready_and_blocked_designs(
        self, tmp_path
    ):
        rng = np.random.default_rng(17)
        pbmc_rows = []
        pbmc_counts = []
        genes = ["ISG15", "IFIT1", "MX1", "OAS1", "STAT1", "ACTB"]
        for donor_index, donor in enumerate(["d1", "d2", "d3"]):
            for condition in ["ctrl", "stim"]:
                t_cells = 4 + donor_index + (3 if condition == "stim" else 0)
                b_cells = 12 - t_cells
                for cell_type, n_cells in (("CD4 T cells", t_cells), ("B cells", b_cells)):
                    for _ in range(n_cells):
                        counts = rng.poisson(4, size=len(genes)).astype(float)
                        if condition == "stim":
                            counts[:5] += rng.poisson(5, size=5)
                        pbmc_counts.append(counts)
                        pbmc_rows.append(
                            {
                                "sample": f"capture_{condition}",
                                "donor": donor,
                                "condition": condition,
                                "cell_type": cell_type,
                                "batch_group": "paired_capture",
                            }
                        )
        pbmc = AnnData(
            X=np.asarray(pbmc_counts),
            obs=pd.DataFrame(pbmc_rows),
            var=pd.DataFrame(index=genes),
        )
        pbmc.layers["counts"] = pbmc.X.copy()
        pbmc.write_h5ad(tmp_path / "kang2018.pbmc.h5ad")

        pdac = AnnData(
            X=np.ones((8, 3)),
            obs=pd.DataFrame(
                {
                    "sampleID": ["tumor_1"] * 4 + ["tumor_2"] * 4,
                    "condition": ["Primary tumor"] * 8,
                    "cell_type": [""] * 8,
                }
            ),
            var=pd.DataFrame(index=["G1", "G2", "G3"]),
        )
        pdac.write_h5ad(tmp_path / "lin2020.pdac.h5ad")

        output_dir = tmp_path / "evidence"
        manifest = run_inference_contract(
            tmp_path,
            output_dir,
            max_genes=6,
            max_cell_types=2,
            min_cells_per_sample=3,
        )

        assert manifest["gate_status"] == "PASS"
        assert manifest["datasets"]["kang2018.pbmc"]["design"]["status"] == "READY"
        assert manifest["datasets"]["lin2020.pdac"]["design"]["status"] == "BLOCKED"
        assert manifest["datasets"]["lin2020.pdac"]["blocked_safely"] is True
        assert set(manifest["artifacts"]) == {
            "metadata_propagation_matrix",
            "real_data_design_audit",
            "pbmc_proportion_estimates",
            "pbmc_proportion_statistics",
            "pbmc_pseudobulk_de",
            "manifest",
        }
        assert all(Path(path).exists() for path in manifest["artifacts"].values())
        de_table = pd.read_csv(output_dir / "pbmc_pseudobulk_de.tsv", sep="\t")
        assert set(de_table["n_experimental_units_condition1"]) == {3}
        assert set(de_table["n_experimental_units_condition2"]) == {3}
