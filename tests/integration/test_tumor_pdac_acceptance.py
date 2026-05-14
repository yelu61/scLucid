"""Tumor real-data acceptance test against lin2020 PDAC.

Runs the full scLucid tumor workflow (QC → preprocess → analysis → tumor
stage) on a subsampled real PDAC dataset and asserts biologically plausible
outcomes. This is the most demanding integration test in the suite —
it exercises the same code path used by the maintained PDAC golden path
script but verifies the *result*, not just the absence of exceptions.

The test is marked ``slow`` + ``integration`` and is skipped if the PDAC
.h5ad file is not present (e.g. on CI where Git LFS is not configured).

Expected runtime: 2-4 minutes on a typical workstation.

Biological assertions are deliberately conservative — they verify the
workflow produces sensible distributions rather than precise numbers, which
would couple the test to the random seed and library versions.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pdac_golden_path.py"
DATA_PATH = REPO_ROOT / "data" / "lin2020.pdac.h5ad"


def _load_script_module():
    spec = importlib.util.spec_from_file_location("run_pdac_golden_path", SCRIPT_PATH)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def pdac_manifest(tmp_path_factory):
    """Run the PDAC golden path once per module; share the manifest."""
    if not DATA_PATH.exists():
        pytest.skip(
            f"PDAC data file not available at {DATA_PATH}; "
            "skipping tumor real-data acceptance test."
        )

    module = _load_script_module()
    output_dir = tmp_path_factory.mktemp("pdac_acceptance")

    manifest = module.run_pdac_golden_path(
        data_path=DATA_PATH,
        output_dir=output_dir,
        n_cells=500,
        n_top_genes=2000,
        n_pcs=20,
        n_neighbors=15,
        random_state=42,
        include_annotation=True,
        overwrite=True,
        show_progress=False,
    )
    return manifest, output_dir


@pytest.mark.slow
@pytest.mark.integration
class TestLin2020PDACAcceptance:
    """Acceptance tests against the lin2020 PDAC dataset (subsampled to 500 cells).

    These tests verify that the documented tumor workflow produces
    biologically plausible outputs, not just that it does not crash.
    Each test reads from the shared module-scoped manifest fixture so the
    full pipeline only runs once per pytest invocation.
    """

    def test_pipeline_runs_to_completion(self, pdac_manifest):
        """The PDAC golden path must produce a manifest of the expected shape."""
        manifest, _ = pdac_manifest
        assert manifest["workflow"] == "pdac_golden_path"
        assert manifest["input_shape"]["n_cells"] == 500
        assert manifest["final_shape"]["n_cells"] > 0
        # HVG selection can drop a small number below the n_top_genes target
        # when low-variance genes are filtered, so the gene count is an upper
        # bound rather than an equality.
        assert 1500 <= manifest["final_shape"]["n_genes"] <= 2000

    def test_qc_retains_substantial_fraction(self, pdac_manifest):
        """QC must not drop more than 70% of cells on this PDAC subset.

        PDAC tissue typically has higher dropout / mitochondrial fraction
        than blood; nevertheless retaining < 30% of input would indicate the
        QC thresholds are too aggressive for tumor data.
        """
        manifest, _ = pdac_manifest
        assert manifest["retention_fraction"] >= 0.3, (
            f"QC retained only {manifest['retention_fraction']:.1%} of cells; "
            "tumor-aware QC may be too aggressive."
        )

    def test_all_workflow_contracts_validate(self, pdac_manifest):
        """Every scLucid stage contract must pass on the output AnnData."""
        manifest, _ = pdac_manifest
        contracts = manifest["contracts"]
        for stage in ("qc", "preprocess", "analysis"):
            assert stage in contracts, f"Missing contract result for stage {stage}"
            assert contracts[stage]["valid"] is True, (
                f"{stage} output contract is invalid: {contracts[stage]}"
            )

    def test_clustering_produces_heterogeneity(self, pdac_manifest):
        """A PDAC sample should cluster into at least 3 transcriptomically
        distinct groups (epithelial / stromal / immune at minimum)."""
        manifest, _ = pdac_manifest
        n_clusters = manifest["obs_summary"]["n_clusters"]
        assert n_clusters is not None and n_clusters >= 3, (
            f"Only {n_clusters} clusters detected; PDAC subset should have "
            "epithelial / stromal / immune at minimum."
        )

    def test_annotation_produces_multiple_cell_types(self, pdac_manifest):
        """Cell-type annotation must call at least 2 distinct types."""
        manifest, _ = pdac_manifest
        n_cell_types = manifest["obs_summary"]["n_cell_types"]
        if n_cell_types is None:
            pytest.skip("Annotation did not run (no marker DB or skipped step).")
        assert n_cell_types >= 2, (
            f"Annotation produced only {n_cell_types} cell type(s); "
            "expected multiple in a heterogeneous tumor sample."
        )

    def test_tumor_stage_runs_at_least_one_step(self, pdac_manifest):
        """The tumor stage should complete at least one of malignancy/CNV/TME."""
        manifest, _ = pdac_manifest
        tumor_steps = manifest["tumor"]["steps_executed"]
        assert tumor_steps, (
            f"Tumor stage executed no steps. Warnings: {manifest['tumor']['warnings']}"
        )

    def test_malignancy_proportion_in_plausible_range(self, pdac_manifest):
        """If is_malignant is computed, its fraction must be in (0%, 95%)."""
        manifest, _ = pdac_manifest
        n_malignant = manifest["obs_summary"].get("n_malignant")
        n_cells = manifest["final_shape"]["n_cells"]
        if n_malignant is None:
            pytest.skip(
                "Malignancy classification did not produce is_malignant column."
            )
        fraction = n_malignant / max(n_cells, 1)
        # PDAC samples vary widely — sometimes mostly stromal, sometimes
        # mostly malignant ductal cells. The hard constraint is just that
        # the classifier did not collapse to all-malignant or none-malignant.
        assert 0.0 < fraction < 0.95, (
            f"Malignant fraction {fraction:.1%} is implausible "
            f"(n_malignant={n_malignant}, n_cells={n_cells})."
        )

    def test_review_summary_artifacts_written(self, pdac_manifest):
        """Each stage's review_summary should be persisted to disk."""
        _, output_dir = pdac_manifest
        for stage in ("qc", "preprocess", "analysis"):
            json_path = output_dir / stage / f"{stage}_review_summary.json"
            assert json_path.exists(), (
                f"{stage} review summary JSON missing at {json_path}; "
                "auditability promise broken."
            )
        assert (output_dir / "validation" / "qc_preprocess_validation.json").exists()
        assert (output_dir / "validation" / "qc_preprocess_validation_table.csv").exists()

    def test_qc_preprocess_validation_ready_for_comparison(self, pdac_manifest):
        """The lightweight scaffold should mark QC/preprocess as comparison-ready."""
        manifest, _ = pdac_manifest
        validation = manifest["validation"]
        assert validation["ready_for_comparative_validation"] is True
        assert validation["readiness_status"] == "ready_for_comparative_validation"
        assert "does not claim scientific superiority" in validation["claim_boundary"]

    def test_final_h5ad_round_trips(self, pdac_manifest):
        """The saved .h5ad must reload without errors."""
        manifest, _ = pdac_manifest
        import anndata as ad

        path = manifest["artifacts"]["final_h5ad"]
        adata = ad.read_h5ad(path)
        assert adata.n_obs == manifest["final_shape"]["n_cells"]
        assert "sclucid" in adata.uns, (
            "Saved AnnData lost the scLucid namespace; audit report would fail."
        )

    def test_export_audit_report_renders(self, pdac_manifest, tmp_path):
        """The audit report exporter must succeed on real PDAC output."""
        import scLucid as scl
        import anndata as ad

        manifest, _ = pdac_manifest
        adata = ad.read_h5ad(manifest["artifacts"]["final_h5ad"])
        report_path = tmp_path / "pdac_audit_report.html"
        result = scl.export_audit_report(
            adata, report_path, title="PDAC Acceptance Audit"
        )
        assert result.exists()
        html = result.read_text()
        assert "Qc" in html or "QC" in html
        assert "Preprocess" in html
        assert "Analysis" in html
        # The title we passed through must be honored.
        assert "PDAC Acceptance Audit" in html
