"""PBMC golden-path integration test."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pbmc_golden_path.py"
DATA_PATH = REPO_ROOT / "data" / "pbmc3k.h5ad"


def _load_script_module():
    spec = importlib.util.spec_from_file_location("run_pbmc_golden_path", SCRIPT_PATH)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.mark.slow
@pytest.mark.integration
def test_pbmc_golden_path_subset_outputs(tmp_path):
    """The PBMC golden path should produce the core acceptance artifacts."""
    module = _load_script_module()
    output_dir = tmp_path / "pbmc_golden"

    manifest = module.run_pbmc_golden_path(
        data_path=DATA_PATH,
        output_dir=output_dir,
        n_cells=300,
        n_top_genes=500,
        n_pcs=20,
        n_neighbors=10,
        random_state=42,
        overwrite=True,
        show_progress=False,
    )

    assert manifest["workflow"] == "pbmc3k_golden_path"
    assert manifest["input_shape"]["n_cells"] == 300
    assert manifest["final_shape"]["n_cells"] > 0
    assert manifest["final_shape"]["n_genes"] == 500
    assert manifest["retention_fraction"] > 0.5
    assert manifest["obs_summary"]["n_clusters"] >= 2
    assert manifest["obs_summary"]["n_cell_types"] is not None

    for stage in ["qc", "preprocess", "analysis"]:
        assert manifest["contracts"][stage]["valid"] is True

    assert (output_dir / "pbmc3k_golden_final.h5ad").exists()
    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "qc" / "qc_review_summary.json").exists()
    assert (output_dir / "preprocess" / "preprocess_review_summary.json").exists()
    assert (output_dir / "analysis" / "analysis_review_summary.json").exists()
    assert manifest["artifacts"]["figures"]
    for figure in manifest["artifacts"]["figures"]:
        assert Path(figure).exists()
