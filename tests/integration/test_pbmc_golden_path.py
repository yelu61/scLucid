"""PBMC golden-path integration test."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

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


@pytest.mark.parametrize(
    "script_path",
    [SCRIPT_PATH, REPO_ROOT / "scripts" / "run_pdac_golden_path.py"],
)
def test_golden_path_hdf5_sanitizer_encodes_empty_metadata_keys(tmp_path, script_path):
    """Empty biological labels must not make compact audit metadata unwritable."""
    spec = importlib.util.spec_from_file_location(f"sanitize_{script_path.stem}", script_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    safe = module.make_hdf5_safe(
        {"retention": {"per_cell_type": {"": {"initial_cells": 3}, "B cells": {}}}}
    )
    assert "__empty_label__" in safe["retention"]["per_cell_type"]
    adata = AnnData(X=np.ones((1, 1), dtype=float))
    adata.uns["sclucid"] = safe
    adata.write_h5ad(tmp_path / f"{script_path.stem}.h5ad")


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

    assert manifest["validation"]["ready_for_comparative_validation"] is True
    assert "does not claim scientific superiority" in manifest["validation"]["claim_boundary"]
    assert (output_dir / "pbmc3k_golden_final.h5ad").exists()
    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "validation" / "qc_preprocess_validation.json").exists()
    assert (output_dir / "validation" / "qc_preprocess_validation_table.csv").exists()
    assert (output_dir / "qc" / "qc_review_summary.json").exists()
    assert (output_dir / "preprocess" / "preprocess_review_summary.json").exists()
    assert (output_dir / "analysis" / "analysis_review_summary.json").exists()
    assert manifest["artifacts"]["figures"]
    for figure in manifest["artifacts"]["figures"]:
        assert Path(figure).exists()
