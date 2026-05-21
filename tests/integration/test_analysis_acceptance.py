"""Analysis acceptance runner integration test."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import scanpy as sc

REPO_ROOT = Path(__file__).parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_analysis_acceptance.py"


def _load_script_module():
    spec = importlib.util.spec_from_file_location("run_analysis_acceptance", SCRIPT_PATH)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _make_preprocessed_h5ad(path: Path) -> Path:
    import anndata

    rng = np.random.default_rng(11)
    genes = [
        "EPCAM",
        "KRT8",
        "KRT18",
        "MUC1",
        "MKI67",
        "TOP2A",
        "PTPRC",
        "CD3D",
        "CD3E",
        "TRAC",
        "COL1A1",
        "DCN",
        "LUM",
        "VWF",
        "PECAM1",
        "CDH5",
    ] + [f"gene_{i}" for i in range(64)]
    X = rng.poisson(2, size=(96, len(genes))).astype(np.float32)
    X[:32, :6] += 8
    X[32:64, 6:10] += 8
    X[64:, 10:16] += 8
    adata = anndata.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.obs["sampleID"] = "synthetic"
    adata.layers["counts"] = X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat")
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=20)
    sc.pp.neighbors(adata)
    adata.write(path)
    return path


@pytest.mark.integration
def test_analysis_acceptance_runner_writes_review_artifacts(tmp_path):
    module = _load_script_module()
    input_path = _make_preprocessed_h5ad(tmp_path / "preprocessed.h5ad")
    output_dir = tmp_path / "analysis_acceptance"

    manifest = module.run_analysis_acceptance(
        input_path=input_path,
        output_dir=output_dir,
        resolutions=(1.0,),
        run_malignancy=True,
        overwrite=True,
        show_progress=False,
        write_h5ad=False,
    )

    assert manifest["workflow"] == "analysis_acceptance"
    assert manifest["acceptance"]["ready_for_real_data_review"] is True
    assert manifest["acceptance"]["metrics"]["n_clusters"] >= 2
    assert manifest["acceptance"]["metrics"]["review_table_rows"] >= 2
    assert manifest["acceptance"]["metrics"]["n_final_labels"] >= 1
    assert manifest["acceptance"]["metrics"]["malignancy_enabled"] is True
    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "validation" / "analysis_acceptance.json").exists()
    assert (output_dir / "analysis" / "annotation_review_table.csv").exists()
    assert (output_dir / "analysis" / "llm_annotation_bundle.json").exists()
    assert (output_dir / "analysis" / "malignancy_interpretation_table.csv").exists()
