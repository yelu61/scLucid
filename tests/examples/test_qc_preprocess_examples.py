from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sparse
from anndata import AnnData, read_h5ad

import scLucid as scl

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXAMPLES = REPO_ROOT / "examples" / "02_simple_api"
ADVANCED_EXAMPLES = REPO_ROOT / "examples" / "03_advanced_notebooks"
USER_DOCS = REPO_ROOT / "docs" / "user"


def _load_example(filename: str):
    path = SIMPLE_EXAMPLES / filename
    spec = importlib.util.spec_from_file_location(f"sclucid_example_{path.stem}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _count_adata(*, sample_key: str) -> AnnData:
    rng = np.random.default_rng(23)
    matrix = sparse.csr_matrix(rng.poisson(1.2, size=(80, 240)).astype(np.int32))
    adata = AnnData(matrix)
    adata.obs_names = [f"cell_{idx}" for idx in range(adata.n_obs)]
    adata.var_names = [f"MT-G{idx}" if idx < 5 else f"G{idx}" for idx in range(adata.n_vars)]
    adata.obs[sample_key] = ["S1"] * 40 + ["S2"] * 40
    adata.layers["counts"] = adata.X.copy()
    return adata


def _notebook_source(filename: str) -> str:
    payload = json.loads((ADVANCED_EXAMPLES / filename).read_text())
    return "\n".join(
        line
        for cell in payload["cells"]
        for line in cell.get("source", [])
    )


def test_canonical_example_runs_public_review_apply_path(tmp_path):
    module = _load_example("qc_preprocess_review.py")
    input_path = tmp_path / "input.h5ad"
    _count_adata(sample_key="sample").write_h5ad(input_path)

    output_path = module.main(data_path=input_path, output_dir=tmp_path / "output")

    assert output_path.exists()
    result = read_h5ad(output_path)
    assert {"counts", "normalized_full"} <= set(result.layers)
    assert result.raw is not None and result.raw.n_vars == result.n_vars
    assert "discovery_feature" in result.var
    assert "X_pca" in result.obsm
    assert "policy_run_evidence" in result.uns["sclucid"]["preprocess"]


def test_manual_preprocess_example_keeps_full_expression_spaces(tmp_path):
    module = _load_example("preprocess_step_by_step.py")
    input_path = tmp_path / "qc_input.h5ad"
    _count_adata(sample_key="sampleID").write_h5ad(input_path)

    output_path = module.main(data_path=input_path, output_dir=tmp_path / "manual")

    result = read_h5ad(output_path)
    assert {"counts", "normalized_full"} <= set(result.layers)
    assert "scaled" not in result.layers
    assert "regressed" not in result.layers
    assert result.raw is not None and result.raw.n_vars == result.n_vars
    assert "discovery_feature" in result.var
    assert "X_pca" in result.obsm
    contract = result.uns["sclucid"]["preprocess"]["representation_contract"]
    assert contract["formal_count_model_source"] == "layers[counts]"


def test_context_examples_return_current_public_decision_types():
    adata = _count_adata(sample_key="sampleID")
    qc_module = _load_example("intelligent_qc.py")
    evaluation_module = _load_example("qc_evaluation.py")
    preprocess_module = _load_example("intelligent_preprocess.py")

    normal, tumor = qc_module.main(adata)
    evaluated_normal, evaluated_tumor = evaluation_module.main(adata)
    preprocess = preprocess_module.main(adata)

    assert all(
        isinstance(card, scl.DecisionCard)
        for card in (normal, tumor, evaluated_normal, evaluated_tumor, preprocess)
    )
    assert isinstance(normal.policy, scl.QCPolicy)
    assert isinstance(preprocess.policy, scl.PreprocessPolicy)


def test_product_examples_do_not_restore_scalar_quality_score():
    checked = list((REPO_ROOT / "examples").rglob("*.py"))
    checked.extend(ADVANCED_EXAMPLES.glob("*.ipynb"))
    for path in checked:
        assert "data_quality_score" not in path.read_text(), path
        assert "readiness_score" not in path.read_text(), path


def test_simple_api_examples_are_syntax_valid():
    for path in SIMPLE_EXAMPLES.glob("*.py"):
        compile(path.read_text(), str(path), "exec")


def test_preprocess_examples_preserve_four_space_semantics():
    manual_source = (SIMPLE_EXAMPLES / "preprocess_step_by_step.py").read_text()
    notebook_source = _notebook_source("Step1B-Preprocessing_Audit.ipynb")

    for source in (manual_source, notebook_source):
        assert "normalized_full" in source
        assert "discovery_feature" in source
        assert "formal_count_model_source" in source
        assert "temporary discovery" in source.lower()

    assert "subset=True" not in notebook_source
    assert "required_layers = {'counts', 'normalized', 'regressed', 'scaled'}" not in notebook_source


def test_user_docs_name_the_policy_api_and_four_space_contract():
    best_practices = (USER_DOCS / "best_practices.md").read_text()
    usage_layers = (USER_DOCS / "usage_layers.md").read_text()
    data_contracts = (USER_DOCS / "data_contracts.md").read_text()

    assert "recommend_qc_policy" in best_practices
    assert "recommend_preprocess_policy" in best_practices
    assert "compatibility" in usage_layers and "sensitivity APIs" in usage_layers
    assert 'adata.layers["normalized_full"]' in data_contracts
    assert 'adata.var["discovery_feature"]' in data_contracts
