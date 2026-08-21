from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from validation.dataset_registry import DATASETS
from validation.qc_preprocess.run_locked_qc_acceptance import run_locked_qc_acceptance
from validation.qc_preprocess.truth_pack import TruthDatasetSpec, build_truth_pack


def _write_fixture(path: Path) -> None:
    rng = np.random.default_rng(17)
    adata = AnnData(rng.poisson(1.0, size=(80, 120)).astype(np.int32))
    adata.obs_names = [f"cell_{idx}" for idx in range(80)]
    adata.var_names = [f"G{idx}" for idx in range(120)]
    adata.obs["sample"] = ["good"] * 40 + ["bad"] * 40
    adata.obs["n_genes_by_counts"] = [1000.0] * 40 + [45.0] * 40
    adata.obs["total_counts"] = [3000.0] * 40 + [80.0] * 40
    adata.obs["pct_counts_mt"] = [5.0] * 40 + [75.0] * 40
    adata.obs["pct_counts_in_top_20_genes"] = [30.0] * 40 + [97.0] * 40
    adata.write_h5ad(path)


def _build_pack(tmp_path: Path, name: str = "pack") -> Path:
    source = tmp_path / "fixture.h5ad"
    if not source.exists():
        _write_fixture(source)
    pack = tmp_path / name
    build_truth_pack(
        [
            TruthDatasetSpec(
                key="fixture",
                path=str(source),
                tissue="test",
                dataset_type="tumor_tissue",
                sample_key="sample",
            )
        ],
        pack,
        seed=31,
        primary_per_library=20,
        challenge_per_axis=5,
    )
    return pack


def _write_fixture_contract(tmp_path: Path) -> Path:
    source = Path(__file__).parents[1] / "validation/qc_preprocess/acceptance_contract.json"
    payload = json.loads(source.read_text())
    payload["development_datasets"] = ["fixture"]
    payload["external_projects"] = []
    output = tmp_path / "acceptance_contract.json"
    output.write_text(json.dumps(payload))
    return output


def test_truth_pack_is_deterministic_and_prediction_blinded(tmp_path):
    first = _build_pack(tmp_path, "pack_a")
    second = _build_pack(tmp_path, "pack_b")

    first_cells = (first / "reviewer" / "cell_evidence.tsv").read_text()
    second_cells = (second / "reviewer" / "cell_evidence.tsv").read_text()
    assert first_cells == second_cells
    assert "original_obs_name" not in first_cells
    assert "sampling_tier" not in first_cells
    assert "selector_call" not in first_cells
    assert "cell_" not in first_cells
    assert "original_obs_name" in (first / "sealed" / "cell_key.tsv").read_text()


def test_public_mixology_registry_keeps_qc_truth_scope_explicit():
    spec = next(item for item in DATASETS if item.key == "public_mixology")

    assert "mixology_identity" in spec.required_obs
    assert "identity" in spec.qc_roles[0]
    assert "not low-quality-cell truth" in spec.benchmark_notes


def test_incomplete_labels_block_acceptance_without_running_predictions(tmp_path):
    pack = _build_pack(tmp_path)
    contract = _write_fixture_contract(tmp_path)
    report = run_locked_qc_acceptance(
        pack,
        tmp_path / "report",
        contract_path=contract,
    )

    assert report["status"] == "BLOCKED"
    assert report["label_gate"]["status"] == "BLOCKED"
    assert report["source_integrity"]["status"] == "NOT_RUN"
    assert (tmp_path / "report" / "locked_qc_acceptance.json").exists()


def test_frozen_labels_unblind_and_reach_locked_gates(tmp_path):
    pack = _build_pack(tmp_path)
    reviewer = pack / "reviewer"
    sealed = pack / "sealed"

    sample_labels = pd.read_csv(reviewer / "sample_labels.tsv", sep="\t", dtype=str).fillna("")
    sample_key = pd.read_csv(sealed / "sample_key.tsv", sep="\t", dtype=str)
    sample_truth = sample_key.set_index("case_id")["original_sample"].map(
        {"good": "KEEP", "bad": "REMOVE"}
    )
    sample_labels["expert_label"] = sample_labels["case_id"].map(sample_truth)
    sample_labels["reviewer_id"] = "test-reviewer"
    sample_labels.to_csv(reviewer / "sample_labels.tsv", sep="\t", index=False)

    cell_labels = pd.read_csv(reviewer / "cell_labels.tsv", sep="\t", dtype=str).fillna("")
    cell_key = pd.read_csv(sealed / "cell_key.tsv", sep="\t", dtype=str)
    cell_truth = cell_key.set_index("case_id")["original_sample"].map(
        {"good": "KEEP", "bad": "REMOVE"}
    )
    cell_labels["expert_label"] = cell_labels["case_id"].map(cell_truth)
    cell_labels["reviewer_id"] = "test-reviewer"
    cell_labels.to_csv(reviewer / "cell_labels.tsv", sep="\t", index=False)

    contract = _write_fixture_contract(tmp_path)
    report = run_locked_qc_acceptance(
        pack,
        tmp_path / "report",
        contract_path=contract,
    )

    assert report["status"] in {"PASS", "FAIL"}
    assert report["label_gate"]["status"] == "PASS"
    assert report["source_integrity"]["status"] == "PASS"
    assert report["cell_endpoint"]["status"] in {"PASS", "FAIL"}
    assert report["sample_endpoint"]["status"] in {"PASS", "FAIL"}
    stored = json.loads((tmp_path / "report" / "locked_qc_acceptance.json").read_text())
    assert stored["status"] == report["status"]
