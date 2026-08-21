from __future__ import annotations

import hashlib
import json
from pathlib import Path

from validation.acquisition.build_qc_truth_acquisition_plan import (
    TARGET_DATASET_IDS,
    build_plan,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO_ROOT / "validation" / "dataset_evidence_registry.json"
CONTRACT_PATH = REPO_ROOT / "validation" / "qc_preprocess" / "acceptance_contract.json"


def _dataset(plan: dict, dataset_id: str) -> dict:
    return next(row for row in plan["datasets"] if row["dataset_id"] == dataset_id)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_plan_is_registry_backed_and_keeps_acquisition_separate_from_science(tmp_path):
    registry = json.loads(REGISTRY_PATH.read_text())
    registry_by_id = {row["dataset_id"]: row for row in registry["datasets"]}

    plan = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=tmp_path / "qc_truth",
    )

    assert tuple(row["dataset_id"] for row in plan["datasets"]) == TARGET_DATASET_IDS
    assert plan["summary"]["scientifically_passed"] == 0
    assert plan["claim_boundary"].startswith("Acquisition readiness is not scientific PASS")
    for row in plan["datasets"]:
        source = registry_by_id[row["dataset_id"]]
        assert row["accessions"] == source["accessions"]
        assert row["license"] == source["license"]
        assert row["redistribution"] == source["license"]["redistribution"]
        assert row["raw_availability"] == source["download"]
        assert row["required_metadata"] == source["required_metadata"]
        assert set(row["qc_endpoint_ids"]) <= set(source["endpoint_ids"])
        assert row["overall_status"] == "NOT_READY"
        assert row["scientific_validation_status"] == "NOT_RUN"
        assert row["next_action"]["action_id"] == "RUN_LOCAL_PREFLIGHT"
        assert len(row["next_action"]) == 2


def test_check_local_is_read_only_and_missing_files_fail_closed(tmp_path):
    target_root = tmp_path / "qc_truth"
    assert not target_root.exists()

    plan = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=target_root,
        check_local=True,
    )

    assert not target_root.exists()
    assert plan["summary"]["blocked"] == len(TARGET_DATASET_IDS)
    for row in plan["datasets"]:
        assert row["source_urls"]
        assert row["overall_status"] == "BLOCKED"
        assert row["preparation_status"] == "NOT_READY"
        assert row["download_status"].startswith("BLOCKED_")
        assert row["next_action"]["action_id"].startswith("ACQUIRE_")


def test_verified_local_bundle_is_ready_for_validation_but_not_scientific_pass(tmp_path):
    dataset_id = "cell_hashing_gse108313"
    target_root = tmp_path / "qc_truth"
    dataset_dir = target_root / dataset_id

    unchecked = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=target_root,
    )
    template = _dataset(unchecked, dataset_id)
    for artifact in template["artifacts"]:
        path = Path(artifact["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"test payload for {artifact['role']}\n")

    checksum_lines = [
        f"{_sha256(Path(artifact['path']))}  {artifact['checksum_key']}"
        for artifact in template["artifacts"]
    ]
    (dataset_dir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n")
    (dataset_dir / "provenance.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset_id,
                "accessions": ["GSE108313"],
                "source_urls": ["https://www.ncbi.nlm.nih.gov/geo/"],
                "retrieved_at": "2026-08-20T00:00:00Z",
                "files": [artifact["checksum_key"] for artifact in template["artifacts"]],
            }
        )
    )
    (dataset_dir / "metadata_manifest.json").write_text(
        json.dumps({"available_fields": template["required_metadata"]})
    )

    checked = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=target_root,
        check_local=True,
    )
    row = _dataset(checked, dataset_id)

    assert row["overall_status"] == "READY_FOR_VALIDATION"
    assert row["download_status"] == "LOCAL_ARTIFACTS_VERIFIED"
    assert row["local_preflight"]["checksum"]["status"] == "VERIFIED"
    assert row["local_preflight"]["provenance"]["status"] == "VERIFIED"
    assert row["local_preflight"]["metadata"]["status"] == "VERIFIED"
    assert row["scientific_validation_status"] == "NOT_RUN"
    assert row["next_action"]["action_id"] == "RUN_REGISTERED_QC_ENDPOINTS"


def test_bad_checksum_blocks_preparation_without_inventing_replacement_hash(tmp_path):
    dataset_id = "tenx_hgmm_6k"
    target_root = tmp_path / "qc_truth"
    dataset_dir = target_root / dataset_id
    unchecked = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=target_root,
    )
    template = _dataset(unchecked, dataset_id)
    for artifact in template["artifacts"]:
        path = Path(artifact["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(artifact["role"])

    checksum_lines = [
        f"{'1' * 64}  {artifact['checksum_key']}" for artifact in template["artifacts"]
    ]
    (dataset_dir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n")
    (dataset_dir / "provenance.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset_id,
                "accessions": ["6k_HGMM_3p_v2"],
                "source_urls": ["https://www.10xgenomics.com/"],
                "retrieved_at": "2026-08-20T00:00:00Z",
                "files": [artifact["checksum_key"] for artifact in template["artifacts"]],
            }
        )
    )
    (dataset_dir / "metadata_manifest.json").write_text(
        json.dumps({"available_fields": template["required_metadata"]})
    )

    checked = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=target_root,
        check_local=True,
    )
    row = _dataset(checked, dataset_id)

    assert row["overall_status"] == "BLOCKED"
    assert row["local_preflight"]["checksum"]["status"] == "INVALID"
    assert "checksum mismatch" in row["local_preflight"]["checksum"]["error"]
    assert row["next_action"]["action_id"] == "VERIFY_CHECKSUMS"
    assert "1" * 64 not in json.dumps(row["next_action"])


def test_unfiltered_and_label_requirements_are_endpoint_specific(tmp_path):
    plan = build_plan(
        REGISTRY_PATH,
        CONTRACT_PATH,
        repo_root=tmp_path,
        target_root=tmp_path / "qc_truth",
    )

    cell_hashing = _dataset(plan, "cell_hashing_gse108313")
    assert not cell_hashing["requirements"]["unfiltered_droplets"]["required"]
    assert (
        "not the complete unfiltered RNA droplet universe"
        in cell_hashing["requirements"]["unfiltered_droplets"]["reason"]
    )
    assert _dataset(plan, "tenx_hgmm_6k")["requirements"]["unfiltered_droplets"]["required"]
    assert not _dataset(plan, "kang2018_pbmc")["requirements"]["unfiltered_droplets"]["required"]
    assert (
        _dataset(plan, "emtab2600_microscopy_quality")["requirements"]["unfiltered_droplets"][
            "required"
        ]
        is None
    )
    assert not _dataset(plan, "tenx_hgmm_6k")["requirements"]["label_recovery"]["required"]
