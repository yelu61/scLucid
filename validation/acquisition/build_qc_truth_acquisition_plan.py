"""Build a fail-closed acquisition plan for public QC truth datasets.

This module deliberately stops at acquisition readiness.  A reachable source
URL or a locally present file is not scientific validation evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCHEMA_NAME = "scLucidQCTruthAcquisitionPlan"
SCHEMA_VERSION = "1.0.0"
TARGET_DATASET_IDS = (
    "emtab2600_microscopy_quality",
    "cell_hashing_gse108313",
    "tenx_hgmm_6k",
    "kang2018_pbmc",
)
_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


@dataclass(frozen=True)
class ArtifactRequirement:
    role: str
    relative_path: str | None = None
    registry_local_path: bool = False
    description: str = ""


@dataclass(frozen=True)
class DatasetAcquisitionRequirement:
    unfiltered_droplets_required: bool | None
    unfiltered_droplets_reason: str
    label_recovery_required: bool
    label_recovery_strategy: str
    artifacts: tuple[ArtifactRequirement, ...]


ACQUISITION_REQUIREMENTS: dict[str, DatasetAcquisitionRequirement] = {
    "emtab2600_microscopy_quality": DatasetAcquisitionRequirement(
        unfiltered_droplets_required=None,
        unfiltered_droplets_reason=(
            "Not a modern droplet experiment; reconstruct the complete set of "
            "capture sites, including microscopy-labelled empty sites."
        ),
        label_recovery_required=True,
        label_recovery_strategy=(
            "Recover microscopy labels and the capture-site to sequenced-library "
            "mapping from the source publication/supplement; preserve ambiguous "
            "sites as uncertain."
        ),
        artifacts=(
            ArtifactRequirement(
                "ena_file_manifest",
                "source/ena_file_manifest.tsv",
                description="ENA run/file manifest used to reconstruct source reads.",
            ),
            ArtifactRequirement(
                "microscopy_quality_labels",
                "labels/microscopy_quality_labels.tsv",
                description="Intact, damaged, empty, and uncertain capture-site labels.",
            ),
            ArtifactRequirement(
                "capture_site_barcode_mapping",
                "labels/capture_site_barcode_mapping.tsv",
                description="Mapping between microscopy capture sites and count columns.",
            ),
            ArtifactRequirement(
                "prepared_truth_object",
                "prepared/emtab2600_qc_truth.h5ad",
                description="Prepared count object with independent microscopy truth.",
            ),
        ),
    ),
    "cell_hashing_gse108313": DatasetAcquisitionRequirement(
        unfiltered_droplets_required=False,
        unfiltered_droplets_reason=(
            "GSE108313 provides candidate RNA barcodes and an unfiltered HTO table, "
            "but not the complete unfiltered RNA droplet universe. Use it for HTO "
            "doublet/sample-hash calibration, not complete RNA cell calling."
        ),
        label_recovery_required=True,
        label_recovery_strategy=(
            "Recover HTO classifications, sample hashes, RNA barcodes, and explicit "
            "negative/unknown classes without treating HTO failure as a singlet."
        ),
        artifacts=(
            ArtifactRequirement(
                "rna_candidate_matrix",
                "source/GSM3501446_MixCellLines-RNA.umi.txt.gz",
                description="GEO RNA UMI matrix for the candidate barcode set.",
            ),
            ArtifactRequirement(
                "unfiltered_hto_matrix",
                "source/GSM3501447_MixCellLines-HTO-counts.csv.gz",
                description="GEO HTO count table retaining negative/background events.",
            ),
            ArtifactRequirement(
                "hto_assignments",
                "labels/hto_assignments.tsv",
                description="Orthogonal HTO singlet, multiplet, negative, and unknown labels.",
            ),
            ArtifactRequirement(
                "rna_hto_barcode_mapping",
                "labels/rna_hto_barcode_mapping.tsv",
                description="Auditable mapping between RNA and HTO barcode identities.",
            ),
            ArtifactRequirement(
                "prepared_truth_object",
                "prepared/cell_hashing_qc_truth.h5ad",
                description="Prepared count object with HTO evidence retained.",
            ),
        ),
    ),
    "tenx_hgmm_6k": DatasetAcquisitionRequirement(
        unfiltered_droplets_required=True,
        unfiltered_droplets_reason=(
            "Cell-calling and ambient-contamination endpoints require droplets below "
            "the vendor filtered-cell boundary."
        ),
        label_recovery_required=False,
        label_recovery_strategy=(
            "Derive species identity and mixed-species evidence from human and mouse "
            "UMI counts using a versioned preparation rule; this is derived truth, "
            "not a recovered manual label."
        ),
        artifacts=(
            ArtifactRequirement(
                "unfiltered_rna_matrix",
                "source/hgmm_6k_raw_gene_bc_matrices.tar.gz",
                description="Vendor unfiltered human-mouse count matrix.",
            ),
            ArtifactRequirement(
                "filtered_cell_matrix",
                "source/hgmm_6k_filtered_gene_bc_matrices.tar.gz",
                description="Vendor filtered matrix retained as a comparator, not truth.",
            ),
            ArtifactRequirement(
                "species_count_table",
                "labels/species_umi_counts.tsv",
                description="Per-barcode human and mouse UMI counts with derivation metadata.",
            ),
            ArtifactRequirement(
                "prepared_truth_object",
                "prepared/tenx_hgmm_qc_truth.h5ad",
                description="Prepared count object with species-derived evidence.",
            ),
        ),
    ),
    "kang2018_pbmc": DatasetAcquisitionRequirement(
        unfiltered_droplets_required=False,
        unfiltered_droplets_reason=(
            "The registered Kang QC endpoint is doublet calibration against demuxlet; "
            "cell calling and ambient correction are not assigned to this dataset."
        ),
        label_recovery_required=True,
        label_recovery_strategy=(
            "Recover and verify demuxlet singlet, doublet, and ambiguous calls together "
            "with donor, condition, and capture identifiers."
        ),
        artifacts=(
            ArtifactRequirement(
                "demuxlet_label_table",
                "labels/demuxlet_labels.tsv",
                description="Barcode-level demuxlet calls, including ambiguous calls.",
            ),
            ArtifactRequirement(
                "prepared_truth_object",
                registry_local_path=True,
                description="Registry-declared prepared Kang object.",
            ),
        ),
    ),
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else root / path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_checksum_manifest(path: Path) -> tuple[dict[str, str], str | None]:
    entries: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        return {}, f"checksum manifest is unreadable: {exc}"

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not _SHA256_RE.fullmatch(parts[0]):
            return {}, f"invalid checksum line {line_number}"
        relative_path = parts[1].lstrip("* ")
        if (
            not relative_path
            or Path(relative_path).is_absolute()
            or ".." in Path(relative_path).parts
        ):
            return {}, f"unsafe checksum path on line {line_number}"
        entries[relative_path] = parts[0].lower()
    if not entries:
        return {}, "checksum manifest contains no entries"
    return entries, None


def _check_checksums(*, manifest_path: Path, artifacts: list[dict[str, Any]]) -> dict[str, Any]:
    if not manifest_path.is_file():
        return {
            "status": "MISSING",
            "path": str(manifest_path),
            "verified_files": [],
            "error": "checksum manifest is required and was not found",
        }

    entries, error = _read_checksum_manifest(manifest_path)
    if error:
        return {
            "status": "INVALID",
            "path": str(manifest_path),
            "verified_files": [],
            "error": error,
        }

    verified: list[str] = []
    errors: list[str] = []
    for artifact in artifacts:
        artifact_path = Path(artifact["path"])
        if not artifact_path.is_file():
            continue
        relative_path = artifact["checksum_key"]
        expected = entries.get(relative_path)
        if expected is None:
            errors.append(f"missing checksum entry: {relative_path}")
            continue
        if _sha256(artifact_path) != expected:
            errors.append(f"checksum mismatch: {relative_path}")
            continue
        verified.append(relative_path)

    return {
        "status": "VERIFIED" if not errors else "INVALID",
        "path": str(manifest_path),
        "verified_files": verified,
        "error": "; ".join(errors) if errors else None,
    }


def _check_provenance(path: Path, dataset_id: str) -> dict[str, Any]:
    if not path.is_file():
        return {
            "status": "MISSING",
            "path": str(path),
            "error": "provenance manifest is required and was not found",
        }
    try:
        provenance = _load_json(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {"status": "INVALID", "path": str(path), "error": str(exc)}

    required = {"dataset_id", "accessions", "source_urls", "retrieved_at", "files"}
    missing = sorted(required - set(provenance))
    if provenance.get("dataset_id") != dataset_id:
        missing.append("dataset_id_match")
    empty = sorted(
        key for key in required & set(provenance) if provenance.get(key) in (None, "", [])
    )
    problems = sorted(set(missing + empty))
    return {
        "status": "VERIFIED" if not problems else "INVALID",
        "path": str(path),
        "error": f"missing or invalid fields: {', '.join(problems)}" if problems else None,
    }


def _check_metadata(path: Path, required_metadata: list[str]) -> dict[str, Any]:
    if not path.is_file():
        return {
            "status": "MISSING",
            "path": str(path),
            "available_fields": [],
            "missing_fields": required_metadata,
        }
    try:
        manifest = _load_json(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {
            "status": "INVALID",
            "path": str(path),
            "available_fields": [],
            "missing_fields": required_metadata,
            "error": str(exc),
        }
    available = sorted(set(manifest.get("available_fields", [])))
    missing = sorted(set(required_metadata) - set(available))
    return {
        "status": "VERIFIED" if not missing else "INCOMPLETE",
        "path": str(path),
        "available_fields": available,
        "missing_fields": missing,
    }


def _artifact_path(
    artifact: ArtifactRequirement,
    *,
    dataset: dict[str, Any],
    dataset_dir: Path,
    repo_root: Path,
) -> Path:
    if artifact.registry_local_path:
        local_path = dataset.get("local_path")
        if not local_path:
            return dataset_dir / "prepared" / "registry_local_path_missing"
        return _resolve(Path(local_path), repo_root)
    if artifact.relative_path is None:
        raise ValueError(f"artifact {artifact.role} has no path")
    return dataset_dir / artifact.relative_path


def _single_next_action(
    *,
    checked: bool,
    artifact_checks: list[dict[str, Any]],
    checksum_check: dict[str, Any],
    provenance_check: dict[str, Any],
    metadata_check: dict[str, Any],
    dataset_id: str,
) -> dict[str, str]:
    if not checked:
        return {
            "action_id": "RUN_LOCAL_PREFLIGHT",
            "instruction": f"Run --check-local for {dataset_id} before acquisition claims.",
        }
    for artifact in artifact_checks:
        if artifact["status"] != "PRESENT":
            return {
                "action_id": f"ACQUIRE_{artifact['role'].upper()}",
                "instruction": f"Acquire or prepare the required artifact at {artifact['path']}.",
            }
    if provenance_check["status"] != "VERIFIED":
        return {
            "action_id": "RECORD_PROVENANCE",
            "instruction": (
                "Create a source-derived provenance.json with real accession, URL, "
                "retrieval time, and file records."
            ),
        }
    if checksum_check["status"] != "VERIFIED":
        return {
            "action_id": "VERIFY_CHECKSUMS",
            "instruction": (
                "Generate checksums.sha256 from the acquired files, then rerun the "
                "preflight; do not use placeholder hashes."
            ),
        }
    if metadata_check["status"] != "VERIFIED":
        return {
            "action_id": "RECOVER_REQUIRED_METADATA",
            "instruction": (
                "Recover the missing registered metadata and record verified field "
                "names in metadata_manifest.json."
            ),
        }
    return {
        "action_id": "RUN_REGISTERED_QC_ENDPOINTS",
        "instruction": (
            "Run the registered QC endpoints and create RunEvidence; acquisition "
            "readiness alone is not a scientific PASS."
        ),
    }


def _dataset_plan(
    *,
    dataset: dict[str, Any],
    qc_heads: set[str],
    repo_root: Path,
    target_root: Path,
    check_local: bool,
) -> dict[str, Any]:
    dataset_id = dataset["dataset_id"]
    requirement = ACQUISITION_REQUIREMENTS[dataset_id]
    dataset_dir = target_root / dataset_id
    artifacts = []
    for artifact in requirement.artifacts:
        artifact_path = _artifact_path(
            artifact,
            dataset=dataset,
            dataset_dir=dataset_dir,
            repo_root=repo_root,
        )
        checksum_key = (
            f"registry_local/{artifact_path.name}"
            if artifact.registry_local_path
            else artifact_path.relative_to(dataset_dir).as_posix()
        )
        artifacts.append(
            {
                "role": artifact.role,
                "path": str(artifact_path),
                "checksum_key": checksum_key,
                "description": artifact.description,
                "required": True,
            }
        )
    artifact_checks = [
        {
            **artifact,
            "status": ("PRESENT" if check_local and Path(artifact["path"]).is_file() else "MISSING")
            if check_local
            else "NOT_CHECKED",
        }
        for artifact in artifacts
    ]

    checksum_path = dataset_dir / "checksums.sha256"
    provenance_path = dataset_dir / "provenance.json"
    metadata_path = dataset_dir / "metadata_manifest.json"
    if check_local:
        checksum_check = _check_checksums(
            manifest_path=checksum_path,
            artifacts=artifact_checks,
        )
        provenance_check = _check_provenance(provenance_path, dataset_id)
        metadata_check = _check_metadata(metadata_path, dataset["required_metadata"])
    else:
        checksum_check = {"status": "NOT_CHECKED", "path": str(checksum_path)}
        provenance_check = {"status": "NOT_CHECKED", "path": str(provenance_path)}
        metadata_check = {
            "status": "NOT_CHECKED",
            "path": str(metadata_path),
            "available_fields": [],
            "missing_fields": dataset["required_metadata"],
        }

    checks_ready = (
        check_local
        and all(item["status"] == "PRESENT" for item in artifact_checks)
        and checksum_check["status"] == "VERIFIED"
        and provenance_check["status"] == "VERIFIED"
        and metadata_check["status"] == "VERIFIED"
    )
    if checks_ready:
        acquisition_status = "LOCAL_ARTIFACTS_VERIFIED"
        preparation_status = "READY_FOR_VALIDATION"
        overall_status = "READY_FOR_VALIDATION"
    elif check_local:
        acquisition_status = "BLOCKED_MISSING_OR_UNVERIFIED_LOCAL_EVIDENCE"
        preparation_status = "NOT_READY"
        overall_status = "BLOCKED"
    else:
        acquisition_status = "NOT_CHECKED"
        preparation_status = "NOT_READY"
        overall_status = "NOT_READY"

    accessions = dataset.get("accessions", [])
    return {
        "dataset_id": dataset_id,
        "display_name": dataset.get("display_name"),
        "priority": dataset.get("priority"),
        "accessions": accessions,
        "source_urls": [row["url"] for row in accessions if row.get("url")],
        "license": dataset.get("license"),
        "redistribution": (dataset.get("license") or {}).get("redistribution"),
        "raw_availability": dataset.get("download"),
        "required_metadata": dataset.get("required_metadata", []),
        "declared_endpoint_ids": dataset.get("endpoint_ids", []),
        "qc_endpoint_ids": [
            endpoint for endpoint in dataset.get("endpoint_ids", []) if endpoint in qc_heads
        ],
        "target_local_directory": str(dataset_dir),
        "registry_local_path": dataset.get("local_path"),
        "requirements": {
            "unfiltered_droplets": {
                "required": requirement.unfiltered_droplets_required,
                "reason": requirement.unfiltered_droplets_reason,
            },
            "label_recovery": {
                "required": requirement.label_recovery_required,
                "strategy": requirement.label_recovery_strategy,
            },
            "checksum": {
                "algorithm": "sha256",
                "manifest": str(checksum_path),
                "source": "must_be_computed_from_acquired_files",
                "placeholder_values_allowed": False,
            },
            "provenance": {
                "manifest": str(provenance_path),
                "required_fields": [
                    "dataset_id",
                    "accessions",
                    "source_urls",
                    "retrieved_at",
                    "files",
                ],
            },
        },
        "artifacts": artifact_checks,
        "local_preflight": {
            "performed": check_local,
            "checksum": checksum_check,
            "provenance": provenance_check,
            "metadata": metadata_check,
        },
        "download_status": acquisition_status,
        "preparation_status": preparation_status,
        "scientific_validation_status": "NOT_RUN",
        "overall_status": overall_status,
        "next_action": _single_next_action(
            checked=check_local,
            artifact_checks=artifact_checks,
            checksum_check=checksum_check,
            provenance_check=provenance_check,
            metadata_check=metadata_check,
            dataset_id=dataset_id,
        ),
        "limitations": dataset.get("limitations", []),
    }


def build_plan(
    registry_path: Path,
    contract_path: Path,
    *,
    repo_root: Path,
    target_root: Path | None = None,
    check_local: bool = False,
) -> dict[str, Any]:
    """Build the acquisition plan without downloading or changing source data."""

    registry = _load_json(registry_path)
    contract = _load_json(contract_path)
    datasets_by_id = {row["dataset_id"]: row for row in registry["datasets"]}
    missing = sorted(set(TARGET_DATASET_IDS) - set(datasets_by_id))
    if missing:
        raise ValueError(f"registry is missing required datasets: {', '.join(missing)}")

    qc_heads = set(contract["qc_validation_design"]["evidence_heads"])
    required_qc_portfolio = set(contract["required_dataset_portfolio"]["qc"])
    outside_portfolio = sorted(set(TARGET_DATASET_IDS) - required_qc_portfolio)
    if outside_portfolio:
        raise ValueError(
            "datasets are not in the locked QC portfolio: " + ", ".join(outside_portfolio)
        )

    target_root = target_root or repo_root / "data" / "external" / "qc_truth"
    target_root = _resolve(target_root, repo_root)
    dataset_plans = [
        _dataset_plan(
            dataset=datasets_by_id[dataset_id],
            qc_heads=qc_heads,
            repo_root=repo_root,
            target_root=target_root,
            check_local=check_local,
        )
        for dataset_id in TARGET_DATASET_IDS
    ]
    ready = sum(row["overall_status"] == "READY_FOR_VALIDATION" for row in dataset_plans)
    blocked = sum(row["overall_status"] == "BLOCKED" for row in dataset_plans)
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "source_registry": str(registry_path),
        "source_registry_schema_version": registry.get("schema_version"),
        "source_acceptance_contract": str(contract_path),
        "source_acceptance_contract_schema_version": contract.get("schema_version"),
        "operating_mode": "AUDIT_THEN_PLAN_EVIDENCE",
        "scope": "public QC truth acquisition and local preflight only",
        "claim_boundary": (
            "Acquisition readiness is not scientific PASS. URL presence, repository "
            "registration, or an unchecked local file cannot satisfy a QC endpoint."
        ),
        "check_local": check_local,
        "datasets": dataset_plans,
        "summary": {
            "dataset_count": len(dataset_plans),
            "ready_for_validation": ready,
            "blocked": blocked,
            "not_checked": len(dataset_plans) - ready - blocked,
            "scientifically_passed": 0,
        },
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("validation/dataset_evidence_registry.json"),
    )
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("validation/qc_preprocess/acceptance_contract.json"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--target-root", type=Path)
    parser.add_argument(
        "--check-local",
        action="store_true",
        help="Read-only check of required local files, provenance, metadata, and hashes.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional report path. Source data are never changed.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    repo_root = args.repo_root.resolve()
    plan = build_plan(
        _resolve(args.registry, repo_root),
        _resolve(args.contract, repo_root),
        repo_root=repo_root,
        target_root=args.target_root,
        check_local=args.check_local,
    )
    rendered = json.dumps(plan, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        output = _resolve(args.output, repo_root)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
