#!/usr/bin/env python3
"""Validate the dataset registry and report acquisition/evidence readiness.

The report deliberately separates a declared benchmark, an accessible dataset,
and a passed scientific endpoint. A dataset is never promoted because a URL or
local file merely exists.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

MODULES = {"qc", "preprocess", "analysis"}
PRIORITIES = {"P0", "P1", "P2", "ENG"}
PASSING_EVIDENCE = {"PASS", "PASS_BASELINE"}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def validate_registry(registry: dict[str, Any]) -> list[str]:
    issues: list[str] = []
    if registry.get("schema_name") != "scLucidDatasetEvidenceRegistry":
        issues.append("Unexpected registry schema_name.")
    endpoints = registry.get("endpoint_definitions", {})
    datasets = registry.get("datasets", [])
    if not endpoints:
        issues.append("No endpoint definitions are registered.")
    if not datasets:
        issues.append("No datasets are registered.")

    seen: set[str] = set()
    for row in datasets:
        dataset_id = str(row.get("dataset_id", ""))
        if not dataset_id:
            issues.append("A dataset row is missing dataset_id.")
            continue
        if dataset_id in seen:
            issues.append(f"Duplicate dataset_id: {dataset_id}")
        seen.add(dataset_id)
        if row.get("priority") not in PRIORITIES:
            issues.append(f"{dataset_id}: invalid priority.")
        modules = set(row.get("modules", []))
        if not modules or not modules <= MODULES:
            issues.append(f"{dataset_id}: invalid or empty modules.")
        if not row.get("accessions"):
            issues.append(f"{dataset_id}: accession/download entry is missing.")
        download = row.get("download", {})
        for key in ("raw_reads", "raw_counts", "processed"):
            if not download.get(key):
                issues.append(f"{dataset_id}: download.{key} is missing.")
        license_record = row.get("license", {})
        for key in ("status", "identifier", "url", "redistribution"):
            if license_record.get(key) in (None, ""):
                issues.append(f"{dataset_id}: license.{key} is missing.")
        if not row.get("required_metadata"):
            issues.append(f"{dataset_id}: required_metadata is empty.")
        for endpoint_id in row.get("endpoint_ids", []):
            if endpoint_id not in endpoints:
                issues.append(f"{dataset_id}: unknown endpoint {endpoint_id}.")

    for endpoint_id, endpoint in endpoints.items():
        if endpoint.get("module") not in MODULES:
            issues.append(f"{endpoint_id}: invalid module.")
        for field in ("estimand", "experimental_unit", "metrics", "acceptance"):
            if not endpoint.get(field):
                issues.append(f"{endpoint_id}: {field} is missing.")
    return issues


def validate_required_endpoint_portfolio(
    registry: dict[str, Any],
    required_endpoint_portfolio: dict[str, dict[str, list[str]]],
) -> list[str]:
    """Validate exact dataset/endpoint bindings used by release gates."""
    issues: list[str] = []
    if not isinstance(required_endpoint_portfolio, dict):
        return ["Required endpoint portfolio must be a mapping."]
    endpoints = registry.get("endpoint_definitions", {})
    datasets = {
        str(row.get("dataset_id")): row for row in registry.get("datasets", [])
    }
    for module, endpoint_map in required_endpoint_portfolio.items():
        if module not in MODULES:
            issues.append(f"Required endpoint portfolio has unknown module: {module}")
            continue
        if not isinstance(endpoint_map, dict):
            issues.append(f"Required endpoint portfolio for {module} must be a mapping.")
            continue
        if not endpoint_map:
            issues.append(f"Required endpoint portfolio for {module} is empty.")
        for endpoint_id, dataset_ids in endpoint_map.items():
            endpoint = endpoints.get(endpoint_id)
            if endpoint is None:
                issues.append(f"Required endpoint portfolio references unknown endpoint: {endpoint_id}")
                continue
            if endpoint.get("module") != module:
                issues.append(
                    f"{endpoint_id}: required under {module}, but registered under "
                    f"{endpoint.get('module')}."
                )
            if not isinstance(dataset_ids, list):
                issues.append(f"{endpoint_id}: required datasets must be a list.")
                continue
            if not dataset_ids:
                issues.append(f"{endpoint_id}: required dataset list is empty.")
            string_dataset_ids = [
                dataset_id for dataset_id in dataset_ids if isinstance(dataset_id, str)
            ]
            if len(string_dataset_ids) != len(set(string_dataset_ids)):
                issues.append(f"{endpoint_id}: required dataset list contains duplicates.")
            for dataset_id in dataset_ids:
                if not isinstance(dataset_id, str) or not dataset_id:
                    issues.append(f"{endpoint_id}: required dataset_id must be a non-empty string.")
                    continue
                dataset = datasets.get(dataset_id)
                if dataset is None:
                    issues.append(
                        f"{endpoint_id}: required dataset is not registered: {dataset_id}"
                    )
                elif endpoint_id not in dataset.get("endpoint_ids", []):
                    issues.append(
                        f"{dataset_id}: does not declare required endpoint {endpoint_id}."
                    )
    return issues


def _downloadable(row: dict[str, Any]) -> bool:
    values = " ".join(str(value) for value in row["download"].values())
    return "AVAILABLE" in values or "GENERATED_LOCAL" in values


def _evidence_lookup(index: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for run in index.get("runs", []):
        key = (str(run["dataset_id"]), str(run["endpoint_id"]))
        if key in lookup:
            raise ValueError(f"Duplicate evidence binding: {key[0]} / {key[1]}")
        lookup[key] = run
    return lookup


def build_report(
    registry_path: Path,
    evidence_index_path: Path,
    *,
    repo_root: Path,
    required_portfolio: dict[str, list[str]] | None = None,
    required_endpoint_portfolio: dict[str, dict[str, list[str]]] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    registry = _load_json(registry_path)
    evidence_index = _load_json(evidence_index_path)
    endpoint_definitions = registry.get("endpoint_definitions", {})
    issues = validate_registry(registry)
    required_endpoint_portfolio = required_endpoint_portfolio or {}
    issues.extend(
        validate_required_endpoint_portfolio(registry, required_endpoint_portfolio)
    )
    evidence = _evidence_lookup(evidence_index)
    declared_by_dataset = {
        str(row["dataset_id"]): set(row.get("endpoint_ids", []))
        for row in registry.get("datasets", [])
    }
    for dataset_id, endpoint_id in evidence:
        if dataset_id not in declared_by_dataset:
            issues.append(f"Evidence index references unknown dataset: {dataset_id}")
        elif endpoint_id not in declared_by_dataset[dataset_id]:
            issues.append(
                f"Evidence index binds undeclared endpoint {endpoint_id} to {dataset_id}."
            )
    rows: list[dict[str, Any]] = []
    binding_status: dict[tuple[str, str], str] = {}

    for dataset in registry.get("datasets", []):
        dataset_id = dataset["dataset_id"]
        local_path = dataset.get("local_path")
        local_exists = bool(local_path and (repo_root / local_path).exists())
        endpoint_statuses: list[str] = []
        endpoint_status_by_id: dict[str, str] = {}
        missing_artifacts: list[str] = []
        for endpoint_id in dataset.get("endpoint_ids", []):
            run = evidence.get((dataset_id, endpoint_id))
            if run is None:
                endpoint_statuses.append("NOT_RUN")
                endpoint_status_by_id[endpoint_id] = "NOT_RUN"
                binding_status[(dataset_id, endpoint_id)] = "NOT_RUN"
                continue
            status = str(run.get("status", "NOT_RUN"))
            artifact = run.get("artifact")
            if status in PASSING_EVIDENCE and not artifact:
                status = "MISSING_ARTIFACT"
                missing_artifacts.append("<not declared>")
            elif artifact and not (repo_root / artifact).exists():
                status = "MISSING_ARTIFACT"
                missing_artifacts.append(str(artifact))
            endpoint_statuses.append(status)
            endpoint_status_by_id[endpoint_id] = status
            binding_status[(dataset_id, endpoint_id)] = status

        module_statuses: dict[str, str] = {}
        for module in MODULES:
            statuses = [
                endpoint_status_by_id[endpoint_id]
                for endpoint_id in dataset.get("endpoint_ids", [])
                if endpoint_definitions.get(endpoint_id, {}).get("module") == module
            ]
            if not statuses:
                module_statuses[module] = "NOT_APPLICABLE"
            elif all(status in PASSING_EVIDENCE for status in statuses):
                module_statuses[module] = "PASS"
            elif any(status == "FAIL" for status in statuses):
                module_statuses[module] = "FAIL"
            elif any(status in {"BLOCKED", "MISSING_ARTIFACT"} for status in statuses):
                module_statuses[module] = "BLOCKED"
            elif any(status == "REVIEW" for status in statuses):
                module_statuses[module] = "REVIEW"
            else:
                module_statuses[module] = "NOT_RUN"

        if not dataset.get("endpoint_ids"):
            scientific_status = "ENGINEERING_ONLY"
        elif endpoint_statuses and all(status in PASSING_EVIDENCE for status in endpoint_statuses):
            scientific_status = "PASS"
        elif any(status == "FAIL" for status in endpoint_statuses):
            scientific_status = "FAIL"
        elif any(status in {"BLOCKED", "MISSING_ARTIFACT"} for status in endpoint_statuses):
            scientific_status = "BLOCKED"
        elif any(status == "REVIEW" for status in endpoint_statuses):
            scientific_status = "REVIEW"
        else:
            scientific_status = "NOT_RUN"

        rows.append(
            {
                "dataset_id": dataset_id,
                "priority": dataset["priority"],
                "release_scope": dataset["release_scope"],
                "modules": ";".join(dataset["modules"]),
                "truth_types": ";".join(dataset["truth_types"]),
                "accession": ";".join(item["id"] for item in dataset["accessions"]),
                "raw_reads": dataset["download"]["raw_reads"],
                "raw_counts": dataset["download"]["raw_counts"],
                "processed": dataset["download"]["processed"],
                "license_status": dataset["license"]["status"],
                "redistribution": dataset["license"]["redistribution"],
                "local_path": local_path or "",
                "local_exists": local_exists,
                "acquisition_status": (
                    "LOCAL_READY"
                    if local_exists
                    else "DOWNLOADABLE"
                    if _downloadable(dataset)
                    else "CONTROLLED_OR_REVIEW"
                ),
                "required_metadata": ";".join(dataset["required_metadata"]),
                "endpoint_ids": ";".join(dataset.get("endpoint_ids", [])),
                "endpoint_statuses": ";".join(endpoint_statuses),
                "scientific_status": scientific_status,
                "qc_status": module_statuses["qc"],
                "preprocess_status": module_statuses["preprocess"],
                "analysis_status": module_statuses["analysis"],
                "missing_artifacts": ";".join(missing_artifacts),
            }
        )

    required_portfolio = required_portfolio or {}
    by_id = {row["dataset_id"]: row for row in rows}
    endpoint_gates: dict[str, dict[str, Any]] = {}
    for module, endpoint_map in sorted(required_endpoint_portfolio.items()):
        endpoint_gates[module] = {}
        for endpoint_id, dataset_ids in endpoint_map.items():
            statuses = {
                dataset_id: binding_status.get((dataset_id, endpoint_id), "NOT_RUN")
                for dataset_id in dataset_ids
            }
            not_passed = [
                dataset_id
                for dataset_id, status in statuses.items()
                if status not in PASSING_EVIDENCE
            ]
            endpoint_gates[module][endpoint_id] = {
                "status": "PASS" if dataset_ids and not not_passed else "BLOCKED",
                "required_datasets": list(dataset_ids),
                "dataset_statuses": statuses,
                "datasets_without_passing_evidence": not_passed,
            }

    module_gates: dict[str, dict[str, Any]] = {}
    for module in sorted(MODULES):
        legacy_required = list(required_portfolio.get(module, []))
        legacy_missing = [
            dataset_id for dataset_id in legacy_required if dataset_id not in by_id
        ]
        legacy_not_passed = [
            dataset_id
            for dataset_id in legacy_required
            if dataset_id in by_id and by_id[dataset_id][f"{module}_status"] != "PASS"
        ]
        module_endpoint_gates = endpoint_gates.get(module, {})
        blocked_endpoints = [
            endpoint_id
            for endpoint_id, gate in module_endpoint_gates.items()
            if gate["status"] != "PASS"
        ]
        if module_endpoint_gates:
            required = list(
                dict.fromkeys(
                    dataset_id
                    for gate in module_endpoint_gates.values()
                    for dataset_id in gate["required_datasets"]
                )
            )
            missing = [dataset_id for dataset_id in required if dataset_id not in by_id]
            not_passed = list(
                dict.fromkeys(
                    dataset_id
                    for gate in module_endpoint_gates.values()
                    for dataset_id in gate["datasets_without_passing_evidence"]
                )
            )
            gate_status = "PASS" if not blocked_endpoints else "BLOCKED"
            gate_basis = "required_endpoint_portfolio"
        else:
            required = legacy_required
            missing = legacy_missing
            not_passed = legacy_not_passed
            gate_status = "PASS" if required and not missing and not_passed == [] else "BLOCKED"
            gate_basis = "required_dataset_portfolio"
        module_gates[module] = {
            "status": gate_status,
            "gate_basis": gate_basis,
            "required_datasets": required,
            "missing_registry_rows": missing,
            "datasets_without_passed_endpoints": not_passed,
            "required_endpoints": list(module_endpoint_gates),
            "endpoints_without_passed_evidence": blocked_endpoints,
        }

    license_counts: dict[str, int] = defaultdict(int)
    scientific_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        license_counts[row["license_status"]] += 1
        scientific_counts[row["scientific_status"]] += 1
    report = {
        "schema_version": "sclucid_dataset_evidence_readiness_v1.1",
        "status": "PASS" if not issues else "BLOCKED",
        "registry_schema_status": "PASS" if not issues else "BLOCKED",
        "issues": issues,
        "n_datasets": len(rows),
        "n_endpoints": len(registry.get("endpoint_definitions", {})),
        "license_status_counts": dict(sorted(license_counts.items())),
        "scientific_status_counts": dict(sorted(scientific_counts.items())),
        "endpoint_gates": endpoint_gates,
        "module_gates": module_gates,
        "claim_boundary": (
            "Registry PASS means the portfolio is specified, not that scientific endpoints passed."
        ),
    }
    return report, rows


def write_report(report: dict[str, Any], rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "dataset_evidence_readiness.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output_dir / "dataset_evidence_readiness.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    lines = [
        "# scLucid dataset evidence readiness",
        "",
        f"Registry schema: **{report['registry_schema_status']}**",
        "",
        "A registered or downloadable dataset is not a passed scientific endpoint.",
        "",
        "## Module gates",
        "",
    ]
    for module, gate in report["module_gates"].items():
        lines.append(f"- {module}: {gate['status']}")
        if gate["endpoints_without_passed_evidence"]:
            lines.append(
                "  - blocked evidence heads: "
                + ", ".join(gate["endpoints_without_passed_evidence"])
            )
        if gate["datasets_without_passed_endpoints"]:
            lines.append(
                "  - awaiting passed evidence: "
                + ", ".join(gate["datasets_without_passed_endpoints"])
            )
    lines.extend(["", "## Dataset rows", ""])
    for row in rows:
        lines.append(
            f"- {row['dataset_id']}: acquisition={row['acquisition_status']}; "
            f"scientific={row['scientific_status']}; "
            f"qc={row['qc_status']}; preprocess={row['preprocess_status']}; "
            f"analysis={row['analysis_status']}; license={row['license_status']}"
        )
    lines.append("")
    (output_dir / "dataset_evidence_readiness.md").write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("validation/dataset_evidence_registry.json"),
    )
    parser.add_argument(
        "--evidence-index",
        type=Path,
        default=Path("validation/evidence_run_index.json"),
    )
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("validation/qc_preprocess/acceptance_contract.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/dataset_registry"),
    )
    args = parser.parse_args()
    contract = _load_json(args.contract)
    report, rows = build_report(
        args.registry,
        args.evidence_index,
        repo_root=Path.cwd(),
        required_portfolio=contract.get("required_dataset_portfolio", {}),
        required_endpoint_portfolio=contract.get("required_endpoint_portfolio", {}),
    )
    write_report(report, rows, args.output_dir)
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
