#!/usr/bin/env python3
"""Unblind frozen expert labels and run the locked scLucid QC acceptance gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scLucid as scl
from validation.qc_preprocess.locked_acceptance import evaluate_qc_policy_acceptance
from validation.qc_preprocess.truth_pack import sha256_file, validate_frozen_labels


def _write_report(report: dict[str, Any], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "locked_qc_acceptance.json"
    md_path = output_dir / "locked_qc_acceptance.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Locked QC acceptance",
        "",
        f"Status: **{report['status']}**",
        "",
        report["claim_boundary"],
        "",
        "## Evidence gates",
        "",
    ]
    for name in (
        "label_gate",
        "dataset_coverage",
        "source_integrity",
        "cell_endpoint",
        "sample_endpoint",
    ):
        payload = report.get(name, {})
        lines.append(f"- {name}: {payload.get('status', 'NOT_RUN')}")
        if payload.get("reason"):
            lines.append(f"  - {payload['reason']}")
    lines.extend(
        [
            "",
            "## Next action",
            "",
            report["next_action"],
            "",
        ]
    )
    md_path.write_text("\n".join(lines))
    return {"json": json_path, "markdown": md_path}


def _blocked_report(reason: str, issues: list[str]) -> dict[str, Any]:
    return {
        "schema_version": "sclucid_locked_qc_acceptance_v1",
        "status": "BLOCKED",
        "dataset_registry_ids": [
            "lin2020_pdac",
            "kang2018_pbmc",
            "real_project_panel",
        ],
        "preregistered_endpoint_ids": [
            "qc_selector_superiority",
            "qc_catastrophic_sample_detection",
            "qc_doublet_calibration",
        ],
        "label_gate": {"status": "BLOCKED", "reason": reason, "issues": issues},
        "dataset_coverage": {"status": "NOT_RUN"},
        "source_integrity": {"status": "NOT_RUN"},
        "cell_endpoint": {"status": "NOT_RUN"},
        "sample_endpoint": {"status": "NOT_RUN"},
        "secondary_endpoints": {
            "doublet_auprc": "NOT_RUN_SEPARATE_EVIDENCE_HEAD",
            "ambient": "NOT_EVALUABLE_WITHOUT_RAW_DROPLETS_OR_TRUSTWORTHY_POPULATIONS",
        },
        "claim_boundary": "Scientific superiority is not established.",
        "next_action": "Complete and freeze all blinded expert labels, then rerun this command.",
    }


def _dataset_coverage(manifest: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    observed = {str(record["key"]) for record in manifest["datasets"]}
    expected_development = set(map(str, contract.get("development_datasets", [])))
    expected_projects = set(map(str, contract.get("external_projects", [])))
    missing_development = sorted(expected_development - observed)
    missing_projects = sorted(
        project
        for project in expected_projects
        if not any(key == project or key.startswith(f"{project}_") for key in observed)
    )
    passed = not missing_development and not missing_projects
    return {
        "status": "PASS" if passed else "BLOCKED",
        "scope": "frozen_truth_pack_inputs_only",
        "observed_dataset_keys": sorted(observed),
        "missing_development_datasets": missing_development,
        "missing_external_projects": missing_projects,
        "reason": (
            "All inputs registered for this frozen truth pack are represented; this is not the external P0 portfolio gate."
            if passed
            else "The frozen truth pack does not cover every registered development dataset/project."
        ),
    }


def _doublet_endpoint(summary_path: Path | None) -> dict[str, Any]:
    if summary_path is None or not summary_path.exists():
        return {"status": "NOT_RUN", "reason": "No demuxlet benchmark summary was attached."}
    payload = json.loads(summary_path.read_text())
    return {
        "status": "EVALUATED",
        "dataset": "kang2018.pbmc",
        "best_method": payload.get("best_method"),
        "best_method_auprc": payload.get("best_method_auprc"),
        "best_method_auc": payload.get("best_method_auc"),
        "best_method_f1": payload.get("best_method_f1"),
        "claim_boundary": (
            "Method calibration evidence for Kang/demuxlet only; no automatic-removal claim."
        ),
    }


def _sample_gate(
    sample_labels: pd.DataFrame,
    sample_key: pd.DataFrame,
    policy_sample_decisions: list[dict[str, Any]],
    contract: dict[str, Any],
) -> dict[str, Any]:
    labels = sample_labels.set_index("case_id")["expert_label"]
    mapping = sample_key.set_index("case_id")
    decisions = {
        (str(row["dataset_key"]), str(row["original_sample"])): str(row["policy_status"])
        for row in policy_sample_decisions
    }
    rows: list[dict[str, str]] = []
    for case_id, label in labels.items():
        if label == "UNCERTAIN":
            continue
        source = mapping.loc[case_id]
        predicted = decisions.get((str(source["dataset_key"]), str(source["original_sample"])))
        rows.append(
            {
                "case_id": case_id,
                "expert_label": label,
                "policy_status": predicted or "MISSING",
            }
        )
    frame = pd.DataFrame(rows)
    if (
        frame.empty
        or not (frame["expert_label"] == "KEEP").any()
        or not (frame["expert_label"] == "REMOVE").any()
    ):
        return {
            "status": "NOT_EVALUABLE",
            "reason": "Both expert KEEP and REMOVE sample/library labels are required.",
            "rows": rows,
        }
    remove = frame["expert_label"] == "REMOVE"
    keep = frame["expert_label"] == "KEEP"
    detection = float((frame.loc[remove, "policy_status"] == "BLOCKED").mean())
    false_block = float((frame.loc[keep, "policy_status"] == "BLOCKED").mean())
    required_detection = float(
        contract["qc_primary_endpoint"]["catastrophic_sample_detection_rate"]
    )
    required_false_block = float(
        contract["qc_primary_endpoint"]["locked_high_quality_false_block_rate"]
    )
    passed = detection >= required_detection and false_block <= required_false_block
    return {
        "status": "PASS" if passed else "FAIL",
        "catastrophic_sample_detection_rate": detection,
        "locked_high_quality_false_block_rate": false_block,
        "thresholds": {
            "minimum_detection_rate": required_detection,
            "maximum_false_block_rate": required_false_block,
        },
        "rows": rows,
    }


def run_locked_qc_acceptance(
    pack_dir: Path,
    output_dir: Path,
    *,
    contract_path: Path,
    sample_labels_path: Path | None = None,
    cell_labels_path: Path | None = None,
    doublet_summary_path: Path | None = None,
) -> dict[str, Any]:
    reviewer_dir = pack_dir / "reviewer"
    sealed_dir = pack_dir / "sealed"
    contract = json.loads(contract_path.read_text())
    manifest = json.loads((sealed_dir / "manifest.json").read_text())
    coverage = _dataset_coverage(manifest, contract)
    sample_labels_path = sample_labels_path or reviewer_dir / "sample_labels.tsv"
    cell_labels_path = cell_labels_path or reviewer_dir / "cell_labels.tsv"
    sample_labels, sample_issues = validate_frozen_labels(
        reviewer_dir / "sample_evidence.tsv", sample_labels_path
    )
    cell_labels, cell_issues = validate_frozen_labels(
        reviewer_dir / "cell_evidence.tsv", cell_labels_path
    )
    issues = [*sample_issues, *cell_issues]
    if issues:
        report = _blocked_report("Blinded labels are incomplete or invalid.", issues)
        report["dataset_coverage"] = coverage
        report["secondary_endpoints"]["doublet_auprc"] = _doublet_endpoint(
            doublet_summary_path
        )
        _write_report(report, output_dir)
        return report

    if coverage["status"] != "PASS":
        coverage_issues = [
            f"Missing development dataset: {key}"
            for key in coverage["missing_development_datasets"]
        ] + [
            f"Missing external project: {key}" for key in coverage["missing_external_projects"]
        ]
        report = _blocked_report("Registered truth-source coverage is incomplete.", coverage_issues)
        report["label_gate"] = {"status": "PASS"}
        report["dataset_coverage"] = coverage
        report["secondary_endpoints"]["doublet_auprc"] = _doublet_endpoint(
            doublet_summary_path
        )
        report["next_action"] = "Add the missing registered source and build a new versioned truth pack."
        _write_report(report, output_dir)
        return report

    source_issues: list[str] = []
    for record in manifest["datasets"]:
        path = Path(record["path"])
        if not path.exists():
            source_issues.append(f"Missing frozen source: {record['dataset_alias']}")
        elif sha256_file(path) != record["source_sha256"]:
            source_issues.append(f"Source fingerprint changed: {record['dataset_alias']}")
    if source_issues:
        report = _blocked_report("Frozen source integrity failed.", source_issues)
        report["secondary_endpoints"]["doublet_auprc"] = _doublet_endpoint(
            doublet_summary_path
        )
        report["label_gate"] = {"status": "PASS"}
        report["dataset_coverage"] = coverage
        report["source_integrity"] = {
            "status": "BLOCKED",
            "reason": "One or more source fingerprints changed.",
            "issues": source_issues,
        }
        report["next_action"] = (
            "Restore the frozen source versions or build a new versioned truth pack."
        )
        _write_report(report, output_dir)
        return report

    sealed_cells = pd.read_csv(sealed_dir / "cell_key.tsv", sep="\t", dtype=str)
    sealed_samples = pd.read_csv(sealed_dir / "sample_key.tsv", sep="\t", dtype=str)
    label_lookup = cell_labels.set_index("case_id")["expert_label"]
    primary_cells = sealed_cells[sealed_cells["sampling_tier"] == "primary_uniform"].copy()
    primary_cells["expert_label"] = primary_cells["case_id"].map(label_lookup)

    aggregate_remove: list[str] = []
    aggregate_baselines: dict[str, list[str]] = {
        name: [] for name in contract["qc_primary_endpoint"]["baselines"]
    }
    policy_sample_decisions: list[dict[str, Any]] = []
    dataset_runs: list[dict[str, Any]] = []
    for record in manifest["datasets"]:
        dataset_key = str(record["key"])
        dataset_alias = str(record["dataset_alias"])
        adata = ad.read_h5ad(record["path"])
        context = scl.ProjectContext(**record["context"])
        card = scl.recommend_qc_policy(adata, context)
        subset = primary_cells[primary_cells["dataset_key"] == dataset_key]
        obs_to_case = subset.set_index("original_obs_name")["case_id"].to_dict()
        aggregate_remove.extend(
            obs_to_case[name] for name in card.policy.remove_obs_names if name in obs_to_case
        )
        candidates = {row["name"]: row for row in card.policy.candidate_policies}
        for baseline in aggregate_baselines:
            calls = candidates.get(baseline, {}).get("flagged_obs_names", [])
            aggregate_baselines[baseline].extend(
                obs_to_case[name] for name in calls if name in obs_to_case
            )
        for row in card.policy.sample_decisions:
            policy_sample_decisions.append(
                {
                    "dataset_key": dataset_key,
                    "original_sample": str(row["sample"]),
                    "policy_status": str(row["status"]),
                }
            )
        dataset_runs.append(
            {
                "dataset_alias": dataset_alias,
                "policy_status": card.status,
                "n_primary_cases": int(len(subset)),
            }
        )
        del adata

    composite_policy = SimpleNamespace(
        remove_obs_names=aggregate_remove,
        candidate_policies=[
            {"name": name, "flagged_obs_names": calls}
            for name, calls in aggregate_baselines.items()
        ],
    )
    truth = primary_cells.set_index("case_id")["expert_label"]
    groups = primary_cells.set_index("case_id")["library_alias"]
    endpoint = contract["qc_primary_endpoint"]
    cell_endpoint = evaluate_qc_policy_acceptance(
        composite_policy,
        truth,
        groups,
        baseline_names=endpoint["baselines"],
        max_keep_false_removal=float(endpoint["max_keep_false_removal_rate"]),
        min_absolute_recall_gain=float(endpoint["min_absolute_low_quality_recall_gain"]),
        n_bootstrap=1000,
        seed=int(manifest["seed"]),
    )
    sample_endpoint = _sample_gate(
        sample_labels,
        sealed_samples,
        policy_sample_decisions,
        contract,
    )
    if cell_endpoint["status"] == "NOT_EVALUABLE" or sample_endpoint["status"] == "NOT_EVALUABLE":
        status = "BLOCKED"
    elif cell_endpoint["status"] == "PASS" and sample_endpoint["status"] == "PASS":
        status = "PASS"
    else:
        status = "FAIL"
    report = {
        "schema_version": "sclucid_locked_qc_acceptance_v1",
        "status": status,
        "dataset_registry_ids": [
            "lin2020_pdac",
            "kang2018_pbmc",
            "real_project_panel",
        ],
        "preregistered_endpoint_ids": [
            "qc_selector_superiority",
            "qc_catastrophic_sample_detection",
            "qc_doublet_calibration",
        ],
        "label_gate": {"status": "PASS"},
        "dataset_coverage": coverage,
        "source_integrity": {"status": "PASS", "n_sources": len(manifest["datasets"])},
        "cell_endpoint": cell_endpoint,
        "sample_endpoint": sample_endpoint,
        "dataset_runs": dataset_runs,
        "secondary_endpoints": {
            "doublet_auprc": _doublet_endpoint(doublet_summary_path),
            "ambient": "SEPARATE_RAW_DROPLET_EVIDENCE_REQUIRED",
            "secondary_metric_challenge_cases": int(
                (sealed_cells["sampling_tier"] == "secondary_metric_challenge").sum()
            ),
        },
        "claim_boundary": (
            "PASS supports superiority only for this frozen truth pack and registered baselines."
            if status == "PASS"
            else "Scientific superiority is not established."
        ),
        "next_action": (
            "Freeze this report and proceed to leave-one-project-out preprocessing validation."
            if status == "PASS"
            else "Inspect the failed or non-evaluable gate without changing the locked thresholds."
        ),
    }
    _write_report(report, output_dir)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pack-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path(__file__).with_name("acceptance_contract.json"),
    )
    parser.add_argument("--sample-labels", type=Path)
    parser.add_argument("--cell-labels", type=Path)
    parser.add_argument("--doublet-summary", type=Path)
    args = parser.parse_args()
    report = run_locked_qc_acceptance(
        args.pack_dir,
        args.output_dir,
        contract_path=args.contract,
        sample_labels_path=args.sample_labels,
        cell_labels_path=args.cell_labels,
        doublet_summary_path=args.doublet_summary,
    )
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
