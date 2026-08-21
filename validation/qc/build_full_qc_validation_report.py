#!/usr/bin/env python3
"""Build the complete QC evidence-head readiness report without scoring.

The dataset readiness report is an execution index, while the acceptance
contract is the release specification.  This report joins the two at the exact
``endpoint x required dataset`` level.  It deliberately does not average head
statuses: every registered QC evidence head must pass for the overall QC gate
to pass.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

DEFAULT_OUTPUT_DIR = Path("validation_outputs/current/qc_full_gate")
PASSING_EVIDENCE_DEFAULT = {"PASS", "PASS_BASELINE"}


HEAD_METADATA: dict[str, dict[str, Any]] = {
    "qc_input_contract": {
        "evidence_class": "ENGINEERING_CONTRACT",
        "supports": "Count/provenance integrity, fail-closed sample identity, and read-only review behavior.",
        "does_not_support": "Biological low-quality-cell detection or superiority over a QC baseline.",
        "next_action": "Run the controlled and real-project input-contract checks and bind their RunEvidence.",
    },
    "qc_profile_selection": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Assay-aware profile selection without a universal mitochondrial threshold.",
        "does_not_support": "A single QC threshold that generalizes across tissues, assays, or protocols.",
        "next_action": "Validate profile selection on locked high-quality, failed, and real-project libraries.",
    },
    "qc_cell_calling": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "True-cell and low-RNA-cell recovery from unfiltered droplets at controlled empty-droplet FDR.",
        "does_not_support": "Cell-calling performance for projects that provide only a filtered matrix.",
        "next_action": "Run cell-calling benchmarks with microscopy, hashing, and cross-species truth.",
    },
    "qc_ambient_correction": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Reduction of known ambient signal while preserving native expression and identity.",
        "does_not_support": "Automatic ambient correction when raw droplets or a credible soup reference are absent.",
        "next_action": "Run ambient-removal and native-expression preservation endpoints on registered raw-droplet truth.",
    },
    "qc_catastrophic_sample_detection": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Detection of known failed libraries without falsely blocking locked high-quality libraries.",
        "does_not_support": "Reliable cell-level filtering within every non-catastrophic sample.",
        "next_action": "Freeze blinded sample labels and run grouped catastrophic-sample acceptance.",
    },
    "qc_damage_classification": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "KEEP/REMOVE discrimination and calibrated uncertainty for damaged or broken cells.",
        "does_not_support": "Treating stress, apoptosis, mitochondrial fraction, or cluster position alone as cell damage truth.",
        "next_action": "Complete independent damage labels and evaluate binary and three-class calibration endpoints.",
    },
    "qc_doublet_calibration": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Dataset-specific doublet ranking and calibration against experimental multiplet labels.",
        "does_not_support": "Universal automatic removal from one unconfirmed doublet method.",
        "next_action": "Run all registered experimental-label datasets and calibrate AUPRC regret and probabilities.",
    },
    "qc_rare_population_preservation": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Rare-lineage and tumor-population retention under the registered QC policy.",
        "does_not_support": "Biological validity of every retained rare cluster or tumor program.",
        "next_action": "Measure false removal and lineage-retention gaps on tumor and expert-reviewed projects.",
    },
    "qc_selector_superiority": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Selector recall gain over both global-threshold and per-sample MAD baselines at the KEEP guardrail.",
        "does_not_support": "Universal superiority outside the registered datasets, baselines, and estimand.",
        "next_action": "Finish blinded labels and run grouped-bootstrap selector comparisons across required datasets.",
    },
    "qc_iterative_review": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": "Evidence-dependent multi-round review that stops or defers when no independent evidence is added.",
        "does_not_support": "Removal based solely on embedding or cluster position, or a fixed universal number of QC rounds.",
        "next_action": "Execute round-by-round counterfactual review and verify each new removal has independent evidence.",
    },
    "qc_policy_execution": {
        "evidence_class": "ENGINEERING_CONTRACT",
        "supports": "Exact, deterministic, idempotent, count-preserving application of an explicitly reviewed policy.",
        "does_not_support": "Scientific correctness of the policy that was applied.",
        "next_action": "Run controlled and real-project policy agreement, idempotence, and counts-preservation checks.",
    },
    "qc_decisioncard_ux": {
        "evidence_class": "UX_VALIDATION",
        "supports": "Correct next-action comprehension and reduced manual configuration in registered real projects.",
        "does_not_support": "Scientific superiority of the recommended QC policy.",
        "next_action": "Complete the three real-project runs without manual workarounds and record RunEvidence completeness.",
    },
    "qc_scalability": {
        "evidence_class": "ENGINEERING_ROBUSTNESS",
        "supports": "Deterministic sparse execution, runtime, and peak-memory behavior on declared reference hardware and sizes.",
        "does_not_support": "A hardware-independent wall-time guarantee or scientific generalization.",
        "next_action": "Run registered scale tiers on reference hardware and record runtime, peak memory, density, and repeatability.",
    },
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _passing_statuses(contract: dict[str, Any]) -> set[str]:
    configured = contract.get("portfolio_gate_policy", {}).get("passing_evidence_statuses", [])
    return set(configured) if configured else set(PASSING_EVIDENCE_DEFAULT)


def _head_metadata(endpoint_id: str) -> dict[str, Any]:
    return HEAD_METADATA.get(
        endpoint_id,
        {
            "evidence_class": "UNCLASSIFIED",
            "supports": "Only the explicitly registered endpoint estimand.",
            "does_not_support": "Any stronger scientific, engineering, or UX claim.",
            "next_action": f"Define the evidence class and execute the required datasets for {endpoint_id}.",
        },
    )


def build_report(
    readiness_path: Path,
    acceptance_contract_path: Path,
) -> dict[str, Any]:
    """Return evidence-head readiness from exact required dataset bindings."""
    readiness = _load_json(readiness_path)
    contract = _load_json(acceptance_contract_path)
    evidence_heads = list(contract.get("qc_validation_design", {}).get("evidence_heads", []))
    required_by_head = contract.get("required_endpoint_portfolio", {}).get("qc", {})
    source_gates = readiness.get("endpoint_gates", {}).get("qc", {})
    passing_statuses = _passing_statuses(contract)

    heads: list[dict[str, Any]] = []
    for endpoint_id in evidence_heads:
        source_gate = source_gates.get(endpoint_id, {})
        required_datasets = list(required_by_head.get(endpoint_id, []))
        source_statuses = source_gate.get("dataset_statuses", {})
        dataset_statuses = {
            dataset_id: str(source_statuses.get(dataset_id, "NOT_RUN"))
            for dataset_id in required_datasets
        }
        nonpassing = [
            dataset_id
            for dataset_id, status in dataset_statuses.items()
            if status not in passing_statuses
        ]
        blockers: list[str] = []
        if endpoint_id not in required_by_head:
            blockers.append(
                "The acceptance contract has no required dataset binding for this head."
            )
        elif not required_datasets:
            blockers.append("The acceptance contract binds no required datasets to this head.")
        if endpoint_id not in source_gates:
            blockers.append("The dataset readiness report has no endpoint gate for this head.")
        blockers.extend(
            f"{dataset_id}: {dataset_statuses[dataset_id]}" for dataset_id in nonpassing
        )

        status = "PASS" if required_datasets and not blockers else "BLOCKED"
        metadata = _head_metadata(endpoint_id)
        heads.append(
            {
                "endpoint_id": endpoint_id,
                "status": status,
                "source_gate_status": source_gate.get("status", "NOT_RUN"),
                "evidence_class": metadata["evidence_class"],
                "required_datasets": required_datasets,
                "dataset_statuses": dataset_statuses,
                "blockers": blockers,
                "next_action": metadata["next_action"],
                "claim_boundary": {
                    "supports": metadata["supports"],
                    "does_not_support": metadata["does_not_support"],
                },
            }
        )

    blocked_heads = [head["endpoint_id"] for head in heads if head["status"] != "PASS"]
    overall_status = "PASS" if heads and not blocked_heads else "BLOCKED"
    next_action = next(
        (head["next_action"] for head in heads if head["status"] != "PASS"),
        "All registered QC evidence heads passed; preserve the evidence lock for release review.",
    )
    return {
        "schema_version": "sclucid_full_qc_validation_readiness_v1",
        "status": overall_status,
        "source_readiness": str(readiness_path),
        "source_acceptance_contract": str(acceptance_contract_path),
        "no_aggregate_quality_score": True,
        "passing_evidence_statuses": sorted(passing_statuses),
        "nonpassing_status_examples": [
            "NOT_EVALUABLE",
            "SIMULATION_PASS_NOT_EXTERNAL",
            "CONTRACT_PASS_NOT_PERFORMANCE",
        ],
        "head_count": len(heads),
        "passed_head_count": len(heads) - len(blocked_heads),
        "blocked_head_count": len(blocked_heads),
        "blocked_heads": blocked_heads,
        "evidence_heads": heads,
        "next_action": next_action,
        "claim_boundary": {
            "supported": (
                "Only heads marked PASS have complete passing evidence for every "
                "dataset required by the locked acceptance contract."
            ),
            "unsupported": (
                "Registry availability, simulated performance, contract conformance, or "
                "a partial set of passed heads does not establish QC scientific readiness."
            ),
        },
    }


def write_report(report: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "full_qc_validation_readiness.json"
    markdown_path = output_dir / "full_qc_validation_readiness.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# scLucid full QC validation readiness",
        "",
        f"Overall status: **{report['status']}**",
        "",
        "No aggregate quality score is used. Every required evidence head must pass.",
        "",
        "| Evidence head | Status | Evidence class | Required datasets |",
        "| --- | --- | --- | --- |",
    ]
    for head in report["evidence_heads"]:
        datasets = ", ".join(head["required_datasets"]) or "none"
        lines.append(
            f"| `{head['endpoint_id']}` | {head['status']} | "
            f"{head['evidence_class']} | {datasets} |"
        )

    lines.extend(["", "## Evidence-head details", ""])
    for head in report["evidence_heads"]:
        lines.extend(
            [
                f"### {head['endpoint_id']}",
                "",
                f"- Status: **{head['status']}**",
                f"- Source gate status: {head['source_gate_status']}",
                f"- Evidence class: {head['evidence_class']}",
                "- Dataset statuses: "
                + (
                    "; ".join(
                        f"{dataset_id}={status}"
                        for dataset_id, status in head["dataset_statuses"].items()
                    )
                    or "none"
                ),
                "- Blockers: " + ("; ".join(head["blockers"]) or "none"),
                f"- Next action: {head['next_action']}",
                f"- Supports: {head['claim_boundary']['supports']}",
                f"- Does not support: {head['claim_boundary']['does_not_support']}",
                "",
            ]
        )

    lines.extend(
        [
            "## Overall claim boundary",
            "",
            f"- Supported: {report['claim_boundary']['supported']}",
            f"- Unsupported: {report['claim_boundary']['unsupported']}",
            "",
            f"Next action: {report['next_action']}",
            "",
        ]
    )
    markdown_path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--readiness",
        type=Path,
        default=Path("validation_outputs/current/dataset_registry/dataset_evidence_readiness.json"),
    )
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("validation/qc_preprocess/acceptance_contract.json"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    report = build_report(args.readiness, args.contract)
    write_report(report, args.output_dir)
    print(
        json.dumps(
            {"status": report["status"], "output_dir": str(args.output_dir)},
            indent=2,
        )
    )
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
