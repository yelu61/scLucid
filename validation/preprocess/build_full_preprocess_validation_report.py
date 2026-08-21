#!/usr/bin/env python3
"""Build the complete Preprocess evidence-head readiness report without scoring.

The dataset readiness report records executed ``endpoint x dataset`` evidence,
while the acceptance contract defines which evidence heads and datasets are
required.  This report joins those two sources without averaging: every
contract-declared Preprocess head must have ``PASS`` or ``PASS_BASELINE``
evidence on every required dataset before the overall gate can pass.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

DEFAULT_OUTPUT_DIR = Path("validation_outputs/current/preprocess_full_gate")
PASSING_EVIDENCE = {"PASS", "PASS_BASELINE"}


HEAD_METADATA: dict[str, dict[str, str]] = {
    "pp_input_representation_contract": {
        "evidence_class": "ENGINEERING_CONTRACT",
        "supports": (
            "Preservation and explicit semantics of counts, normalized-full, discovery, "
            "and integrated representations for their declared consumers."
        ),
        "does_not_support": (
            "Biological superiority of any normalization, feature-selection, or integration method."
        ),
        "next_action": (
            "Run representation-contract checks on controlled fixtures and all registered real projects."
        ),
    },
    "pp_normalization_selection": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Consumer-aware selection among eligible normalization candidates on the registered tasks."
        ),
        "does_not_support": (
            "A universally best normalization method across assays, tissues, or downstream consumers."
        ),
        "next_action": (
            "Compare the simple log baseline and eligible alternatives on controlled and external datasets."
        ),
    },
    "pp_feature_selection": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Stable unsupervised feature selection that preserves registered identities and rare populations."
        ),
        "does_not_support": (
            "Treating selected features as biological truth or forcing prior markers into discovery features."
        ),
        "next_action": (
            "Evaluate batch-aware HVG and deviance candidates across controlled and tumor datasets."
        ),
    },
    "pp_selector_regret": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Held-out task utility of the selected policy within the preregistered candidate family."
        ),
        "does_not_support": (
            "Optimality outside the registered datasets, tasks, candidates, and regret definition."
        ),
        "next_action": (
            "Run leave-one-dataset-or-protocol-out selection and bind regret evidence for every required dataset."
        ),
    },
    "pp_graph_stability": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Neighborhood and discovery-structure stability across registered seeds and parameter variants."
        ),
        "does_not_support": (
            "Biological validity of every cluster or quantitative interpretation of an embedding."
        ),
        "next_action": (
            "Run graph-stability sensitivity analyses across the registered parameter and dataset panel."
        ),
    },
    "pp_integration_need_confounding": {
        "evidence_class": "SCIENTIFIC_SAFETY_GATE",
        "supports": (
            "Fail-closed integration decisions when correction is unnecessary or batch is confounded with biology."
        ),
        "does_not_support": (
            "Recovery of biological contrasts that are not identifiable from the study design."
        ),
        "next_action": (
            "Test no-integration selection and confounding blockers on controlled and real-project designs."
        ),
    },
    "pp_integration_pareto": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Integration choice only when it improves registered batch objectives without breaching biology guardrails."
        ),
        "does_not_support": (
            "A universal integration ranking or use of latent or graph output as corrected expression."
        ),
        "next_action": (
            "Compare every eligible integration candidate with the mandatory unintegrated baseline."
        ),
    },
    "pp_identity_preservation": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Preservation of externally known identities across registered protocol, technology, and depth shifts."
        ),
        "does_not_support": (
            "Accuracy of identities lacking independent truth or generalization to unregistered tissues."
        ),
        "next_action": (
            "Measure identity and rare-class preservation on all required external-truth datasets."
        ),
    },
    "pp_tumor_structure_preservation": {
        "evidence_class": "SCIENTIFIC_PERFORMANCE",
        "supports": (
            "Preservation of registered lineage, patient, rare-population, and tumor-program structure."
        ),
        "does_not_support": (
            "Mechanistic validity of tumor programs or universal benefit in tumor scRNA-seq."
        ),
        "next_action": (
            "Run locked tumor-structure endpoints on public cohorts and the real-project panel."
        ),
    },
    "pp_policy_execution": {
        "evidence_class": "ENGINEERING_CONTRACT",
        "supports": (
            "Deterministic, idempotent, count-preserving application of an explicitly reviewed policy."
        ),
        "does_not_support": "Scientific correctness of the policy that was applied.",
        "next_action": (
            "Run exact policy-agreement, idempotence, and count-preservation checks on required datasets."
        ),
    },
    "pp_decisioncard_ux": {
        "evidence_class": "UX_VALIDATION",
        "supports": (
            "Understandable next actions and reduced manual configuration in registered real projects."
        ),
        "does_not_support": "Scientific superiority of the recommended preprocessing policy.",
        "next_action": (
            "Complete real-project runs without manual schema bypasses or project-specific patches."
        ),
    },
    "pp_scalability": {
        "evidence_class": "ENGINEERING_ROBUSTNESS",
        "supports": (
            "Deterministic sparse execution, runtime, and peak-memory behavior on declared reference hardware."
        ),
        "does_not_support": (
            "A hardware-independent wall-time guarantee or scientific generalization from scale alone."
        ),
        "next_action": (
            "Run registered scale tiers and record runtime, peak memory, density, and repeatability."
        ),
    },
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _head_ids(design: Mapping[str, Any]) -> tuple[list[str], list[str]]:
    """Return ordered endpoint IDs and fail-closed contract issues."""
    raw_heads = design.get("evidence_heads", [])
    if not isinstance(raw_heads, list):
        return [], ["preprocess_validation_design.evidence_heads must be a list."]

    head_ids: list[str] = []
    issues: list[str] = []
    for index, raw_head in enumerate(raw_heads):
        if isinstance(raw_head, str):
            endpoint_id = raw_head
        elif isinstance(raw_head, Mapping):
            endpoint_id = str(raw_head.get("endpoint_id", ""))
        else:
            endpoint_id = ""
        if not endpoint_id:
            issues.append(f"Preprocess evidence head at index {index} has no endpoint_id.")
            continue
        if endpoint_id in head_ids:
            issues.append(f"Duplicate preprocess evidence head: {endpoint_id}.")
            continue
        head_ids.append(endpoint_id)
    if not head_ids:
        issues.append("No preprocess evidence heads are declared by the acceptance contract.")
    return head_ids, issues


def _contract_metadata(design: Mapping[str, Any], endpoint_id: str) -> dict[str, Any]:
    """Allow contract-owned metadata while retaining conservative defaults."""
    metadata: dict[str, Any] = {}
    raw_heads = design.get("evidence_heads", [])
    if isinstance(raw_heads, list):
        for raw_head in raw_heads:
            if isinstance(raw_head, Mapping) and raw_head.get("endpoint_id") == endpoint_id:
                metadata.update(raw_head)
                break
    for field in ("evidence_head_metadata", "evidence_head_definitions"):
        by_head = design.get(field, {})
        if isinstance(by_head, Mapping) and isinstance(by_head.get(endpoint_id), Mapping):
            metadata.update(by_head[endpoint_id])
    return metadata


def _metadata(design: Mapping[str, Any], endpoint_id: str) -> dict[str, str]:
    contract_metadata = _contract_metadata(design, endpoint_id)
    built_in = HEAD_METADATA.get(endpoint_id, {})
    claim_boundary = contract_metadata.get("claim_boundary", {})
    if not isinstance(claim_boundary, Mapping):
        claim_boundary = {}
    return {
        "evidence_class": str(
            contract_metadata.get("evidence_class", built_in.get("evidence_class", "UNCLASSIFIED"))
        ),
        "supports": str(
            claim_boundary.get(
                "supports",
                contract_metadata.get(
                    "supports",
                    built_in.get(
                        "supports", "Only the explicitly registered Preprocess endpoint estimand."
                    ),
                ),
            )
        ),
        "does_not_support": str(
            claim_boundary.get(
                "does_not_support",
                contract_metadata.get(
                    "does_not_support",
                    built_in.get(
                        "does_not_support",
                        "Any broader scientific, engineering, UX, or external-generalization claim.",
                    ),
                ),
            )
        ),
        "next_action": str(
            contract_metadata.get(
                "next_action",
                built_in.get(
                    "next_action",
                    (
                        "Classify this evidence head and generate PASS or PASS_BASELINE "
                        f"RunEvidence for {endpoint_id}."
                    ),
                ),
            )
        ),
    }


def build_report(
    readiness_path: Path,
    acceptance_contract_path: Path,
) -> dict[str, Any]:
    """Return head-level Preprocess readiness from exact required bindings."""
    readiness = _load_json(readiness_path)
    contract = _load_json(acceptance_contract_path)
    design = contract.get("preprocess_validation_design", {})
    if not isinstance(design, Mapping):
        design = {}
    head_ids, contract_issues = _head_ids(design)

    required_by_head = contract.get("required_endpoint_portfolio", {}).get("preprocess", {})
    if not isinstance(required_by_head, Mapping):
        required_by_head = {}
        contract_issues.append(
            "required_endpoint_portfolio.preprocess must be an endpoint-to-datasets mapping."
        )
    extra_required_heads = sorted(set(required_by_head) - set(head_ids))
    contract_issues.extend(
        f"Required preprocess endpoint is absent from evidence_heads: {endpoint_id}."
        for endpoint_id in extra_required_heads
    )

    source_gates = readiness.get("endpoint_gates", {}).get("preprocess", {})
    if not isinstance(source_gates, Mapping):
        source_gates = {}

    heads: list[dict[str, Any]] = []
    for endpoint_id in head_ids:
        source_gate = source_gates.get(endpoint_id, {})
        if not isinstance(source_gate, Mapping):
            source_gate = {}
        raw_required_datasets = required_by_head.get(endpoint_id, [])
        required_datasets = (
            list(raw_required_datasets) if isinstance(raw_required_datasets, list) else []
        )
        raw_statuses = source_gate.get("dataset_statuses", {})
        if not isinstance(raw_statuses, Mapping):
            raw_statuses = {}
        dataset_statuses = {
            dataset_id: str(raw_statuses.get(dataset_id, "NOT_RUN"))
            for dataset_id in required_datasets
        }

        blockers: list[str] = []
        if endpoint_id not in required_by_head:
            blockers.append(
                "The acceptance contract has no required dataset binding for this head."
            )
        elif not isinstance(raw_required_datasets, list):
            blockers.append("The required dataset binding for this head is not a list.")
        elif not required_datasets:
            blockers.append("The acceptance contract binds no required datasets to this head.")
        elif any(
            not isinstance(dataset_id, str) or not dataset_id for dataset_id in required_datasets
        ):
            blockers.append("Every required dataset ID must be a non-empty string.")
        elif len(required_datasets) != len(set(required_datasets)):
            blockers.append("The required dataset binding contains duplicates.")

        if endpoint_id not in source_gates:
            blockers.append("The dataset readiness report has no endpoint gate for this head.")
        else:
            source_required = source_gate.get("required_datasets")
            if isinstance(source_required, list) and source_required != required_datasets:
                blockers.append(
                    "The readiness endpoint gate was built from a different required dataset binding."
                )
            if source_gate.get("status") not in PASSING_EVIDENCE:
                blockers.append(
                    "The readiness endpoint gate is not PASS or PASS_BASELINE: "
                    f"{source_gate.get('status', 'NOT_RUN')}."
                )

        blockers.extend(
            f"{dataset_id}: {status}"
            for dataset_id, status in dataset_statuses.items()
            if status not in PASSING_EVIDENCE
        )
        metadata = _metadata(design, endpoint_id)
        heads.append(
            {
                "endpoint_id": endpoint_id,
                "status": "PASS" if required_datasets and not blockers else "BLOCKED",
                "source_gate_status": str(source_gate.get("status", "NOT_RUN")),
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
    overall_status = "PASS" if heads and not blocked_heads and not contract_issues else "BLOCKED"
    next_action = next(
        (head["next_action"] for head in heads if head["status"] != "PASS"),
        (
            "Resolve the acceptance-contract issues and rebuild readiness."
            if contract_issues
            else "All registered Preprocess evidence heads passed; preserve the evidence lock."
        ),
    )
    return {
        "schema_version": "sclucid_full_preprocess_validation_readiness_v1",
        "status": overall_status,
        "source_readiness": str(readiness_path),
        "source_acceptance_contract": str(acceptance_contract_path),
        "no_aggregate_quality_score": True,
        "passing_evidence_statuses": sorted(PASSING_EVIDENCE),
        "nonpassing_status_examples": [
            "NOT_EVALUABLE",
            "SIMULATION_PASS_NOT_EXTERNAL",
            "CONTRACT_PASS_NOT_PERFORMANCE",
            "REVIEW",
        ],
        "contract_issues": contract_issues,
        "head_count": len(heads),
        "passed_head_count": len(heads) - len(blocked_heads),
        "blocked_head_count": len(blocked_heads),
        "blocked_heads": blocked_heads,
        "evidence_heads": heads,
        "next_action": next_action,
        "claim_boundary": {
            "supported": (
                "Only heads marked PASS have PASS or PASS_BASELINE evidence for every "
                "dataset required by the locked acceptance contract."
            ),
            "unsupported": (
                "A controlled-dataset result, simulation, contract conformance, registry "
                "availability, or partial head coverage does not establish complete "
                "Preprocess scientific readiness or universal superiority."
            ),
        },
    }


def write_report(report: dict[str, Any], output_dir: Path) -> None:
    """Write the machine-readable gate and a compact human review."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "full_preprocess_validation_readiness.json"
    markdown_path = output_dir / "full_preprocess_validation_readiness.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# scLucid full Preprocess validation readiness",
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

    if report["contract_issues"]:
        lines.extend(["", "## Contract issues", ""])
        lines.extend(f"- {issue}" for issue in report["contract_issues"])

    lines.extend(["", "## Evidence-head details", ""])
    for head in report["evidence_heads"]:
        statuses = "; ".join(
            f"{dataset_id}={status}" for dataset_id, status in head["dataset_statuses"].items()
        )
        lines.extend(
            [
                f"### {head['endpoint_id']}",
                "",
                f"- Status: **{head['status']}**",
                f"- Source gate status: {head['source_gate_status']}",
                f"- Evidence class: {head['evidence_class']}",
                f"- Dataset statuses: {statuses or 'none'}",
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
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
