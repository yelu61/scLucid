#!/usr/bin/env python3
"""Build a use-tiered QC readiness gap report without an aggregate score."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

PASSING = {"PASS", "PASS_BASELINE"}
FILTERED_COUNT_HEADS = (
    "qc_input_contract",
    "qc_profile_selection",
    "qc_catastrophic_sample_detection",
    "qc_damage_classification",
    "qc_doublet_calibration",
    "qc_rare_population_preservation",
    "qc_selector_superiority",
    "qc_iterative_review",
    "qc_policy_execution",
    "qc_decisioncard_ux",
    "qc_scalability",
)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def summarize_bindings(
    evidence_heads: list[dict[str, Any]],
    *,
    endpoint_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Count exact endpoint-by-dataset statuses; this is coverage, not quality."""
    statuses: list[str] = []
    selected = []
    for head in evidence_heads:
        if endpoint_ids is not None and head["endpoint_id"] not in endpoint_ids:
            continue
        selected.append(head)
        statuses.extend(map(str, head.get("dataset_statuses", {}).values()))
    counts = Counter(statuses)
    return {
        "head_count": len(selected),
        "passed_head_count": sum(head.get("status") == "PASS" for head in selected),
        "binding_count": len(statuses),
        "passing_binding_count": sum(counts[status] for status in PASSING),
        "status_counts": dict(sorted(counts.items())),
        "interpretation": "Coverage count only; it is not an aggregate quality score.",
    }


def _label_completion(path: Path) -> dict[str, int]:
    labels = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    complete = int((labels["expert_label"].str.strip() != "").sum())
    return {
        "total": int(len(labels)),
        "completed": complete,
        "remaining": int(len(labels) - complete),
    }


def _head_map(full_qc: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {head["endpoint_id"]: head for head in full_qc["evidence_heads"]}


def _tier(
    tier_id: str,
    label: str,
    required_heads: tuple[str, ...],
    heads: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    statuses = {head: heads.get(head, {}).get("status", "NOT_RUN") for head in required_heads}
    blocked = [head for head, status in statuses.items() if status != "PASS"]
    return {
        "tier_id": tier_id,
        "label": label,
        "status": "PASS" if not blocked else "BLOCKED",
        "required_heads": list(required_heads),
        "head_statuses": statuses,
        "passed_heads": len(required_heads) - len(blocked),
        "required_head_count": len(required_heads),
        "blocked_heads": blocked,
    }


def build_report(
    full_qc_path: Path,
    concordance_path: Path,
    ux_path: Path,
    sample_labels_path: Path,
    cell_labels_path: Path,
    ambient_path: Path,
    damage_path: Path,
    doublet_path: Path,
) -> dict[str, Any]:
    """Return the exact gap between current QC evidence and declared use tiers."""
    full_qc = _load_json(full_qc_path)
    concordance = _load_json(concordance_path)
    ux = _load_json(ux_path)
    ambient = _load_json(ambient_path)
    damage = _load_json(damage_path)
    doublet = _load_json(doublet_path)
    heads = _head_map(full_qc)

    input_contract_ready = (
        heads["qc_input_contract"]["dataset_statuses"].get("controlled_qc_truth_suite") == "PASS"
    )
    policy_execution_ready = (
        heads["qc_policy_execution"]["dataset_statuses"].get("controlled_qc_truth_suite") == "PASS"
    )
    controlled_ready = input_contract_ready and policy_execution_ready
    projects_ran = concordance.get("datasets_run") == concordance.get(
        "datasets_expected"
    ) and not concordance.get("blockers")
    review_only_status = (
        "AVAILABLE_WITH_GUARDRAILS" if controlled_ready and projects_ran else "BLOCKED"
    )

    filtered_tier = _tier(
        "trusted_filtered_count_qc",
        "Filtered-count QC recommendations trusted after registered validation",
        FILTERED_COUNT_HEADS,
        heads,
    )
    full_tier = _tier(
        "full_qc_core",
        "Full QC CORE including raw-droplet cell calling and ambient correction",
        tuple(heads),
        heads,
    )

    historical_rows = concordance.get("results", [])
    historical_coverages = [
        float(row["historical_removed_flagged_fraction"]) for row in historical_rows
    ]
    damage_endpoint = damage["endpoints"]["qc_damage_classification"]
    damage_metrics = damage_endpoint["aggregate_metrics"]
    report = {
        "schema_version": "sclucid_qc_real_project_readiness_gap_v1",
        "status": "BLOCKED" if full_qc.get("status") != "PASS" else "PASS",
        "current_recommendation": (
            "FULL_QC_CORE"
            if full_qc.get("status") == "PASS"
            else "SUPERVISED_REVIEW_ONLY"
            if controlled_ready and projects_ran
            else "NOT_READY"
        ),
        "no_aggregate_quality_score": True,
        "distance": {
            "all_qc": summarize_bindings(full_qc["evidence_heads"]),
            "filtered_count_qc": summarize_bindings(
                full_qc["evidence_heads"], endpoint_ids=set(FILTERED_COUNT_HEADS)
            ),
            "sample_truth_labels": _label_completion(sample_labels_path),
            "cell_truth_labels": _label_completion(cell_labels_path),
        },
        "use_tiers": [
            {
                "tier_id": "engineering_execution",
                "label": "Read-only review and explicit policy execution contract",
                "status": "PASS_LIMITED" if controlled_ready else "BLOCKED",
                "criteria": {
                    "controlled_input_contract": input_contract_ready,
                    "controlled_policy_execution": policy_execution_ready,
                },
                "claim_boundary": "Engineering behavior only; policy correctness is not established.",
            },
            {
                "tier_id": "supervised_real_project_review",
                "label": "Real-project use with expert review and no unattended deletion",
                "status": review_only_status,
                "criteria": {
                    "all_registered_objects_reviewed": projects_ran,
                    "input_is_not_mutated_by_review": controlled_ready,
                    "scientific_truth_complete": False,
                    "real_project_ux_gate": ux.get("status") == "PASS",
                },
                "claim_boundary": (
                    "Permits a reviewer to inspect and explicitly approve decisions; "
                    "does not permit unattended filtering or superiority claims."
                ),
            },
            filtered_tier,
            full_tier,
        ],
        "decisive_scientific_results": [
            {
                "endpoint": "qc_ambient_correction",
                "dataset": "tenx_hgmm_6k",
                "status": ambient["status"],
                "metrics": ambient["metrics"],
                "interpretation": "The current linear fallback removed no measured cross-species ambient signal.",
            },
            {
                "endpoint": "qc_damage_classification",
                "dataset": "emtab2600_microscopy_quality",
                "status": damage_endpoint["status"],
                "metrics": {
                    "n_damaged_evaluated": damage_metrics["n_damaged_evaluated"],
                    "damaged_cell_recall": damage_metrics["damaged_cell_recall"],
                    "keep_false_removal_rate": damage_metrics["keep_false_removal_rate"],
                },
                "interpretation": "The conservative policy protected KEEP cells but detected none of the evaluated damaged cells.",
            },
            {
                "endpoint": "qc_doublet_calibration",
                "dataset": "tenx_hgmm_6k",
                "status": doublet["status"],
                "metrics": {
                    "selected_auprc": doublet["selected_auprc"],
                    "best_auprc": doublet["best_auprc"],
                    "auprc_regret": doublet["auprc_regret"],
                },
                "interpretation": "A limited cross-species doublet result passed, without supporting automatic deletion.",
            },
        ],
        "real_project_operational_reference": {
            "status": concordance.get("status"),
            "datasets_run": concordance.get("datasets_run"),
            "historical_removed_flagged_fraction_range": (
                [min(historical_coverages), max(historical_coverages)]
                if historical_coverages
                else []
            ),
            "results": historical_rows,
            "interpretation": (
                "The current reviewer is substantially more conservative than historical Step1 membership. "
                "Only blinded labels can determine whether this is beneficial protection or missed low-quality cells."
            ),
        },
        "safe_now": {
            "allowed": [
                "Run the read-only DecisionCard on a filtered count matrix.",
                "Inspect per-sample summaries and all REVIEW/REMOVE cells with an expert.",
                "Apply an explicitly approved policy while retaining counts and RunEvidence.",
                "Use standard project-specific QC as the decision authority when scLucid evidence conflicts.",
            ],
            "not_allowed": [
                "Unattended automatic cell deletion based on the current selector.",
                "Automatic ambient correction with the current linear fallback.",
                "Automatic doublet deletion from one method or score.",
                "Claims that scLucid QC is superior to global thresholds or per-sample MAD.",
                "Treating filtered matrices as evidence for cell-calling or ambient performance.",
            ],
        },
        "shortest_path": {
            "filtered_count_qc": [
                "Complete and lock the 46 sample and 4,195 cell expert labels without viewing method predictions.",
                "Fix damage/low-quality recall, then pass the <=2% KEEP false-removal and >=5-point gain gates.",
                "Adjudicate real-project conflicts and validate rare-lineage/tumor-population preservation.",
                "Demonstrate evidence-dependent multi-round review and selector superiority with grouped bootstrap.",
                "Complete real-project application, idempotence, UX, and production-scale evidence.",
            ],
            "full_raw_droplet_qc_additions": [
                "Supply valid unfiltered-droplet cell-calling truth and pass recall/FDR gates.",
                "Replace or disable the failed ambient fallback and validate a model-based backend on external truth.",
            ],
        },
        "claim_boundary": (
            "scLucid QC is currently suitable as a conservative, supervised review aid. "
            "It is not yet validated as an autonomous or superior QC decision system."
        ),
    }
    return report


def _write_markdown(report: dict[str, Any], path: Path) -> None:
    distance = report["distance"]
    lines = [
        "# QC real-project readiness gap",
        "",
        f"Status: **{report['status']}**",
        "",
        report["claim_boundary"],
        "",
        "## Exact distance (coverage, not a score)",
        "",
        (
            f"- Full QC: {distance['all_qc']['passed_head_count']}/"
            f"{distance['all_qc']['head_count']} complete heads; "
            f"{distance['all_qc']['passing_binding_count']}/"
            f"{distance['all_qc']['binding_count']} exact bindings pass."
        ),
        (
            f"- Filtered-count QC: {distance['filtered_count_qc']['passed_head_count']}/"
            f"{distance['filtered_count_qc']['head_count']} complete heads; "
            f"{distance['filtered_count_qc']['passing_binding_count']}/"
            f"{distance['filtered_count_qc']['binding_count']} exact bindings pass."
        ),
        (
            f"- Blinded labels remaining: {distance['sample_truth_labels']['remaining']} sample, "
            f"{distance['cell_truth_labels']['remaining']} cell."
        ),
        "",
        "## Use tiers",
        "",
    ]
    lines.extend(f"- {tier['label']}: **{tier['status']}**" for tier in report["use_tiers"])
    lines.extend(["", "## Safe current boundary", ""])
    lines.extend(f"- Allowed: {item}" for item in report["safe_now"]["allowed"])
    lines.extend(f"- Not allowed: {item}" for item in report["safe_now"]["not_allowed"])
    lines.extend(["", "## Shortest path", ""])
    lines.extend(
        f"{idx}. {item}"
        for idx, item in enumerate(report["shortest_path"]["filtered_count_qc"], start=1)
    )
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-qc",
        type=Path,
        default=Path("validation_outputs/current/qc_full_gate/full_qc_validation_readiness.json"),
    )
    parser.add_argument(
        "--concordance",
        type=Path,
        default=Path(
            "validation_outputs/current/qc_real_projects/real_project_qc_concordance.json"
        ),
    )
    parser.add_argument(
        "--ux",
        type=Path,
        default=Path("validation_outputs/current/real_project_ux/real_project_ux_acceptance.json"),
    )
    parser.add_argument(
        "--sample-labels",
        type=Path,
        default=Path("validation_outputs/current/qc_truth_pack/reviewer/sample_labels.tsv"),
    )
    parser.add_argument(
        "--cell-labels",
        type=Path,
        default=Path("validation_outputs/current/qc_truth_pack/reviewer/cell_labels.tsv"),
    )
    parser.add_argument(
        "--ambient",
        type=Path,
        default=Path("validation_outputs/current/qc_hgmm/hgmm_ambient_truth_benchmark.json"),
    )
    parser.add_argument(
        "--damage",
        type=Path,
        default=Path(
            "validation_outputs/current/qc_emtab2600/cell_calling_damage_truth_evaluation.json"
        ),
    )
    parser.add_argument(
        "--doublet",
        type=Path,
        default=Path("validation_outputs/current/qc_hgmm/hgmm_doublet_truth_benchmark.json"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/current/qc_real_project_gap")
    )
    args = parser.parse_args()
    report = build_report(
        args.full_qc,
        args.concordance,
        args.ux,
        args.sample_labels,
        args.cell_labels,
        args.ambient,
        args.damage,
        args.doublet,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "qc_real_project_readiness_gap.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    _write_markdown(report, args.output_dir / "qc_real_project_readiness_gap.md")
    print(
        json.dumps(
            {
                "status": report["status"],
                "current_recommendation": report["current_recommendation"],
            },
            indent=2,
        )
    )
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
