#!/usr/bin/env python3
"""Evaluate or initialize the locked real-project QC/Preprocess UX gate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.qc_preprocess.locked_acceptance import evaluate_real_project_ux_acceptance

TEMPLATE_COLUMNS = [
    "project",
    "legacy_config_fields",
    "current_config_fields",
    "manual_predicted_doublet_deletion",
    "manual_review_summary_edit",
    "schema_bypass",
    "project_specific_patch_count",
    "run_evidence_status",
    "run_evidence_path",
    "reviewer",
    "notes",
]


def _template(projects: list[str]) -> pd.DataFrame:
    return pd.DataFrame([{column: project if column == "project" else "" for column in TEMPLATE_COLUMNS} for project in projects])


def run(input_path: Path | None, output_dir: Path, contract_path: Path) -> dict:
    contract = json.loads(contract_path.read_text())
    projects = list(map(str, contract["external_projects"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    template_path = output_dir / "real_project_ux_records.tsv"
    if input_path is None:
        records = _template(projects)
        records.to_csv(template_path, sep="\t", index=False)
    else:
        records = pd.read_csv(input_path, sep="\t", dtype=str).fillna("")
        records.to_csv(template_path, sep="\t", index=False)

    result = evaluate_real_project_ux_acceptance(
        records,
        expected_projects=projects,
        min_config_reduction=float(contract["ux_primary_endpoint"]["min_config_field_reduction"]),
    )
    report = {
        "schema_version": "sclucid_real_project_ux_acceptance_v1",
        **result,
        "input_status": "RECORDED" if input_path is not None else "TEMPLATE_ONLY",
        "next_action": (
            "No action; the recorded projects passed the UX gate."
            if result["status"] == "PASS"
            else "Execute the maintained four-action path in each project and complete the evidence record."
        ),
    }
    (output_dir / "real_project_ux_acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Real-project QC/Preprocess UX acceptance",
        "",
        f"Status: **{report['status']}**",
        "",
        report["claim_boundary"],
        "",
        "## Projects",
        "",
    ]
    lines.extend(
        f"- {row['project']}: {row['status']} — {row.get('reason', 'record evaluated')}"
        for row in report["projects"]
    )
    lines.extend(["", "## Next action", "", report["next_action"], ""])
    (output_dir / "real_project_ux_acceptance.md").write_text("\n".join(lines))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("validation/qc_preprocess/acceptance_contract.json"),
    )
    args = parser.parse_args()
    report = run(args.input, args.output_dir, args.contract)
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
