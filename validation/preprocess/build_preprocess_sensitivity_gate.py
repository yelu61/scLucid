#!/usr/bin/env python3
"""Aggregate preregistered preprocess sensitivity runs into one locked gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def build_gate(
    report_paths: list[Path],
    *,
    min_variants: int = 3,
    max_regret: float = 0.05,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for path in report_paths:
        payload = json.loads(path.read_text())
        candidate = payload.get("candidate_acceptance", {})
        dataset_rows = candidate.get("datasets", [])
        product_policy = payload.get("product_policy") or {}
        selected = product_policy.get("selected_candidate")
        if selected is None and len(dataset_rows) == 1:
            selected = dataset_rows[0].get("selected")
        regret = dataset_rows[0].get("regret") if len(dataset_rows) == 1 else None
        rows.append(
            {
                "source": str(path),
                "status": candidate.get("status", "NOT_RUN"),
                "selected_candidate": selected,
                "regret": regret,
                "integration_review": payload.get("integration_review", {}).get(
                    "status", "NOT_RUN"
                ),
            }
        )
    selected = {row["selected_candidate"] for row in rows if row["selected_candidate"]}
    enough_variants = len(rows) >= min_variants
    all_candidate_pass = bool(rows) and all(row["status"] == "PASS" for row in rows)
    stable_selection = len(selected) == 1 and all(
        row["selected_candidate"] for row in rows
    )
    all_regret_pass = bool(rows) and all(
        row["regret"] is not None and float(row["regret"]) <= max_regret for row in rows
    )
    passed = bool(
        enough_variants
        and all_candidate_pass
        and stable_selection
        and all_regret_pass
    )
    return {
        "schema_version": "sclucid_preprocess_sensitivity_gate_v1",
        "status": "PASS" if passed else "BLOCKED",
        "n_variants": len(rows),
        "selected_candidates": sorted(selected),
        "checks": {
            "minimum_variants": "PASS" if enough_variants else "BLOCKED",
            "candidate_acceptance": "PASS" if all_candidate_pass else "BLOCKED",
            "selection_stability": "PASS" if stable_selection else "BLOCKED",
            "regret_guardrail": "PASS" if all_regret_pass else "BLOCKED",
            "integration_guardrail": "SEPARATE_ENDPOINT",
        },
        "thresholds": {
            "min_variants": min_variants,
            "max_regret": max_regret,
            "selected_candidate_consistency_required": 1.0,
        },
        "runs": rows,
        "claim_boundary": (
            "Sensitivity PASS supports stability only across the registered mixology parameter grid."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-variants", type=int, default=3)
    parser.add_argument("--max-regret", type=float, default=0.05)
    args = parser.parse_args()
    gate = build_gate(
        args.report,
        min_variants=args.min_variants,
        max_regret=args.max_regret,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(gate, indent=2) + "\n")
    print(json.dumps({"status": gate["status"], "output": str(args.output)}, indent=2))
    return 0 if gate["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
