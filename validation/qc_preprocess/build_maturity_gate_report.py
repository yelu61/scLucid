#!/usr/bin/env python3
"""Combine QC, controlled preprocess, and real-project UX gates without scoring."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def build_report(
    qc_path: Path,
    preprocess_path: Path,
    ux_path: Path,
    portfolio_path: Path | None = None,
    full_qc_path: Path | None = None,
    full_preprocess_path: Path | None = None,
) -> dict:
    qc = json.loads(qc_path.read_text())
    preprocess = json.loads(preprocess_path.read_text())
    ux = json.loads(ux_path.read_text())
    portfolio = json.loads(portfolio_path.read_text()) if portfolio_path else {}
    full_qc = json.loads(full_qc_path.read_text()) if full_qc_path else {}
    full_preprocess = (
        json.loads(full_preprocess_path.read_text()) if full_preprocess_path else {}
    )
    portfolio_gates = portfolio.get("module_gates", {})
    qc_portfolio_ready = portfolio_gates.get("qc", {}).get("status") == "PASS"
    preprocess_portfolio_ready = portfolio_gates.get("preprocess", {}).get("status") == "PASS"
    analysis_portfolio_ready = portfolio_gates.get("analysis", {}).get("status") == "PASS"
    qc_ready = qc.get("status") == "PASS"
    full_qc_ready = full_qc.get("status") == "PASS" if full_qc_path else True
    full_preprocess_ready = (
        full_preprocess.get("status") == "PASS" if full_preprocess_path else True
    )
    preprocess_controlled_ready = preprocess.get("candidate_acceptance", {}).get(
        "status"
    ) == "PASS" and preprocess.get("integration_review", {}).get("status") in {
        "PASS_BASELINE",
        "PASS",
    }
    preprocess_release_ready = (
        preprocess_controlled_ready
        and preprocess.get("release_gate", {}).get("status") == "PASS"
        and preprocess_portfolio_ready
        and full_preprocess_ready
    )
    ux_ready = ux.get("status") == "PASS"
    qc_release_ready = qc_ready and qc_portfolio_ready and full_qc_ready
    overall = "READY" if qc_release_ready and preprocess_release_ready and ux_ready else "BLOCKED"
    return {
        "schema_version": "sclucid_qc_preprocess_maturity_gate_v1",
        "status": overall,
        "module_status": {
            "qc": "CORE" if overall == "READY" else "REVIEW",
            "preprocess": "CORE" if overall == "READY" else "REVIEW",
            "analysis": ("CORE" if overall == "READY" and analysis_portfolio_ready else "REVIEW"),
        },
        "downstream_feature_development": ("UNFROZEN" if overall == "READY" else "FROZEN"),
        "gates": {
            "qc_locked_scientific_acceptance": {
                "status": "PASS" if qc_release_ready else "BLOCKED",
                "locked_endpoint_status": qc.get("status", "NOT_RUN"),
                "dataset_coverage": qc.get("dataset_coverage", {}).get("status", "NOT_RUN"),
                "label_gate": qc.get("label_gate", {}).get("status", "NOT_RUN"),
                "source": str(qc_path),
            },
            "qc_full_head_acceptance": {
                "status": "PASS" if full_qc_ready else "BLOCKED",
                "blocked_heads": full_qc.get("blocked_heads", []),
                "source": str(full_qc_path) if full_qc_path else None,
                "compatibility_mode": full_qc_path is None,
            },
            "preprocess_controlled_acceptance": {
                "status": "PASS" if preprocess_controlled_ready else "BLOCKED",
                "candidate_acceptance": preprocess.get("candidate_acceptance", {}).get(
                    "status", "NOT_RUN"
                ),
                "integration_review": preprocess.get("integration_review", {}).get(
                    "status", "NOT_RUN"
                ),
                "source": str(preprocess_path),
            },
            "preprocess_full_head_acceptance": {
                "status": "PASS" if full_preprocess_ready else "BLOCKED",
                "blocked_heads": full_preprocess.get("blocked_heads", []),
                "source": str(full_preprocess_path) if full_preprocess_path else None,
                "compatibility_mode": full_preprocess_path is None,
            },
            "preprocess_external_release": {
                "status": "PASS" if preprocess_release_ready else "BLOCKED",
                "source": str(preprocess_path),
            },
            "real_project_ux": {
                "status": ux.get("status", "NOT_RUN"),
                "source": str(ux_path),
            },
            "analysis_scientific_acceptance": {
                "status": "PASS" if analysis_portfolio_ready else "BLOCKED",
                "source": str(portfolio_path) if portfolio_path else None,
            },
            "dataset_portfolio": {
                "status": (
                    "PASS" if qc_portfolio_ready and preprocess_portfolio_ready else "BLOCKED"
                ),
                "qc": portfolio_gates.get("qc", {}).get("status", "NOT_RUN"),
                "preprocess": portfolio_gates.get("preprocess", {}).get("status", "NOT_RUN"),
                "analysis": portfolio_gates.get("analysis", {}).get("status", "NOT_RUN"),
                "source": str(portfolio_path) if portfolio_path else None,
            },
        },
        "claim_boundary": {
            "supported": [
                "The controlled mixology selector and representation contracts were evaluated."
            ],
            "unsupported": [
                "QC superiority, external tumor preprocessing benefit, Analysis performance, and CORE release."
            ],
        },
        "prioritized_actions": [
            "Complete and freeze blinded QC expert labels.",
            "Execute the maintained four-action path in all three real projects and record UX evidence.",
            "Repeat preprocessing regret and Pareto review on external tumor projects.",
            "Run every required P0 portfolio endpoint; registration or downloadability is not a pass.",
            "Unfreeze Analysis/Tumor only if every gate passes.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--preprocess", type=Path, required=True)
    parser.add_argument("--ux", type=Path, required=True)
    parser.add_argument("--portfolio", type=Path, required=True)
    parser.add_argument(
        "--full-qc",
        type=Path,
        default=Path("validation_outputs/current/qc_full_gate/full_qc_validation_readiness.json"),
    )
    parser.add_argument(
        "--full-preprocess",
        type=Path,
        default=Path(
            "validation_outputs/current/preprocess_full_gate/"
            "full_preprocess_validation_readiness.json"
        ),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(
        args.qc,
        args.preprocess,
        args.ux,
        args.portfolio,
        args.full_qc,
        args.full_preprocess,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "maturity_gate.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# scLucid QC/Preprocess maturity gate",
        "",
        f"Status: **{report['status']}**",
        "",
        "No aggregate quality score is used.",
        "",
        "## Gates",
        "",
    ]
    lines.extend(f"- {name}: {payload['status']}" for name, payload in report["gates"].items())
    lines.extend(["", "## Prioritized actions", ""])
    lines.extend(
        f"{index}. {action}" for index, action in enumerate(report["prioritized_actions"], 1)
    )
    lines.append("")
    (args.output_dir / "maturity_gate.md").write_text("\n".join(lines))
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0 if report["status"] == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
