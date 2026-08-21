from __future__ import annotations

import json

from validation.preprocess.build_preprocess_sensitivity_gate import build_gate


def _report(path, *, selected="standard_unintegrated", regret=0.01, status="PASS"):
    path.write_text(
        json.dumps(
            {
                "candidate_acceptance": {
                    "status": status,
                    "datasets": [{"regret": regret}],
                },
                "product_policy": {"selected_candidate": selected},
                "integration_review": {"status": "PASS_BASELINE"},
            }
        )
    )
    return path


def test_sensitivity_gate_requires_three_consistent_passing_variants(tmp_path):
    reports = [
        _report(tmp_path / f"run_{index}.json", regret=0.01 * index)
        for index in range(1, 4)
    ]

    gate = build_gate(reports)

    assert gate["status"] == "PASS"
    assert gate["selected_candidates"] == ["standard_unintegrated"]


def test_sensitivity_gate_blocks_selector_instability(tmp_path):
    reports = [
        _report(tmp_path / "a.json"),
        _report(tmp_path / "b.json"),
        _report(tmp_path / "c.json", selected="complex"),
    ]

    gate = build_gate(reports)

    assert gate["status"] == "BLOCKED"
    assert gate["checks"]["selection_stability"] == "BLOCKED"
