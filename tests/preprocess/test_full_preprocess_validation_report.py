from __future__ import annotations

import json

import pytest

from validation.preprocess.build_full_preprocess_validation_report import (
    build_report,
    write_report,
)

PREPROCESS_HEADS = [
    "pp_input_representation_contract",
    "pp_normalization_selection",
    "pp_feature_selection",
    "pp_selector_regret",
    "pp_graph_stability",
    "pp_integration_need_confounding",
    "pp_integration_pareto",
    "pp_identity_preservation",
    "pp_tumor_structure_preservation",
    "pp_policy_execution",
    "pp_decisioncard_ux",
    "pp_scalability",
]


def _write(path, payload):
    path.write_text(json.dumps(payload))
    return path


def _inputs(tmp_path, *, changed_head=None, changed_status="PASS"):
    required = {head: [f"dataset_{index}"] for index, head in enumerate(PREPROCESS_HEADS)}
    gates = {}
    for head, datasets in required.items():
        status = changed_status if head == changed_head else "PASS"
        gates[head] = {
            "status": "PASS" if status in {"PASS", "PASS_BASELINE"} else "BLOCKED",
            "required_datasets": datasets,
            "dataset_statuses": {datasets[0]: status},
        }
    contract = {
        "preprocess_validation_design": {"evidence_heads": PREPROCESS_HEADS},
        "required_endpoint_portfolio": {"preprocess": required},
        "portfolio_gate_policy": {
            "passing_evidence_statuses": [
                "PASS",
                "PASS_BASELINE",
                "CONTRACT_PASS_NOT_PERFORMANCE",
            ]
        },
    }
    readiness = {"endpoint_gates": {"preprocess": gates}}
    return (
        _write(tmp_path / "readiness.json", readiness),
        _write(tmp_path / "contract.json", contract),
    )


def test_all_contract_declared_heads_must_pass_without_averaging(tmp_path):
    readiness, contract = _inputs(tmp_path)

    report = build_report(readiness, contract)

    assert report["status"] == "PASS"
    assert report["head_count"] == len(PREPROCESS_HEADS)
    assert report["passed_head_count"] == len(PREPROCESS_HEADS)
    assert report["blocked_heads"] == []
    assert report["no_aggregate_quality_score"] is True
    assert report["passing_evidence_statuses"] == ["PASS", "PASS_BASELINE"]
    assert {head["evidence_class"] for head in report["evidence_heads"]} >= {
        "SCIENTIFIC_PERFORMANCE",
        "ENGINEERING_CONTRACT",
        "UX_VALIDATION",
    }


@pytest.mark.parametrize(
    "nonpassing_status",
    [
        "NOT_EVALUABLE",
        "SIMULATION_PASS_NOT_EXTERNAL",
        "CONTRACT_PASS_NOT_PERFORMANCE",
        "REVIEW",
    ],
)
def test_nonperformance_or_incomplete_evidence_cannot_pass(tmp_path, nonpassing_status):
    readiness, contract = _inputs(
        tmp_path,
        changed_head="pp_normalization_selection",
        changed_status=nonpassing_status,
    )

    report = build_report(readiness, contract)
    normalization = next(
        head
        for head in report["evidence_heads"]
        if head["endpoint_id"] == "pp_normalization_selection"
    )

    assert report["status"] == "BLOCKED"
    assert report["passed_head_count"] == len(PREPROCESS_HEADS) - 1
    assert report["blocked_heads"] == ["pp_normalization_selection"]
    assert normalization["status"] == "BLOCKED"
    assert normalization["blockers"] == [
        "The readiness endpoint gate is not PASS or PASS_BASELINE: BLOCKED.",
        f"dataset_1: {nonpassing_status}",
    ]


def test_acceptance_contract_binding_is_authoritative(tmp_path):
    readiness, contract = _inputs(tmp_path)
    readiness_payload = json.loads(readiness.read_text())
    gate = readiness_payload["endpoint_gates"]["preprocess"]["pp_selector_regret"]
    gate["required_datasets"] = ["stale_dataset"]
    readiness.write_text(json.dumps(readiness_payload))

    report = build_report(readiness, contract)
    selector = next(
        head for head in report["evidence_heads"] if head["endpoint_id"] == "pp_selector_regret"
    )

    assert report["status"] == "BLOCKED"
    assert selector["dataset_statuses"] == {"dataset_3": "PASS"}
    assert selector["blockers"] == [
        "The readiness endpoint gate was built from a different required dataset binding."
    ]


def test_contract_metadata_and_unknown_head_are_supported_dynamically(tmp_path):
    endpoint_id = "pp_future_consumer_contract"
    contract = {
        "preprocess_validation_design": {
            "evidence_heads": [endpoint_id],
            "evidence_head_metadata": {
                endpoint_id: {
                    "evidence_class": "ENGINEERING_CONTRACT",
                    "next_action": "Run the future consumer contract.",
                    "claim_boundary": {
                        "supports": "The future consumer contract only.",
                        "does_not_support": "Scientific performance.",
                    },
                }
            },
        },
        "required_endpoint_portfolio": {"preprocess": {endpoint_id: ["fixture"]}},
    }
    readiness = {
        "endpoint_gates": {
            "preprocess": {
                endpoint_id: {
                    "status": "PASS_BASELINE",
                    "required_datasets": ["fixture"],
                    "dataset_statuses": {"fixture": "PASS_BASELINE"},
                }
            }
        }
    }

    report = build_report(
        _write(tmp_path / "readiness.json", readiness),
        _write(tmp_path / "contract.json", contract),
    )
    head = report["evidence_heads"][0]

    assert report["status"] == "PASS"
    assert head["evidence_class"] == "ENGINEERING_CONTRACT"
    assert head["next_action"] == "Run the future consumer contract."
    assert head["claim_boundary"] == {
        "supports": "The future consumer contract only.",
        "does_not_support": "Scientific performance.",
    }


def test_missing_head_binding_or_extra_required_endpoint_fails_closed(tmp_path):
    readiness, contract = _inputs(tmp_path)
    contract_payload = json.loads(contract.read_text())
    del contract_payload["required_endpoint_portfolio"]["preprocess"]["pp_graph_stability"]
    contract_payload["required_endpoint_portfolio"]["preprocess"]["pp_unlisted"] = ["dataset_x"]
    contract.write_text(json.dumps(contract_payload))

    report = build_report(readiness, contract)
    graph = next(
        head for head in report["evidence_heads"] if head["endpoint_id"] == "pp_graph_stability"
    )

    assert report["status"] == "BLOCKED"
    assert graph["status"] == "BLOCKED"
    assert graph["blockers"][0] == (
        "The acceptance contract has no required dataset binding for this head."
    )
    assert report["contract_issues"] == [
        "Required preprocess endpoint is absent from evidence_heads: pp_unlisted."
    ]


def test_write_report_uses_locked_names_and_exposes_boundaries(tmp_path):
    readiness, contract = _inputs(
        tmp_path,
        changed_head="pp_tumor_structure_preservation",
        changed_status="REVIEW",
    )
    report = build_report(readiness, contract)
    output_dir = tmp_path / "out"

    write_report(report, output_dir)

    json_path = output_dir / "full_preprocess_validation_readiness.json"
    markdown_path = output_dir / "full_preprocess_validation_readiness.md"
    assert json_path.exists()
    assert markdown_path.exists()
    markdown = markdown_path.read_text()
    assert "No aggregate quality score is used" in markdown
    assert "pp_tumor_structure_preservation" in markdown
    assert "Does not support" in markdown
    assert "dataset_8=REVIEW" in markdown
