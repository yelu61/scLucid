from __future__ import annotations

import json

import pytest

from validation.qc.build_full_qc_validation_report import build_report, write_report

QC_HEADS = [
    "qc_input_contract",
    "qc_profile_selection",
    "qc_cell_calling",
    "qc_ambient_correction",
    "qc_catastrophic_sample_detection",
    "qc_damage_classification",
    "qc_doublet_calibration",
    "qc_rare_population_preservation",
    "qc_selector_superiority",
    "qc_iterative_review",
    "qc_policy_execution",
    "qc_decisioncard_ux",
    "qc_scalability",
]


def _write(path, payload):
    path.write_text(json.dumps(payload))
    return path


def _inputs(tmp_path, *, changed_head=None, changed_status="PASS"):
    required = {head: [f"dataset_{index}"] for index, head in enumerate(QC_HEADS)}
    gates = {}
    for head, datasets in required.items():
        status = changed_status if head == changed_head else "PASS"
        gates[head] = {
            "status": "PASS" if status in {"PASS", "PASS_BASELINE"} else "BLOCKED",
            "required_datasets": datasets,
            "dataset_statuses": {datasets[0]: status},
        }
    contract = {
        "qc_validation_design": {"evidence_heads": QC_HEADS},
        "required_endpoint_portfolio": {"qc": required},
        "portfolio_gate_policy": {"passing_evidence_statuses": ["PASS", "PASS_BASELINE"]},
    }
    readiness = {"endpoint_gates": {"qc": gates}}
    return (
        _write(tmp_path / "readiness.json", readiness),
        _write(tmp_path / "contract.json", contract),
    )


def test_all_thirteen_heads_must_pass_without_averaging(tmp_path):
    readiness, contract = _inputs(tmp_path)
    report = build_report(readiness, contract)

    assert report["status"] == "PASS"
    assert report["head_count"] == 13
    assert report["passed_head_count"] == 13
    assert report["blocked_heads"] == []
    assert report["no_aggregate_quality_score"] is True
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
    ],
)
def test_partial_or_nonexternal_evidence_cannot_pass_a_head(tmp_path, nonpassing_status):
    readiness, contract = _inputs(
        tmp_path,
        changed_head="qc_ambient_correction",
        changed_status=nonpassing_status,
    )

    report = build_report(readiness, contract)
    ambient = next(
        head for head in report["evidence_heads"] if head["endpoint_id"] == "qc_ambient_correction"
    )

    assert report["status"] == "BLOCKED"
    assert report["passed_head_count"] == 12
    assert report["blocked_heads"] == ["qc_ambient_correction"]
    assert ambient["status"] == "BLOCKED"
    assert ambient["blockers"] == [f"dataset_3: {nonpassing_status}"]


def test_contract_binding_is_authoritative_and_missing_source_gate_blocks(tmp_path):
    readiness, contract = _inputs(tmp_path)
    payload = json.loads(readiness.read_text())
    del payload["endpoint_gates"]["qc"]["qc_policy_execution"]
    readiness.write_text(json.dumps(payload))

    report = build_report(readiness, contract)
    policy = next(
        head for head in report["evidence_heads"] if head["endpoint_id"] == "qc_policy_execution"
    )

    assert report["status"] == "BLOCKED"
    assert policy["dataset_statuses"] == {"dataset_10": "NOT_RUN"}
    assert policy["blockers"] == [
        "The dataset readiness report has no endpoint gate for this head.",
        "dataset_10: NOT_RUN",
    ]


def test_write_report_uses_locked_output_names_and_exposes_boundaries(tmp_path):
    readiness, contract = _inputs(
        tmp_path,
        changed_head="qc_doublet_calibration",
        changed_status="REVIEW",
    )
    report = build_report(readiness, contract)
    output_dir = tmp_path / "out"

    write_report(report, output_dir)

    json_path = output_dir / "full_qc_validation_readiness.json"
    markdown_path = output_dir / "full_qc_validation_readiness.md"
    assert json_path.exists()
    assert markdown_path.exists()
    markdown = markdown_path.read_text()
    assert "No aggregate quality score is used" in markdown
    assert "qc_doublet_calibration" in markdown
    assert "Does not support" in markdown
    assert "dataset_6=REVIEW" in markdown
