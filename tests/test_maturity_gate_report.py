from __future__ import annotations

import json

from validation.qc_preprocess.build_maturity_gate_report import build_report


def _write(path, payload):
    path.write_text(json.dumps(payload))
    return path


def test_maturity_gate_keeps_downstream_frozen_until_all_inputs_pass(tmp_path):
    qc = _write(
        tmp_path / "qc.json",
        {
            "status": "BLOCKED",
            "dataset_coverage": {"status": "PASS"},
            "label_gate": {"status": "BLOCKED"},
        },
    )
    preprocess = _write(
        tmp_path / "preprocess.json",
        {
            "candidate_acceptance": {"status": "PASS"},
            "integration_review": {"status": "PASS_BASELINE"},
            "release_gate": {"status": "BLOCKED"},
        },
    )
    ux = _write(tmp_path / "ux.json", {"status": "BLOCKED"})

    report = build_report(qc, preprocess, ux)

    assert report["status"] == "BLOCKED"
    assert report["gates"]["preprocess_controlled_acceptance"]["status"] == "PASS"
    assert report["gates"]["preprocess_external_release"]["status"] == "BLOCKED"
    assert report["downstream_feature_development"] == "FROZEN"


def test_maturity_gate_requires_executed_dataset_portfolio(tmp_path):
    qc = _write(
        tmp_path / "qc.json",
        {
            "status": "PASS",
            "dataset_coverage": {"status": "PASS"},
            "label_gate": {"status": "PASS"},
        },
    )
    preprocess = _write(
        tmp_path / "preprocess.json",
        {
            "candidate_acceptance": {"status": "PASS"},
            "integration_review": {"status": "PASS_BASELINE"},
            "release_gate": {"status": "PASS"},
        },
    )
    ux = _write(tmp_path / "ux.json", {"status": "PASS"})
    portfolio = _write(
        tmp_path / "portfolio.json",
        {
            "module_gates": {
                "qc": {"status": "PASS"},
                "preprocess": {"status": "BLOCKED"},
                "analysis": {"status": "BLOCKED"},
            }
        },
    )

    report = build_report(qc, preprocess, ux, portfolio)

    assert report["status"] == "BLOCKED"
    assert report["gates"]["dataset_portfolio"]["preprocess"] == "BLOCKED"
    assert report["gates"]["preprocess_external_release"]["status"] == "BLOCKED"
    assert report["module_status"]["analysis"] == "REVIEW"
    assert report["gates"]["analysis_scientific_acceptance"]["status"] == "BLOCKED"


def test_maturity_gate_requires_every_qc_evidence_head_when_report_is_provided(tmp_path):
    qc = _write(
        tmp_path / "qc.json",
        {
            "status": "PASS",
            "dataset_coverage": {"status": "PASS"},
            "label_gate": {"status": "PASS"},
        },
    )
    preprocess = _write(
        tmp_path / "preprocess.json",
        {
            "candidate_acceptance": {"status": "PASS"},
            "integration_review": {"status": "PASS_BASELINE"},
            "release_gate": {"status": "PASS"},
        },
    )
    ux = _write(tmp_path / "ux.json", {"status": "PASS"})
    portfolio = _write(
        tmp_path / "portfolio.json",
        {
            "module_gates": {
                "qc": {"status": "PASS"},
                "preprocess": {"status": "PASS"},
                "analysis": {"status": "PASS"},
            }
        },
    )
    full_qc = _write(
        tmp_path / "full_qc.json",
        {"status": "BLOCKED", "blocked_heads": ["qc_ambient_correction"]},
    )

    report = build_report(qc, preprocess, ux, portfolio, full_qc)

    assert report["status"] == "BLOCKED"
    assert report["gates"]["qc_full_head_acceptance"]["status"] == "BLOCKED"
    assert report["gates"]["qc_full_head_acceptance"]["blocked_heads"] == [
        "qc_ambient_correction"
    ]
    assert report["module_status"]["qc"] == "REVIEW"


def test_maturity_gate_requires_every_preprocess_head_when_report_is_provided(tmp_path):
    qc = _write(
        tmp_path / "qc.json",
        {
            "status": "PASS",
            "dataset_coverage": {"status": "PASS"},
            "label_gate": {"status": "PASS"},
        },
    )
    preprocess = _write(
        tmp_path / "preprocess.json",
        {
            "candidate_acceptance": {"status": "PASS"},
            "integration_review": {"status": "PASS_BASELINE"},
            "release_gate": {"status": "PASS"},
        },
    )
    ux = _write(tmp_path / "ux.json", {"status": "PASS"})
    portfolio = _write(
        tmp_path / "portfolio.json",
        {
            "module_gates": {
                "qc": {"status": "PASS"},
                "preprocess": {"status": "PASS"},
                "analysis": {"status": "PASS"},
            }
        },
    )
    full_preprocess = _write(
        tmp_path / "full_preprocess.json",
        {"status": "BLOCKED", "blocked_heads": ["pp_tumor_structure_preservation"]},
    )

    report = build_report(
        qc,
        preprocess,
        ux,
        portfolio,
        None,
        full_preprocess,
    )

    assert report["status"] == "BLOCKED"
    assert report["gates"]["preprocess_full_head_acceptance"]["status"] == "BLOCKED"
    assert report["gates"]["preprocess_full_head_acceptance"]["blocked_heads"] == [
        "pp_tumor_structure_preservation"
    ]
    assert report["module_status"]["preprocess"] == "REVIEW"
