from __future__ import annotations

import json
from pathlib import Path

import pytest

from validation.build_dataset_evidence_registry_report import (
    build_report,
    validate_registry,
    validate_required_endpoint_portfolio,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO_ROOT / "validation" / "dataset_evidence_registry.json"
CONTRACT_PATH = REPO_ROOT / "validation" / "qc_preprocess" / "acceptance_contract.json"


def test_dataset_registry_has_complete_machine_readable_contract():
    registry = json.loads(REGISTRY_PATH.read_text())
    contract = json.loads(CONTRACT_PATH.read_text())

    assert validate_registry(registry) == []
    assert validate_required_endpoint_portfolio(
        registry, contract["required_endpoint_portfolio"]
    ) == []
    dataset_ids = {row["dataset_id"] for row in registry["datasets"]}
    for required in contract["required_dataset_portfolio"].values():
        assert set(required) <= dataset_ids

    future = {
        row["dataset_id"]
        for row in registry["datasets"]
        if row["release_scope"].startswith("future_")
    }
    assert future
    assert not future.intersection(contract["required_dataset_portfolio"]["analysis"])


def test_readiness_report_does_not_confuse_downloadable_with_scientific_pass(tmp_path):
    registry = {
        "schema_name": "scLucidDatasetEvidenceRegistry",
        "schema_version": "1.0.0",
        "endpoint_definitions": {
            "e1": {
                "module": "qc",
                "estimand": "truth",
                "experimental_unit": "library",
                "metrics": ["recall"],
                "acceptance": {"recall_min": 1.0},
            }
        },
        "datasets": [
            {
                "dataset_id": "downloadable",
                "priority": "P0",
                "release_scope": "current_core",
                "modules": ["qc"],
                "truth_types": ["GT-A"],
                "accessions": [{"id": "A1", "url": "https://example.org"}],
                "download": {
                    "raw_reads": "AVAILABLE",
                    "raw_counts": "AVAILABLE",
                    "processed": "AVAILABLE",
                },
                "license": {
                    "status": "REVIEW",
                    "identifier": "review",
                    "url": "https://example.org",
                    "redistribution": "BLOCKED",
                },
                "local_path": None,
                "required_metadata": ["library"],
                "endpoint_ids": ["e1"],
            }
        ],
    }
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(registry))
    evidence_path = tmp_path / "evidence.json"
    evidence_path.write_text(json.dumps({"runs": []}))

    report, rows = build_report(
        registry_path,
        evidence_path,
        repo_root=tmp_path,
        required_portfolio={"qc": ["downloadable"]},
    )

    assert report["registry_schema_status"] == "PASS"
    assert rows[0]["acquisition_status"] == "DOWNLOADABLE"
    assert rows[0]["scientific_status"] == "NOT_RUN"
    assert report["module_gates"]["qc"]["status"] == "BLOCKED"


def test_module_gate_ignores_unrun_endpoints_from_other_modules(tmp_path):
    registry = {
        "schema_name": "scLucidDatasetEvidenceRegistry",
        "endpoint_definitions": {
            "pp": {
                "module": "preprocess",
                "estimand": "truth",
                "experimental_unit": "sample",
                "metrics": ["utility"],
                "acceptance": {"regret_max": 0.05},
            },
            "an": {
                "module": "analysis",
                "estimand": "truth",
                "experimental_unit": "sample",
                "metrics": ["accuracy"],
                "acceptance": {"regret_max": 0.05},
            },
        },
        "datasets": [
            {
                "dataset_id": "multi",
                "priority": "P0",
                "release_scope": "test",
                "modules": ["preprocess", "analysis"],
                "truth_types": ["GT-B"],
                "accessions": [{"id": "X"}],
                "download": {
                    "raw_reads": "AVAILABLE",
                    "raw_counts": "AVAILABLE",
                    "processed": "AVAILABLE",
                },
                "license": {
                    "status": "REVIEW",
                    "identifier": "test",
                    "url": "https://example.org",
                    "redistribution": "NO",
                },
                "required_metadata": ["sample"],
                "endpoint_ids": ["pp", "an"],
            }
        ],
    }
    (tmp_path / "registry.json").write_text(json.dumps(registry))
    (tmp_path / "index.json").write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "dataset_id": "multi",
                        "endpoint_id": "pp",
                        "status": "PASS",
                        "artifact": "artifact.json",
                    },
                    {
                        "dataset_id": "multi",
                        "endpoint_id": "an",
                        "status": "REVIEW",
                        "artifact": "artifact.json",
                    },
                ]
            }
        )
    )
    (tmp_path / "artifact.json").write_text("{}")

    report, rows = build_report(
        tmp_path / "registry.json",
        tmp_path / "index.json",
        repo_root=tmp_path,
        required_portfolio={"preprocess": ["multi"]},
    )

    assert rows[0]["preprocess_status"] == "PASS"
    assert rows[0]["analysis_status"] == "REVIEW"
    assert report["module_gates"]["preprocess"]["status"] == "PASS"


def test_exact_endpoint_gate_does_not_require_every_declared_dataset_endpoint(tmp_path):
    registry = {
        "schema_name": "scLucidDatasetEvidenceRegistry",
        "endpoint_definitions": {
            endpoint: {
                "module": "qc",
                "estimand": endpoint,
                "experimental_unit": "library",
                "metrics": ["metric"],
                "acceptance": {"metric_min": 1.0},
            }
            for endpoint in ("input", "scientific")
        },
        "datasets": [
            {
                "dataset_id": "fixture",
                "priority": "P0",
                "release_scope": "test",
                "modules": ["qc"],
                "truth_types": ["ENG"],
                "accessions": [{"id": "generated"}],
                "download": {
                    "raw_reads": "NOT_APPLICABLE",
                    "raw_counts": "GENERATED_LOCAL",
                    "processed": "GENERATED_LOCAL",
                },
                "license": {
                    "status": "GENERATED_LOCAL",
                    "identifier": "test",
                    "url": "LOCAL",
                    "redistribution": "YES",
                },
                "required_metadata": ["library"],
                "endpoint_ids": ["input", "scientific"],
            }
        ],
    }
    (tmp_path / "registry.json").write_text(json.dumps(registry))
    (tmp_path / "index.json").write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "dataset_id": "fixture",
                        "endpoint_id": "input",
                        "status": "PASS",
                        "artifact": "artifact.json",
                    }
                ]
            }
        )
    )
    (tmp_path / "artifact.json").write_text("{}")

    report, rows = build_report(
        tmp_path / "registry.json",
        tmp_path / "index.json",
        repo_root=tmp_path,
        required_portfolio={"qc": ["fixture"]},
        required_endpoint_portfolio={"qc": {"input": ["fixture"]}},
    )

    assert rows[0]["qc_status"] == "NOT_RUN"
    assert report["endpoint_gates"]["qc"]["input"]["status"] == "PASS"
    assert report["module_gates"]["qc"]["status"] == "PASS"
    assert report["module_gates"]["qc"]["gate_basis"] == "required_endpoint_portfolio"
    assert report["module_gates"]["qc"]["datasets_without_passed_endpoints"] == []


@pytest.mark.parametrize(
    "non_passing_status",
    [
        "NOT_EVALUABLE",
        "CONTRACT_PASS_NOT_PERFORMANCE",
        "SIMULATION_PASS_NOT_EXTERNAL",
    ],
)
def test_non_performance_status_does_not_pass_a_required_endpoint(
    tmp_path, non_passing_status
):
    registry = {
        "schema_name": "scLucidDatasetEvidenceRegistry",
        "endpoint_definitions": {
            "ambient": {
                "module": "qc",
                "estimand": "ambient",
                "experimental_unit": "library",
                "metrics": ["metric"],
                "acceptance": {"metric_min": 1.0},
            }
        },
        "datasets": [
            {
                "dataset_id": "filtered_only",
                "priority": "P0",
                "release_scope": "test",
                "modules": ["qc"],
                "truth_types": ["RV"],
                "accessions": [{"id": "local"}],
                "download": {
                    "raw_reads": "NOT_AVAILABLE",
                    "raw_counts": "NOT_AVAILABLE",
                    "processed": "AVAILABLE",
                },
                "license": {
                    "status": "REVIEW",
                    "identifier": "test",
                    "url": "LOCAL",
                    "redistribution": "NO",
                },
                "required_metadata": ["library"],
                "endpoint_ids": ["ambient"],
            }
        ],
    }
    (tmp_path / "registry.json").write_text(json.dumps(registry))
    (tmp_path / "index.json").write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "dataset_id": "filtered_only",
                        "endpoint_id": "ambient",
                        "status": non_passing_status,
                        "artifact": "artifact.json",
                    }
                ]
            }
        )
    )
    (tmp_path / "artifact.json").write_text("{}")

    report, _ = build_report(
        tmp_path / "registry.json",
        tmp_path / "index.json",
        repo_root=tmp_path,
        required_endpoint_portfolio={"qc": {"ambient": ["filtered_only"]}},
    )

    assert report["endpoint_gates"]["qc"]["ambient"]["status"] == "BLOCKED"
    assert report["module_gates"]["qc"]["status"] == "BLOCKED"


def test_passing_status_without_artifact_fails_closed(tmp_path):
    registry = {
        "schema_name": "scLucidDatasetEvidenceRegistry",
        "endpoint_definitions": {
            "input": {
                "module": "qc",
                "estimand": "input integrity",
                "experimental_unit": "library",
                "metrics": ["integrity"],
                "acceptance": {"integrity_min": 1.0},
            }
        },
        "datasets": [
            {
                "dataset_id": "fixture",
                "priority": "ENG",
                "release_scope": "test",
                "modules": ["qc"],
                "truth_types": ["ENG"],
                "accessions": [{"id": "generated"}],
                "download": {
                    "raw_reads": "NOT_APPLICABLE",
                    "raw_counts": "GENERATED_LOCAL",
                    "processed": "GENERATED_LOCAL",
                },
                "license": {
                    "status": "GENERATED_LOCAL",
                    "identifier": "test",
                    "url": "LOCAL",
                    "redistribution": "YES",
                },
                "required_metadata": ["library"],
                "endpoint_ids": ["input"],
            }
        ],
    }
    (tmp_path / "registry.json").write_text(json.dumps(registry))
    (tmp_path / "index.json").write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "dataset_id": "fixture",
                        "endpoint_id": "input",
                        "status": "PASS",
                    }
                ]
            }
        )
    )

    report, rows = build_report(
        tmp_path / "registry.json",
        tmp_path / "index.json",
        repo_root=tmp_path,
        required_endpoint_portfolio={"qc": {"input": ["fixture"]}},
    )

    assert rows[0]["endpoint_statuses"] == "MISSING_ARTIFACT"
    assert report["endpoint_gates"]["qc"]["input"]["status"] == "BLOCKED"


def test_required_endpoint_portfolio_validation_rejects_invalid_bindings():
    registry = json.loads(REGISTRY_PATH.read_text())
    issues = validate_required_endpoint_portfolio(
        registry,
        {
            "qc": {
                "pp_selector_regret": ["scmixology_gse118767"],
                "qc_input_contract": ["missing_dataset"],
                "qc_cell_calling": ["lin2020_pdac", "lin2020_pdac"],
            }
        },
    )

    assert any("registered under preprocess" in issue for issue in issues)
    assert any("required dataset is not registered" in issue for issue in issues)
    assert any("contains duplicates" in issue for issue in issues)
    assert any("does not declare required endpoint qc_cell_calling" in issue for issue in issues)


def test_required_endpoint_portfolio_validation_fails_closed_on_malformed_shapes():
    registry = json.loads(REGISTRY_PATH.read_text())

    assert validate_required_endpoint_portfolio(registry, []) == [
        "Required endpoint portfolio must be a mapping."
    ]
    assert validate_required_endpoint_portfolio(registry, {"qc": []}) == [
        "Required endpoint portfolio for qc must be a mapping."
    ]
    issues = validate_required_endpoint_portfolio(
        registry, {"qc": {"qc_input_contract": [{"not": "a dataset id"}]}}
    )
    assert any("must be a non-empty string" in issue for issue in issues)


def test_evidence_index_cannot_bind_an_undeclared_dataset_endpoint(tmp_path):
    evidence_path = tmp_path / "index.json"
    evidence_path.write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "dataset_id": "pbmc3k_engineering_fixture",
                        "endpoint_id": "qc_input_contract",
                        "status": "PASS",
                        "artifact": "artifact.json",
                    }
                ]
            }
        )
    )
    (tmp_path / "artifact.json").write_text("{}")

    report, _ = build_report(
        REGISTRY_PATH,
        evidence_path,
        repo_root=tmp_path,
    )

    assert report["registry_schema_status"] == "BLOCKED"
    assert any(
        "binds undeclared endpoint qc_input_contract" in issue
        for issue in report["issues"]
    )
