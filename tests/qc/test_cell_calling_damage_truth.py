from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from validation.qc.evaluate_cell_calling_damage_truth import (
    evaluate_truth_table,
    write_report,
)

REGISTRY = Path("validation/dataset_evidence_registry.json")


def _truth_table() -> pd.DataFrame:
    rows = []
    for library in ("lib_a", "lib_b"):
        for index in range(20):
            rows.append(
                {
                    "library": library,
                    "truth_cell_class": "intact",
                    "truth_low_rna": index < 4,
                    "predicted_cell_call": True,
                    "predicted_qc_decision": "REVIEW" if index == 0 else "KEEP",
                    "predicted_damage_probability": 0.02,
                    "baseline_expert_global_qc_decision": "KEEP",
                    "baseline_per_sample_mad_qc_decision": "KEEP",
                }
            )
        for index in range(10):
            rows.append(
                {
                    "library": library,
                    "truth_cell_class": "damaged",
                    "truth_low_rna": False,
                    "predicted_cell_call": True,
                    "predicted_qc_decision": "REVIEW" if index == 0 else "REMOVE",
                    "predicted_damage_probability": 0.95,
                    "baseline_expert_global_qc_decision": ("REMOVE" if index < 4 else "KEEP"),
                    "baseline_per_sample_mad_qc_decision": ("REMOVE" if index < 5 else "KEEP"),
                }
            )
        for _ in range(10):
            rows.append(
                {
                    "library": library,
                    "truth_cell_class": "empty",
                    "truth_low_rna": False,
                    "predicted_cell_call": False,
                    "predicted_qc_decision": "REVIEW",
                    "predicted_damage_probability": None,
                    "baseline_expert_global_qc_decision": "KEEP",
                    "baseline_per_sample_mad_qc_decision": "KEEP",
                }
            )
    rows.append(
        {
            "library": "lib_b",
            "truth_cell_class": "uncertain",
            "truth_low_rna": False,
            "predicted_cell_call": True,
            "predicted_qc_decision": "REMOVE",
            "predicted_damage_probability": 0.5,
            "baseline_expert_global_qc_decision": "REMOVE",
            "baseline_per_sample_mad_qc_decision": "REMOVE",
        }
    )
    return pd.DataFrame(rows)


def test_endpoints_pass_with_per_library_metrics_and_fixed_bootstrap():
    frame = _truth_table()

    first = evaluate_truth_table(
        frame,
        registry_path=REGISTRY,
        n_bootstrap=100,
        seed=17,
    )
    second = evaluate_truth_table(
        frame,
        registry_path=REGISTRY,
        n_bootstrap=100,
        seed=17,
    )

    assert first["status"] == "PASS"
    assert first["endpoint_status"] == {
        "qc_cell_calling": "PASS",
        "qc_damage_classification": "PASS",
    }
    cell_call = first["endpoints"]["qc_cell_calling"]
    assert cell_call["aggregate_metrics"]["true_cell_recall"] == 1.0
    assert cell_call["aggregate_metrics"]["low_rna_cell_recall"] == 1.0
    assert cell_call["aggregate_metrics"]["empty_droplet_fdr"] == 0.0
    assert cell_call["aggregate_metrics"]["n_truth_uncertain_excluded"] == 1
    assert len(cell_call["by_library"]) == 2
    assert (
        cell_call["grouped_bootstrap"]
        == second["endpoints"]["qc_cell_calling"]["grouped_bootstrap"]
    )
    damage = first["endpoints"]["qc_damage_classification"]
    assert damage["aggregate_metrics"]["keep_false_removal_rate"] == 0.0
    assert damage["aggregate_metrics"]["damaged_cell_recall"] == 1.0
    assert damage["coverage"]["prediction_review_or_uncertain_rows"] == 4
    assert damage["coverage"]["truth_uncertain_rows_excluded"] == 1
    assert damage["coverage"]["empty_rows_excluded"] == 20
    assert {row["baseline"] for row in damage["comparisons"]} == {
        "expert_global",
        "per_sample_mad",
    }
    assert all(row["absolute_recall_gain"] >= 0.05 for row in damage["comparisons"])
    assert damage["calibration"]["status"] == "EVALUATED"
    assert first["truth_independence"]["status"] == "NOT_VERIFIED_FROM_TABLE"
    assert "Predictions are never used" in first["claim_boundary"]["unsupported"][2]


def test_missing_baseline_blocks_damage_instead_of_skipping_comparison():
    frame = _truth_table().drop(columns="baseline_per_sample_mad_qc_decision")

    report = evaluate_truth_table(
        frame,
        registry_path=REGISTRY,
        n_bootstrap=20,
    )

    assert report["endpoint_status"]["qc_cell_calling"] == "PASS"
    assert report["endpoint_status"]["qc_damage_classification"] == "BLOCKED"
    assert report["status"] == "BLOCKED"
    assert report["endpoints"]["qc_damage_classification"]["blockers"] == [
        "Missing required baseline column for per_sample_mad: baseline_per_sample_mad_qc_decision"
    ]


def test_cell_calling_fails_registered_threshold_without_changing_damage_truth():
    frame = _truth_table()
    intact = frame["truth_cell_class"].eq("intact")
    frame.loc[frame.index[intact][:8], "predicted_cell_call"] = False

    report = evaluate_truth_table(
        frame,
        registry_path=REGISTRY,
        n_bootstrap=20,
    )

    cell_call = report["endpoints"]["qc_cell_calling"]
    assert cell_call["status"] == "FAIL"
    assert (
        cell_call["aggregate_metrics"]["true_cell_recall"]
        < cell_call["thresholds"]["true_cell_recall_min"]
    )
    assert report["endpoints"]["qc_damage_classification"]["status"] == "PASS"
    assert report["status"] == "FAIL"


def test_thresholds_are_loaded_from_registry(tmp_path):
    registry = json.loads(REGISTRY.read_text())
    registry["endpoint_definitions"]["qc_cell_calling"]["acceptance"]["true_cell_recall_min"] = 1.01
    registry["endpoint_definitions"]["qc_damage_classification"]["acceptance"][
        "absolute_recall_gain_vs_each_baseline_min"
    ] = 1.01
    custom_registry = tmp_path / "registry.json"
    custom_registry.write_text(json.dumps(registry))

    report = evaluate_truth_table(
        _truth_table(),
        registry_path=custom_registry,
        n_bootstrap=20,
    )

    assert report["endpoint_status"]["qc_cell_calling"] == "FAIL"
    assert report["endpoint_status"]["qc_damage_classification"] == "FAIL"
    assert report["threshold_source"] == str(custom_registry)


def test_tsv_cli_writes_json_and_markdown(tmp_path):
    table_path = tmp_path / "truth.tsv"
    output_dir = tmp_path / "out"
    _truth_table().to_csv(table_path, sep="\t", index=False)

    completed = subprocess.run(
        [
            sys.executable,
            "validation/qc/evaluate_cell_calling_damage_truth.py",
            str(table_path),
            "--registry",
            str(REGISTRY),
            "--output-dir",
            str(output_dir),
            "--n-bootstrap",
            "20",
            "--seed",
            "9",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    summary = json.loads(completed.stdout)
    payload = json.loads((output_dir / "cell_calling_damage_truth_evaluation.json").read_text())

    assert summary["status"] == "PASS"
    assert payload["bootstrap"] == {"unit": "library", "seed": 9, "n_bootstrap": 20}
    assert (output_dir / "cell_calling_damage_truth_evaluation.md").exists()


def test_write_report_refuses_nonfinite_json_by_contract(tmp_path):
    report = evaluate_truth_table(
        _truth_table(),
        registry_path=REGISTRY,
        n_bootstrap=20,
    )

    json_path, markdown_path = write_report(report, tmp_path)

    assert "NaN" not in json_path.read_text()
    assert "Does not" not in markdown_path.read_text()
    assert "does not prove" in markdown_path.read_text()
