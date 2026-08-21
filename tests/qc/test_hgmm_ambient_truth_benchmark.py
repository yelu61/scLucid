from __future__ import annotations

import numpy as np

from validation.qc.run_hgmm_ambient_truth_benchmark import (
    _classify,
    _species_metrics,
)

THRESHOLDS = {
    "contamination_reduction_min": 0.5,
    "native_marker_loss_max": 0.02,
    "identity_accuracy_loss_max": 0.02,
}


def test_species_metrics_and_gate_pass() -> None:
    before = np.array([[100, 10], [10, 100]], dtype=float)
    after = np.array([[100, 4], [4, 100]], dtype=float)
    classes = np.array(["HUMAN_SINGLET", "MOUSE_SINGLET"])

    metrics = _species_metrics(before, after, classes)

    assert metrics["contamination_reduction"] == 0.6
    assert metrics["native_marker_loss"] == 0.0
    assert metrics["identity_accuracy_loss"] == 0.0
    assert _classify(metrics, THRESHOLDS) == "PASS"


def test_species_gate_fails_without_contamination_reduction() -> None:
    metrics = {
        "contamination_reduction": 0.0,
        "native_marker_loss": 0.0,
        "identity_accuracy_loss": 0.0,
    }

    assert _classify(metrics, THRESHOLDS) == "FAIL"
