from validation.qc.run_hgmm_doublet_truth_benchmark import _classify


def test_doublet_regret_gate_passes_at_threshold() -> None:
    assert _classify(0.55, 0.60, 0.05) == "PASS"


def test_doublet_regret_gate_fails_above_threshold() -> None:
    assert _classify(0.54, 0.60, 0.05) == "FAIL"
