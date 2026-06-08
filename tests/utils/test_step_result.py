"""Tests for StepResult semantics and utilities."""

from scLucid.utils.step_result import (
    StepResult,
    rollup_step_status,
    step_results_from_storage,
    step_results_to_storage,
    summarize_step_results,
)


def test_step_result_default_fields():
    r = StepResult(name="test_step")
    assert r.name == "test_step"
    assert r.status == "completed"
    assert r.evidence_level == "heuristic"
    assert r.outputs == {}
    assert r.warnings == []
    assert r.error is None


def test_step_result_from_exception_failed():
    exc = ValueError("something broke")
    r = StepResult.from_exception(name="bad_step", exc=exc)
    assert r.status == "failed"
    assert r.evidence_level == "unavailable"
    assert "ValueError" in r.error
    assert "something broke" in r.error


def test_step_result_from_exception_degraded():
    exc = RuntimeError("optional dependency missing")
    r = StepResult.from_exception(name="ok_step", exc=exc, degraded=True)
    assert r.status == "degraded"
    assert r.evidence_level == "unavailable"


def test_step_result_skipped():
    r = StepResult.skipped(name="skip_step", reason="missing input")
    assert r.status == "skipped"
    assert r.warnings == ["missing input"]


def test_step_result_storage_roundtrip():
    results = [
        StepResult(name="a", status="completed", evidence_level="validated_core"),
        StepResult.from_exception(name="b", exc=KeyError("x")),
        StepResult.skipped(name="c", reason="no data"),
    ]
    storage = step_results_to_storage(results)
    assert isinstance(storage, dict)
    recovered = step_results_from_storage(storage)
    assert len(recovered) == 3
    assert recovered[0].name == "a"
    assert recovered[0].status == "completed"
    assert recovered[1].status == "failed"
    assert recovered[2].status == "skipped"


def test_rollup_step_status():
    assert rollup_step_status([]) == "completed"
    assert (
        rollup_step_status(
            [
                StepResult(name="a", status="completed"),
                StepResult(name="b", status="completed"),
            ]
        )
        == "completed"
    )
    assert (
        rollup_step_status(
            [
                StepResult(name="a", status="completed"),
                StepResult(name="b", status="degraded"),
            ]
        )
        == "degraded"
    )
    assert (
        rollup_step_status(
            [
                StepResult(name="a", status="completed"),
                StepResult(name="b", status="skipped"),
                StepResult(name="c", status="failed"),
            ]
        )
        == "failed"
    )


def test_summarize_step_results():
    results = [
        StepResult(name="a", status="completed", evidence_level="validated_core"),
        StepResult(
            name="b", status="degraded", evidence_level="heuristic", warnings=["w1"]
        ),
        StepResult.from_exception(name="c", exc=RuntimeError("x")),
    ]
    summary = summarize_step_results(results)
    assert summary["n_steps"] == 3
    assert summary["overall_status"] == "failed"
    assert summary["by_status"]["completed"] == ["a"]
    assert "w1" in summary["warnings"]
    assert any("RuntimeError" in e for e in summary["errors"])


def test_step_result_to_storage_is_hdf5_safe():
    r = StepResult(name="x", outputs={"n": 5})
    d = r.to_storage_dict()
    assert isinstance(d, dict)
    assert d["error"] == ""  # None should not leak into HDF5


def test_step_results_from_storage_accepts_legacy_list():
    storage = [
        StepResult(name="a").to_storage_dict(),
        StepResult(name="b", status="skipped").to_storage_dict(),
    ]
    recovered = step_results_from_storage(storage)
    assert [r.name for r in recovered] == ["a", "b"]
