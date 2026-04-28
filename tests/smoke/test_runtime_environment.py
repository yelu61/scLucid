"""Smoke tests for runtime safety defaults."""

from pathlib import Path

import pytest


@pytest.mark.smoke
def test_runtime_environment_creates_writable_cache_dirs(monkeypatch, tmp_path):
    from scLucid.runtime import setup_runtime_environment

    monkeypatch.setenv("SCLUCID_RUNTIME_DIR", str(tmp_path))
    for key in ("XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR", "JOBLIB_TEMP_FOLDER"):
        monkeypatch.delenv(key, raising=False)

    env = setup_runtime_environment(force=True)

    assert env["MPLBACKEND"] == "Agg"
    for key in ("XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR", "JOBLIB_TEMP_FOLDER"):
        path = Path(env[key])
        assert path.exists()
        assert path.is_dir()
        assert str(path).startswith(str(tmp_path))


@pytest.mark.smoke
def test_runtime_parallel_defaults_to_sequential_in_ci(monkeypatch):
    from scLucid.runtime import effective_n_jobs, run_joblib_or_sequential

    monkeypatch.setenv("CI", "true")
    monkeypatch.delenv("SCLUCID_SAFE_PARALLEL", raising=False)

    assert effective_n_jobs(-1) == 1
    assert run_joblib_or_sequential(lambda x: x + 1, [1, 2, 3], n_jobs=-1) == [2, 3, 4]


@pytest.mark.smoke
def test_joblib_failure_falls_back_to_sequential(monkeypatch):
    from scLucid.runtime import run_joblib_or_sequential

    for key in ("CI", "GITHUB_ACTIONS", "GITLAB_CI", "BUILDKITE", "JENKINS_URL"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.delenv("SCLUCID_SAFE_PARALLEL", raising=False)

    assert run_joblib_or_sequential(lambda x: x * 2, [1, 2, 3], n_jobs=2, backend="bad") == [
        2,
        4,
        6,
    ]
