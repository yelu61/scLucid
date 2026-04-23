"""Runtime behavior tests for global configuration APIs."""

from pathlib import Path

import pytest

from scLucid.config import (
    GlobalConfig,
    config_context,
    get_config,
    reset_config,
    set_config,
)


@pytest.mark.config
@pytest.mark.unit
def test_get_config_returns_singleton():
    first = get_config()
    second = get_config()
    assert first is second


@pytest.mark.config
@pytest.mark.unit
def test_set_config_updates_known_fields_and_rejects_unknown():
    reset_config()
    set_config(n_jobs=2, verbosity=1)
    cfg = get_config()
    assert cfg.n_jobs == 2
    assert cfg.verbosity == 1

    with pytest.raises(ValueError):
        set_config(not_a_real_key=1)


@pytest.mark.config
@pytest.mark.unit
def test_config_context_restores_values_on_exit_and_exception():
    reset_config()
    baseline = get_config().n_jobs

    with config_context(n_jobs=3):
        assert get_config().n_jobs == 3

    assert get_config().n_jobs == baseline

    with pytest.raises(RuntimeError), config_context(n_jobs=5):
        assert get_config().n_jobs == 5
        raise RuntimeError("boom")

    assert get_config().n_jobs == baseline


@pytest.mark.config
@pytest.mark.unit
def test_reset_config_restores_defaults():
    reset_config()
    set_config(n_jobs=4, verbosity=2)
    assert get_config().n_jobs == 4

    reset_config()
    cfg = get_config()
    assert cfg.n_jobs == -1
    assert cfg.verbosity == 1


@pytest.mark.config
@pytest.mark.unit
def test_global_config_converts_string_paths(tmp_path):
    log_path = tmp_path / "logs" / "runtime.log"
    cfg = GlobalConfig(log_file=str(log_path), cache_dir=str(tmp_path / "cache"))

    assert isinstance(cfg.log_file, Path)
    assert isinstance(cfg.cache_dir, Path)
    assert cfg.log_file.parent.exists()
