"""Runtime safeguards for fresh installs, tests, and headless environments."""

from __future__ import annotations

import logging
import os
import tempfile
import warnings
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence, TypeVar

log = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")


def is_ci_environment() -> bool:
    """Return whether scLucid appears to be running in CI."""
    return any(
        os.environ.get(key)
        for key in (
            "CI",
            "GITHUB_ACTIONS",
            "GITLAB_CI",
            "BUILDKITE",
            "JENKINS_URL",
        )
    )


def _runtime_root() -> Path:
    return Path(os.environ.get("SCLUCID_RUNTIME_DIR", tempfile.gettempdir())) / "sclucid"


def _ensure_dir(path: str | os.PathLike[str]) -> str:
    resolved = Path(path).expanduser()
    resolved.mkdir(parents=True, exist_ok=True)
    return str(resolved)


def setup_runtime_environment(*, force: bool = False) -> dict[str, str]:
    """Set writable runtime paths before heavy scientific libraries import."""
    root = _runtime_root()
    defaults = {
        "MPLBACKEND": "Agg",
        "XDG_CACHE_HOME": str(root / "cache"),
        "MPLCONFIGDIR": str(root / "matplotlib"),
        "NUMBA_CACHE_DIR": str(root / "numba"),
        "JOBLIB_TEMP_FOLDER": str(root / "joblib"),
    }

    for key, value in defaults.items():
        if force or not os.environ.get(key):
            os.environ[key] = value

    for key in ("XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR", "JOBLIB_TEMP_FOLDER"):
        os.environ[key] = _ensure_dir(os.environ[key])

    if is_ci_environment() or os.environ.get("SCLUCID_SAFE_PARALLEL"):
        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            os.environ.setdefault(key, "1")
        os.environ.setdefault("SCLUCID_SAFE_PARALLEL", "1")

    return {
        "MPLBACKEND": os.environ["MPLBACKEND"],
        "XDG_CACHE_HOME": os.environ["XDG_CACHE_HOME"],
        "MPLCONFIGDIR": os.environ["MPLCONFIGDIR"],
        "NUMBA_CACHE_DIR": os.environ["NUMBA_CACHE_DIR"],
        "JOBLIB_TEMP_FOLDER": os.environ["JOBLIB_TEMP_FOLDER"],
    }


def effective_n_jobs(n_jobs: int | None = None, *, max_jobs: int | None = None) -> int:
    """Normalize job counts for stable local and CI execution."""
    requested = -1 if n_jobs is None else int(n_jobs)
    if requested == 0 or requested < -1:
        warnings.warn(
            f"Invalid n_jobs={requested}; falling back to sequential execution.",
            RuntimeWarning,
            stacklevel=2,
        )
        return 1

    if os.environ.get("SCLUCID_SAFE_PARALLEL") or is_ci_environment():
        return 1

    cpu_count = os.cpu_count() or 1
    resolved = cpu_count if requested == -1 else requested
    if max_jobs is not None:
        resolved = min(resolved, max(1, int(max_jobs)))
    return max(1, resolved)


def run_joblib_or_sequential(
    func: Callable[[T], R],
    items: Iterable[T],
    *,
    n_jobs: int | None = None,
    backend: str | None = None,
    prefer: str | None = None,
    description: str = "parallel task",
) -> list[R]:
    """Run a simple map with joblib and fall back to sequential execution."""
    materialized: Sequence[T] = list(items)
    resolved_jobs = effective_n_jobs(n_jobs, max_jobs=len(materialized) or 1)
    if resolved_jobs == 1 or len(materialized) <= 1:
        return [func(item) for item in materialized]

    try:
        from joblib import Parallel, delayed

        kwargs: dict[str, Any] = {"n_jobs": resolved_jobs}
        if backend:
            kwargs["backend"] = backend
        if prefer:
            kwargs["prefer"] = prefer
        return Parallel(**kwargs)(delayed(func)(item) for item in materialized)
    except Exception as exc:
        log.warning(
            "joblib execution failed for %s (%s). Falling back to sequential execution.",
            description,
            exc,
        )
        return [func(item) for item in materialized]


__all__ = [
    "effective_n_jobs",
    "is_ci_environment",
    "run_joblib_or_sequential",
    "setup_runtime_environment",
]
