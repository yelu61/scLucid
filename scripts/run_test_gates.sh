#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-src}"
export MPLBACKEND="${MPLBACKEND:-Agg}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/sclucid/cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/sclucid/matplotlib}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/sclucid/numba}"
export JOBLIB_TEMP_FOLDER="${JOBLIB_TEMP_FOLDER:-/tmp/sclucid/joblib}"
export SCLUCID_SAFE_PARALLEL="${SCLUCID_SAFE_PARALLEL:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$JOBLIB_TEMP_FOLDER"

if [[ "${SCLUCID_TEST_ENV_ACTIVE:-0}" != "1" ]]; then
  env_args=()
  if [[ -n "${SCLUCID_TEST_ENV_PATH:-}" ]]; then
    env_args=(-p "$SCLUCID_TEST_ENV_PATH")
  elif [[ -n "${SCLUCID_TEST_ENV:-}" ]]; then
    env_args=(-n "$SCLUCID_TEST_ENV")
  fi

  if [[ "${#env_args[@]}" -gt 0 ]]; then
    runner="${MAMBA_EXE:-}"
    if [[ -z "$runner" ]]; then
      runner="$(command -v micromamba || command -v mamba || true)"
    fi
    if [[ -z "$runner" ]]; then
      echo "SCLUCID_TEST_ENV was set, but neither micromamba nor mamba is on PATH." >&2
      echo "Set MAMBA_EXE=/path/to/mamba or activate the environment manually." >&2
      exit 2
    fi
    export SCLUCID_TEST_ENV_ACTIVE=1
    exec "$runner" run "${env_args[@]}" "$0" "$@"
  fi
fi

echo "==> Python"
python --version

echo "==> Syntax check"
python -m compileall -q src/scLucid tests

if ! python - <<'PY'
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("pytest") else 1)
PY
then
  cat >&2 <<'EOF'
pytest is not installed in this environment.
Install development dependencies, then rerun:

  python -m pip install -e ".[dev]"

EOF
  exit 2
fi

echo "==> Smoke tests"
python -m pytest -q tests/smoke

echo "==> Core contract/config tests"
python -m pytest -q \
  tests/test_contracts.py \
  tests/test_configs.py \
  tests/config/test_global_config_runtime.py \
  tests/utils/test_public_api.py

echo "==> Lightweight integration tests"
python -m pytest -q tests/integration/test_full_pipeline.py

echo "All lightweight gates passed."
