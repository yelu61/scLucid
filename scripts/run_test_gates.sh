#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export PYTHONPATH="${PYTHONPATH:-src}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mplcfg}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba-cache}"

# Syntax check for changed Python files; fall back to package/tests if no git diff is available.
changed_py=()
while IFS= read -r line; do
  [[ -n "$line" && -f "$line" ]] && changed_py+=("$line")
done < <(git diff --name-only -- '*.py' 2>/dev/null || true)

if [[ ${#changed_py[@]} -eq 0 ]]; then
  while IFS= read -r line; do
    [[ -n "$line" && -f "$line" ]] && changed_py+=("$line")
  done < <(find src/scLucid tests -type f -name '*.py' | sort)
fi

if [[ ${#changed_py[@]} -gt 0 ]]; then
  python - <<'PY' "${changed_py[@]}"
import py_compile
import sys

files = sys.argv[1:]
for path in files:
    py_compile.compile(path, doraise=True)
print(f"py_compile ok for {len(files)} files")
PY
fi

# Import/public surface smoke checks.
pytest -q tests/smoke/test_imports.py

# Baseline module and integration checks.
pytest -q \
  tests/config/test_global_config_runtime.py \
  tests/qc/test_qc_public_api.py \
  tests/preprocess/test_preprocess_public_api.py \
  tests/analysis/test_analysis_public_api.py \
  tests/tools/test_tools_public_api.py \
  tests/utils/test_public_api.py \
  tests/integration/test_full_pipeline.py
