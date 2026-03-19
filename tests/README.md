# Testing Guide

## Baseline test gates

Run the standard pre-merge checks:

```bash
scripts/run_test_gates.sh
```

This runs:
- `py_compile` syntax checks on changed Python files (fallback: all package/tests files)
- Import/public-surface smoke tests
- Core config/public API tests
- Lightweight integration coverage for QC metrics, preprocessing normalization, and clustering entry points

## Useful direct commands

```bash
# Full test collection sanity
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest --collect-only -q

# Smoke-only
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest -m smoke -q

# Integration-only
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest -m integration -q
```
