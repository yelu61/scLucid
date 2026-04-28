# Testing Guide

## Baseline test gates

Run the standard pre-merge checks:

```bash
scripts/run_test_gates.sh
```

To run the gates in the local single-cell environment without manually
activating it:

```bash
MAMBA_EXE=/opt/homebrew/bin/mamba \
SCLUCID_TEST_ENV_PATH=/Users/luye/micromamba/envs/scrna-env \
scripts/run_test_gates.sh
```

This runs:
- `compileall` syntax checks on package and test files
- Import/public-surface smoke tests
- Core config/public API tests
- Lightweight integration coverage for QC metrics, preprocessing normalization, and clustering entry points

If the script reports that `pytest` is missing, install the development extras:

```bash
python -m pip install -e ".[dev]"
```

## Useful direct commands

```bash
# Full test collection sanity
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest --collect-only -q

# Smoke-only
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest -m smoke -q

# Integration-only
PYTHONPATH=src MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba-cache pytest -m integration -q
```
