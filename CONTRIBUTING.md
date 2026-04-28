# Contributing to scLucid

Thanks for helping make scLucid more useful and reliable. This project is still
moving quickly, so the most valuable contributions are changes that make real
single-cell workflows easier to run, inspect, and reproduce.

## Development Setup

```bash
git clone https://github.com/yelu61/scLucid.git
cd scLucid
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

For modules that rely on optional scientific stacks, install the relevant extra:

```bash
python -m pip install -e ".[all,dev]"
```

On this development machine, the maintained single-cell environment is available
at ``/Users/luye/micromamba/envs/scrna-env`` and can be used directly by the
test gates:

```bash
MAMBA_EXE=/opt/homebrew/bin/mamba \
SCLUCID_TEST_ENV_PATH=/Users/luye/micromamba/envs/scrna-env \
scripts/run_test_gates.sh
```

On macOS, if SciPy/Scanpy fails with a missing Fortran runtime, prefer a clean
conda-forge environment:

```bash
conda create -n sclucid-dev -c conda-forge python=3.11 scipy scanpy
conda activate sclucid-dev
python -m pip install -e ".[dev]"
```

## Test Gates

Run the lightweight pre-merge gates before opening a PR:

```bash
scripts/run_test_gates.sh
```

The script checks syntax, import smoke tests, contract/config tests, and a small
integration subset when `pytest` is available. If `pytest` is not installed, it
prints the install command and exits with a clear error.

For a fuller local check:

```bash
python -m pytest -m "not slow" tests
```

## Contribution Guidelines

- Preserve the canonical AnnData contract documented in
  `docs/source/data_contracts.rst`.
- Store workflow outputs under `adata.uns["sclucid"][module]`.
- Keep public workflow decisions reviewable through `review_summary` records.
- Add or update tests when changing workflow behavior, storage keys, or user
  facing defaults.
- Keep wrapper modules honest: when wrapping mature tools, document the supported
  subset and add parity or smoke tests where practical.

## Pull Request Checklist

- [ ] The change is scoped to one workflow, module, or contract surface.
- [ ] New public behavior is documented.
- [ ] Tests or smoke coverage were added for regressions.
- [ ] `scripts/run_test_gates.sh` passes locally, or the reason it cannot run is
      included in the PR notes.
