# Contributing to scLucid

Thanks for considering a contribution! scLucid is a single-cell RNA-seq analysis
toolkit focused on tumor-aware workflows, traceable parameter selection, and
honest Python ports of mature R bioinformatics packages.

This document covers the practical mechanics. The architectural philosophy
lives in `docs/source/usage_layers.rst`, `docs/source/data_contracts.rst`, and
`docs/source/qc_preprocess_maturity.rst` — read those before proposing a major
change.

---

## Quick start

```bash
git clone https://github.com/yelu61/scLucid.git
cd scLucid
pip install -e ".[dev,docs]"
pre-commit install   # optional but recommended
```

If you maintain a dedicated single-cell Python environment, point the test
gates at it:

```bash
MAMBA_EXE=/opt/homebrew/bin/mamba \
SCLUCID_TEST_ENV_PATH=/path/to/your/scrna-env \
scripts/run_test_gates.sh
```

---

## Where things live

```
src/scLucid/
  qc/             quality control: doublets, filters, intelligent QC, benchmarks
  preprocess/     normalize, HVG, scale, integration, neighbors, PCA, UMAP
  analysis/       clustering, annotation, differential expression, proportion
  tumor/          CNV, malignancy, microenvironment, evolution, therapy
  recommendation/ cross-stage parameter recommendation engine
  tools/          R-package ports (pyBayesPrism, pyMonocle3, pyCellChat, pyDWLS)
                  + lightweight wrappers (pySCENIC, cellphonedb, infercnv, spatial)
  plotting/       publication themes + embedding/marker/feature plots
  utils/          contracts, context, validation, storage, profiling
```

The three user-facing API layers are:

- **Workflow**: `scl.run_pipeline()`, `scl.run_standard_qc()`, etc.
- **Simple API**: `scl.pp.normalize_data()`, `scl.qc.calculate_qc_metric()`, etc.
- **Advanced notebooks**: `examples/03_advanced_notebooks/`

A new feature must work cleanly in **all three layers** (or come with a clear
note about which layer it intentionally skips).

---

## Required reading before a substantial change

- `docs/source/data_contracts.rst` — the stable AnnData and review-summary
  conventions. Most contributions should preserve them; if you change them,
  bump `SCHEMA_VERSION` in `src/scLucid/utils/contracts.py`.
- `docs/source/workflow_hardening.rst` — how real-data validation works
  (PBMC + PDAC golden paths).
- `docs/PLUGIN_DEVELOPMENT_GUIDE.md` — extension points for custom
  scoring/annotation/filter methods.

---

## Coding conventions

- **Formatting**: `black` (line-length 100).
- **Linting**: `ruff check src/scLucid tests` must pass.
- **Type checks**: `mypy src/scLucid` should not regress.
- **Docstrings**: NumPy-style with `Parameters / Returns: / Raises: / Examples:`
  sections (note the trailing colons — that's the project idiom matched by
  Sphinx Napoleon).
- **Logging**: `log = logging.getLogger(__name__)` at module top. Use
  `log.info` for milestones, `log.debug` for internals, `log.warning` for
  recoverable degraded paths, no `print`.
- **Configs**: Inherit `SclucidBaseConfig` from `src/scLucid/base_config.py`.
  Don't introduce plain dataclasses for configuration.
- **Fitted attributes**: trailing underscore (`signature_matrix_`, `results_`).

---

## Tests

Tests are the spec — many functions in scLucid have a pre-existing test file
that pins the public surface. Always check `tests/` before designing a new API.

```bash
# Smoke tests (fast, runs on every CI matrix entry)
pytest tests/smoke

# Full suite (~4 min)
pytest

# Single module
pytest tests/qc -v --no-cov

# With coverage report for a specific file
pytest tests/utils/test_validation.py \
       --cov=src/scLucid/utils/validation --cov-report=term-missing
```

Markers:

- `@pytest.mark.unit` — fast isolated test (default)
- `@pytest.mark.integration` — may use real data
- `@pytest.mark.slow` — skip in quick runs
- `@pytest.mark.optional` — depends on an extra dependency

A new test class without a marker defaults to `unit`. Add an explicit marker
for anything slow or environment-dependent.

`pytest.skip(...)` is reserved for *real* environment differences (missing
optional dep, dataset not available). Do not use it to mask missing
implementations — that creates "green CI, broken code" debt.

---

## Coverage policy

- The CI floor is set in `pyproject.toml` (`--cov-fail-under`). Don't lower it
  without discussion.
- Aim to keep new modules above 60% coverage and any change-touched file from
  regressing.
- Files in `src/scLucid/utils/` should target 90%+ (they are core contracts).

---

## Pull request checklist

Before opening a PR:

- [ ] `pre-commit run --all-files` clean
- [ ] `pytest` clean
- [ ] Added or updated tests for new code paths
- [ ] Updated relevant documentation (`docs/source/*.rst`, examples, README)
- [ ] Added an entry under "Unreleased" in the changelog if user-facing
- [ ] Confirmed `import scLucid` produces zero `ImportWarning`

When the PR touches a workflow contract, the AnnData layout, or a public API:

- [ ] Bumped `SCHEMA_VERSION` if needed
- [ ] Updated `docs/source/data_contracts.rst`
- [ ] Verified PBMC golden path still runs:
      `scripts/run_pbmc_golden_path.py --n-cells 300 --output-dir /tmp/pbmc_check --overwrite`

---

## Reporting bugs

Use the GitHub issue template at `.github/ISSUE_TEMPLATE/bug_report.md`.
Real-data bug reports should include:

- the dataset shape (`adata.shape`) and species
- the exact `scl.run_*` call or notebook cell
- the full traceback
- `scLucid.__version__` and a `pip freeze | grep -E 'scanpy|anndata|sclucid'`

For tumor-specific reports, please include `cancer_type` and whether the data
is single-sample, multi-sample, primary, metastatic, or treated/untreated.

---

## License

scLucid is MIT licensed. By contributing, you agree your contribution is also
MIT licensed.
