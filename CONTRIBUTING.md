# Contributing to scLucid

Thanks for considering a contribution! scLucid is a single-cell RNA-seq analysis
framework focused on tumor-aware interpretation, traceable parameter selection,
explicit inference semantics, and audit-ready workflow outputs.

This document covers the practical mechanics. The architectural philosophy
lives in `README.md`, `docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`,
`docs/source/usage_layers.rst`, `docs/source/data_contracts.rst`, and
`docs/source/qc_preprocess_maturity.rst`. Read those before proposing a major
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
  tools/          evidence modules and ecosystem adapters:
                  bulk, spatial, pyBayesPrism, pyDWLS, pyMonocle3, pyCellChat,
                  pySCENIC, cellphonedb, infercnv
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

## Product boundaries

scLucid is not a broad multi-omics toolbox and it is not a claim that every
result is automatically better than Scanpy, Seurat, or specialist tools. New
contributions should strengthen one of these routes:

- the core `qc -> preprocess -> analysis` single-cell workflow;
- tumor interpretation (`malignancy`, CNV evidence, TME, therapy, heterogeneity,
  ecosystem/ecotype summaries);
- evidence modules under `scLucid.tools` such as bulk/spatial utilities;
- selective R/Python parity with validation and dependency isolation;
- engineering quality: lightweight imports, sparse-aware execution, reports,
  benchmarks, and clear error messages.

Avoid adding a method merely because it exists upstream. It must fit a real
tumor single-cell workflow need and have a validation story.

---

## Required reading before a substantial change

- `docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md` — current five-direction
  strategic plan and implementation milestones.
- `docs/source/data_contracts.rst` — the stable AnnData and review-summary
  conventions. Most contributions should preserve them; if you change them,
  bump `SCHEMA_VERSION` in `src/scLucid/utils/contracts.py`.
- `docs/source/workflow_hardening.rst` — how real-data validation works
  (PBMC + PDAC golden paths).
- `docs/BULK_SPATIAL_DESIGN.md` — namespace and storage contract for
  `scLucid.tools.bulk` and `scLucid.tools.spatial`.
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

Generated coverage artifacts (`coverage.xml`, `htmlcov/`, `.coverage*`) should
not be committed.

---

## Pull request checklist

Before opening a PR:

- [ ] `pre-commit run --all-files` clean
- [ ] `pytest` clean
- [ ] Added or updated tests for new code paths
- [ ] Updated relevant documentation (`docs/source/*.rst`, examples, README)
- [ ] Added an entry under "Unreleased" in the changelog if user-facing
- [ ] Confirmed `import scLucid` produces zero `ImportWarning`
- [ ] Confirmed no generated artifacts are included (`htmlcov/`, `coverage.xml`,
      ad-hoc output directories)

When the PR touches a workflow contract, the AnnData layout, or a public API:

- [ ] Bumped `SCHEMA_VERSION` if needed
- [ ] Updated `docs/source/data_contracts.rst`
- [ ] Verified PBMC golden path still runs:
      `scripts/run_pbmc_golden_path.py --n-cells 300 --output-dir /tmp/pbmc_check --overwrite`

When the PR adds or wraps an external method:

- [ ] Dependency is optional or already core by design
- [ ] Missing dependency path has a clear error message
- [ ] Method version/parameters are recorded in `adata.uns["sclucid"]`
- [ ] Output includes inference-level semantics where relevant
- [ ] Validation/parity limitations are documented

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
