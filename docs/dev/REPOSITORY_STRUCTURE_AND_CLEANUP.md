# Repository Structure And Cleanup Policy

This note defines where active code, examples, tests, and validation assets
belong. It is intentionally practical: when a file no longer fits one of these
roles, archive or remove it instead of leaving another parallel workflow in the
root tree.

## Top-Level Roles

| Path | Role | Cleanup rule |
|---|---|---|
| `src/` | Package source code and public/private implementation modules. | Keep only maintained implementation code. Deprecated compatibility shims need tests and a removal note. |
| `data/` | Local real-data fixtures and benchmark dataset manifests. | Keep manifests and small checked-in fixtures. Large downloaded/generated files stay ignored. |
| `docs/` | User docs, API docs, design governance, roadmap, and archived provenance. | Current docs should point to one source of truth. Historical drafts belong in `docs/archive/`. |
| `examples/` | User-facing usage examples and notebooks. | Keep examples aligned with public APIs and current workflow contracts. |
| `scripts/` | Maintained developer/user entrypoints that are not package APIs. | Keep only reusable entrypoints. One-off validation scripts should move into `validation/` or be removed. |
| `tests/` | Unit, contract, integration, and smoke tests for current behavior. | Do not delete by age alone. Remove tests only when the protected behavior or compatibility promise is removed. |
| `validation/` | Executable scientific benchmark runners and dataset registries. | This is the formal evidence layer. Keep runners current; generated outputs do not live here. |
| `validation_outputs/` | Generated benchmark outputs. | Ignored by git. Keep only intentionally retained evidence packages; clean older rerun directories when superseded. |

## Scripts

`scripts/` should stay small. A script belongs here only if it is a maintained
entrypoint that a developer or user may reasonably run again:

- public API audit generation;
- benchmark data ingestion or fixture preparation;
- legacy dataset contract normalization;
- golden-path or acceptance runners;
- test-gate orchestration.

Ad-hoc scientific validation scripts should not remain in `scripts/`. Promote
them to `validation/<module>/` if they are part of formal evidence, or remove
them when superseded.

## Tests

The test tree intentionally has several layers:

| Test type | Examples | Keep when |
|---|---|---|
| Unit and contract tests | `tests/qc/`, `tests/preprocess/`, `tests/analysis/` | They protect public functions, review summaries, layer contracts, or decision schemas. |
| Benchmark-runner tests | `tests/qc/*benchmark.py`, `tests/preprocess/*benchmark.py` | They exercise formal `validation/` runner internals with tiny fixtures. |
| Integration and golden paths | `tests/integration/` | They verify maintained vertical slices and expected artifacts. |
| Legacy compatibility tests | files with `legacy` in the name | The compatibility alias or missing-metadata degradation path is still promised. |

A test is a cleanup candidate when it only checks a removed API, references a
deleted script, duplicates a stronger contract test without adding a distinct
failure mode, or writes hard-coded outputs into ignored project directories
instead of `tmp_path`.

## Validation

`validation/` is not a dumping ground for old experiments. It should contain
the executable benchmark families that support scLucid's scientific claims:

- QC threshold, tumor-fidelity, doublet, and ambient evidence;
- preprocess layer-contract, HVG preservation, batch diagnostic, and graph
  stability evidence;
- analysis annotation, pseudobulk DE, and proportion consistency evidence.

Generated result directories such as `validation_outputs/qc_threshold_benchmark_v2/`
should be periodically consolidated. Keep the latest evidence package and any
named manuscript/source-data package; remove superseded rerun directories after
their results have been summarized in docs or source-data tables.

## Generated Files

The following are local products and should not be committed:

- `.DS_Store`, `__pycache__/`, `.pytest_cache/`, `.ruff_cache/`, coverage files;
- `results/`, `override/`, `test_results/`, and similar temporary output dirs;
- `validation_outputs/` unless a future release explicitly promotes a frozen
  source-data package.

Use the maintained cleanup entrypoint before reviewing status or preparing a
commit:

```bash
python scripts/clean_workspace.py
```

The default cleanup removes deterministic cache files and `override/`. To also
remove ignored local output directories such as `results/`, run:

```bash
python scripts/clean_workspace.py --include-outputs
```
