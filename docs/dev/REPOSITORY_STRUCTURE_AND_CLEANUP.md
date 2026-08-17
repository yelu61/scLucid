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

## P0.5 Repository Consolidation Audit (2026-08-14)

This section is the evidence log for the bounded P0.5 cleanup. It does not
replace the strategic plan or create a parallel backlog.

### Pre-cleanup baseline

| Measure | Count / status |
|---|---:|
| Tracked files | 574 |
| Files present in the worktree, excluding `.git` | 1,180 |
| Package Python files | 227 |
| Test Python files | 134 |
| Documentation Markdown files | 80 |
| Present tracked example notebooks | 5 |
| Script and validation files | 53 |
| P0 focused/contract test gate | 136 passed |
| P0 real PBMC/PDAC acceptance gate | 13 passed |
| P0 safe core gate with local real-data fixtures | 1,321 passed, 6 skipped; four resource tests separately blocked by a stale path contract |

The working tree already contained decision/context/audit/documentation work
and the P0 fixes. Those changes were treated as protected inputs. The deleted
status of `Step1-QC_and_Preprocessing_v2.ipynb` was present before P0.5 and was
not changed.

### Candidate classification

| Candidate | Classification | Evidence and action |
|---|---|---|
| `src/scLucid/{qc,preprocess,analysis,tumor}` and their contract tests | **Keep** | Canonical QC → Preprocess → Analysis → Tumor spine; protected by unit, contract, integration, and acceptance tests. |
| `tools/bulk`, `tools/spatial`, pyMonocle3 and other optional evidence modules | **Keep** | Optional support evidence, not the main product spine. Low coverage alone is not deletion evidence. |
| `docs/marker_resources/*.jsonl` | **Consolidate** | These tracked files are the live resource sources. Code/tests still used the pre-migration `docs/*.jsonl` paths; consumers were corrected instead of copying a second queue. |
| `docs/roadmap/README.md` and `docs/roadmap/index.md` | **Consolidate** | `index.md` is canonical. The former long README was archived for provenance and replaced by a short pointer. |
| Point-in-time QC/API audits in `docs/dev/` | **Keep / Archive marker** | Existing superseded banners distinguish them from current contracts. The real-project execution note received an explicit point-in-time/external-input banner. |
| `analysis.bulk`, bulk legacy aliases, analysis malignancy wrapper, `plotting.main` | **Deprecate / compatibility** | They retain potential external callers and compatibility tests. No callable API was removed in P0.5; decisions remain in `API_TRIAGE.md`. |
| TCGA and CNV-reference placeholders plus reserved bulk limma method | **Experimental / blocked** | All are explicitly triaged; none may be treated as executed evidence. Behavior remains unchanged pending a versioned maintainer decision. |
| Caches, coverage HTML/XML, `override/`, bytecode | **Delete** | Ignored deterministic products identified by `scripts/clean_workspace.py`; 65 targets removed. |
| `build_review_summary_evidence_tables.py.bak` | **Delete** | Untracked older subset of the tracked script. It contained no unique project-mode logic; the formal script also adds benchmark mode. `*.bak` is now ignored. |
| `Step1-QC_and_Preprocessing_v2.ipynb` deletion | **Keep current state / decision required** | The pre-existing deletion is ambiguous. P0.5 neither restored the file nor converted the deletion into a permanent cleanup decision. |

### Test and validation boundaries

| Layer | Canonical location | Contract |
|---|---|---|
| Unit | `tests/qc`, `tests/preprocess`, `tests/analysis`, `tests/tumor`, `tests/tools` | Small behavioral and failure-mode checks. |
| Cross-module contract | `tests/test_context.py`, `tests/test_decision.py`, `tests/test_contracts.py`, `tests/test_semantics_contracts.py` | Context, planning, stage storage, evidence semantics, and inference boundaries. |
| Integration | `tests/integration/` | Executable vertical slices and artifact handoffs. |
| Real-data acceptance | PBMC golden path and Lin2020 PDAC acceptance tests | Maintained real-data outcomes and round-trip artifacts; data remain local and ignored. |
| Scientific validation | `validation/<module>/` with thin runner tests under `tests/` | Benchmark evidence, provenance, and claim calibration. Generated outputs belong in ignored `validation_outputs/`. |

Golden-path scripts and acceptance tests are not duplicates: scripts own
artifact generation, while tests own executable assertions. No runner or
fixture was removed without stronger equivalence evidence.

### Source-of-truth order

For this repository, implementation disputes are resolved in this order:

1. code plus passing contract/acceptance tests;
2. tracked resource sources and executable validation runners;
3. maintained API and user documentation;
4. roadmap and developer audit notes;
5. archived design or curation provenance.

The four authoritative navigation entries remain `README.md`,
`docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md`,
`docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`, and `docs/roadmap/index.md`.

### Post-cleanup outcome

| Measure | Before | After | Interpretation |
|---|---:|---:|---|
| Tracked files | 574 | 574 | No tracked source, resource, test, or user artifact was deleted. |
| Files present in the worktree, excluding `.git` | 1,180 | 582 | Deterministic ignored caches and generated products were removed. |
| Directories present in the worktree, excluding `.git` | not recorded | 93 | Reported for the post-cleanup baseline. |
| Package Python files | 227 | 227 | No package implementation file was removed. |
| Test Python files | 134 | 134 | No test layer was collapsed without equivalence evidence. |
| Documentation Markdown files | 80 | 82 | Archive governance and the preserved roadmap provenance added two documents. |
| Present example notebooks | 5 | 5 | Notebook content was not changed by P0.5. |
| Script and validation files | 53 | 35 | The decrease is cache/bytecode removal; maintained runners and source data were retained. |

The stale marker-resource consumers now use the canonical tracked files under
`docs/marker_resources/`; no resource content was fabricated or duplicated.
The completed gates were 188 focused tests, all four formerly blocked resource
tests, 13 PBMC/PDAC real-data acceptance tests, and the safe core suite with
1,326 passed and 6 skipped (100 explicitly deselected out-of-scope tests).
Core-import smoke checks, strict MkDocs, generated-API consistency, Python
compilation, targeted Ruff checks, and `git diff --check` also passed. These
checks preserve the P0 scientific behavior while tightening repository and
optional-dependency boundaries.
