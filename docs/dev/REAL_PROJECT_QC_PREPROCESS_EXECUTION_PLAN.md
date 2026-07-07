# Real-project QC + preprocess + analysis execution plan

This note records how to use the two local Step1 notebooks as evidence for
hardening scLucid's QC, preprocessing, and downstream analysis review modules.
The notebooks are intentionally not tracked by git; they are project-specific
working references.

## Current evidence sources

Local notebooks:

- `Step1-QC_and_Preprocessing.ipynb`: LPJ/BBB-style real project.
- `Step1-QC_and_Preprocessing_CT26.ipynb`: CT26 tumor project.

Repository evidence:

- QC validation scripts in `validation/qc/`.
- Preprocess validation scripts in `validation/preprocess_analysis/`.
- Real data fixtures in `data/`.
- Existing validation outputs under `validation_outputs/`.

## Initial notebook readout

Both notebooks already use the canonical QC path:

1. Build `QCWorkflowConfig`.
2. Generate expected doublet rates with `scl.qc.generate_doublet_rates`.
3. Run reviewer-first QC with `scl.qc.run_iterative_qc`.
4. Inspect `qc_review_summary`, retention, doublets, benchmark summaries, and
   action items.
5. Save the QC-filtered object.

This means the QC module is close to the real-project path. Remaining QC work
should focus on evidence depth and validation, not another entrypoint rewrite.

The preprocessing path now has two evidence tracks:

1. The project-specific audited path remains notebook-driven after bootstrap:
   `scl.pp.run_preprocessing(..., steps=["gene_filtering", "normalization", "set_raw"])`
   owns normalization and raw/layer setup, followed by manual HVG, regression,
   integration, and embedding review.
2. The canonical module-owned comparison path is now available through
   `scl.pp.run_iterative_preprocessing(...)` and can be enabled in each notebook
   with `RUN_CANONICAL_ITERATIVE_PREPROCESSING_COMPARISON=True`.
3. The recommendation/advisor path is now under `scl.recommendation` and writes
   `recommendation_review/workflow_recommendations.json` for parameter
   comparison.
4. HVG choice and stability can still be reviewed manually with `suggest_hvg_choice`,
   `evaluate_hvg_stability`, `select_and_audit_hvgs`.
5. Final object audit and manual review-summary finalization remain
   notebook-local, but now include recommendation and canonical comparison
   summaries when those checkpoints are run.

The remaining proof gap is no longer the absence of a preprocess entrypoint. It
is whether executed real-project outputs contain enough comparable evidence to
explain QC filtering, doublet handling, HVG strategy, integration decisions, and
final representation choice.

## Execution order

### Phase 1: Notebook-to-module gap matrix

Create a table mapping each notebook block to one of:

- covered by stable API,
- covered by low-level API but not workflow-owned,
- notebook-only logic that should become a module feature,
- project-specific logic that should stay local.

Priority mappings:

| Notebook block | Current API coverage | Product decision |
|---|---|---|
| QC execution | `scl.qc.run_iterative_qc` | Keep as canonical real-project QC path. |
| QC retention/doublet review | QC review summary + notebook audit | Move missing high-signal tables into QC reports if still notebook-only. |
| Preprocess bootstrap | `scl.pp.run_preprocessing` partial steps | Keep for project-specific stepwise audit and resume. |
| Canonical preprocess comparison | `scl.pp.run_iterative_preprocessing` | Execute on a copy and compare against notebook-selected decisions. |
| Recommendation advisor | `scl.recommendation.recommend_analysis_parameters` | Compare suggested parameters with project choices. |
| Dual HVG strategy | `find_hvgs`, `suggest_hvg_choice`, `select_and_audit_hvgs` | Compare manual dual-HVG decisions with canonical workflow evidence. |
| HVG stability | `evaluate_hvg_stability` | Verify stability metrics appear in validation evidence tables. |
| Cell-cycle regression diagnostic | `diagnose_cell_cycle_regression` | Verify regression decisions are explainable in review summaries. |
| Diagnostic pre-integration embedding | `run_embedding_pipeline` | Compare with iterative preprocessing diagnostic embedding summary. |
| Integration decision | `decide_integration`, `batch_correction` | Validate key-decision evidence rows for run/skip and overcorrection risk. |
| Final graph/UMAP optimization | `run_embedding_pipeline` | Validate final representation and graph evidence rows. |
| Final object audit | notebook-local checks + review finalizer | Feed review-summary completeness and decision evidence tables. |

### Phase 2: Execute updated real-project notebooks

Each local Step1 notebook now contains:

- a `scl.recommendation` advisor checkpoint,
- the existing project-specific audited preprocessing path,
- an optional `scl.pp.run_iterative_preprocessing` comparison path controlled by
  `RUN_CANONICAL_ITERATIVE_PREPROCESSING_COMPARISON`,
- final manual review summaries that include recommendation and canonical
  comparison payloads when available.

Run the notebooks inside the real project environments, then keep the final
`.h5ad` objects and sidecar JSON/TSV files under each project's result
directory.

### Phase 3: Preprocess structure cleanup

Current cleanup decisions:

- `src/scLucid/preprocess/trace.py` remains a single file because it is not yet
  structurally complex enough to justify a subpackage.
- `src/scLucid/preprocess/workflow/` is a subpackage because it now contains the
  standard step-control workflow and the iterative reviewer-first workflow.
- `src/scLucid/preprocess/intelligent/` moved to
  `src/scLucid/recommendation/preprocess/`; recommendation logic no longer lives
  under `scl.pp`.
- Transitional aliases have been removed from the preprocess public API.

### Phase 4: Validation proof layer

Use the existing validation inventory as the scientific proof surface.

QC proof should aggregate:

- threshold decision auditability,
- tumor-aware marker/program preservation,
- doublet calibration and demuxlet overlap,
- ambient diagnostic contract and residual/correction evidence,
- retention fairness by sample/cell type/group,
- benchmark scorecard and reviewer action items.

Preprocess proof should aggregate:

- layer contract correctness,
- HVG marker/program preservation,
- HVG strategy overlap and stability,
- batch correction diagnostic recommendations,
- overcorrection risk,
- graph/cluster stability,
- rare population preservation,
- final object completeness.

Analysis proof should aggregate:

- clustering and annotation review evidence,
- pseudobulk-first DE routing and publication-inference validity,
- cell-type proportion/composition method choice,
- inference-level labels for exploratory versus sample-level claims,
- review-summary action items and downstream readiness.

The next validation deliverable is a joint QC+preprocess+analysis evidence package:

```text
validation_outputs/qc_preprocess_real_project/
  real_project_gap_matrix.tsv
  qc_preprocess_evidence_report.md

validation_outputs/qc_preprocess_review_evidence/
  review_summary_completeness.tsv
  key_decision_evidence.tsv
  qc_preprocess_analysis_review_evidence_report.md
```

### Phase 5: Notebook simplification

After executing the comparison path on real data, decide whether the local Step1
notebooks can shrink to:

```python
adata = scl.qc.run_iterative_qc(...)
adata = scl.pp.run_iterative_preprocessing(...)
adata = scl.al.run_standard_analysis(..., run_proportion=True, pseudobulk_first=True)
scl.recommendation.recommend_analysis_parameters(...)
validation/qc_preprocess/build_review_summary_evidence_tables.py ...
```

The notebooks should remain useful for interpretation, visualization, and manual
review, but core logic should live in scLucid modules.

## Immediate next task

Run the updated notebooks in the real project environments, then build evidence
tables from their outputs:

```bash
python validation/qc_preprocess/build_real_project_gap_matrix.py
python validation/qc_preprocess/build_review_summary_evidence_tables.py \
  /path/to/project/results \
  --out-dir validation_outputs/qc_preprocess_review_evidence
```

The first script validates notebook structure and module coverage. The second
script validates executed review summaries and key decisions across QC,
preprocess, recommendation, and analysis.
