# scLucid QC + Preprocess Validation

This directory contains benchmark scaffolds for turning local real datasets into
reviewable evidence. The goal is not to prove that every default is optimal.
The goal is to make each QC and preprocessing decision inspectable, comparable,
and biologically defensible.

## Dataset Evidence Map

| Dataset | QC evidence | Preprocess evidence | Why it matters |
|---|---|---|---|
| `pbmc3k` | Fixed-threshold baseline, fast smoke path | Layer contract, Scanpy-style baseline | Shows the default path behaves on a familiar non-tumor dataset. |
| `lin2020.pdac` | Tumor-aware QC, sample retention bias | Tumor marker preservation, PDAC golden path | First tumor case for high-mito/high-stress warnings and malignant/TME signal preservation. |
| `schlesinger2020.pdac` | Tumor-aware generalization, single-sample behavior | Tumor marker preservation | Prevents the PDAC claim from depending on one fixture. |
| `zilionis2019.nsclc` | Tumor/blood retention, cell-type retention | TME marker/state preservation | Second cancer type with paired blood control and author labels. |
| `lee2020.crc` | Tumor/normal retention, patient-level retention bias | Patient-aware integration diagnostics, marker preservation | Large patient-diverse tumor dataset for retention and overcorrection evidence. |
| `baron2016.pancreas` | Normal reference | Donor batch diagnostic, integration risk | Tests whether scLucid recommends correction conservatively instead of flipping a default switch. |
| `kang2018.pbmc` | Demuxlet singlet/doublet/ambiguous evidence | Donor/stimulation-aware preprocessing | Makes doublet claims measurable and perturbation structure visible. |
| `cellbender_tiny` | Ambient RNA and empty-droplet diagnostic contract | Not a preprocess benchmark | Tiny fixture for ambient diagnostic plumbing, not biological performance claims. |

## Phase 2 QC Benchmark

Primary outputs should live under `validation_outputs/qc_*` and eventually feed
the QC evidence package:

- `qc_evidence_package/qc_source_data.tsv`: unified QC source-data
  table with harmonized panels (`2A` ambient contract, `2B` threshold decision
  quality, `2C` tumor program fidelity, `2D` doublet evidence, `2E`
  reviewer-facing tumor narrative).
- `qc_evidence_package/qc_claim_scorecard.tsv`: claim-level evidence status for
  QC auditability, tumor-aware biological fidelity, doublet calibration,
  ambient diagnostic contract, and dataset coverage.
- `qc_evidence_package/qc_evidence_report.md`: compact reviewer-oriented report
  that states which claims are supported, partial, or contract-only.
- `adata.uns["sclucid"]["qc"]["review_summary"]["data"]["qc_handoff_readiness"]`:
  QC-to-preprocess contract with the recommended counts layer,
  review/sensitivity cell fractions, decision-column availability,
  safe-to-continue flags, blockers, and required downstream handling.
- `qc_workflow_decision_table.tsv`: reviewer table with recommended, applied,
  source, confidence, evidence, review_required, affected_cells,
  biological_guardrail, strategy_rank, recommended_policy,
  strategy_composite_score, decision_narrative, and risk_note.
- `retention_marker_fidelity.tsv`: retention rate, sample retention bias,
  cell-type retention bias, marker fidelity before/after QC.
- `qc_strategy_scorecard.tsv`: policy-level scorecard for overall retention,
  stratified retention fairness, marker fidelity, tumor mt-removal safety, and
  composite rank.
- `tumor_program_retention.tsv`: malignant/TME marker and program retention in
  PDAC, NSCLC, and CRC.
- `tumor_qc_strategy_scorecard.tsv`: tumor-aware biological fidelity scorecard
  that flags high-mt removed cells carrying tumor/stress/proliferation signal.
- `tumor_qc_biological_fidelity_narrative.tsv`: reviewer-facing dataset/strategy
  narrative with recommended-policy status, strategy rank, high-mt signal,
  worst retained group, review-required status, and decision rationale.
- `qc_doublet_evidence/doublet_evidence.tsv`: Scrublet, pyscdblfinder,
  algorithm-plus-heuristic fusion, transparent heuristic fallback, and demuxlet
  overlap metrics.
- `qc_ambient_evidence/ambient_evidence.tsv`: ambient/empty-droplet diagnostic
  contract summary.
- `adata.uns["sclucid"]["qc"]["review_summary"]["data"]["doublet_evidence_summary"]`:
  report-layer doublet summary with prediction rates, score ranges,
  heterotypic/homotypic risk metadata, and external-evidence notes.
- `scdblfinder_python_vs_r_reference.tsv`,
  `scdblfinder_python_vs_r_reference_by_group.tsv`, and
  `scdblfinder_python_vs_r_disagreement_cells.tsv`: same-cell parity between
  scLucid's `pyscdblfinder` path and Bioconductor scDblFinder reference output.
- `scdblfinder_python_vs_r_seed_stability.tsv`,
  `scdblfinder_python_vs_r_group_stability.tsv`, and
  `scdblfinder_python_vs_r_full_group_concentration.tsv`: random-seed stability
  and donor/sample/cell-type concentration diagnostics for Python/R
  scDblFinder disagreements.
- `doublet_threshold_calibration.tsv`: demuxlet-grounded score-threshold scan
  for methods with continuous scores, used to diagnose high-AUC/low-recall
  behavior such as conservative Scrublet defaults.
- `doublet_algorithm_weight_recommendation.tsv`: demuxlet-grounded comparison
  of algorithm-only versus algorithm-plus-heuristic fusion rows across candidate
  `algorithm_weight` values. This table is the evidence surface for whether
  scLucid should keep the default 0.7 weight, favor algorithm-only behavior, or
  require manual review for a dataset.
- `doublet_benchmark_report_summary.json`: compact best-method, Python/R
  parity, disagreement-group, threshold-calibration, and algorithm-weight
  payload that can be attached to
  `adata.uns["sclucid"]["qc"]["doublet_benchmark_evidence"]` for normal QC
  reports.

Minimum credible Phase 2 comparisons:

- Scanpy fixed threshold.
- Seurat-style fixed threshold.
- scLucid adaptive recommendation.
- scLucid tumor-aware recommendation for tumor datasets.
- Doublet baselines on `kang2018.pbmc`, with `ambs` reported separately from
  singlet/doublet metrics.

Current claim-scorecard interpretation:

- QC decision auditability: supported by threshold decision tables across the
  local benchmark inventory.
- Tumor-aware biological fidelity: supported as a marker/program-retention
  proxy across PDAC, NSCLC, and CRC.
- Doublet calibration: supported on Kang demuxlet with review required because
  demuxlet mainly validates genotype-detectable donor doublets.
- Ambient diagnostics: contract-only until a full raw 10x ambient benchmark is
  added.

## Phase 3 Preprocess Benchmark

Primary outputs should live under `validation_outputs/preprocess_*` and
eventually feed the preprocess evidence package:

- `layer_contract_report.tsv`: layer transition table for counts, normalized,
  log, scaled, PCA, neighbors, UMAP, and `.raw` semantics.
- `hvg_marker_preservation.tsv`: standard HVG, custom marker/program masks,
  direct/union/intersection/auto choices, semantic protected-auto, and
  budget-preserving retained marker/program genes.
- `hvg_strategy_summary.tsv` and `hvg_set_overlap.tsv`: strategy-level marker
  fidelity, program retention, HVG budget, and overlap diagnostics.
- `batch_clustering_stability.tsv`: no-correction vs opt-in correction
  diagnostics, biological conservation, marker fidelity, and overcorrection
  warnings.
- `batch_method_comparison.tsv`: no correction, Harmony, BBKNN, and scVI
  method status plus batch mixing and biological conservation metrics.
- `inference_semantics_guardrails.tsv`: handoff checks that downstream
  DE/proportion outputs are marked exploratory/descriptive/inferential.

Minimum credible Phase 3 comparisons:

- scLucid standard preprocessing with no integration.
- scLucid diagnostic-only integration recommendation.
- Opt-in Harmony on batch-heavy datasets when available.
- Scanpy standard workflow baseline.

## Combined Evidence + Real-Project Feedback Loop

The next evidence milestone is a combined QC/preprocess evidence package. It
should consume the stable review-summary contracts rather than re-derive meaning
from raw workflow internals:

- QC contributes threshold decisions, retention/tumor/doublet/ambient evidence,
  `qc_handoff_readiness`, and claim-level status.
- Preprocess contributes layer contracts, HVG preservation, batch/graph
  diagnostics, `analysis_handoff_readiness`, and downstream-use safety flags.
- Real-project acceptance runs should record where the contract was clear,
  where analysts still needed manual interpretation, and which defaults caused
  avoidable review churn.

## First Scripts

Run these lightweight manifest builders before writing heavier benchmarks:

```bash
python validation/qc/build_qc_benchmark_manifest.py
python validation/preprocess/build_preprocess_benchmark_manifest.py
python validation/qc/run_threshold_benchmark.py
python validation/qc/run_tumor_biological_fidelity_benchmark.py
python validation/qc/run_doublet_evidence_benchmark.py
python validation/qc/run_ambient_evidence_benchmark.py
python validation/preprocess/run_layer_contract_benchmark.py
python validation/preprocess/run_hvg_marker_preservation_benchmark.py
python validation/preprocess/run_batch_correction_diagnostic_benchmark.py
python validation/preprocess/run_graph_stability_benchmark.py
python validation/qc/build_qc_evidence_package.py
python validation/preprocess/build_preprocess_evidence_package.py
python validation/qc_preprocess/build_qc_preprocess_evidence_package.py
```

They verify metadata readiness and write dataset/strategy/Figure-panel plans
under `validation_outputs/`. The executable runners then add the first
real-data evidence tables for QC threshold decisions, tumor biological
fidelity, demuxlet-grounded doublet evidence, ambient/empty-droplet contracts,
preprocessing layer contracts, HVG marker/program preservation, and
batch-correction diagnostic recommendations, plus PCA/neighbor/cluster handoff
stability evidence.
