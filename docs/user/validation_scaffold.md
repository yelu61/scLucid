# Validation Scaffold

scLucid uses three coordinated validation tracks. The current lightweight
scaffold locks down QC/preprocess workflow maturity without claiming that
scLucid is scientifically superior to Scanpy, Seurat, scran, or other
standard workflows. Real-project product acceptance runs in parallel on
active projects, while formal comparative scientific validation starts once
the affected Analysis/Tumor slice has a stable evidence contract.

## Current Scope

The current scaffold validates whether a golden-path run is:

- auditable: review summaries and warning counts are present
- reproducible: input/final shapes and retention are recorded
- preprocessing-ready: count, normalized, raw, HVG, PCA, graph, and UMAP
  state can be inspected
- ready for later comparative validation

For QC, the lightweight scaffold now also has a QC evidence
package under `validation_outputs/qc_evidence_package/`. This package
consolidates existing threshold, tumor-aware, doublet, and ambient
validation outputs into a reviewable source-data table and a claim
scorecard.

It does **not** validate:

- superiority over standard workflows
- optimal biological filtering thresholds
- cross-dataset scientific accuracy
- publication-level benchmark conclusions

The claim scorecard should be used to keep these limitations explicit.
For example, Kang demuxlet labels support heterotypic donor-doublet
evidence but do not fully validate homotypic doublets, while CellBender
tiny validates ambient diagnostic plumbing rather than
ambient-correction performance.

## Artifacts

Golden paths write validation outputs under `<output_dir>/validation/`:

- `qc_preprocess_validation.json`
- `qc_preprocess_validation_table.csv`

The JSON includes the full scaffold manifest. The CSV is a compact
review table with one row per metric, including status and
interpretation.

QC evidence runners write current Phase 2 outputs under
`validation_outputs/qc_*`:

- `qc_evidence_package/qc_source_data.tsv`: harmonized QC
  source data for QC threshold decisions, tumor biological fidelity,
  doublet evidence, and ambient contract checks
- `qc_evidence_package/qc_claim_scorecard.tsv`: claim-level status table
  for QC auditability, tumor-aware biological fidelity, doublet
  calibration, ambient diagnostic contract, and dataset coverage
- `qc_evidence_package/qc_dataset_coverage.tsv`: dataset role and QC
  panel coverage
- `qc_evidence_package/qc_evidence_report.md`: compact reviewer-facing
  summary

The Analysis inference-contract runner writes local, generated evidence under
`validation_outputs/analysis_inference_contract/`:

- `metadata_propagation_matrix.tsv`: resolved metadata names and statistical
  roles across context, proportion, and pseudobulk
- `real_data_design_audit.tsv`: dataset-level READY/BLOCKED design facts
- `pbmc_proportion_statistics.tsv`: paired donor-level CLR results
- `pbmc_pseudobulk_de.tsv`: paired donor-blocked logCPM model results
- `analysis_inference_evidence_manifest.json`: input hashes, selected genes and
  cell types, design audits, artifact paths, and the gate status

Run it against the local real-data inventory:

``` bash
python validation/analysis/run_inference_contract_benchmark.py \
  --data-dir /path/to/scLucid/data
```

Kang2018 PBMC is executable only after defining one aggregation observation per
`donor × condition`; `donor` remains both the independent experimental unit and
paired block. Lin2020 PDAC is intentionally recorded as BLOCKED for condition
inference because it contains one observed condition and no usable cell-type
labels. The runner does not invent a control group or annotation.

## Programmatic Use

``` python
import scLucid as scl

validation = scl.ut.build_qc_preprocess_validation(
    adata,
    run_manifest=manifest,
    dataset_role="pbmc_baseline",
    workflow_name="pbmc3k_golden_path",
)
scl.ut.write_validation_outputs(validation, "results/golden/pbmc3k/validation")
```

## Recommended Timing

Run the three tracks in parallel rather than waiting for every module to be
fully redesigned:

1. contract validation: keep QC/preprocess and Analysis review schemas
   executable in CI;
2. product acceptance: use AK112/LPJ-like active projects to record parameter
   friction, misunderstood outputs, analyst overrides, and rerun scope;
3. scientific validation: extend `qc_preprocess_analysis_validation` through
   PBMC, PDAC, and a second tumor cohort, then add optional external workflow
   comparisons.

When a real-project or benchmark run exposes a blocker, improve that specific
QC, preprocess, analysis, or tumor boundary and rerun the same vertical slice.
Do not postpone usability validation until all four modules are considered
complete.
