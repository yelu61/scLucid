# QC Module Benchmark Gap Audit

This note tracks the remaining gap between the current QC module and a
benchmark-grade reference module.

## What Is Now Evidence-Backed

- Threshold policies are compared against Scanpy fixed thresholds, Seurat-style
  fixed thresholds, scLucid adaptive thresholds, and scLucid tumor-aware
  thresholds in `validation/qc/run_threshold_benchmark.py`.
- The threshold benchmark now emits a strategy scorecard with overall retention,
  stratified retention fairness, marker fidelity, tumor mt-removal safety, and a
  composite rank.
- Tumor-aware QC is separately evaluated in
  `validation/qc/run_tumor_biological_fidelity_benchmark.py`, including marker
  retention, tumor/stress/proliferation program retention, high-mt removed-cell
  biology, and sample/cell-type retention bias.
- Doublet evidence is grounded against Kang 2018 demuxlet labels and compares
  transparent heuristic fallbacks with Scrublet and scDblFinder/pyscdblfinder
  when dependencies are available.
- Ambient/empty-droplet support has contract evidence on CellBender tiny, but
  not yet a full ambient correction performance benchmark.
- The Figure 2 evidence package builder
  (`validation/qc/build_figure2_qc_evidence_package.py`) now consolidates
  threshold, tumor-aware, doublet, and ambient outputs into
  `validation_outputs/qc_figure2_package/figure2_qc_source_data.tsv`,
  `qc_claim_scorecard.tsv`, `qc_dataset_coverage.tsv`, and
  `qc_evidence_report.md`.

## Current Pilot Lessons

- scLucid should not claim universal QC superiority on normal/PBMC datasets.
  Scanpy fixed thresholds can be equivalent or stronger there.
- scLucid's real differentiator is tumor-aware biological preservation: avoid
  mechanically deleting high-mt cells when they preserve malignant, stress,
  cell-cycle, or TME programs.
- The original adaptive and tumor-aware policies were too similar. The validation
  runner now tests a clearer tumor-aware policy: more permissive low-quality
  bounds and a high-mt guardrail that does not fall below 20% in tumor contexts.
- Scrublet score AUC can be useful even when thresholded recall is low. The next
  doublet improvement should calibrate thresholds by expected rate and external
  evidence, not treat default binary calls as final truth.

## Remaining QC Gaps

- The refined tumor-aware policy currently lives in validation logic. It should
  be promoted into the canonical QC threshold decision path so `run_standard_qc`
  and notebook/manual APIs use the same policy that the benchmark tests.
- Threshold decision output should expose a single reviewer table per strategy:
  recommended, applied, source, confidence, evidence, affected cells,
  biological-risk note, and final action. The pieces exist, but the user-facing
  narrative is still split between suggestion, sample thresholds, filtering, and
  tumor flags.
- Sample-aware thresholds need clearer boundaries: pooled, hierarchical, and
  independent modes should be reported as policy choices, not mixed with
  threshold estimation details.
- Doublet reporting should distinguish algorithm score quality, thresholded
  call quality, heuristic fallback, heterotypic risk, homotypic risk, and
  external hashing/genotype evidence in one compact table.
- Ambient RNA still needs a real raw 10x dataset with known ambient burden or
  external CellBender/SoupX outputs before scLucid can make correction
  performance claims.
- Figure 2 source data is now stabilized, but plotting code and visual QA still
  need to be built around `validation_outputs/qc_figure2_package/`.

## Recommended Next Code Changes

- Move the validation-tested tumor-aware threshold guardrail into
  `scLucid.qc.filtering` / `scLucid.qc.workflow` as a named policy layer.
- Add a canonical QC policy scorecard helper to `scLucid.qc.benchmark` and reuse
  it from validation scripts and reports.
- Extend `generate_qc_report` to include the strategy scorecard and tumor
  biological fidelity summary when present.
- Add threshold calibration for Scrublet using demuxlet/hashing-aware expected
  rate evidence, with binary calls clearly marked as calibrated or default.
- Add plotting and visual QA for the stable Figure 2 source-data package.
- Add full raw 10x ambient evidence and homotypic/HTO/synthetic doublet
  evidence before expanding accuracy claims.
