# QC + Preprocess Evidence Pilot

This report summarizes the first real-data evidence runners for QC Phase 2 and
Preprocess Phase 3. The current outputs are pilot tables, not final manuscript
figures.

## QC Evidence

### Threshold Benchmark

Script:

```bash
python validation/qc/run_threshold_benchmark.py
```

Outputs:

- `validation_outputs/qc_threshold_benchmark/qc_threshold_decision_table.tsv`
- `validation_outputs/qc_threshold_benchmark/qc_retention_summary.tsv`
- `validation_outputs/qc_threshold_benchmark/qc_marker_fidelity.tsv`
- `validation_outputs/qc_threshold_benchmark/qc_strategy_scorecard.tsv`
- `validation_outputs/qc_threshold_benchmark/figure2_threshold_data.tsv`

Evidence role:

- Compare Scanpy fixed threshold, Seurat fixed threshold, scLucid adaptive, and
  scLucid tumor-aware policies.
- Produce reviewer-facing threshold rows with recommended/applied/source,
  confidence, evidence, review-required status, affected cells, biological
  guardrail, and risk note.
- Score each strategy by overall retention, stratified retention fairness,
  marker fidelity, tumor mt-removal safety, and a composite rank.
- Attach strategy-level evidence back to every threshold decision row:
  `strategy_rank`, `recommended_policy`, `strategy_composite_score`,
  `strategy_risk_note`, and `decision_narrative`.

Pilot observations:

- On normal/PBMC datasets, Scanpy fixed thresholds can remain a strong or
  equivalent baseline. scLucid should not claim universal superiority there.
- On tumor datasets, the refined tumor-aware policy now behaves differently
  from generic adaptive QC: it uses a more permissive high-mt guardrail and
  emphasizes biological review rather than mechanical deletion.
- In the current pilot, scLucid tumor-aware ranks first for PDAC and NSCLC in
  the threshold scorecard, ties Scanpy on CRC, and Seurat 5% mt is consistently
  flagged as high-risk on tumor datasets.
- The threshold decision table now reads as a reviewer table rather than a raw
  threshold dump: each threshold records whether its parent policy is the
  benchmark-selected recommendation and why that policy ranked where it did.

### Tumor Biological Fidelity

Script:

```bash
python validation/qc/run_tumor_biological_fidelity_benchmark.py
```

Outputs:

- `validation_outputs/qc_tumor_fidelity/tumor_marker_retention.tsv`
- `validation_outputs/qc_tumor_fidelity/tumor_program_retention.tsv`
- `validation_outputs/qc_tumor_fidelity/sample_celltype_retention_bias.tsv`
- `validation_outputs/qc_tumor_fidelity/tumor_qc_strategy_scorecard.tsv`
- `validation_outputs/qc_tumor_fidelity/tumor_qc_biological_fidelity_narrative.tsv`
- `validation_outputs/qc_tumor_fidelity/figure2_tumor_fidelity_data.tsv`

Pilot observations:

- In the PDAC pilot subset, Seurat-style 5% mitochondrial filtering removed
  many more high-mt cells than Scanpy 20% or scLucid adaptive/tumor-aware
  thresholds.
- The PDAC sample-retention table flags sample-specific loss under fixed
  thresholds, which is exactly the sample-bias failure mode tumor-aware QC must
  make visible.
- The tumor program table records whether high-mt removed cells retain
  epithelial, cell-cycle, hypoxia/stress, or stromal programs, so high-mt cells
  are reviewed biologically rather than discarded mechanically.
- The tumor scorecard ranks strategies by group retention, program retention,
  high-mt removed-cell signal, and biological harm risk. In the pilot,
  tumor-aware QC is strongest in NSCLC and comparable to Scanpy 20% in PDAC/CRC;
  Seurat 5% mt is repeatedly flagged for biological harm risk.
- The biological-fidelity narrative table translates the scorecard into one
  row per dataset/strategy with recommended-policy status, strategy rank,
  high-mt program/marker signal, worst retained sample or cell type, review
  requirement, and a reviewer-facing decision narrative. Figure 2 data now
  includes a `2E` panel for tumor QC review-required/narrative context.

### Doublet Evidence

Script:

```bash
python validation/qc/run_doublet_evidence_benchmark.py
```

Outputs:

- `validation_outputs/qc_doublet_evidence/doublet_evidence.tsv`
- `validation_outputs/qc_doublet_evidence/doublet_method_overlap.tsv`
- `validation_outputs/qc_doublet_evidence/doublet_ground_truth_label_counts.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_reference.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_reference_by_group.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_disagreement_cells.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_seed_stability.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_group_stability.tsv`
- `validation_outputs/qc_doublet_evidence/scdblfinder_python_vs_r_full_group_concentration.tsv`
- `validation_outputs/qc_doublet_evidence/doublet_threshold_calibration.tsv`
- `validation_outputs/qc_doublet_evidence/doublet_algorithm_weight_recommendation.tsv`
- `validation_outputs/qc_doublet_evidence/doublet_benchmark_report_summary.json`

Evidence role:

- Uses Kang 2018 demuxlet singlet/doublet labels as external evidence.
- Reports ambiguous `ambs` cells separately from singlet/doublet metrics.
- Compares transparent heuristic fallbacks, Scrublet, and
  the scLucid `pyscdblfinder` path when dependencies are available. Each method
  row records availability, runtime, predicted rate, precision, recall, F1, AUC,
  fallback status, and review-required state.
- Adds explicit algorithm-plus-heuristic fusion rows, e.g.
  `scdblfinder_python_pyscdblfinder_plus_heuristic_w0.70`, so scLucid's
  production strategy is benchmarked separately from algorithm-only and
  heuristic-only baselines.
- Scans candidate `algorithm_weight` values for weighted-average fusion and
  writes a recommendation table with F1/precision/recall/AUC deltas versus the
  algorithm-only baseline.
- Bioconductor scDblFinder is treated as an external reference, not as another
  in-module method. Provide a CSV from the R package with
  `--r-scdblfinder-reference` to generate
  `scdblfinder_python_vs_r_reference.tsv`, which reports score correlation,
  prediction agreement, Jaccard overlap, AUC/F1 delta, and review status.

### Figure 2 Evidence Package

Script:

```bash
python validation/qc/build_figure2_qc_evidence_package.py
```

Outputs:

- `validation_outputs/qc_figure2_package/figure2_qc_source_data.tsv`
- `validation_outputs/qc_figure2_package/qc_claim_scorecard.tsv`
- `validation_outputs/qc_figure2_package/qc_dataset_coverage.tsv`
- `validation_outputs/qc_figure2_package/qc_evidence_report.md`

Evidence role:

- Consolidates threshold, tumor-aware, doublet, and ambient validation outputs
  into a single Figure 2 source-data table.
- Harmonizes panels so `2A` is workflow/ambient contract, `2B` is threshold
  decision quality, `2C` is tumor program fidelity, `2D` is doublet evidence,
  and `2E` is reviewer-facing tumor narrative.
- Emits a claim-level scorecard that distinguishes supported claims from
  partial or contract-only claims, preventing overstatement of current QC
  accuracy/advancedness.

Pilot observations:

- On the Kang 2018 stratified subset, the scLucid `pyscdblfinder` path
  outperformed the fallback heuristics and Scrublet by F1 and AUC. Scrublet had
  good score AUC but very low recall under the current thresholding, so
  threshold calibration is the next improvement target.
- Bioconductor scDblFinder 1.24.10 has been installed locally and same-cell
  Kang 2018 references have been generated for both 6,000-cell stratified
  subsets and the full Kang dataset. On full Kang, parity passes the current
  review target: score Spearman 0.814, prediction agreement 96.44%, prediction
  Jaccard 0.748, Python AUC 0.909 vs R AUC 0.881, and Python F1 0.668 vs R F1
  0.628.
- The 6,000-cell seed-stability run uses seeds 11, 23, and 37. Prediction
  agreement is stable at 95.98-96.50%, Jaccard is 0.720-0.749, Python AUC is
  0.885-0.888, and Python F1 is 0.613-0.633. Spearman is 0.790-0.810, so two
  6,000-cell subsets remain review-required by the strict 0.8 score-correlation
  threshold even though classification agreement is stable. This identifies
  score-scale calibration as the remaining parity gap rather than a broad
  classification failure.
- The parity runner now also writes group-level and cell-level disagreement
  tables. In full Kang, disagreement is not uniform: the highest rates occur in
  cell subtype 2, donors 1085/1249/1154, sample A/B/C, and dendritic cells.
  Across the three 6,000-cell seeds, dendritic cells are consistently the
  highest-disagreement cell type, while sample and donor groups explain more
  of the low-correlation structure than condition alone. This means the higher
  Python AUC should be interpreted as a reproducible demuxlet-ranking advantage,
  with remaining implementation review focused on donor/sample-specific score
  calibration.
- Standard QC review summaries now include `doublet_evidence_summary`, so
  doublet predictions, score ranges, heterotypic/homotypic risk metadata, and
  external-evidence notes surface in the normal QC report path instead of living
  only in validation outputs.
- The threshold-calibration table shows why Scrublet needs calibration rather
  than replacement: in the 6,000-cell Kang pilot, Scrublet has useful ranking
  signal (AUC 0.818) but its default call set has recall 0.023. A
  benchmark-calibrated threshold for recall >=0.5 raises recall to 0.525 and F1
  from 0.044 to 0.496. The scLucid pyscdblfinder path is already close to its
  default F1 optimum, with additional recall only available by trading away
  precision.
- `doublet_benchmark_report_summary.json` compresses best method, Python/R
  parity, top disagreement groups, and threshold-calibration warnings into a
  payload that can be attached to
  `adata.uns["sclucid"]["qc"]["doublet_benchmark_evidence"]` and surfaced in
  ordinary QC reports.
- Heuristic methods remain useful as transparent fallback and risk
  decomposition evidence, not as the primary doublet benchmark claim.

### Ambient / Empty Droplet Evidence

Script:

```bash
python validation/qc/run_ambient_evidence_benchmark.py
```

Outputs:

- `validation_outputs/qc_ambient_evidence/ambient_evidence.tsv`
- `validation_outputs/qc_ambient_evidence/ambient_barcode_class_summary.tsv`

Evidence role:

- Validates empty-droplet diagnostic contracts on CellBender tiny.
- Does not claim full ambient RNA correction performance; a larger raw 10x
  matrix is still needed before making performance claims.

## Preprocess Evidence

### Layer Contract

Script:

```bash
python validation/preprocess_analysis/run_layer_contract_benchmark.py
```

Outputs:

- `validation_outputs/preprocess_layer_contract/preprocess_input_contract.tsv`
- `validation_outputs/preprocess_layer_contract/layer_contract_report.tsv`
- `validation_outputs/preprocess_layer_contract/marker_panel_coverage.tsv`
- `validation_outputs/preprocess_layer_contract/batch_diagnostic_inputs.tsv`

Evidence role:

- Checks count-layer availability, `.X` expectations, `.raw`, `obsm`, dataset
  metadata, candidate batch keys, and marker-panel coverage before preprocessing.

### HVG Marker / Program Preservation

Script:

```bash
python validation/preprocess_analysis/run_hvg_marker_preservation_benchmark.py
```

Outputs:

- `validation_outputs/preprocess_hvg_preservation/hvg_marker_preservation.tsv`
- `validation_outputs/preprocess_hvg_preservation/program_gene_retention.tsv`
- `validation_outputs/preprocess_hvg_preservation/hvg_strategy_summary.tsv`
- `validation_outputs/preprocess_hvg_preservation/hvg_set_overlap.tsv`
- `validation_outputs/preprocess_hvg_preservation/figure3_hvg_data.tsv`

Pilot observations:

- Standard variance HVG can drop lineage markers or tumor/stress/cell-cycle
  program genes in several datasets.
- The benchmark now compares standard HVG, custom marker/program masks, direct
  standard selection, union, intersection, overlap-only auto, semantic
  protected-auto, and scLucid budget-preserving retained HVGs.
- In the six-dataset pilot, standard variance HVG retained about 45% of lineage
  marker panels and 49% of tumor/program panels on average. The
  marker/program-retained strategy retained 100% of present marker/program
  genes while keeping the 2,000-gene HVG budget.
- Union and semantic protected-auto also retained 100% of marker/program genes
  but expanded the HVG set to about 2,024 genes on average. This is useful when
  maximum biological coverage matters more than a fixed feature budget.
- The overlap-only auto rule selected intersection for all six datasets because
  curated marker/program sets have very low Jaccard overlap with 2,000-gene
  variance HVGs. That result directly motivated the API addition of
  `set_roles`, so protected marker/program masks are treated as biology to
  retain rather than as an ordinary disagreeing HVG method.

### Batch Correction Diagnostic

Script:

```bash
python validation/preprocess_analysis/run_batch_correction_diagnostic_benchmark.py
```

Outputs:

- `validation_outputs/preprocess_batch_diagnostic/batch_diagnostic_summary.tsv`
- `validation_outputs/preprocess_batch_diagnostic/batch_mixing_vs_biology.tsv`
- `validation_outputs/preprocess_batch_diagnostic/overcorrection_risk.tsv`
- `validation_outputs/preprocess_batch_diagnostic/batch_method_comparison.tsv`
- `validation_outputs/preprocess_batch_diagnostic/figure3_batch_data.tsv`

Pilot observations:

- Pancreas donor/sample keys are correction candidates with modest biology
  association.
- NSCLC/CRC sample keys can be strongly confounded with condition, so they are
  diagnostic-only unless a downstream embedding comparison shows acceptable
  biological conservation.
- Kang sample is confounded with stimulation, so correction should protect
  condition structure rather than defaulting to integration.
- The method-comparison table records no correction, Harmony, BBKNN, and scVI
  status. In the current environment BBKNN is dependency-missing; Harmony and
  scVI run on pilot subsets. scVI is flagged for overcorrection risk on Kang
  when label silhouette becomes negative.

### PCA / Graph Handoff Stability

Script:

```bash
python validation/preprocess_analysis/run_graph_stability_benchmark.py
```

Outputs:

- `validation_outputs/preprocess_graph_stability/pca_neighbors_stability.tsv`
- `validation_outputs/preprocess_graph_stability/clustering_seed_stability.tsv`
- `validation_outputs/preprocess_graph_stability/rare_population_preservation.tsv`
- `validation_outputs/preprocess_graph_stability/figure3_graph_data.tsv`

Evidence role:

- Records PCA variance captured by candidate `n_pcs` settings.
- Compares clustering stability across random seeds.
- Compares neighbor overlap between low/high PC settings.
- Flags rare annotated populations whose cluster concentration proxy is low.

## Remaining Work

- Add homotypic/synthetic or HTO doublet evidence beyond Kang demuxlet's
  genotype-detectable donor-doublet reference.
- Promote graph/PCA/neighbors/clustering stability from pilot proxy metrics to
  workflow-level evidence generated by `run_preprocessing`.
- Add real full-size ambient RNA validation data before claiming ambient
  correction performance.
- Build plotting and visual QA around the stabilized Figure 2 QC source-data
  package and the existing Figure 3 preprocess source tables.
