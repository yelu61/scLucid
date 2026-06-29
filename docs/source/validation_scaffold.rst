Validation Scaffold
===================

scLucid uses a two-stage validation strategy. The current lightweight scaffold
locks down QC/preprocess workflow maturity without claiming that scLucid is
scientifically superior to Scanpy, Seurat, scran, or other standard workflows.
Formal comparative validation should happen after the analysis module reaches
the same auditability level as QC and preprocessing.

Current Scope
-------------

The current scaffold validates whether a golden-path run is:

- auditable: review summaries and warning counts are present
- reproducible: input/final shapes and retention are recorded
- preprocessing-ready: count, normalized, raw, HVG, PCA, graph, and UMAP state
  can be inspected
- ready for later comparative validation

For QC, the lightweight scaffold now also has a Figure 2 evidence package under
``validation_outputs/qc_figure2_package/``. This package consolidates existing
threshold, tumor-aware, doublet, and ambient validation outputs into a
reviewable source-data table and a claim scorecard.

It does **not** validate:

- superiority over standard workflows
- optimal biological filtering thresholds
- cross-dataset scientific accuracy
- publication-level benchmark conclusions

The claim scorecard should be used to keep these limitations explicit. For
example, Kang demuxlet labels support heterotypic donor-doublet evidence but do
not fully validate homotypic doublets, while CellBender tiny validates ambient
diagnostic plumbing rather than ambient-correction performance.

Artifacts
---------

Golden paths write validation outputs under ``<output_dir>/validation/``:

- ``qc_preprocess_validation.json``
- ``qc_preprocess_validation_table.csv``

The JSON includes the full scaffold manifest. The CSV is a compact review table
with one row per metric, including status and interpretation.

QC evidence runners write current Phase 2 outputs under
``validation_outputs/qc_*``:

- ``qc_figure2_package/figure2_qc_source_data.tsv``: harmonized Figure 2 source
  data for QC threshold decisions, tumor biological fidelity, doublet evidence,
  and ambient contract checks
- ``qc_figure2_package/qc_claim_scorecard.tsv``: claim-level status table for
  QC auditability, tumor-aware biological fidelity, doublet calibration, ambient
  diagnostic contract, and dataset coverage
- ``qc_figure2_package/qc_dataset_coverage.tsv``: dataset role and Figure 2
  panel coverage
- ``qc_figure2_package/qc_evidence_report.md``: compact reviewer-facing summary

Programmatic Use
----------------

.. code-block:: python

   import scLucid as scl

   validation = scl.ut.build_qc_preprocess_validation(
       adata,
       run_manifest=manifest,
       dataset_role="pbmc_baseline",
       workflow_name="pbmc3k_golden_path",
   )
   scl.ut.write_validation_outputs(validation, "results/golden/pbmc3k/validation")

Recommended Timing
------------------

Use this scaffold now to stabilize QC/preprocess maturity claims. After
analysis has comparable review-summary and evidence contracts, extend the
validation layer into ``qc_preprocess_analysis_validation`` with PBMC, PDAC,
cross-dataset tumor validation, and optional external workflow comparisons.
