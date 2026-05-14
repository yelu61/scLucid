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

It does **not** validate:

- superiority over standard workflows
- optimal biological filtering thresholds
- cross-dataset scientific accuracy
- publication-level benchmark conclusions

Artifacts
---------

Golden paths write validation outputs under ``<output_dir>/validation/``:

- ``qc_preprocess_validation.json``
- ``qc_preprocess_validation_table.csv``

The JSON includes the full scaffold manifest. The CSV is a compact review table
with one row per metric, including status and interpretation.

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
