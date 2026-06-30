Quick Start
===========

This page shows the recommended minimal path for using scLucid on a new dataset.
For longer runnable scripts, see :doc:`examples`. For full narrative analyses,
see :doc:`notebooks`.

Recommended Learning Order
--------------------------

1. Start with ``scLucid.run_pipeline()`` for the supported QC -> preprocessing -> analysis path
2. Inspect the review summaries stored under ``adata.uns["sclucid"]``
3. Drop down to stage-specific functions when you need explicit control
4. Export reviewer-facing summaries before making biological claims

Minimal End-To-End Example
--------------------------

.. code-block:: python

    import scanpy as sc
    import scLucid as scl

    adata = sc.read_h5ad("data/pbmc3k.h5ad")
    adata.layers["counts"] = adata.X.copy()

    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        dataset_type="pbmc_or_blood",
        species="human",
        qc_save_dir="results/qc",
        preprocess_save_dir="results/preprocess",
        show_progress=True,
    )

    adata.write("results/final_annotated.h5ad")

What This Path Gives You
------------------------

- QC trace and ``qc_reviewer_table`` under ``adata.uns["sclucid"]["qc"]``
- QC review sidecars when ``save_dir`` is set
- standard preprocessing outputs such as normalized layers, HVG metadata, PCA,
  and neighbors/UMAP
- preprocessing layer contract and reviewer table describing
  ``counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph``
- clustering labels in ``adata.obs``
- annotation evidence, analysis output contract, inference policy, and
  analysis reviewer table in ``adata.obs`` and ``adata.uns``

Light Default, Optional Enhancements
------------------------------------

The recommended path is intentionally light by default. It avoids mandatory R
dependencies and avoids aggressive correction steps unless the data provide a
reason to use them.

Default preprocessing:

- filter genes detected in too few cells after QC and before normalization
- normalize raw counts with library-size normalization and ``log1p``
- select HVGs with ``flavor="auto"``, which resolves to dependency-light
  ``seurat`` on log-normalized inputs
- scale, run PCA, build neighbors, and compute UMAP
- skip regression and batch correction unless explicitly enabled

Optional enhancements:

- ``normalization.method="scran"`` keeps Scanpy's
  ``scanpy.external.pp.scran_normalize`` path for users who already have a
  working R/scran environment
- ``hvg.flavor="seurat_v3"`` can be used on raw-count inputs when
  ``scikit-misc`` is installed
- ``run_integration=True`` with Harmony/scVI/scANVI/BBKNN/ComBat should be used
  only after inspecting batch effects and over-correction risk

Choosing Between Default And Intelligent Preprocessing
------------------------------------------------------

Use `PreprocessingWorkflowConfig.default()` when:

- you want the canonical light-dependency package path
- your dataset is standard scRNA-seq with familiar batch structure
- you value stability, signal preservation, and simplicity over parameter search

Use `run_intelligent_preprocessing()` when:

- you want data-driven parameter suggestions
- you want a reviewer-facing summary before applying recommendations
- you need help choosing HVG / PCA / neighbors / integration settings

Related Repository Entry Points
-------------------------------

- ``examples/01_workflow/basic_pipeline.py``: shortest maintained workflow-layer script
- ``examples/02_simple_api/qc_step_by_step.py``: composable QC inspection path
- ``examples/02_simple_api/preprocess_step_by_step.py``: composable preprocessing path
- ``examples/03_advanced_notebooks/``: full notebook analyses with richer outputs
- ``scripts/run_pbmc_golden_path.py``: real-data acceptance script with manifest output

When To Use Stage-Specific Functions
------------------------------------

Use ``run_standard_qc()``, ``run_preprocessing()``, ``cluster_cells()``, and
``run_annotation()`` directly when you are building a manuscript workflow,
testing a single module, or overriding a specific parameter family. The unified
pipeline is the recommended first screen; stage-specific functions are the
expert path.

See :doc:`usage_layers` for the full product-layer model and :doc:`qc_preprocess_maturity`
for the QC/preprocessing hardening standard.
