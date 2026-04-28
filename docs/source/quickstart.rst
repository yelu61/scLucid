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

- QC trace under ``adata.uns["sclucid"]["qc"]``
- QC review sidecars when ``save_dir`` is set
- standard preprocessing outputs such as normalized layers, HVG metadata, PCA, and neighbors/UMAP
- clustering labels in ``adata.obs``
- annotation evidence and annotation outputs in ``adata.obs`` and ``adata.uns``

Choosing Between Default And Intelligent Preprocessing
------------------------------------------------------

Use `PreprocessingWorkflowConfig.default()` when:

- you want the canonical package path
- your dataset is standard scRNA-seq with familiar batch structure
- you value stability and simplicity over parameter search

Use `run_intelligent_preprocessing()` when:

- you want data-driven parameter suggestions
- you want a reviewer-facing summary before applying recommendations
- you need help choosing HVG / PCA / neighbors / integration settings

Related Repository Entry Points
-------------------------------

- ``examples/quickstart.py``: shortest runnable script
- ``examples/preprocessing.py``: standard path + intelligent path + manual path
- ``examples/annotation_report.py``: reviewer-facing annotation export
- ``notebooks/``: full notebook analyses with richer outputs

When To Use Stage-Specific Functions
------------------------------------

Use ``run_standard_qc()``, ``run_preprocessing()``, ``cluster_cells()``, and
``run_annotation()`` directly when you are building a manuscript workflow,
testing a single module, or overriding a specific parameter family. The unified
pipeline is the recommended first screen; stage-specific functions are the
expert path.
