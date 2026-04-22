Quick Start
===========

This page shows the recommended minimal path for using scLucid on a new dataset.
For longer runnable scripts, see :doc:`examples`. For full narrative analyses,
see :doc:`notebooks`.

Recommended Learning Order
--------------------------

1. Run QC with `run_standard_qc()`
2. Run preprocessing with `PreprocessingWorkflowConfig.default()` or intelligent preprocessing
3. Run clustering and annotation
4. Export reviewer-facing summaries before making biological claims

Minimal End-To-End Example
--------------------------

.. code-block:: python

    import scanpy as sc
    import scLucid as scl
    from scLucid.qc import run_standard_qc, QCWorkflowConfig
    from scLucid.preprocess import run_preprocessing, PreprocessingWorkflowConfig
    from scLucid.analysis import cluster_cells, run_annotation, ClusteringConfig, AnnotationConfig

    adata = sc.read_h5ad("data/pbmc3k.h5ad")
    adata.layers["counts"] = adata.X.copy()

    qc_config = QCWorkflowConfig(
        save_dir="results/qc",
        use_recommendations=True,
        threshold_mode="hierarchical",
    )
    adata = run_standard_qc(adata, config=qc_config)

    preprocess_config = PreprocessingWorkflowConfig.default(save_dir="results/preprocess")
    adata = run_preprocessing(adata, config=preprocess_config)

    adata = cluster_cells(
        adata,
        ClusteringConfig(method="leiden", resolution=1.0),
    )

    adata = run_annotation(
        adata,
        config=AnnotationConfig(
            cluster_key="leiden_clusters",
            marker_species="human",
            final_method="combined",
        ),
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
