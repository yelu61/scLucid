Quick Start Guide
=================

This guide will help you get started with scLucid in 5 minutes.

Basic Workflow
--------------

scLucid provides a complete pipeline for single-cell RNA-seq analysis:

1. **Quality Control (QC)** - Filter low-quality cells and detect doublets
2. **Preprocessing** - Normalize, find HVGs, scale, and run PCA
3. **Analysis** - Clustering, annotation, and differential expression

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import scLucid

    # 1. Quality Control
    from scLucid.qc import QCWorkflowConfig, run_standard_qc

    qc_config = QCWorkflowConfig(
        min_genes=200,
        pc_mt=20.0,
        doublet_method="scrublet"
    )
    adata = run_standard_qc(adata, config=qc_config)

    # 2. Preprocessing
    from scLucid.preprocess import WorkflowConfig, run_preprocessing

    preprocess_config = WorkflowConfig(
        target_sum=1e4,
        n_top_genes=2000,
        n_pcs=50
    )
    adata = run_preprocessing(adata, config=preprocess_config)

    # 3. Clustering and Annotation
    from scLucid.analysis import ClusteringConfig, AnnotationConfig

    clustering_config = ClusteringConfig(method="leiden", resolution=1.0)
    from scLucid.analysis import cluster_cells, annotate_clusters

    adata = cluster_cells(adata, config=clustering_config)
    adata = annotate_clusters(adata, marker_manager)

Configuration System
---------------------

scLucid uses **Pydantic** for automatic configuration validation:

.. code-block:: python

    from scLucid.qc import QCThresholds

    # Valid config - validates automatically
    thresholds = QCThresholds(min_genes=200, pc_mt=20.0)

    # Invalid config - raises ValidationError
    try:
        invalid = QCThresholds(min_genes=-1)  # Error: must be >= 0
    except Exception as e:
        print(f"Validation error: {e}")

Save and Load Configs
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Save config to file
    thresholds.to_json_file("qc_config.json")

    # Load config from file
    loaded = QCThresholds.from_json_file("qc_config.json")

Global Settings
---------------

Configure runtime behavior:

.. code-block:: python

    from scLucid.config import set_config, config_context

    # Permanent change
    set_config(
        n_jobs=4,           # Use 4 CPU cores
        verbosity=2,        # DEBUG level logging
        figure_dpi=300      # High-res figures
    )

    # Temporary change (within context)
    with config_context(n_jobs=8):
        # Uses 8 cores temporarily
        run_preprocessing(adata, config)

    # Original settings restored after context

Common Tasks
------------

Calculate QC Metrics
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from scLucid.qc import calculate_qc_metric

    adata = calculate_qc_metric(
        adata,
        config=qc_config.reporting,
        species="human"
    )

    # Metrics stored in adata.obs:
    # - n_genes_by_counts
    # - total_counts
    # - pct_counts_mt
    # - pct_counts_ribo

Filter Cells
~~~~~~~~~~~~

.. code-block:: python

    from scLucid.qc import filter_cells, QCThresholds

    thresholds = QCThresholds(
        min_genes=200,
        max_genes=5000,
        pc_mt=20.0
    )

    adata_filtered = filter_cells(adata, thresholds)

Detect Doublets
~~~~~~~~~~~~~~~

.. code-block:: python

    from scLucid.qc import predict_doublets, DoubletConfig

    doublet_config = DoubletConfig(
        method="scrublet",
        expected_doublet_rate=0.06
    )

    adata = predict_doublets(adata, config=doublet_config)

Normalize Data
~~~~~~~~~~~~~~

.. code-block:: python

    from scLucid.preprocess import normalize_data, NormalizationConfig

    norm_config = NormalizationConfig(
        target_sum=1e4,
        output_layer="normalized"
    )

    adata = normalize_data(adata, config=norm_config)

    # Result in adata.layers["normalized"]

Find Highly Variable Genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from scLucid.preprocess import find_hvgs, HVGConfig

    hvg_config = HVGConfig(
        method="seurat",
        n_top_genes=2000
    )

    adata = find_hvgs(adata, config=hvg_config)

    # HVGs marked in adata.var["highly_variable"]

Batch Correction
~~~~~~~~~~~~~~~

.. code-block:: python

    from scLucid.preprocess import batch_correction, IntegrationConfig

    integration_config = IntegrationConfig(
        method="harmony",
        batch_key="sample"
    )

    adata = batch_correction(adata, config=integration_config)

What's Next
-----------

- :doc:`../notebooks/01_quality_control` - Interactive QC tutorial
- :doc:`../notebooks/02_preprocessing` - Interactive preprocessing tutorial
- :doc:`../notebooks/03_clustering_annotation` - Clustering and annotation tutorial
- :doc:`../notebooks/04_advanced_topics` - Advanced features and optimization
- :doc:`best_practices` - Recommended practices
- :doc:`api/qc` - QC API reference
- :doc:`api/preprocess` - Preprocessing API reference
- :doc:`api/analysis` - Analysis API reference

Getting Help
------------

- **GitHub Issues**: https://github.com/yourusername/scLucid/issues
- **Documentation**: https://sclucid.readthedocs.io/
- **Examples**: See ``examples/`` directory in the repository
