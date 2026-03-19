Best Practices
==============

This guide covers recommended practices for using scLucid effectively.

Configuration Management
------------------------

**Use Pydantic Configs for Reproducibility**

.. code-block:: python

    # ✅ GOOD: Save and load configs
    from scLucid.qc import QCWorkflowConfig

    config = QCWorkflowConfig(min_genes=200, pc_mt=20.0)
    config.to_json_file("analysis_config.json")

    # Later, reproduce exact settings
    config = QCWorkflowConfig.from_json_file("analysis_config.json")

    # ❌ BAD: Hard-coded parameters scattered in code
    adata = run_standard_qc(adata, min_genes=200, pc_mt=20.0)

**Use Global Config for Runtime Settings**

.. code-block:: python

    from scLucid.config import set_config

    # Set at start of analysis
    set_config(
        n_jobs=8,           # Parallel processing
        verbosity=1,        # INFO level logging
        figure_dpi=300,     # Publication-quality figures
        log_file="analysis.log"
    )

Layer Convention
-----------------

scLucid follows a strict layer naming convention:

- ``counts`` - Raw UMI counts (never modify)
- ``normalized`` - Log1p(CPM) normalized data
- ``scaled`` - Z-score scaled data
- ``regressed`` - After covariate regression

**Always preserve raw counts**:

.. code-block:: python

    # ✅ GOOD: Store raw counts before processing
    adata.layers["counts"] = adata.X.copy()

    # ❌ BAD: Overwrite .X without backup
    adata.X = normalized_data  # Raw counts lost!

**Use appropriate layers for analysis**:

.. code-block:: python

    # Normalization: use raw counts
    adata = normalize_data(adata, input_layer="counts")

    # HVG selection: use normalized data
    adata = find_hvgs(adata, input_layer="normalized")

    # Clustering: use scaled data
    adata = cluster_cells(adata, use_rep="X_pca")

Quality Control
---------------

**Start with QC, don't skip it**

.. code-block:: python

    # Always run QC first
    from scLucid.qc import run_standard_qc

    adata_qc = run_standard_qc(adata)

    # Check filtering statistics
    print(f"Cells before: {adata.n_obs}")
    print(f"Cells after: {adata_qc.n_obs}")
    print(f"Filtered: {(1 - adata_qc.n_obs/adata.n_obs)*100:.1f}%")

**Adapt thresholds to your data**

.. code-block:: python

    # Use adaptive thresholds for automatic QC
    from scLucid.qc import AdaptiveThresholdLearner

    learner = AdaptiveThresholdLearner(method="gmm")
    threshold = learner.fit(adata.obs["pct_counts_mt"])
    print(f"Suggested mt threshold: {threshold}%")

**Visualize QC metrics**

.. code-block:: python

    # Always inspect QC plots
    sc.pl.violin(adata, ["n_genes_by_counts", "pct_counts_mt"], groupby=False)

    # Check doublet scores
    sc.pl.violin(adata, ["doublet_score"])

Preprocessing
-------------

**Choose appropriate normalization**

.. code-block:: python

    # Standard: CPM + log1p
    from scLucid.preprocess import NormalizationConfig

    config = NormalizationConfig(target_sum=1e4)

    # For sparse data: use SCTransform-style
    config = NormalizationConfig(
        target_sum=1e4,
        log_transform=True,
        normalize_per_cell=True
    )

**Validate HVG selection**

.. code-block:: python

    # Always plot HVGs
    sc.pl.highly_variable_genes(adata)

    # Check number of HVGs
    n_hvgs = adata.var["highly_variable"].sum()
    print(f"Found {n_hvgs} highly variable genes")

    # Expected: 2000-3000 for most datasets

**Batch correction only when needed**

.. code-block:: python

    # Diagnose batch effects first
    from scLucid.preprocess import diagnose_batch_effects

    diagnose_batch_effects(adata, batch_key="sample")

    # Only correct if strong batch effects present
    if needs_correction:
        from scLucid.preprocess import batch_correction
        adata = batch_correction(adata, method="harmony", batch_key="sample")

Clustering
----------

**Optimize resolution parameter**

.. code-block:: python

    # Don't use default resolution=1.0 blindly
    from scLucid.analysis import find_resolution

    optimal_res = find_resolution(
        adata,
        resolution_range=(0.2, 2.0, 10)
    )

    adata = cluster_cells(adata, resolution=optimal_res)

**Validate cluster quality**

.. code-block:: python

    # Check silhouette score
    from scikit_metrics import silhouette_score

    score = silhouette_score(adata.obsm["X_pca"], adata.obs["leiden"])

    # Good clusters: score > 0.25
    if score < 0.1:
        print("Warning: Poor clustering quality")

Annotation
-----------

**Combine multiple annotation methods**

.. code-block:: python

    from scLucid.analysis import AnnotationConfig, run_annotation

    # Use ensemble for robust annotation
    config = AnnotationConfig(
        method="ensemble",
        methods=["markers", "celltypist"],
        min_confidence=0.5,
        voting_strategy="majority"
    )

    adata = run_annotation(adata, config)

**Validate annotations manually**

.. code-block:: python

    # Always inspect marker expression
    marker_genes = {
        "CD3D": "T cells",
        "CD19": "B cells",
        "CD14": "Monocytes"
    }

    sc.pl.umap(adata, color=list(marker_genes.keys()))

Performance Optimization
-------------------------

**Enable parallel processing**

.. code-block:: python

    from scLucid.config import set_config

    # Use all available cores
    set_config(n_jobs=-1)

    # Or specify number
    set_config(n_jobs=8)

**Use caching for repeated analyses**

.. code-block:: python

    from scLucid.qc import CacheConfig, enable_cache

    # Enable caching
    enable_cache(cache_dir="~/.sclucid/qc_cache")

    # Subsequent runs will be faster
    adata = run_standard_qc(adata)

**Low-memory mode for large datasets**

.. code-block:: python

    from scLucid.config import set_config

    set_config(
        low_memory_mode=True,
        chunk_size=1000  # Process 1000 cells at a time
    )

Reproducibility
---------------

**Set random seed**

.. code-block:: python

    from scLucid.config import set_config

    set_config(random_state=42)  # For reproducible results

**Version your configs**

.. code-block:: python

    import json
    from pathlib import Path

    # Save config with metadata
    config_dict = {
        "config": config.to_dict(),
        "version": scLucid.__version__,
        "date": "2025-02-08",
        "data": "pbmc_batch1"
    }

    Path("config_v1.json").write_text(json.dumps(config_dict, indent=2))

**Store provenance in AnnData**

.. code-block:: python

    # Config is automatically stored
    adata.uns["sclucid"]["qc"]["config_used"]  # QC config
    adata.uns["sclucid"]["preprocess"]["workflow_config"]  # Preprocess config

Common Pitfalls
---------------

❌ **Don't**: Skip QC preprocessing
✅ **Do**: Always run QC first

❌ **Don't**: Use default parameters blindly
✅ **Do**: Adapt to your data characteristics

❌ **Don't**: Overwrite .X without backup
✅ **Do**: Preserve raw counts in .layers["counts"]

❌ **Don't**: Trust automated annotations blindly
✅ **Do**: Validate with marker gene inspection

❌ **Don't**: Ignore batch effects
✅ **Do**: Diagnose and correct if needed

Further Reading
---------------

- :doc:`quickstart` - Get started quickly
- :doc:`../notebooks/01_quality_control` - Interactive QC tutorial
- :doc:`../notebooks/02_preprocessing` - Interactive preprocessing tutorial
- :doc:`../notebooks/03_clustering_annotation` - Clustering and annotation tutorial
- :doc:`../notebooks/04_advanced_topics` - Advanced features and optimization
- :doc:`../api/qc` - QC API documentation
