.. scLucid documentation master file, created by sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scLucid!
===================

**Evidence-driven single-cell RNA-seq analysis toolkit**

scLucid is a comprehensive Python toolkit for single-cell RNA-seq analysis that emphasizes **traceable, configurable, and biologically interpretable** analysis.

.. image:: https://img.shields.io/badge/python-3.9+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python Version

.. image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

Key Features
------------

**🔬 Complete Analysis Pipeline**

- Quality control with adaptive thresholds
- Preprocessing with checkpoint/resume support
- Clustering, annotation, and differential expression
- Batch correction and integration

**🎯 Evidence-Driven Analysis**

- Configurable workflows with Pydantic validation
- Automatic provenance tracking
- Reproduducible analyses with version-controlled configs

**🧬 Biology-Aware**

- Hierarchical marker gene system
- Cell type annotation with multiple methods
- Lineage-aware doublet detection

**⚡ Performance**

- Parallel processing support
- Incremental QC for large datasets
- Optional caching for faster re-runs

Quick Example
-------------

.. code-block:: python

    import scLucid

    # Complete pipeline in 3 steps
    from scLucid.qc import run_standard_qc
    from scLucid.preprocess import run_preprocessing
    from scLucid.analysis import cluster_cells, annotate_clusters

    # 1. Quality control
    adata = run_standard_qc(adata)

    # 2. Preprocessing
    adata = run_preprocessing(adata)

    # 3. Analysis
    adata = cluster_cells(adata)
    adata = annotate_clusters(adata, marker_manager)

Getting Started
---------------

.. toctree::
   :maxdepth: 2

   installation
   quickstart

Interactive Tutorials
---------------------

.. toctree::
   :maxdepth: 2

   notebooks/01_quality_control
   notebooks/02_preprocessing
   notebooks/03_clustering_annotation
   notebooks/04_advanced_topics

Best Practices
--------------

.. toctree::
   :maxdepth: 2

   best_practices
   faq

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api/qc
   api/preprocess
   api/analysis
   api/tools
   api/plotting
   api/utils

Code Examples
------------

Runnable scripts for quick reference:

- **quickstart.py** - Complete pipeline in 5 minutes
- **qc_pipeline.py** - Quality control workflow
- **preprocessing.py** - Preprocessing pipeline

See the :doc:`examples` directory.

Contributing
------------

Contributions are welcome! Please see our GitHub repository for guidelines:

- **Bug Reports**: https://github.com/yourusername/scLucid/issues
- **Pull Requests**: https://github.com/yourusername/scLucid/pulls
- **Documentation**: https://github.com/yourusername/scLucid/tree/main/docs

License
-------

scLucid is released under the MIT License. See LICENSE file for details.

Citation
--------

If you use scLucid in your research, please cite::

    TODO: Add citation information

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
