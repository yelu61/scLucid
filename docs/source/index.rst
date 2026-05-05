scLucid Documentation
=====================

scLucid is a configurable single-cell RNA-seq analysis toolkit centered on
three ideas:

- **Traceable workflows** for QC, preprocessing, annotation, and differential analysis
- **Biology-aware decisions** such as marker-driven QC and multi-evidence annotation
- **Runnable entrypoints** for users who want either a one-liner workflow or a reviewable script

How To Use This Repository
--------------------------

The repository is organized into three user-facing layers:

- **docs/**: stable explanations, recommended workflows, API reference, and best practices
- **examples/**: short runnable scripts showing supported usage patterns
- **notebooks/**: longer end-to-end analyses on real or realistic datasets, including outputs worth reviewing

If you are new to scLucid, start with:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   usage_layers
   data_contracts
   best_practices
   workflow_hardening
   qc_preprocess_maturity

Workflow Guides
---------------

.. toctree::
   :maxdepth: 2

   examples
   notebooks
   proportion_methods_guide
   faq

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api/qc
   api/preprocess
   api/analysis
   api/recommendation
   api/tools
   api/plotting
   api/utils

What Goes Where
---------------

Use **docs** when you need:

- the official recommended path
- parameter semantics
- method-selection guidance
- stable references for package users

Use **examples** when you need:

- a short runnable script
- a template to adapt to your own project
- a minimal demonstration of one workflow or one reporting surface

Use **notebooks** when you need:

- a full analysis narrative
- real-data or project-style execution
- richer plots, intermediate decisions, and reviewer-facing outputs

Indices And Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
