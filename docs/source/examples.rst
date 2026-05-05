Examples
========

The ``examples/`` directory contains short runnable scripts.
These scripts are templates, not the authoritative source of package policy.
The recommended defaults are defined in :doc:`quickstart` and :doc:`best_practices`.

Current Example Set
-------------------

``examples/01_workflow/basic_pipeline.py``
    Minimal end-to-end QC -> preprocessing -> clustering -> annotation flow.

``examples/01_workflow/prepare_data.py``
    Loading 10x-style inputs, attaching metadata, and ensuring
    ``adata.layers["counts"]`` is present.

``examples/01_workflow/plugin_development.py``
    Minimal plugin extension example for users who need custom steps.

``examples/02_simple_api/qc_step_by_step.py``
    Composable QC path: metrics, recommendations, doublets, marking, filtering,
    and report export.

``examples/02_simple_api/preprocess_step_by_step.py``
    Composable preprocessing path:

    - standard default path
    - intelligent preprocessing with review summary
    - manual step-by-step path

``examples/02_simple_api/intelligent_qc.py``
    Focused demonstration of data-driven QC recommendation concepts.

``examples/02_simple_api/intelligent_preprocess.py``
    Focused demonstration of recommendation generation and review summary usage.

``examples/02_simple_api/annotation_report.py``
    Export reviewer-facing annotation reports and sidecars.

``examples/02_simple_api/annotation_workflow.py``
    Show how to replace long notebook cells with reusable API calls for marker
    review, manual mapping, scoring, and composition plots.

``examples/03_advanced_notebooks/``
    Complete notebook workflows for full analysis narratives and reviewable
    intermediate outputs. The recommended product-facing sequence is:

    - ``Step1A-QC_Audit.ipynb`` - QC benchmark path and
      ``Step1-sce_cleaned.h5ad``
    - ``Step1B-Preprocessing_Audit.ipynb`` - QC handoff, layer audit,
      preprocessing step evidence, and ``Step2-sce_preprocessed.h5ad``
    - ``Step2-Annotation_and_Malignancy.ipynb`` - clustering, annotation,
      malignancy review, CNV-aware interpretation, and
      ``Step3-sce_annotated.h5ad``
    - ``Step3-Standard_Downstream.ipynb`` - composition, proportion,
      differential expression, and enrichment
    - ``Step4-Signature_and_Target_Analysis.ipynb`` - project-specific
      signatures, focused cell-state analysis, and target-oriented exports

``scripts/run_pbmc_golden_path.py``
    Maintained real-data acceptance script. This is longer than an example
    because it writes a manifest, figures, review artifacts, and final ``.h5ad``.

Layer Mapping
-------------

Use these files to demonstrate the three-layer product design:

- Workflow layer: ``examples/01_workflow/basic_pipeline.py``
- Simple API layer: ``examples/02_simple_api/qc_step_by_step.py`` and
  ``examples/02_simple_api/preprocess_step_by_step.py``. For the QC/preprocess
  benchmark path, use ``examples/02_simple_api/qc_preprocess_review.py``.
- Advanced layer: ``examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb``,
  ``examples/03_advanced_notebooks/Step1B-Preprocessing_Audit.ipynb``,
  ``examples/03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb``,
  ``examples/03_advanced_notebooks/Step3-Standard_Downstream.ipynb``,
  ``examples/03_advanced_notebooks/Step4-Signature_and_Target_Analysis.ipynb``,
  and ``scripts/run_pbmc_golden_path.py``

How To Maintain Examples
------------------------

- each example should remain runnable as a standalone script
- each example should demonstrate one supported usage pattern clearly
- examples should link back to the recommended path, not redefine it
- if an example demonstrates an advanced or experimental route, say so explicitly
