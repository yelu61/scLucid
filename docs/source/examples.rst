Examples
========

The ``examples/`` directory contains short runnable scripts.
These scripts are templates, not the authoritative source of package policy.
The recommended defaults are defined in :doc:`quickstart` and :doc:`best_practices`.

Current Example Set
-------------------

``examples/quickstart.py``
    Minimal end-to-end QC → preprocessing → clustering → annotation flow.

``examples/preprocessing.py``
    Three supported preprocessing modes:

    - standard default path
    - intelligent preprocessing with review summary
    - manual step-by-step path

``examples/intelligent_qc_example.py``
    Focused demonstration of data-driven QC recommendation concepts.

``examples/intelligent_preprocessing_example.py``
    Focused demonstration of recommendation generation and review summary usage.

``examples/annotation_report.py``
    Export reviewer-facing annotation reports and sidecars.

``examples/curated_annotation_workflow.py``
    Show how to replace long notebook cells with reusable API calls for marker review,
    manual mapping, scoring, and composition plots.

How To Maintain Examples
------------------------

- each example should remain runnable as a standalone script
- each example should demonstrate one supported usage pattern clearly
- examples should link back to the recommended path, not redefine it
- if an example demonstrates an advanced or experimental route, say so explicitly
