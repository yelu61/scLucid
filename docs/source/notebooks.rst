Notebooks
=========

Project notebooks are intentionally kept outside the documentation source tree.
In the current repository layout they live under
``examples/03_advanced_notebooks/``.

Why Notebooks Are Separate
--------------------------

Notebooks serve a different purpose from docs and examples:

- **docs** explain the supported package interface
- **examples** provide short runnable scripts
- **notebooks** show complete analysis stories, often with real data, richer plots,
  intermediate reasoning, and reviewer-facing outputs

For a package that is being prepared for publication, this separation helps keep
the stable documentation concise while still preserving full analysis narratives.

Recommended Notebook Scope
--------------------------

Use notebooks for:

- real-data walkthroughs
- manuscript-style analysis flows
- detailed parameter exploration
- outputs that are valuable to inspect but too long for docs or examples

Do not use notebooks as the only place where the recommended package path is defined.
That guidance should remain in :doc:`quickstart` and :doc:`best_practices`.

Current Notebook Directory
--------------------------

See ``examples/03_advanced_notebooks/`` for:

- ``Step1A-QC_Audit.ipynb`` - benchmark-grade QC audit notebook. It writes
  ``data/processed/Step1-sce_cleaned.h5ad``.
- ``Step1B-Preprocessing_Audit.ipynb`` - benchmark-grade preprocessing audit
  notebook with QC handoff context, layer transition evidence, and
  preprocessing ``step_evidence_summary``. It writes
  ``data/processed/Step2-sce_preprocessed.h5ad``.
- ``Step2-Annotation_and_Malignancy.ipynb`` - clustering, annotation,
  malignancy review, and CNV-aware interpretation. It writes
  ``data/processed/Step3-sce_annotated.h5ad``.
- ``Step3-Standard_Downstream.ipynb`` - composition, proportion,
  differential expression, and enrichment.
- ``Step4-Signature_and_Target_Analysis.ipynb`` - project-specific signatures,
  focused cell-state analysis, and target-oriented exports.
- ``Step1-QC_and_Preprocessing.ipynb`` and
  ``Step2-Celltype_annotation.ipynb`` - legacy unsplit project references.
- ``04_advanced_topics.ipynb``
- ``04_differential_expression.ipynb``
- ``05_trajectory_inference.ipynb``
- ``06_advanced_tools.ipynb``

If these notebooks evolve into a publication-facing tutorial set, keep them in
``examples/03_advanced_notebooks/`` and maintain a short README that explains
their scope and expected inputs.
