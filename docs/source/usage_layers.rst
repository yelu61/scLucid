Usage Layers
============

scLucid is designed as a three-layer product rather than a single API style.
The layers target different users and different moments in the same research
project.

Layer 1: Workflow
-----------------

The workflow layer is the first screen for most users. It runs the supported
QC -> preprocessing -> analysis path with conservative defaults and records the
decisions under ``adata.uns["sclucid"]``.

Use this layer when:

- a user wants a complete first pass with minimal code
- the project needs a reproducible baseline before manual tuning
- the goal is to compare datasets with the same supported path

Primary entrypoints:

- ``scLucid.run_pipeline``
- ``scLucid.qc.run_standard_qc``
- ``scLucid.preprocess.run_preprocessing``
- ``scLucid.analysis.run_standard_analysis``

Expected output:

- filtered and processed ``AnnData``
- stage review summaries
- contract validation records
- optional sidecar reports and figures when ``save_dir`` is set

Layer 2: Simple API
-------------------

The simple API layer exposes each workflow stage as composable steps. It is for
analysts who want to inspect a decision, replace one method, or run a stage
again with a different config while keeping the rest of the package conventions.

Use this layer when:

- QC thresholds need explicit review before filtering
- preprocessing parameters need to be compared
- annotation evidence needs manual inspection
- a user wants more control without rewriting the full workflow

Primary entrypoints:

- ``scLucid.qc.calculate_qc_metric``
- ``scLucid.qc.recommend_intelligent_qc``
- ``scLucid.qc.mark_low_quality_cell``
- ``scLucid.qc.filter_cells``
- ``scLucid.preprocess.normalize_data``
- ``scLucid.preprocess.find_hvgs``
- ``scLucid.preprocess.scale_data``
- ``scLucid.preprocess.batch_correction``

Expected output:

- the same AnnData conventions as the workflow layer
- inspectable intermediate objects and tables
- reviewer-facing reports for the decisions the user chose manually

Layer 3: Advanced
-----------------

The advanced layer is for real exploratory analysis where every decision should
be visible. It is usually a notebook or project script that uses configs,
review summaries, sidecar artifacts, and custom checkpoints together.

Use this layer when:

- a manuscript workflow needs a complete audit trail
- tumor-specific assumptions need manual review
- multiple parameter choices must be compared before finalizing
- the user needs custom hooks, checkpoints, or project-specific metadata

Primary artifacts:

- ``examples/03_advanced_notebooks/``
- project notebooks based on the same contracts
- golden-path scripts such as ``scripts/run_pbmc_golden_path.py``

Recommended product-facing notebook sequence:

1. ``Step1A-QC_Audit.ipynb`` -> ``Step1-sce_cleaned.h5ad``
2. ``Step1B-Preprocessing_Audit.ipynb`` -> ``Step2-sce_preprocessed.h5ad``
3. ``Step2-Annotation_and_Malignancy.ipynb`` -> ``Step3-sce_annotated.h5ad``
4. ``Step3-Standard_Downstream.ipynb``
5. ``Step4-Signature_and_Target_Analysis.ipynb``

Expected output:

- final ``.h5ad``
- review summaries for each stage
- figures for inspection
- machine-readable manifest
- explicit notes for user overrides and biological assumptions

How The Layers Work Together
----------------------------

The three layers should not become three separate products. They should share
the same data contracts, config names, review summary envelope, and output
locations.

The frozen layer contract is available in code:

.. code-block:: python

   from scLucid.utils import get_api_layer_spec

   print(get_api_layer_spec("workflow"))

A common project flow is:

1. Run the workflow layer to get a baseline.
2. Inspect QC and preprocessing review summaries.
3. Drop into the simple API layer for the stage that needs adjustment.
4. Promote the final decisions into an advanced notebook or golden-path script.

Documentation Responsibilities
------------------------------

``docs/`` should explain policy:

- which layer to start with
- which defaults are recommended
- how review summaries should be interpreted
- which features are stable versus experimental

``examples/`` should show runnable usage:

- one short script per supported layer or scenario
- minimal assumptions near the top of the file
- no hidden package policy that is absent from docs

``examples/03_advanced_notebooks/`` should show full analysis narratives:

- richer intermediate plots
- parameter review sections
- real-data or project-style execution
- final outputs that can be inspected by a reviewer
- the same module maturity, compact summary, QC handoff, layer transition, and
  step-evidence contracts used by the workflow/simple API layers when a
  notebook implements QC or preprocessing manually

Product Acceptance Bar
----------------------

A layer is considered product-ready only when:

- the code path is covered by at least one runnable example
- the docs describe when to use it and when not to use it
- the output follows the scLucid AnnData contract
- failures produce actionable errors or review warnings
- a lightweight or golden-path test protects the expected behavior
