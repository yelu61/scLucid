Best Practices
==============

This guide defines the recommended division of labor between the major scLucid
entrypoints and the repository artifacts around them.

QC
--

Recommended default:

- use `run_standard_qc()` as the primary entrypoint
- keep `use_recommendations=True` unless you have a strong reason to lock thresholds manually
- prefer `threshold_mode="hierarchical"` for multi-sample datasets
- review the stored QC trace and summary outputs before finalizing thresholds

When to override defaults:

- use pooled thresholds only when samples are intentionally treated as one shared population
- set explicit thresholds when project constraints require reproducible fixed cutoffs
- treat tumor-aware behavior as a cautionary layer, not a substitute for human review

Preprocessing
-------------

Recommended default:

- use `PreprocessingWorkflowConfig.default()` for the standard path
- reserve `run_intelligent_preprocessing()` for datasets where parameter choice is uncertain
- keep the default path as the package's canonical preprocessing route in manuscripts and examples

When to use intelligent preprocessing:

- when batch correction is uncertain
- when HVG / PCA / neighbors settings need reviewable justification
- when you want the exported preprocessing review summary before applying a workflow

Examples Vs Docs Vs Notebooks
-----------------------------

The repository should be read in this order:

- **docs**: stable explanation of the recommended path
- **examples**: short runnable scripts for each supported usage pattern
- **notebooks**: complete analyses with intermediate decisions and richer plots

Do not use examples as the authoritative source for package policy.
The authoritative recommended path should live in docs, and examples should
implement that policy rather than redefine it.

Suggested Usage Patterns
------------------------

Use `examples/` for:

- quick copy-and-adapt scripts
- testing one workflow in isolation
- minimal reporting examples

Use `notebooks/` for:

- full project walkthroughs
- real-data or publication-style analyses
- exploration with rich intermediate figures

Use `docs/` for:

- installation and onboarding
- recommended defaults
- method selection guidance
- public API reference

Reproducibility
---------------

- store raw counts in ``adata.layers["counts"]`` before preprocessing
- prefer config objects over scattered keyword arguments
- keep workflow outputs under stable ``save_dir`` locations
- review sidecar outputs before treating automated decisions as final

Workflow Hardening
------------------

When improving scLucid, prefer vertical workflow slices over broad module-by-module
expansion. Start with a small reproducible dataset, run the full supported path,
inspect the review summaries, and polish the module boundary that fails first.

Recommended validation tiers:

- ``data/pbmc3k.h5ad`` for the fast normal-tissue baseline
- ``data/lin2020.pdac.h5ad`` for the first tumor golden path
- ``data/schlesinger2020.pdac.h5ad`` for tumor generalization
- active project data for final product acceptance and biological plausibility

See :doc:`workflow_hardening` for the detailed execution plan.

QC And Preprocessing First
--------------------------

For productization, QC and preprocessing should be treated as the first two
benchmark modules. They are the foundation for annotation, differential
expression, tumor analysis, and visualization.

Recommended order:

- stabilize the shared workflow contract
- harden QC on PBMC and tumor data
- harden preprocessing on the QC output
- then move to analysis and annotation polish

See :doc:`qc_preprocess_maturity` for the module-level maturity checklist.
