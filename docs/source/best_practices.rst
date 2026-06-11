Best Practices
==============

This guide defines the recommended division of labor between the major scLucid
entrypoints and the repository artifacts around them.

Recommended Pipeline Shape
--------------------------

scLucid's maintained path is a light-dependency default with optional
enhancements. The default should be easy to install, reproducible on ordinary
Python single-cell environments, and conservative about removing biological
signal. Heavier methods remain available when the project has enough evidence
to justify them.

Light default:

- QC: metric calculation, adaptive/hierarchical threshold suggestions,
  conservative multi-criterion filtering, and Scrublet/heuristic doublet support
- Preprocessing: count preservation, library-size normalization + log1p,
  low-detection gene filtering, dependency-light HVG selection, scaling, PCA,
  neighbors, and UMAP
- Analysis: clustering review, marker/evidence tables, annotation consensus,
  differential expression, and proportion summaries

Optional enhancements:

- ``normalization.method="scran"`` through ``scanpy.external.pp.scran_normalize``
  when an R/scran environment is already available
- ``hvg.flavor="seurat_v3"`` for raw-count HVG selection when ``scikit-misc``
  is installed
- Harmony, Scanorama, scVI/scANVI, BBKNN, or ComBat when batch effects are
  documented and biological groups are not fully confounded with batch
- heavier doublet methods such as SOLO or DoubletDetection when their optional
  dependencies are installed and the dataset size justifies them

Review evidence, not automatic deletion:

- ambient RNA and empty-droplet diagnostics are heuristic review prompts; they
  do not replace CellBender/scAR/EmptyDrops-style project-level decisions
- cell-cycle regression diagnostics report associations and group imbalance;
  regression should be enabled only when the biology/technical tradeoff is
  explicit
- HVG biological protection can rescue marker/pathway/tumor genes, and any cap
  on rescued genes is reported with a deterministic truncation policy
- integration diagnostics, post-hoc QC cluster review, and doublet evidence
  should guide manual review before removing cells or collapsing biology

Removed from the recommended path:

- ScDblFinder wrapper execution paths
- project-level ambient RNA correction as a default preprocess stage
- custom rpy2/Bioconductor execution paths outside Scanpy's optional scran bridge

QC
--

Recommended default:

- use `run_standard_qc()` as the primary entrypoint
- keep `use_recommendations=True` unless you have a strong reason to lock thresholds manually
- prefer `threshold_mode="hierarchical"` for multi-sample datasets
- keep filtering conservative by requiring multiple independent low-quality
  criteria before removing cells
- review the stored QC trace and summary outputs before finalizing thresholds

When to override defaults:

- use pooled thresholds only when samples are intentionally treated as one shared population
- set explicit thresholds when project constraints require reproducible fixed cutoffs
- treat tumor-aware behavior as a cautionary layer, not a substitute for human review

Preprocessing
-------------

Recommended default:

- use `PreprocessingWorkflowConfig.default()` for the standard light-dependency path
- reserve `run_intelligent_preprocessing()` for datasets where parameter choice is uncertain
- keep the default path as the package's canonical preprocessing route in manuscripts and examples
- do not regress out `total_counts` or `pct_counts_mt` by default; enable
  regression only when diagnostics show a technical covariate dominates the
  biological signal
- do not run batch correction by default; enable integration only after
  inspecting batch mixing, sample structure, and the risk of over-correction

When to use intelligent preprocessing:

- when batch correction is uncertain
- when HVG / PCA / neighbors settings need reviewable justification
- when you want the exported preprocessing review summary before applying a workflow

Recommended stage handoff:

- preserve raw counts in ``adata.layers["counts"]``
- filter genes detected in too few cells at the start of preprocessing
  (default ``min_cells_per_gene=3``), after cell-level QC and before
  normalization/HVG selection
- store log-normalized expression in ``adata.layers["normalized"]`` and
  ``adata.raw`` before optional regression or HVG subsetting
- use HVGs for PCA/neighbors, but keep ``adata.raw`` available for marker,
  annotation, and differential-expression review
- use integrated embeddings for visualization/clustering only when justified;
  keep unintegrated normalized expression for marker and DE interpretation

Analysis defaults:

- treat cell-level marker tests as discovery/exploratory evidence
- use sample-level pseudobulk for condition DE when biological replicates exist
- use CLR sample-level tests for cell-type proportions; scLucid closes count,
  proportion, percentage, and sub-composition inputs before CLR transformation
- when using ``linear_model_logcpm`` for pseudobulk DE, the default covariance
  estimator is HC3 robust standard error; set ``robust_cov_type="nonrobust"``
  only for compatibility with ordinary OLS behavior

Stage-Wise Uns Compaction
-------------------------

``adata.uns["sclucid"]`` is useful during review, but it can grow after QC,
preprocessing, and analysis because it may contain benchmark summaries,
diagnostic tables, differential-expression tables, integration evidence, and
other intermediate objects. After exporting sidecar reports, compact the object
before handing it to the next stage or writing a portable ``.h5ad``:

.. code-block:: python

   # Inspect first.
   scl.compact_sclucid_uns(adata, dry_run=True)

   # Keep review summaries, configs, artifacts, context, warnings, and contract
   # metadata; remove heavier intermediate result payloads.
   scl.compact_sclucid_uns(adata)

Use ``keep_keys=[...]`` when a specific top-level module payload is still needed
inside ``uns``. For example, keep ``"de"`` during active marker exploration, but
drop it after exporting DE tables to CSV/Parquet sidecars.

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
