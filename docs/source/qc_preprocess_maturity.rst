QC And Preprocessing Maturity Plan
==================================

QC and preprocessing should become the first benchmark modules for scLucid.
They define the user's trust in the rest of the workflow: if filtering,
normalization, feature selection, and graph construction are not explainable,
the later annotation and tumor modules cannot be trusted.

Target State
------------

The target is not simply "more functions". A mature QC or preprocessing module
should satisfy five standards:

- stable public entrypoints for workflow, simple API, and advanced usage
- explicit AnnData input and output contracts
- review summaries that explain recommendations, applied parameters, and user overrides
- real-data golden-path validation on PBMC and tumor datasets
- clear docs and examples that match the maintained code path

The maintained code path is a light-dependency default with optional
enhancements. Defaults should not require R, scVI-tools, Scanorama,
DoubletDetection, or other heavy extras. Those methods remain supported as
explicit opt-in extensions when the dataset and environment justify them.

QC As The First Benchmark Module
--------------------------------

QC should answer four user questions:

1. What thresholds were recommended?
2. Why were they recommended?
3. What was actually applied?
4. What biological or technical risk remains after filtering?

Required QC outputs:

- ``adata.obs`` metrics such as ``n_genes_by_counts``, ``total_counts``, and
  ``pct_counts_mt`` when available
- low-quality and doublet flags
- ``adata.uns["sclucid"]["qc"]["workflow_config"]``
- ``adata.uns["sclucid"]["qc"]["review_summary"]``
- optional report sidecars under the configured ``save_dir``

QC entries that can usually be compacted after review sidecars are exported:

- detailed benchmark payloads duplicated inside ``review_summary``
- doublet diagnostic parameter/evidence payloads when exported tables exist
- raw recommendation objects when applied/original configs and review summaries
  are retained

QC hardening tasks:

- keep ``run_standard_qc`` as the canonical workflow entrypoint
- keep ``recommend_intelligent_qc`` executable as a standalone simple API tool
- keep ScDblFinder and ambient RNA correction outside the default QC/preprocess
  path; removed wrappers should not reappear unless there is a clear
  dependency and maintenance plan
- make user overrides explicit in the review summary
- test tumor-aware behavior on PDAC data where high mitochondrial content may
  be a warning rather than an automatic removal criterion
- test edge cases: missing mitochondrial genes, single-sample data, small cell
  counts, sparse matrices, and absent ``sampleID``

Preprocessing As The Second Benchmark Module
--------------------------------------------

Preprocessing should answer four user questions:

1. Which expression layer was used at each step?
2. Why were HVG, PCA, neighbors, and batch-correction parameters chosen?
3. Was biological signal protected from over-correction?
4. Is the output ready for clustering and annotation?

Required preprocessing outputs:

- ``adata.layers["normalized"]`` when normalization runs
- ``adata.var["highly_variable"]`` or the configured HVG key
- ``adata.obsm["X_pca"]``
- neighbors graph and ``adata.obsm["X_umap"]`` when graph steps run
- ``adata.uns["sclucid"]["preprocess"]["workflow_config"]``
- ``adata.uns["sclucid"]["preprocess"]["review_summary"]``

Preprocessing entries that can usually be compacted after the handoff to
analysis:

- normalization/HVG/scaling/integration diagnostic dictionaries once their
  summaries are in ``review_summary``
- large integration evaluation payloads and temporary intelligent-preprocess
  recommendation details after sidecar export
- redundant layer statistics when counts, normalized layer, PCA, neighbors, and
  review summary are retained

Analysis entries that can usually be compacted before long-term storage:

- full DE, enrichment, proportion, and annotation-review tables after exporting
  them to sidecar files
- intermediate clustering-resolution grids after the chosen cluster labels and
  review summary are retained

Preprocessing hardening tasks:

- keep ``run_preprocessing`` as the canonical workflow entrypoint
- make layer transitions explicit in the review summary
- make HVG selection evidence inspectable
- keep regression and batch correction opt-in by default, and document when to
  enable them
- document when to skip regression or HVG subsetting
- keep ``scanpy.external.pp.scran_normalize`` as the only supported optional
  scran path; avoid custom rpy2/Bioconductor execution branches
- warn when tumor data are batch-corrected in a way that may remove malignant,
  clone, patient, or microenvironment signal
- test small datasets where PCA components and neighbors must be clipped safely

Recommended Implementation Order
--------------------------------

1. Freeze the minimal shared contract.
   Keep ``adata.uns["sclucid"]``, review summary envelopes, and canonical keys
   stable before refactoring internals.

2. Stabilize QC.
   Run PBMC and PDAC slices, then fix every unclear threshold, warning, or
   report field discovered by those runs.

3. Stabilize preprocessing.
   Use the QC output as input, then harden normalization, HVG, PCA, neighbors,
   and optional integration until the handoff to analysis is predictable.

4. Only then polish analysis.
   Analysis quality depends on the first two modules. Annotation and tumor
   interpretation should be improved after QC and preprocessing are reliable.

Definition Of Done For A Module
-------------------------------

A module can be treated as a scLucid benchmark module when all of the following
are true:

- the workflow entrypoint passes lightweight tests
- at least one real-data golden path exercises the module
- public examples cover workflow and simple API use
- docs explain default behavior, override behavior, and review artifacts
- review summaries contain enough evidence for another analyst to audit the run
- output can be serialized to ``.h5ad`` after compacting heavy artifacts
- known limitations are documented rather than hidden

Suggested Test Tiers
--------------------

- smoke tests: imports, public API availability, config construction
- unit tests: threshold logic, layer selection, HVG behavior, edge cases
- lightweight integration: synthetic or small AnnData through the full module
- golden path: PBMC and PDAC subsets with saved artifacts
- project acceptance: active real-world datasets with biological plausibility review
