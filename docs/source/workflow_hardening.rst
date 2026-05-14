Workflow Hardening Plan
=======================

After the core contracts are stable, scLucid should be improved through vertical
workflow slices rather than by expanding every module at once. A vertical slice
starts from raw or near-raw AnnData, runs the supported workflow, inspects the
stored evidence, and fixes the weakest module boundary discovered by that run.

Why Vertical Slices
-------------------

scLucid's value is not simply that it wraps many tools. The package is meant to
make real single-cell decisions reviewable: QC thresholds, preprocessing
parameters, annotation evidence, tumor-aware warnings, and publication-ready
outputs should line up in one traceable workflow.

Vertical slices test that promise better than isolated unit tests because they
exercise:

- data loading and count-layer assumptions
- QC recommendations and filtering effects
- preprocessing parameter choices
- clustering and annotation handoff
- tumor-specific analysis when relevant
- review-summary and contract outputs
- figure/report generation

Recommended Dataset Roles
-------------------------

Use repository datasets for regression and benchmark stability:

- ``data/pbmc3k.h5ad``: normal immune baseline; fast public regression target
- ``data/lin2020.pdac.h5ad``: tumor tissue workflow target
- ``data/schlesinger2020.pdac.h5ad``: second tumor dataset for generalization

Use real project datasets for product discovery and acceptance:

- identify workflow friction that synthetic or public datasets miss
- test sample metadata conventions, tumor context, and reviewer-facing outputs
- decide whether automated recommendations match biological expectations
- define the final acceptance bar for a module

This split keeps the package reproducible while still letting real research
needs drive the roadmap.

Execution Order
---------------

1. PBMC golden path
   Run QC, preprocessing, clustering, annotation, and figure export on
   ``data/pbmc3k.h5ad``. The goal is a fast, stable non-tumor baseline. The
   maintained entrypoint is ``scripts/run_pbmc_golden_path.py``.

2. PDAC tumor golden path
   Run QC, preprocessing, standard analysis, and the first tumor-specific slice
   on ``data/lin2020.pdac.h5ad``. The initial tumor slice should focus on
   malignancy/CNV/TME rather than every advanced tumor module.

3. Cross-dataset tumor validation
   Repeat the same workflow on ``data/schlesinger2020.pdac.h5ad`` and compare
   retained cells, QC decisions, annotation coverage, tumor warnings, and
   runtime.

4. Real-project acceptance pass
   Apply the current golden workflow to the user's active project data. Gaps
   discovered here should become targeted issues rather than
   broad rewrites.

5. Module polishing
   Polish the module that blocks the vertical workflow most often. In practice
   the likely order is QC, preprocessing, annotation, tumor, then plotting.

Acceptance Artifacts
--------------------

Each golden workflow should produce:

- final ``.h5ad`` output
- QC review summary and benchmark report
- preprocessing review summary
- analysis or annotation review outputs
- tumor execution trace when tumor stages are used
- a small figure set suitable for visual inspection
- a machine-readable run manifest with versions, runtime, dataset shape, and
  contract validation results

Suggested Local Command Pattern
-------------------------------

Use the maintained single-cell environment for local gates:

.. code-block:: bash

   MAMBA_EXE=/opt/homebrew/bin/mamba \
   SCLUCID_TEST_ENV_PATH=/Users/luye/micromamba/envs/scrna-env \
   scripts/run_test_gates.sh

Longer real-data workflow checks should be separate from pre-merge gates so the
fast test loop remains usable.

PBMC Golden Path
----------------

The PBMC golden path is the first vertical-slice regression target. It should
remain small enough for local iteration while still exercising the user-facing
workflow across QC, preprocessing, analysis, annotation, plotting, contracts,
and final artifact writing.

Run a quick subset check:

.. code-block:: bash

   /Users/luye/micromamba/envs/scrna-env/bin/python \
     scripts/run_pbmc_golden_path.py \
     --n-cells 300 \
     --n-top-genes 500 \
     --n-pcs 20 \
     --n-neighbors 10 \
     --output-dir results/golden/pbmc3k_subset \
     --overwrite

Run the integration acceptance test:

.. code-block:: bash

   /Users/luye/micromamba/envs/scrna-env/bin/python \
     -m pytest -q tests/integration/test_pbmc_golden_path.py \
     -m "slow and integration"

Expected outputs include:

- ``manifest.json`` with versions, shapes, runtime, contract validation, and
  artifact paths
- ``validation/qc_preprocess_validation.json`` and
  ``validation/qc_preprocess_validation_table.csv`` with lightweight
  QC/preprocess maturity metrics. These artifacts mark runs as ready for later
  comparative validation; they do not claim superiority over external
  workflows.
- ``pbmc3k_golden_final.h5ad`` with compact scLucid review summaries retained
  under ``adata.uns["sclucid"]``
- QC, preprocessing, and analysis output directories
- UMAP figures for Leiden clusters, automatic cell-type labels, and sample IDs

The final ``.h5ad`` intentionally stores compact review summaries instead of
embedding every intermediate evidence object. Full evidence should remain in
the run manifest and stage output directories so the saved AnnData stays
portable and does not fail HDF5 serialization.
