Data Contracts
==============

scLucid workflows exchange data through an AnnData object plus a stable result
namespace under ``adata.uns["sclucid"]``. These contracts are intentionally
small: modules may store richer details, but the keys below are the shared
surface that downstream stages and tests can rely on.

Canonical AnnData Keys
----------------------

Layers:

- ``adata.layers["counts"]``: raw counts used as the preferred starting point
- ``adata.layers["normalized"]``: normalized expression from preprocessing
- ``adata.layers["scaled"]``: scaled expression when generated

Observation columns:

- ``sampleID``: canonical sample/batch-like sample identifier
- ``n_genes_by_counts``: QC gene-count metric
- ``total_counts``: QC library-size metric
- ``pct_counts_mt``: mitochondrial percentage when available
- ``low_quality``: QC low-quality flag
- ``leiden_clusters``: canonical clustering label
- ``cell_type_auto``: canonical automated annotation label

Embeddings:

- ``adata.obsm["X_pca"]``: PCA representation
- ``adata.obsm["X_umap"]``: UMAP representation
- ``adata.obsm["spatial"]``: spatial coordinates when present

Workflow Namespace
------------------

Every workflow stage stores its stable outputs under:

.. code-block:: python

   adata.uns["sclucid"][module]

The core module names are:

- ``qc``
- ``preprocess``
- ``analysis``
- ``tumor``
- ``tools``

Each stage should expose these standard keys when applicable:

- ``workflow_config``: serialized effective configuration
- ``steps_executed``: ordered workflow steps
- ``review_summary``: reviewer-facing execution and decision summary
- ``contract``: input/output validation records

Review Summary Envelope
-----------------------

Review summaries are normalized by ``scLucid.utils.normalize_review_summary``.
The stable envelope includes:

- ``schema_version``
- ``module``
- ``workflow_name``
- ``generated_at``
- ``steps_executed``
- ``data_shape``
- ``warnings``
- ``config``
- ``contract``

Module-specific summaries may add richer fields such as QC decision tables,
benchmark summaries, evidence bundles, or downstream recommendations. They
should preserve the envelope so automated checks can validate them.

Stage Boundaries
----------------

The canonical stage contracts are available at ``scLucid.utils.STAGE_CONTRACTS``
and can be inspected as dictionaries with:

.. code-block:: python

   from scLucid.utils import get_stage_contract, get_contract_spec

   print(get_stage_contract("preprocess"))
   print(get_contract_spec())

To validate a processed AnnData object:

.. code-block:: python

   from scLucid.utils import validate_stage_contract, validate_all_stage_contracts

   qc_result = validate_stage_contract(adata, "qc", when="output")
   all_results = validate_all_stage_contracts(adata, when="output")

These helpers are used by ``scLucid.run_pipeline`` to record validation results
under ``adata.uns["sclucid"][module]["contract"]``.
