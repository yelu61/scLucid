Data Contracts
==============

scLucid workflows exchange data through an AnnData object plus a stable result
namespace under ``adata.uns["sclucid"]``. These contracts are intentionally
small: modules may store richer details, but the keys below are the shared
surface that downstream stages and tests can rely on.

Contract Version
----------------

The current frozen contract schema is ``1.0``. The machine-readable contract is
available with:

.. code-block:: python

   from scLucid.utils import get_contract_spec

   spec = get_contract_spec()

The contract freezes three surfaces:

- the three user-facing API layers
- the canonical AnnData keys
- the minimal workflow boundary for QC -> preprocessing -> analysis

API Layer Contract
------------------

The public API is organized into three frozen layers. They can be inspected with:

.. code-block:: python

   from scLucid.utils import get_api_layer_spec

   print(get_api_layer_spec())
   print(get_api_layer_spec("workflow"))

The layer order is:

- ``workflow``: one-call or stage-level workflow entrypoints for the supported baseline path
- ``simple_api``: composable QC and preprocessing functions for inspection and overrides
- ``advanced``: notebooks and golden-path scripts for complete audit trails

The workflow and simple API layers are protected by tests that resolve the
documented entrypoints. The advanced layer is artifact-based: notebooks and
golden-path scripts are treated as maintained user-facing assets.

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

The root namespace and each module namespace include ``_metadata`` with the
contract schema version and timestamps. This is deliberately separate from
module outputs so readers can distinguish data products from framework metadata.

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
- ``config_lineage``: global, inherited, explicit, and effective config record
- ``contract``: input/output validation records
- ``artifacts``: paths to saved JSON, Markdown, figures, manifests, or other sidecars
- ``errors``: structured workflow errors when a stage records recoverable failure state

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
- ``config_lineage``
- ``artifacts``

Module-specific summaries may add richer fields such as QC decision tables,
benchmark summaries, evidence bundles, or downstream recommendations. They
should preserve the envelope so automated checks can validate them.

For backward compatibility during the schema transition, normalized review
summaries also expose a shallow ``data`` view containing the same fields as the
flat envelope. New code should read the flat keys directly, while older
notebooks can continue to read ``review_summary["data"]`` until they are
migrated.

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
under ``adata.uns["sclucid"][module]["contract"]``. The pipeline records both
input and output validations for each stage it executes.

Minimal Workflow Contract
-------------------------

The frozen minimal workflow contract is available with:

.. code-block:: python

   from scLucid.utils import get_minimal_workflow_contract

   print(get_minimal_workflow_contract())

The current minimal workflow is:

.. code-block:: text

   qc -> preprocess -> analysis

Each executed stage should expose at least:

- ``workflow_config``
- ``steps_executed``
- ``review_summary``

The unified pipeline additionally records:

- ``adata.uns["sclucid"]["pipeline_context"]``
- ``adata.uns["sclucid"]["analysis_context"]``
- input and output contract validation under each executed stage namespace

This is the stable framework boundary. Individual modules may add richer
details, but downstream code should depend only on this minimal surface unless
it explicitly opts into a module-specific contract.

Config Lineage
--------------

Each pipeline stage records how its effective configuration was derived:

.. code-block:: python

   adata.uns["sclucid"]["qc"]["config_lineage"]

The lineage contains:

- ``global``: package-level defaults
- ``inherited``: pipeline context inferred from data and user hints
- ``stage``: explicit stage config supplied by the caller
- ``effective``: serialized config actually stored by the stage
- ``precedence``: the intended override order

Saved Artifacts
---------------

Review-facing files should be registered under ``artifacts`` whenever a stage
writes sidecars to disk. For example:

.. code-block:: python

   adata.uns["sclucid"]["preprocess"]["artifacts"]["preprocess_review_summary_json"]

This keeps the final AnnData object as the index of a run while allowing large
or format-specific outputs to live beside it on disk.
