Marker Resources
================

scLucid routes marker-dependent workflows through
``scLucid.utils.Manager`` and ``get_marker_manager()``. Marker curation lives in
resource files rather than analysis code, so annotation, QC, scoring, doublet
heuristics, and tumor interpretation can share the same evidence base.

Resource Layout
---------------

Core marker registries:

- ``marker_registry_human.toml``
- ``marker_registry_mouse.toml``

These registries use explicit sections:

- ``compartment``: broad compartments such as immune, stromal, epithelial, neural
- ``lineage``: broad biological lineages
- ``subtype``: reusable subtypes
- ``state``: lineage-aware cell states
- ``artifact``: QC, stress, ribosomal, mitochondrial, hemoglobin, ambient RNA signals
- ``functional_program``: module-scoring programs, not direct global labels

Tissue and tumor resources:

- ``marker_tissue_human.toml``: tissue-specific normal parenchymal and local subtype markers
- ``marker_tumor_human.toml``: epithelial support, malignancy programs,
  tumor type hints, cancer states, and diploid reference anchors

Gene-set resources such as ``genesets_cancer_signatures.json`` and
``genesets_cancer_hallmarks.json`` are intentionally kept separate. They may
overlap TOML resources at the gene level, but they are broad scoring/enrichment
resources rather than concise annotation marker registries.

Manager Views
-------------

Use views to keep biological evidence layers separate:

- ``compartment_annotation``: broad compartments only
- ``lineage_annotation``: compartments plus broad lineages
- ``subtype_annotation``: subtype and tissue-context annotation
- ``state_annotation``: cell states and cancer states
- ``artifact_annotation`` or ``qc_artifact``: QC/artifact signatures
- ``program_scoring``: functional programs and gene-set managers
- ``tumor_interpretation``: tumor-context evidence and malignancy interpretation

Example:

.. code-block:: python

   from scLucid.utils import get_marker_manager

   lineage_mgr = get_marker_manager("human", view="lineage_annotation")
   subtype_mgr = get_marker_manager("human", tissue="Lung", view="subtype_annotation")
   state_mgr = get_marker_manager("human", view="state_annotation")
   program_mgr = get_marker_manager("human", view="program_scoring")
   tumor_mgr = get_marker_manager(
       "human",
       cancer_type="Lung Cancer",
       view="tumor_interpretation",
   )

Curation Contract
-----------------

Each marker entry should include or inherit metadata such as ``kind``,
``granularity``, ``compartment``, ``lineage``, ``scope``, ``applies_to``,
``evidence_tier``, ``source_type``, ``review_status``, and routing flags such as
``use_for_global_annotation``, ``use_for_state_annotation``, and
``use_for_malignancy_interpretation``.

For the detailed curation rules and the LLM extraction prompt used for review
articles and pan-cancer atlas papers, see
``docs/MARKER_RESOURCE_CURATION.md`` in the repository root.
