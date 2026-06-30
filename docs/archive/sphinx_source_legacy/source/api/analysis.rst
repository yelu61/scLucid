Analysis Module API Reference
==============================

The analysis module is the active hardening target after QC and preprocessing.
Its maintained workflow is evidence-first: clustering is reviewed before marker
discovery, annotation is built from multiple evidence sources, final labels are
stored with review status, and post-hoc QC risks are surfaced without automatic
deletion.

.. automodule:: scLucid.analysis
   :members:
   :undoc-members:
   :show-inheritance:

Analysis Review Contract
------------------------

``run_standard_analysis`` writes a structured review summary under
``adata.uns["sclucid"]["analysis"]["review_summary"]``. The current contract
contains:

* ``preprocess_input_context``: PCA, neighbors, UMAP, normalized layer, and HVG
  handoff status from preprocessing.
* ``analysis_inference_policy``: conservative claim boundaries for clustering
  review, marker discovery, condition DE, and exploratory cell-level
  comparisons. Condition DE is directed to sample-level pseudobulk; cell-level
  comparisons are explicitly marked exploratory.
* ``analysis_output_contract``: stable output slots and inference semantics for
  preprocess handoff, clustering, marker discovery, annotation, condition DE,
  and post-hoc QC review.
* ``analysis_decision_summary`` and ``analysis_reviewer_table``: the single
  reviewer-facing decision layer with recommended value, applied value, source,
  confidence, affected output, inference level, biological risk note, and
  manual-review status.
* ``clustering_evidence_summary``: cluster counts and optional resolution-review
  evidence.
* ``annotation_evidence_summary``: marker/reference/LLM evidence availability
  and annotation-review table health.
* ``annotation_consensus_summary``: final label column, confidence, and
  needs-review cells.
* ``posthoc_qc_review_summary``: cluster-level doublet-heavy,
  high-mitochondrial, or stress-high signals detected after clustering. These
  are review prompts, not automatic filtering rules.
* ``malignancy_interpretation_summary``: optional tumor-context call counts,
  malignant/suspect fraction, evidence sources, and low-purity warning.
* ``analysis_readiness`` and ``review_action_items``: downstream readiness and
  human-review tasks.
* ``evidence_bundle`` and ``module_maturity``: shared scLucid evidence and
  contract-completeness views.

The reviewer table is the preferred notebook/API/report surface. Use it before
making biological claims from clusters, labels, condition DE, proportions, or
tumor-context interpretation.

Inspect the compact view:

.. code-block:: python

   compact = scLucid.analysis.summarize_analysis_review_summary(
       adata.uns["sclucid"]["analysis"]["review_summary"]
   )

Validate a result:

.. code-block:: python

   scLucid.analysis.validate_analysis_module_completeness(adata)

Configuration Classes
---------------------

ClusteringConfig
~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.analysis.ClusteringConfig
   :members:
   :undoc-members:

AnnotationConfig
~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.analysis.AnnotationConfig
   :members:
   :undoc-members:

DifferentialConfig
~~~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.analysis.DifferentialConfig
   :members:
   :undoc-members:

EnrichmentConfig
~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.analysis.EnrichmentConfig
   :members:
   :undoc-members:

ProportionConfig
~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.analysis.ProportionConfig
   :members:
   :undoc-members:

Core Functions
--------------

run_standard_analysis
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.run_standard_analysis

cluster_cells
~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.cluster_cells

annotate_clusters
~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.annotate_clusters

run_annotation
~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.run_annotation

run_annotation_evidence
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.run_annotation_evidence

build_annotation_consensus
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.build_annotation_consensus

build_posthoc_qc_review_summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.build_posthoc_qc_review_summary

find_markers
~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.find_markers

run_pseudobulk_de
~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.run_pseudobulk_de

run_enrichment
~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.run_enrichment

Proportion Helpers
------------------

.. autofunction:: scLucid.analysis.compute_celltype_proportion

.. autofunction:: scLucid.analysis.run_statistical_test

.. autofunction:: scLucid.analysis.plot_composition

.. autofunction:: scLucid.analysis.plot_proportion_with_ci

Review Helpers
--------------

.. autofunction:: scLucid.analysis.get_analysis_module_contract

.. autofunction:: scLucid.analysis.build_analysis_output_contract

.. autofunction:: scLucid.analysis.build_analysis_decision_summary

.. autofunction:: scLucid.analysis.build_analysis_reviewer_table

.. autofunction:: scLucid.analysis.summarize_analysis_review_summary

.. autofunction:: scLucid.analysis.validate_analysis_review_summary

.. autofunction:: scLucid.analysis.validate_analysis_module_completeness

characterize_clusters
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.analysis.characterize_clusters
