Preprocessing Module API Reference
===================================

.. automodule:: scLucid.preprocess
   :members:
   :undoc-members:
   :show-inheritance:

Configuration Classes
---------------------

WorkflowConfig
~~~~~~~~~~~~~~

.. autoclass:: scLucid.preprocess.WorkflowConfig
   :members:
   :undoc-members:

NormalizationConfig
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.preprocess.NormalizationConfig
   :members:
   :undoc-members:

HVGConfig
~~~~~~~~~

.. autoclass:: scLucid.preprocess.HVGConfig
   :members:
   :undoc-members:

ScalingConfig
~~~~~~~~~~~~~

.. autoclass:: scLucid.preprocess.ScalingConfig
   :members:
   :undoc-members:

IntegrationConfig
~~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.preprocess.IntegrationConfig
   :members:
   :undoc-members:

Core Functions
--------------

run_preprocessing
~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.preprocess.run_preprocessing

Preprocessing Review Contract
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``run_preprocessing`` stores an auditable review bundle at
``adata.uns["sclucid"]["preprocess"]["review_summary"]``. The bundle includes:

* ``applied_parameter_summary``: effective normalization, HVG, regression,
  scaling, PCA, batch-correction, and graph parameters used by the run.
* ``layer_transition_summary``: expression-layer and embedding transitions from
  counts to normalized/scaled layers, PCA, and optional integrated embeddings.
* ``hvg_selection_evidence_summary``: HVG method, input layer, selected count,
  selected fraction, input statistics, and excluded gene-type evidence.
* ``tumor_aware_batch_correction_warnings``: tumor-context warnings when batch
  correction could over-correct malignant-state, clone, patient, or TME signal.
* ``downstream_analysis_recommendations``: preprocessing-to-analysis handoff
  guidance, including which representation should be considered downstream.
* ``preprocess_readiness`` and ``review_action_items``: machine-readable status,
  score, blockers, review reasons, and prioritized human review actions.
* ``evidence_bundle``: shared scLucid evidence schema view of preprocessing
  evidence, action items, context, and reproducibility-critical parameters.

.. autofunction:: scLucid.preprocess.validate_preprocessing_review_summary

.. autofunction:: scLucid.preprocess.enrich_preprocessing_review_summary

normalize_data
~~~~~~~~~~~~~~

.. autofunction:: scLucid.preprocess.normalize_data

find_hvgs
~~~~~~~~~

.. autofunction:: scLucid.preprocess.find_hvgs

scale_data
~~~~~~~~~~

.. autofunction:: scLucid.preprocess.scale_data

batch_correction
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.preprocess.batch_correction
