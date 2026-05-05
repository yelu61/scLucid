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

* ``qc_input_context``: QC handoff evidence consumed by preprocessing,
  including QC readiness, retained cells, required QC metrics, and whether a
  ``counts`` layer is available.
* ``applied_parameter_summary``: effective normalization, HVG, regression,
  scaling, PCA, batch-correction, and graph parameters used by the run.
* ``layer_transition_summary``: expression-layer and embedding transitions from
  counts to normalized/scaled layers, PCA, and optional integrated embeddings.
* ``step_evidence_summary``: per-step audit records with status, inputs,
  outputs, applied parameters, linked review fields, and review flags.
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
* ``module_maturity``: module-level completeness check against the frozen
  preprocessing contract.

Preprocessing Module Maturity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The preprocessing module is designed to be the second benchmark-grade module
after QC. Its stable contract connects QC output to preprocessing choices and
records the parameter and layer lineage needed for reviewable downstream
analysis. ``step_evidence_summary`` is the most direct way to inspect the
workflow as a series of auditable decisions: each step reports whether it ran,
where its inputs came from, which output key or layer was produced, which
parameters were applied, and which review-summary fields support the step.

.. code-block:: python

   adata = scl.qc.run_standard_qc(adata, config=qc_config)
   adata = scl.pp.run_preprocessing(adata, config=preprocess_config)

   validation = scl.pp.validate_preprocess_module_completeness(adata)
   compact = scl.pp.summarize_preprocess_review_summary(
       adata.uns["sclucid"]["preprocess"]["review_summary"]
   )

   print(validation["valid"])
   print(compact["maturity_status"])
   print(compact["qc_input_available"])
   print(compact["layers_present"])
   print(compact["step_status_counts"])
   print(compact["review_required_steps"])

``validate_preprocess_module_completeness`` checks that the AnnData object
contains the expected preprocessing namespace, review summary, normalized layer,
PCA embedding, and HVG evidence. By default, review-required states are reported
as warnings rather than hard failures, so exploratory analysis can continue
while still making uncertainty explicit.

``summarize_preprocess_review_summary`` returns a compact product-facing view
for notebooks, CLIs, and reports. It includes preprocessing maturity, readiness,
QC handoff status, layers and embeddings present, HVG selection, PCA and UMAP
state, per-step status counts, review-required steps, and downstream analysis
readiness.

.. autofunction:: scLucid.preprocess.get_preprocess_module_contract

.. autofunction:: scLucid.preprocess.validate_preprocess_module_completeness

.. autofunction:: scLucid.preprocess.summarize_preprocess_review_summary

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
