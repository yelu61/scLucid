QC Module API Reference
=======================

.. automodule:: scLucid.qc
   :members:
   :undoc-members:
   :show-inheritance:

Configuration Classes
---------------------

QCThresholds
~~~~~~~~~~~~~

.. autoclass:: scLucid.qc.QCThresholds
   :members:
   :undoc-members:

DoubletConfig
~~~~~~~~~~~~~

.. autoclass:: scLucid.qc.DoubletConfig
   :members:
   :undoc-members:

QCWorkflowConfig
~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.qc.QCWorkflowConfig
   :members:
   :undoc-members:

Core Functions
--------------

run_qc
~~~~~~

.. autofunction:: scLucid.qc.run_qc

recommend_qc_policy
~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.recommend_qc_policy

apply_qc_policy
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.apply_qc_policy

calculate_qc_metric
~~~~~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.calculate_qc_metric

run_standard_qc
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.run_standard_qc

QC Review Contract
~~~~~~~~~~~~~~~~~~

``run_standard_qc`` stores an auditable review bundle at
``adata.uns["sclucid"]["qc"]["review_summary"]["data"]``. The bundle includes:

* ``policy_flow``: the canonical QC decision narrative:
  profile dataset -> propose candidate thresholds -> score biological risk ->
  choose/recommend policy -> emit reviewer table -> optionally apply.
* ``decision_table``: reviewer-facing per-threshold table with recommended value,
  applied value, metric, source, confidence, evidence, active filtering flag,
  affected cells, affected fraction, review-required status, biological guardrail,
  and risk note.
* ``recommended_threshold_summary``: compact map of recommended thresholds, applied
  thresholds, recommendation method, confidence interval, evidence, and final source.
* ``applied_threshold_summary`` and ``user_override_summary``: explicit record of
  final QC thresholds and any user-specified values that overrode recommendations.
* ``evidence_chain``: ordered summary of recommendation, threshold application,
  sample-level thresholds, filtering result, and output health.
* ``execution_trace``: schema version, executed steps, sample context, threshold mode,
  recommendation availability, and tumor-aware status.
* ``sample_threshold_summary`` and ``tumor_aware_summary``: per-sample adaptive
  thresholds plus tumor-aware warnings such as mitochondrial filtering being disabled.
* ``doublet_evidence_summary``: doublet prediction rates, score ranges,
  heterotypic/homotypic risk metadata, external-evidence notes, and optional
  ``benchmark_decision`` fields when validation evidence is attached. The compact
  benchmark decision records the recommended default doublet mode, primary method,
  candidate ``algorithm_weight`` for algorithm-plus-heuristic fusion, benchmark
  deltas versus algorithm-only behavior, and whether manual review is required.
* ``output_health``: downstream-safety checks including retained cells and missing QC
  metrics.
* ``downstream_preprocess_recommendations``: next-step preprocessing guidance derived
  from QC retention, sample structure, tumor context, and available layers.
* ``qc_readiness``: machine-readable verdict, score, blockers, and review reasons for
  deciding whether to proceed to preprocessing.
* ``review_action_items``: prioritized human review actions generated from output
  health, tumor-aware warnings, user overrides, and benchmark results.
* ``reproducibility_manifest``: compact record of executed steps, data shape, layer
  availability, applied thresholds, threshold sources, and config snapshots.
* ``evidence_bundle``: shared scLucid evidence schema view of QC decisions, evidence,
  action items, context, and reproducibility metadata for cross-stage reporting.
* ``module_maturity``: benchmark-module completeness assessment for the QC result,
  including whether the review bundle satisfies the frozen QC module contract.
* ``benchmark_summary``: profile-aware retention, stratified-retention, marker-fidelity,
  risk-level, reason, and action-item checks for PBMC, tissue, tumor, and cell-line style
  datasets when pre/post-filtering data are available.

When ``save_dir`` is set, the same contract is exported as
``qc_review_summary.json`` and summarized in ``qc_review_summary.md``. QC benchmark
results are additionally exported as ``qc_benchmark.json`` and ``qc_benchmark.md``.

Phase 2 validation outputs are consolidated by
``validation/qc/build_figure2_qc_evidence_package.py`` into
``validation_outputs/qc_figure2_package/``. The package contains harmonized
Figure 2 source data plus a claim-level scorecard. This is the recommended
place to inspect whether a QC claim is currently supported, partial, or
contract-only.

.. autofunction:: scLucid.qc.build_qc_decision_table

.. autofunction:: scLucid.qc.enrich_qc_decision_table_for_review

.. autofunction:: scLucid.qc.validate_qc_review_summary

QC Module Maturity
~~~~~~~~~~~~~~~~~~

QC is the first scLucid benchmark module. A benchmark-grade QC result should be
auditable from the AnnData object alone: metrics, decisions, recommendations,
overrides, output health, downstream guidance, action items, and reproducibility
metadata should all be present under ``adata.uns["sclucid"]["qc"]``.

Current evidence status is intentionally scoped: threshold-decision
auditability and tumor-aware biological-fidelity proxies are supported by the
local Figure 2 package; Kang demuxlet supports doublet calibration with review
required; ambient diagnostics remain contract-only until a full raw 10x
benchmark is added.

Inspect the frozen QC module contract:

.. code-block:: python

   import scLucid as scl

   contract = scl.qc.get_qc_module_contract()
   print(contract["stable_entrypoints"])

Validate a QC result:

.. code-block:: python

   validation = scl.qc.validate_qc_module_completeness(adata)
   print(validation["valid"])
   print(validation["summary"])

Create a compact display summary:

.. code-block:: python

   review = adata.uns["sclucid"]["qc"]["review_summary"]
   compact = scl.qc.summarize_qc_review_summary(review)
   print(compact)

.. autofunction:: scLucid.qc.get_qc_module_contract

.. autofunction:: scLucid.qc.validate_qc_module_completeness

.. autofunction:: scLucid.qc.summarize_qc_review_summary

.. autofunction:: scLucid.qc.build_qc_module_maturity_assessment

QC Benchmarking
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.evaluate_qc_benchmark

.. autofunction:: scLucid.qc.build_qc_benchmark_assessment

.. autofunction:: scLucid.qc.compute_retention_metrics

.. autofunction:: scLucid.qc.compute_marker_fidelity

.. autofunction:: scLucid.qc.export_qc_benchmark_report

.. autofunction:: scLucid.qc.render_qc_benchmark_compact_markdown

run_advanced_qc
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.run_advanced_qc

filter_cells
~~~~~~~~~~~~

.. autofunction:: scLucid.qc.filter_cells

predict_doublets
~~~~~~~~~~~~~~~

.. autofunction:: scLucid.qc.predict_doublets

Advanced Features
-----------------

Adaptive Threshold Learning
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: scLucid.qc.AdaptiveThresholdLearner
   :members:
   :undoc-members:

Incremental QC
~~~~~~~~~~~~~

.. autoclass:: scLucid.qc.IncrementalQC
   :members:
   :undoc-members:

Caching
~~~~~~~

.. autoclass:: scLucid.qc.CacheConfig
   :members:
   :undoc-members:
