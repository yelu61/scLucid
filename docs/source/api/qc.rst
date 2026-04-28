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

* ``decision_table``: per-threshold recommendation, applied value, source, confidence,
  evidence, and whether the matching filtering flag was active.
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
* ``benchmark_summary``: profile-aware retention, stratified-retention, marker-fidelity,
  risk-level, reason, and action-item checks for PBMC, tissue, tumor, and cell-line style
  datasets when pre/post-filtering data are available.

When ``save_dir`` is set, the same contract is exported as
``qc_review_summary.json`` and summarized in ``qc_review_summary.md``. QC benchmark
results are additionally exported as ``qc_benchmark.json`` and ``qc_benchmark.md``.

.. autofunction:: scLucid.qc.build_qc_decision_table

.. autofunction:: scLucid.qc.validate_qc_review_summary

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
