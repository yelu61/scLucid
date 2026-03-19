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
