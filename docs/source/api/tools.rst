Tools Module
============

.. automodule:: scLucid.tools
   :members:
   :undoc-members:
   :show-inheritance:

Bulk Deconvolution
------------------

.. autofunction:: scLucid.tools.deconvolve_bulk

.. autofunction:: scLucid.tools.run_deconvolution

.. autofunction:: scLucid.tools.differential_abundance

.. autofunction:: scLucid.tools.correlate_abundance_with_clinical

CNV Analysis
------------

.. autofunction:: scLucid.tools.run_cnv_analysis

.. autofunction:: scLucid.tools.find_tumor

Trajectory Analysis (Monocle3)
-------------------------------

.. autoclass:: scLucid.tools.CellDataSet
   :members:
   :undoc-members:

.. autofunction:: scLucid.tools.new_cell_data_set

.. autofunction:: scLucid.tools.create_cds_from_scanpy

.. autofunction:: scLucid.tools.export_to_scanpy

.. autofunction:: scLucid.tools.preprocess_cds

.. autofunction:: scLucid.tools.reduce_dimension

.. autofunction:: scLucid.tools.cluster_cells

.. autofunction:: scLucid.tools.learn_graph

.. autofunction:: scLucid.tools.order_cells

.. autofunction:: scLucid.tools.graph_test

.. autofunction:: scLucid.tools.top_markers

.. autofunction:: scLucid.tools.plot_cells

Cell-Cell Communication (CellChat)
-----------------------------------

.. autoclass:: scLucid.tools.CellChat
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.CellChatDB
   :members:
   :undoc-members:

.. autofunction:: scLucid.tools.get_default_database

.. autofunction:: scLucid.tools.create_cellchat_from_scanpy

.. autofunction:: scLucid.tools.plot_heatmap

Bulk Deconvolution (BayesPrism)
--------------------------------

.. autoclass:: scLucid.tools.PrismConfig
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.BayesPrismReference
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.BayesPrism
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.BayesPrismEmbedding
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.GibbsSampler
   :members:
   :undoc-members:

.. autofunction:: scLucid.tools.plot_fraction

.. autofunction:: scLucid.tools.plot_correlation

.. autofunction:: scLucid.tools.cleanup_genes

.. autofunction:: scLucid.tools.compute_correlation

.. autofunction:: scLucid.tools.compute_rmse

Bulk Deconvolution (DWLS)
--------------------------

.. autoclass:: scLucid.tools.DWLS
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.SignatureBuilder
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.DampenedWLS
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.MarkerSelector
   :members:
   :undoc-members:

.. autoclass:: scLucid.tools.CrossValidator
   :members:
   :undoc-members:

.. autofunction:: scLucid.tools.solve_nnls

.. autofunction:: scLucid.tools.normalize_data

.. autofunction:: scLucid.tools.filter_genes

.. autofunction:: scLucid.tools.create_pseudo_bulk

Gene Regulatory Network Analysis (pySCENIC)
--------------------------------------------

.. autofunction:: scLucid.tools.run_scenic

.. autofunction:: scLucid.tools.run_scenic_batch

.. autofunction:: scLucid.tools.run_scenic_by_group

.. autofunction:: scLucid.tools.analyze_scenic_results

.. autofunction:: scLucid.tools.export_scenic_report

Cell-Cell Communication (CellPhoneDB)
--------------------------------------

.. autofunction:: scLucid.tools.run_cellphonedb

.. autofunction:: scLucid.tools.run_cellphonedb_batch

.. autofunction:: scLucid.tools.run_cellphonedb_by_group

.. autofunction:: scLucid.tools.summarize_cellphonedb

Compositional Analysis (scCODA)
--------------------------------

.. autofunction:: scLucid.tools.run_sccoda

.. autofunction:: scLucid.tools.run_sccoda_batch

.. autofunction:: scLucid.tools.summarize_sccoda

.. autofunction:: scLucid.tools.plot_sccoda_proportion_with_significance
