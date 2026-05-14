R Package Parity Matrix
========================

scLucid ports four mature R bioinformatics packages to pure Python so that
single-cell workflows that mix R and Python tools can run end-to-end in
Python — without rpy2, without R installation, without serialization round
trips. This page documents what is in parity with the R original, what is
partial, what is intentionally simplified, and what is out of scope.

How To Read This Matrix
-----------------------

Each table compares the R-package API on the left with the scLucid Python
equivalent on the right.

Status legend:

- ✅ **Done** — feature-complete equivalent; algorithm matches the
  published method.
- 🟡 **Partial** — works for the common case but a subset of the original
  options/methods is implemented. Use with awareness of the limitations
  noted.
- 🔴 **Missing** — not yet implemented; planned. Pull requests welcome.
- ⚪ **Out of scope** — intentionally not planned; usually because the R
  feature depends on heavyweight infrastructure (e.g. another R package,
  GPU-only kernels) that does not fit the scLucid philosophy.

The tables target the most-used surface of each R package. If a function you
need is missing from this page entirely, open an issue with the request and
the R reference.

.. contents:: Contents
   :local:
   :depth: 2

pyBayesPrism
------------

**Target**: `BayesPrism <https://github.com/Danko-Lab/BayesPrism>`_ — Chu et
al., *Nature Cancer* 2022. Bayesian cell-type-fraction inference from bulk
RNA-seq using a single-cell reference, with a Gibbs sampler over a Dirichlet
prior on the proportion vector.

scLucid module: ``scLucid.tools.pyBayesPrism``.
Top-level re-exports under ``scLucid.tools`` for the most-used classes.

.. list-table::
   :header-rows: 1
   :widths: 30 30 10 30

   * - R API
     - scLucid equivalent
     - Status
     - Notes
   * - ``new.prism()``
     - ``BayesPrism(reference, mixture, config=PrismConfig(...))``
     - ✅
     - Pydantic ``PrismConfig`` captures all sampler hyperparameters.
   * - ``run.prism()``
     - ``BayesPrism.run_deconvolution()``
     - 🟡
     - Initial NNLS estimate + Gibbs updates implemented; multi-chain
       parallelism uses ``ProcessPoolExecutor`` rather than R's foreach
       backend. Burn-in / chain length match.
   * - ``cleanup.genes()``
     - ``pyBayesPrism.cleanup_genes(...)`` and
       ``BayesPrism.cleanup_genes()``
     - ✅
     - Ribosomal / mitochondrial / sex-chromosome filters plus
       user-supplied patterns.
   * - ``find.outlier.genes()``
     - ``pyBayesPrism.find_outlier_genes(...)``
     - ✅
     - Outlier detection against reference-vs-mixture mismatch.
   * - ``get.fraction(prism)``
     - ``BayesPrism.get_fraction()``
     - ✅
     - Posterior mean fractions; identical schema.
   * - ``get.exp(prism, cell.type=...)``
     - ``BayesPrism.get_expression(cell_type=...)``
     - ✅
     - Cell-type-specific expression Z.
   * - ``compute.cv()``
     - ``BayesPrism.compute_cv()``
     - ✅
     - Coefficient of variation across chains for posterior uncertainty.
   * - ``BayesPrismReference`` (implicit S4 reference holder)
     - ``BayesPrismReference(reference, cell_type_labels, ...)``
     - ✅
     - Standalone class wrapping ``phi`` construction so the reference can
       be reused across multiple bulk mixtures.
   * - Embedding / NMF on Z (``extract.gene.programs``)
     - ``BayesPrismEmbedding(prism).run_nmf()``,
       ``get_gene_programs()``, ``get_program_usage()``,
       ``get_top_genes()``
     - ✅
     - NMF via scikit-learn instead of R's NMF package.
   * - ``GibbsSampler`` low-level access
     - ``pyBayesPrism.GibbsSampler``
     - ✅
     - Exposed for advanced custom workflows; Numba-accelerated kernel
       with Python fallback.
   * - Visualization (``plot.heatmap``, etc.)
     - ``plot_fraction``, ``plot_correlation``, ``plot_stacked_bar``,
       ``plot_gene_programs``, ``plot_program_usage``,
       ``plot_validation_scatter``, ``plot_cv``
     - 🟡
     - Matplotlib-based; covers the most-requested heatmap and bar/scatter
       figures. The R package's interactive Shiny dashboard is out of
       scope.
   * - Pre-built tumor reference catalogs (TCGA-derived)
     - —
     - ⚪
     - Out of scope: ship your own reference. The R package's curated
       references are large and licensed separately.

**Known limitations**

- The Gibbs sampler exposed via ``scLucid.tools.bulk._bayesprism_gibbs_sample``
  is a *simplified* sampler intended for the unified ``deconvolve_bulk()``
  convenience path. For full multi-chain inference, use the
  ``pyBayesPrism.BayesPrism`` class directly.
- No batch-effect updating on the reference matrix (the R "updated theta /
  Z" workflow is implemented but advanced two-step refinement on cell-type-
  specific Z is partial).

pyMonocle3
----------

**Target**: `Monocle3 <https://cole-trapnell-lab.github.io/monocle3/>`_ —
Cao et al., *Nature* 2019. Trajectory inference via principal graph
embedding, Leiden clustering, pseudotime by shortest-path distance, and
Moran's I differential expression along the graph.

scLucid module: ``scLucid.tools.pyMonocle3``.

.. list-table::
   :header-rows: 1
   :widths: 30 30 10 30

   * - R API
     - scLucid equivalent
     - Status
     - Notes
   * - ``new_cell_data_set()`` / ``CellDataSet`` (S4)
     - ``CellDataSet(expression_data, cell_metadata, gene_metadata)`` /
       ``new_cell_data_set(...)``
     - ✅
     - Constructor signature matches; ``save()`` / ``load()`` round-trip
       via pickle.
   * - Conversion to/from Seurat / SingleCellExperiment
     - ``create_cds_from_scanpy(adata)`` /
       ``export_to_scanpy(cds)``
     - ✅
     - The Python ecosystem analogue is AnnData; both directions are
       implemented.
   * - ``preprocess_cds()``
     - ``preprocess_cds(cds, num_dim=50, ...)``
     - ✅
     - Normalize + log + PCA. Methods: ``"PCA"`` and ``"LSI"`` (LSI for
       ATAC-style data).
   * - ``estimate_size_factors()``
     - ``estimate_size_factors(cds)``
     - ✅
     -
   * - ``detect_genes()``
     - ``detect_genes(cds, min_expr=...)``
     - ✅
     -
   * - ``align_cds()`` (batch correction)
     - ``align_cds(cds, alignment_group=...)``
     - 🟡
     - Mean-centering per batch only; R's MNN-based mutual nearest neighbor
       alignment is not implemented (use scanpy's harmony / BBKNN through
       ``scLucid.preprocess.integrate`` instead).
   * - ``reduce_dimension()``
     - ``reduce_dimension(cds, reduction_method="UMAP")``
     - ✅
     - Supports ``"UMAP"``, ``"tSNE"``, ``"PCA"``.
   * - ``cluster_cells()``
     - ``cluster_cells(cds, resolution=...)``
     - ✅
     - Leiden via leidenalg/igraph; resolution parameter equivalent.
   * - ``partition_cells()`` (Louvain partitions)
     - ``partition_cells(cds)``
     - 🟡
     - Connected-component partitions on the kNN graph; R uses Louvain
       directly. Results are usually equivalent for well-separated
       partitions.
   * - ``learn_graph()`` (principal graph via reverse graph embedding)
     - ``learn_graph(cds, ...)``
     - ✅
     - Minimum spanning tree on cluster centroids per partition, with
       optional loop closing and branch pruning. The R package's full
       SimplePPT/DDRTree reverse graph embedding is approximated with a
       cluster-centroid backbone, which is faster but slightly less smooth
       on dense trajectories.
   * - ``order_cells()``
     - ``order_cells(cds, root_cells=...)``
     - ✅
     - String selectors for cell name / cluster id / partition id all
       supported; integer indices and lists also work. Pseudotime via
       Dijkstra shortest path on the principal graph.
   * - ``graph_test()`` (Moran's I along graph)
     - ``graph_test(cds, neighbor_graph="principal_graph", k=25)``
     - ✅
     - Moran's I per gene with permutation p-values.
   * - ``find_gene_modules()`` (UMAP + Louvain on TF-IDF gene matrix)
     - ``calculate_gene_modules(cds, ...)``
     - ✅
     -
   * - ``top_markers()``
     - ``top_markers(cds, group_cells_by=...)``
     - ✅
     -
   * - ``compare_genes()``
     - ``compare_genes(cds, ...)``
     - ✅
     -
   * - ``fit_models()`` / ``coefficient_table()`` regression workflow
     - ``pseudotime_de(cds)`` (limited)
     - 🟡
     - The full quasi-Poisson / negative-binomial regression with
       arbitrary formula support from R is not yet implemented. Use
       ``pyMonocle3.pseudotime_de`` for the standard pseudotime regression
       and ``scLucid.analysis.differential_expression`` for general DE.
   * - ``choose_graph_segments()`` interactive subsetting
     - ``choose_graph_segments(cds, ...)``
     - 🟡
     - Programmatic subsetting works; the R Shiny interactive picker is
       out of scope.
   * - ``plot_cells()``, ``plot_genes_by_partition()``,
       ``plot_pc_variance_explained()``
     - ``plot_cells()``, ``plot_genes_by_group()``,
       ``plot_pseudotime_heatmap()``, ``plot_trajectory()``
     - 🟡
     - Common figure types covered; the R ggplot2-style faceting is
       achievable via matplotlib but requires more user code.

**Known limitations**

- ``learn_graph`` uses a cluster-centroid MST backbone rather than the full
  reverse graph embedding (DDRTree / SimplePPT). Trajectories with very
  fine internal branches may need ``minimal_branch_len`` tuning.
- ``partition_cells`` uses connected components of the kNN graph rather
  than Louvain modularity. Behavior diverges on highly connected datasets.

pyCellChat
----------

**Target**: `CellChat v2 <https://github.com/jinworks/CellChat>`_ — Jin et
al., *Nature Communications* 2021 + 2024 v2 update. Ligand-receptor
inference of cell-cell communication, with multi-condition comparison and
pathway-level summarization.

scLucid module: ``scLucid.tools.pyCellChat``.

.. list-table::
   :header-rows: 1
   :widths: 30 30 10 30

   * - R API
     - scLucid equivalent
     - Status
     - Notes
   * - ``createCellChat(object, ...)``
     - ``CellChat(data, labels, config=CellChatConfig(...))`` /
       ``create_cellchat_from_scanpy(adata, group_by=...)``
     - ✅
     - Adapter for AnnData input is included.
   * - ``CellChatDB.human`` / ``CellChatDB.mouse``
     - ``CellChatDB(species="human")`` / ``get_default_database()``
     - ✅
     - Bundled ligand-receptor pairs across signaling categories. Versions
       documented in the database docstring.
   * - ``subsetData()``
     - ``CellChat.preprocess_data()``
     - ✅
     -
   * - ``identifyOverExpressedGenes()`` /
       ``identifyOverExpressedInteractions()``
     - ``CellChat.identify_overexpressed_genes()``
     - ✅
     - Combines both R functions into one pass.
   * - ``projectData()`` (smoothing on PPI network)
     - —
     - 🔴
     - Missing. The R version smooths expression on the PPI graph before
       communication inference; scLucid currently uses raw subset
       expression. Track in issue #TODO.
   * - ``computeCommunProb()``
     - ``CellChat.compute_communication_prob()``
     - ✅
     - Default mass-action / triMean aggregation; Hill function for L-R
       interaction strength.
   * - ``computeCommunProbPathway()``
     - Available as part of
       ``CellChat.compute_communication_prob(pathway=True)``
     - 🟡
     - Per-pathway aggregation works; signaling-pathway-only object that
       mirrors the R ``netP`` slot is implicit.
   * - ``filterCommunication(min.cells=10)``
     - ``CellChat.filter_communication(min_cells=10)``
     - ✅
     -
   * - ``aggregateNet()`` (sender-receiver matrices)
     - Available via internal ``compute_communication_prob`` output
     - 🟡
     - Aggregated nets are computed; expose as a stable top-level helper
       is planned (issue #TODO).
   * - ``netAnalysis_computeCentrality()``
     - ``compute_centrality(cellchat, ...)`` /
       ``CellChat.compute_network_centrality()``
     - ✅
     - Out-/in-degree, betweenness, eigenvector centrality per cell type.
   * - ``netAnalysis_signalingRole_*`` (dominant senders/receivers)
     - ``identify_roles(cellchat)``
     - ✅
     -
   * - ``mergeCellChat()`` (multi-condition)
     - ``merge_cellchat_objects(objects, names=[...])``
     - ✅
     -
   * - ``compareInteractions()`` / ``rankNet()``
     - ``compare_cellchat_objects(merged)`` /
       ``CellChat.compare_interactions()``
     - ✅
     - Number-of-interactions and interaction-strength comparison across
       conditions.
   * - ``identifyConservedAndContextDependentLR()``
     - ``identify_conserved_pathways(merged)`` /
       ``identify_differential_pathways(merged)``
     - ✅
     -
   * - ``netSimilarity()``
     - ``compute_network_similarity(cellchat_A, cellchat_B)``
     - ✅
     - Network similarity for cross-condition comparison.
   * - Visualization: ``netVisual_circle``, ``netVisual_chord_*``,
       ``netVisual_heatmap``, ``netVisual_bubble``
     - ``plot_circle_network``, ``plot_chord_diagram``, ``plot_heatmap``,
       ``plot_bubble``, ``plot_contribution``,
       ``plot_signaling_gene_expression``
     - 🟡
     - Major figure types covered. The R package's spatial overlay plot
       (``netVisual_spatial``) is missing.
   * - ``identifyCommunicationPatterns()`` (NMF-based)
     - —
     - 🔴
     - Missing. Planned.
   * - Spatial CellChat (v2 feature)
     - —
     - 🔴
     - Not yet ported.

**Known limitations**

- PPI-network smoothing of expression (R ``projectData()``) is not yet
  applied; communication probability is computed on raw subset expression.
  Differences are usually small for high-cell-count datasets and large for
  sparse data.
- ``identifyCommunicationPatterns`` (NMF over outgoing/incoming signaling)
  is on the roadmap but not implemented.

pyDWLS
------

**Target**: `DWLS <https://github.com/dtsoyuzu/DWLS>`_ — Tsoucas et al.,
*Cell Reports* 2019. Dampened Weighted Least Squares for cell-type
proportion inference from bulk RNA-seq using a single-cell reference.

scLucid module: ``scLucid.tools.pyDWLS``.

.. list-table::
   :header-rows: 1
   :widths: 30 30 10 30

   * - R API
     - scLucid equivalent
     - Status
     - Notes
   * - ``solveDampenedWLS(S, B)``
     - ``DampenedWLS(dampen_factor=...).solve(S, b)``
     - ✅
     - Iterative reweighting with weights = ``1/μ²`` capped at a
       dampening-controlled quantile, matching the published algorithm.
   * - ``buildSignatureMatrixUsingSeurat()`` (Seurat differential
       expression based)
     - ``SignatureBuilder(method="mean")`` and
       ``MarkerSelector.select(method=...)``
     - 🟡
     - We provide ratio / difference / fold-change marker selection but do
       not call out to Seurat. Use ``scLucid.analysis.differential_expression``
       upstream if you want test-based marker ranking.
   * - ``buildSignatureMatrixMAST()`` (MAST-based)
     - —
     - ⚪
     - Out of scope: requires the MAST R package; use ``pydeseq2`` or
       ``scanpy.tl.rank_genes_groups`` instead and feed results into
       ``SignatureBuilder``.
   * - ``trimData()`` / shared-gene alignment
     - ``align_data(sig_df, bulk_df, min_common_genes=...)``
     - ✅
     - Case-insensitive gene-name intersection.
   * - ``DWLS`` orchestrator (sequence of build → select → solve)
     - ``DWLS(dampen_factor=...).build_signature_matrix(...).deconvolve(bulk)``
     - ✅
     - Mirrors the standard R workflow as a single class for convenience.
   * - Cross-validation harness
     - ``CrossValidator(n_folds=5).run(sc_data, labels)``
     - ✅
     - K-fold split → held-out pseudo-bulk → recover proportions → Pearson
       + RMSE. No direct R analogue; this is a scLucid extension.
   * - Pseudo-bulk generation
     - ``create_pseudo_bulk(sc_data, labels, n_cells=...)``
     - ✅
     - Reusable utility for CV and other validation experiments.
   * - Normalization helpers
     - ``normalize_data(method="cpm")``, ``filter_genes(...)``,
       ``solve_nnls(A, b)``
     - ✅
     -

**Known limitations**

- No MAST-based signature builder (would require MAST R port or significant
  rewrite; out of scope per project philosophy).
- The ``DWLS`` orchestrator's defaults are tuned for the
  ``scLucid.tools.bulk.deconvolve_bulk(method="DWLS")`` path; advanced users
  should pass ``dampen_factor``, ``n_markers``, and ``signature_method``
  explicitly.

Roadmap And Contributions
-------------------------

Most 🟡 entries above are "good enough for the typical case but missing an
edge". Most 🔴 entries are planned. ⚪ entries are intentional design
choices — these are typically the R features that depend on heavy R
infrastructure that does not fit Python-first workflows.

To suggest or contribute a missing API:

1. Open an issue using the **Feature request** template, noting the
   R reference (function name + R package + version).
2. Specify which user layer (workflow / simple API / advanced notebook)
   benefits.
3. If you have a real-data dataset where the missing feature matters,
   include the dataset shape and the analysis intent.

PRs that close a 🔴 entry must include a numerical-fidelity test
(simulated truth or comparison against the R reference output) — not just a
shape check.
