Frequently Asked Questions
==========================

This page collects common questions about installing, using, and contributing
to scLucid. Sections will grow as questions accumulate; see the
:doc:`quickstart` and :doc:`best_practices` pages for orientation.

Installation
------------

**Q: Which install extra should I use?**

For a typical analysis: ``pip install "sclucid[analysis,de]"``. The full
``[all]`` extra pulls in heavy optional backends (scVI, scVelo, infercnvpy,
squidpy, etc.) and is only needed if you intend to use those tools.

**Q: Why do I see "Could not import optional tools backend"?**

scLucid degrades gracefully when an optional backend's dependencies are
missing. The warning identifies which backend is unavailable — install the
corresponding extra (e.g. ``pip install "sclucid[tools]"`` for squidpy,
scVelo, infercnvpy) to enable it. As of v0.1, ``squidpy`` is required for
the ``spatial`` backend; missing optional dependencies are now skipped
silently rather than raising an :class:`ImportWarning`.

Workflow
--------

**Q: When should I use ``run_pipeline`` versus the per-stage API?**

Use ``run_pipeline(stages=[...])`` for standard analysis: it propagates
``AnalysisContext`` (species, tissue, cancer type) through QC, preprocess,
and analysis, and records contract validation results under
``adata.uns["sclucid"]``. Drop to the per-stage API
(``run_standard_qc``, ``run_preprocessing``, ``run_standard_analysis``)
when you need to inspect or override intermediate results.

**Q: How do I tell scLucid that my data is tumor tissue?**

Pass ``dataset_type="tumor_tissue"`` (and optionally ``cancer_type=...``)
to ``run_pipeline``. The QC stage will then preserve high-mitochondrial
cells as a tumor-aware warning rather than removing them outright, and
later stages can adapt their assumptions accordingly. See
:doc:`workflow_hardening` for the validated tumor path on PDAC.

Reproducibility
---------------

**Q: How do I save the full configuration of a run?**

Every workflow stage writes its effective configuration to
``adata.uns["sclucid"][stage]["workflow_config"]`` and the inheritance
chain to ``adata.uns["sclucid"][stage]["config_lineage"]``. Save the
``.h5ad`` and the stage save directories together; both round-trip.

**Q: Are seeds set automatically?**

Yes. Each workflow accepts a ``random_state`` parameter (default 42) that
is threaded through HVG selection, PCA, neighbor graphs, UMAP, Leiden, and
recommendation samplers.
