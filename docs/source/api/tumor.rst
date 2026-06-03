Tumor Module API Reference
==========================

The tumor namespace contains cancer-specific interpretation tools. The current
recommended path is conservative: QC and preprocessing remain light by default,
analysis produces reviewable annotations and cluster-level QC evidence, and the
tumor module consumes those outputs for malignancy and CNV interpretation.

Do not treat tumor calls as a single ground-truth label. Store malignant,
suspect, non-malignant, and unresolved calls with confidence, evidence basis,
and review requirements.

.. automodule:: scLucid.tumor
   :members:
   :undoc-members:
   :show-inheritance:

Malignancy Interpretation
-------------------------

``run_malignancy_interpretation`` combines annotation priors, tumor marker
context, optional CNV score, and optional malignancy-signature score. It writes
cell-level calls to ``adata.obs`` and review summaries under
``adata.uns["sclucid"]["analysis"]["malignancy"]`` so the analysis review
contract can consume the result.

.. autofunction:: scLucid.tumor.malignancy.run_malignancy_interpretation

.. autofunction:: scLucid.tumor.malignancy.score_malignancy

.. autofunction:: scLucid.tumor.malignancy.classify_malignant_cells

CNV Analysis
------------

The lightweight expression-based CNV tools live in ``scLucid.tumor.cnv``.
The optional ``infercnvpy`` bridge is exposed only when that dependency is
installed.

.. autofunction:: scLucid.tumor.cnv.infer_cnv

.. autofunction:: scLucid.tumor.cnv.calculate_cnv_score

.. autofunction:: scLucid.tumor.cnv.find_tumor_cells

.. autofunction:: scLucid.tumor.cnv.identify_clones

Optional infercnvpy Bridge
~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``infercnvpy`` is installed, ``scLucid.tumor.cnv`` also exposes
``run_cnv_analysis`` and ``find_tumor`` as an optional bridge. These names may
be unavailable in light-dependency environments, so production workflows should
check availability before calling them.

Tumor Microenvironment
----------------------

.. autofunction:: scLucid.tumor.deconvolve_tme

.. autofunction:: scLucid.tumor.estimate_stromal_content

.. autofunction:: scLucid.tumor.analyze_immune_infiltration

Therapy And Heterogeneity
-------------------------

These APIs are useful exploratory surfaces, but they still need project-level
validation before being considered benchmark-grade:

.. autofunction:: scLucid.tumor.predict_therapy_response

.. autofunction:: scLucid.tumor.identify_resistance_mechanisms

.. autofunction:: scLucid.tumor.estimate_intratumoral_heterogeneity
