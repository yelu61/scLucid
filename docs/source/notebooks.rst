Notebooks
=========

Project notebooks are intentionally kept outside the documentation source tree in
the top-level ``notebooks/`` directory.

Why Notebooks Are Separate
--------------------------

Notebooks serve a different purpose from docs and examples:

- **docs** explain the supported package interface
- **examples** provide short runnable scripts
- **notebooks** show complete analysis stories, often with real data, richer plots,
  intermediate reasoning, and reviewer-facing outputs

For a package that is being prepared for publication, this separation helps keep
the stable documentation concise while still preserving full analysis narratives.

Recommended Notebook Scope
--------------------------

Use notebooks for:

- real-data walkthroughs
- manuscript-style analysis flows
- detailed parameter exploration
- outputs that are valuable to inspect but too long for docs or examples

Do not use notebooks as the only place where the recommended package path is defined.
That guidance should remain in :doc:`quickstart` and :doc:`best_practices`.

Current Notebook Directory
--------------------------

See the repository-level ``notebooks/`` directory for:

- ``01_quality_control.ipynb``
- ``02_preprocessing.ipynb``
- ``03_clustering_annotation.ipynb``
- ``04_advanced_topics.ipynb``
- ``04_differential_expression.ipynb``
- ``05_trajectory_inference.ipynb``
- ``06_advanced_tools.ipynb``

If these notebooks evolve into a publication-facing tutorial set, keep them in
``notebooks/`` and add a short README that explains their scope and expected inputs.
