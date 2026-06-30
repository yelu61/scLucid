"""Sphinx configuration for scLucid documentation.

Build with:
    cd docs
    sphinx-build -b html source build/html

Or via the maintained env:
    /Users/luye/micromamba/envs/scrna-env/bin/sphinx-build \
        -b html docs/source docs/build/html
"""

from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

# -- Path setup ---------------------------------------------------------------
# Make the source package importable for autodoc.
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "src"))

# -- Project information -----------------------------------------------------
project = "scLucid"
author = "Ye Lu"
copyright = f"{datetime.now():%Y}, {author}"

try:
    from importlib.metadata import PackageNotFoundError, version as _pkg_version

    release = _pkg_version("sclucid")
except Exception:  # pragma: no cover - fall back during local dev installs
    release = "0.1.0"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
]

# Surface both Google and NumPy style docstring sections; the codebase uses
# a NumPy-style "Parameters / Returns: / Examples:" pattern.
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "inherited-members": False,
}
autodoc_typehints = "description"
autodoc_member_order = "bysource"

# Mock heavy / optional dependencies so autodoc does not fail when they are
# not installed in the docs-build environment. Keep this list narrow — only
# include packages whose top-level import is unconditionally required by a
# scLucid submodule.
autodoc_mock_imports = [
    "scanpy",
    "anndata",
    "scrublet",
    "celltypist",
    "cosg",
    "hdbscan",
    "harmonypy",
    "scanorama",
    "scvi",
    "bbknn",
    "combat",
    "pydeseq2",
    "gseapy",
    "sccoda",
    "scvelo",
    "infercnvpy",
    "pyscenic",
    "loompy",
    "cellphonedb",
    "triku",
    "squidpy",
    "leidenalg",
    "louvain",
    "umap",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {".rst": "restructuredtext"}
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = f"{project} {release}"
html_show_sourcelink = False

html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

# Don't fail the build on unresolved references during early-stage docs work,
# but do surface them as warnings.
nitpicky = False
todo_include_todos = True
