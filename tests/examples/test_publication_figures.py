"""Smoke tests for the publication-figure example scripts.

Each script under ``examples/04_publication_figures/`` should be runnable
in isolation. These tests import each script's ``main()`` and execute it
in a temporary working directory, verifying:

- the script imports cleanly
- ``main()`` runs without error
- the expected PDF artifact is produced
- the PDF starts with the ``%PDF`` magic and is non-empty
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest


EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples" / "04_publication_figures"

SCRIPTS = [
    ("01_umap_annotation.py", "fig01_umap_annotation.pdf"),
    ("02_marker_heatmap.py", "fig02_marker_heatmap.pdf"),
    ("03_volcano_de.py", "fig03_volcano_de.pdf"),
    ("04_cnv_heatmap.py", "fig04_cnv_heatmap.pdf"),
]


def _load_module(script_name: str):
    """Load a publication-figure script as a module."""
    script_path = EXAMPLES_DIR / script_name
    spec = importlib.util.spec_from_file_location(
        f"_pub_fig_test_{script_path.stem}", script_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("script,artifact", SCRIPTS)
def test_publication_figure_script(script, artifact, tmp_path, monkeypatch):
    """Each script must produce its named PDF artifact when ``main()`` runs."""
    # The scripts write to ``results/publication_figures/`` relative to cwd.
    # Move the cwd into the tmp dir so artifacts land there and don't pollute
    # the repo's results/ directory during testing.
    monkeypatch.chdir(tmp_path)
    sys.path.insert(0, str(EXAMPLES_DIR))
    try:
        module = _load_module(script)
        out_path = module.main()
    finally:
        sys.path.pop(0)

    assert isinstance(out_path, Path)
    assert out_path.exists(), f"{script} did not produce {out_path}"
    assert out_path.name == artifact
    # Smallest valid PDF is < 1KB but our figures should be at least a few KB.
    assert out_path.stat().st_size > 2000, (
        f"PDF from {script} is suspiciously small ({out_path.stat().st_size} bytes)"
    )
    # Confirm the file looks like a real PDF.
    with open(out_path, "rb") as fh:
        head = fh.read(8)
    assert head.startswith(b"%PDF"), f"{out_path} does not begin with %PDF magic"


@pytest.mark.parametrize("script,_artifact", SCRIPTS)
def test_publication_figure_script_uses_truetype_fonts(script, _artifact, tmp_path, monkeypatch):
    """Each script must call apply_theme() so PDFs ship editable TrueType fonts.

    apply_theme() sets ``rcParams['pdf.fonttype'] = 42`` which is the
    Illustrator-editable TrueType setting; the alternative ``3`` (Type 3)
    rasterises text. Importing the script must not flip this rcParam to 3.
    """
    monkeypatch.chdir(tmp_path)
    sys.path.insert(0, str(EXAMPLES_DIR))
    try:
        module = _load_module(script)
        module.main()
    finally:
        sys.path.pop(0)

    import matplotlib

    assert matplotlib.rcParams["pdf.fonttype"] == 42
