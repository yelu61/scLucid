#!/usr/bin/env python3
"""Download or generate the CellBender tiny ambient-RNA fixture.

The CellBender quick-start tutorial provides a script that trims the 10x
Genomics ``heart10k`` dataset down to a tiny demo file:

- 500 high-UMI barcodes (likely cells)
- 50,000 low-UMI barcodes (likely empty droplets)
- Top 100 most highly expressed genes

This script clones the CellBender examples repository (shallow clone),
runs the tutorial generator, and copies the resulting
``tiny_raw_feature_bc_matrix.h5ad`` into ``data/cellbender_tiny.h5ad``.

Usage
-----
    python scripts/download_cellbender_tiny.py --output data/cellbender_tiny.h5ad

Requirements
------------
- ``git``
- ``cellbender`` Python package (the generator imports it)
- Internet access to download the full heart10k dataset (~150 MB)

Notes
-----
The tiny fixture is intentionally small and should only be used for testing
ambient-RNA diagnostics. CellBender's documentation warns against trimming or
pre-filtering real datasets before running ``remove-background``.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

CELLBENDER_EXAMPLES_URL = "https://github.com/broadinstitute/CellBender.git"
GENERATOR_SCRIPT = "examples/remove_background/generate_tiny_10x_dataset.py"


def _run_in_env(cmd: list[str], cwd: Path, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command in the given working directory."""
    logger.info("Running: %s", " ".join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=check, text=True, capture_output=True)


def download_cellbender_tiny(output_path: Path, branch: str | None = None) -> Path:
    """Clone CellBender examples, run the tiny-dataset generator, copy result."""
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    clone_cmd = ["git", "clone", "--depth", "1"]
    if branch:
        clone_cmd.extend(["--branch", branch])
    clone_cmd.extend([CELLBENDER_EXAMPLES_URL, "cellbender"])

    with tempfile.TemporaryDirectory(prefix="cellbender_tiny_") as tmp:
        tmp_path = Path(tmp)
        logger.info("Cloning CellBender examples into %s", tmp_path)
        _run_in_env(clone_cmd, cwd=tmp_path)

        repo = tmp_path / "cellbender"
        generator = repo / GENERATOR_SCRIPT
        if not generator.exists():
            raise FileNotFoundError(f"Generator script not found in cloned repo: {generator}")

        logger.info("Generating tiny 10x dataset (downloads heart10k; may take a few minutes)")
        result = _run_in_env(
            [sys.executable, str(generator)],
            cwd=generator.parent,
            check=False,
        )
        if result.returncode != 0:
            logger.error("Generator stdout:\n%s", result.stdout)
            logger.error("Generator stderr:\n%s", result.stderr)
            raise RuntimeError("CellBender tiny dataset generator failed")

        # The generator writes to its own directory.
        tiny_file = generator.parent / "tiny_raw_feature_bc_matrix.h5ad"
        if not tiny_file.exists():
            raise FileNotFoundError(f"Expected output not found: {tiny_file}")

        # Make the output conform to the scLucid benchmark contract.
        adata = ad.read_h5ad(tiny_file)
        counts = adata.X.copy() if adata.X is not None else adata.layers["counts"].copy()
        if sp.issparse(counts):
            counts = counts.astype(np.int32)
        else:
            counts = np.asarray(counts, dtype=np.int32)
        adata.X = counts.copy()
        adata.layers["counts"] = counts.copy()
        barcode_totals = np.asarray(adata.X.sum(axis=1)).ravel()
        # The CellBender tutorial generator selects high-UMI barcodes first,
        # followed by low-UMI barcodes. Ranking by total counts makes that
        # fixture contract explicit and easy to validate.
        barcode_rank = np.argsort(np.argsort(-barcode_totals)) + 1
        likely_cell = barcode_rank <= 500
        adata.obs["sample"] = "heart10k_tiny"
        adata.obs["condition"] = "ambient_fixture"
        adata.obs["cell_type"] = ""
        adata.obs["cell_subtype"] = ""
        adata.obs["barcode_total_counts"] = barcode_totals.astype(np.int64)
        adata.obs["barcode_rank"] = barcode_rank.astype(np.int64)
        adata.obs["likely_cell"] = likely_cell
        adata.obs["likely_empty_droplet"] = ~likely_cell
        adata.obs["dataset"] = "cellbender_tiny"
        adata.var["gene_symbol"] = adata.var_names.astype(str)
        adata.uns.setdefault("sclucid", {})["dataset"] = {
            "dataset_key": "cellbender_tiny",
            "title": "CellBender tiny ambient-RNA fixture (10x heart10k subset)",
            "geo_accession": None,
            "pmid": None,
            "species": "mouse",
            "tissue": "heart",
            "citation": "Derived from the CellBender quick-start tutorial using 10x Genomics heart_10k_v3.",
            "processing": {
                "tool": "scripts/download_cellbender_tiny.py",
                "raw_counts_layer": "counts",
                "normalized_layer": None,
            },
        }
        adata.write_h5ad(tiny_file, compression="gzip")

        shutil.copy2(tiny_file, output_path)
        logger.info("Copied %s -> %s", tiny_file, output_path)

    return output_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/cellbender_tiny.h5ad"),
        help="Destination path for the tiny fixture (default: data/cellbender_tiny.h5ad).",
    )
    parser.add_argument(
        "--branch",
        default=None,
        help="CellBender git branch to clone. If omitted, the repository default branch is used.",
    )
    args = parser.parse_args(argv)

    try:
        download_cellbender_tiny(args.output, branch=args.branch)
        logger.info("Done.")
        return 0
    except Exception as exc:
        logger.exception("Failed to download CellBender tiny fixture: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
