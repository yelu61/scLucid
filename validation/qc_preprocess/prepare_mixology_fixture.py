#!/usr/bin/env python3
"""Convert the official sc_mixology RData release into a locked H5AD fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

SOURCE_URL = "https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class.RData"
REFERENCE_DOI = "10.1038/s41592-019-0425-8"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    output = args.output.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing fixture: {output}")

    helper = Path(__file__).with_name("export_mixology_rdata.R")
    with tempfile.TemporaryDirectory(prefix="sclucid_mixology_") as temp_dir:
        subprocess.run(
            ["Rscript", str(helper), str(source), temp_dir],
            check=True,
        )
        temp = Path(temp_dir)
        counts = mmread(temp / "counts.mtx").tocsr().T
        genes = pd.read_csv(temp / "genes.tsv", sep="\t", dtype=str)
        obs = pd.read_csv(temp / "cells.tsv", sep="\t", dtype=str).fillna("")

    if counts.shape != (len(obs), len(genes)):
        raise RuntimeError(
            f"Shape mismatch: counts={counts.shape}, obs={len(obs)}, genes={len(genes)}"
        )
    identity = obs["cell_line_demuxlet"].str.strip()
    if identity.eq("").any() or identity.nunique() != 3:
        raise RuntimeError("Expected three non-empty demuxlet cell-line identities")
    if np.any(counts.data < 0) or not np.allclose(counts.data, np.rint(counts.data)):
        raise RuntimeError("Counts assay must contain non-negative integer-like values")

    obs = obs.set_index("cell_id")
    obs.index.name = None
    numeric_columns = [
        "mapped_to_MT",
        "number_of_genes",
        "total_count_per_cell",
        "non_mt_percent",
    ]
    for column in numeric_columns:
        if column in obs:
            obs[column] = pd.to_numeric(obs[column], errors="raise")
    non_mt = obs["non_mt_percent"].astype(float)
    if non_mt.between(0.0, 1.0).all():
        obs["pct_counts_mt"] = 100.0 * (1.0 - non_mt)
    elif non_mt.between(0.0, 100.0).all():
        obs["pct_counts_mt"] = 100.0 - non_mt
    else:
        raise RuntimeError("non_mt_percent must be expressed on a 0-1 or 0-100 scale")
    obs["sample"] = obs["protocol"]
    obs["condition"] = "controlled_mixture"
    obs["dataset"] = "public_mixology"
    obs["mixology_identity"] = obs["cell_line_demuxlet"]
    obs["cell_type"] = obs["cell_line_demuxlet"]
    obs["cell_subtype"] = obs["cell_line_demuxlet"]
    var = pd.DataFrame(index=pd.Index(genes["gene"].astype(str), name=None))
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.uns["source"] = {
        "name": "sc_mixology",
        "source_url": SOURCE_URL,
        "source_sha256": _sha256(source),
        "reference_doi": REFERENCE_DOI,
        "source_object": "post-sample-QC SingleCellExperiment objects",
        "gene_merge": "intersection across 10x, CEL-seq2, and Drop-seq",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")

    report = {
        "status": "READY",
        "output": str(output),
        "output_sha256": _sha256(output),
        "source": str(source),
        "source_url": SOURCE_URL,
        "source_sha256": _sha256(source),
        "reference_doi": REFERENCE_DOI,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "protocol_counts": obs["protocol"].value_counts().sort_index().to_dict(),
        "identity_counts": identity.value_counts().sort_index().to_dict(),
        "qc_truth_scope": "not a low-quality-cell truth set; source objects were post sample QC",
    }
    report_path = output.with_suffix(".provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
