#!/usr/bin/env python3
"""Normalize legacy local h5ad fixtures to the scLucid dataset contract.

This script is intentionally conservative: it does not normalize expression,
filter cells, or infer cell types. It only adds missing contract fields needed
for QC/preprocess benchmark manifests and review summaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scipy.sparse as sp


DATASET_META: dict[str, dict[str, Any]] = {
    "pbmc3k": {
        "path": Path("data/pbmc3k.h5ad"),
        "dataset_key": "pbmc3k",
        "title": "10x Genomics PBMC 3k",
        "geo_accession": None,
        "pmid": None,
        "species": "human",
        "tissue": "PBMC",
        "condition": "normal",
        "sample_source": "sampleID",
        "citation": "10x Genomics PBMC 3k fixture distributed through the Scanpy ecosystem.",
    },
    "lin2020.pdac": {
        "path": Path("data/lin2020.pdac.h5ad"),
        "dataset_key": "lin2020.pdac",
        "title": "Lin et al. 2020 PDAC",
        "geo_accession": "GSE154778",
        "pmid": "32988401",
        "species": "human",
        "tissue": "PDAC",
        "condition": "tumor",
        "sample_source": "sampleID",
        "citation": "Lin W, Noel P, Borazanci EH, et al. Single-cell transcriptome analysis of tumor and stromal compartments of pancreatic ductal adenocarcinoma primary tumors and metastatic lesions. Genome Med. 2020;12:80.",
        "limitations": [
            "Development failure-control cohort; author cell labels are missing or unreliable in the local fixture."
        ],
    },
    "schlesinger2020.pdac": {
        "path": Path("data/schlesinger2020.pdac.h5ad"),
        "dataset_key": "schlesinger2020.pdac",
        "title": "Schlesinger et al. 2020 human PDAC sample",
        "geo_accession": "GSE141017",
        "pmid": "32908137",
        "species": "human",
        "tissue": "PDAC",
        "condition": "tumor",
        "sample_source": "sampleID",
        "citation": "Schlesinger Y, Yosefov-Levi O, Kolodkin-Gal D, et al. Single-cell transcriptomes of pancreatic preinvasive lesions and cancer reveal acinar metaplastic cells' heterogeneity. Nat Commun. 2020;11:4516.",
        "limitations": [
            "The local human fixture is the single sample GSM4293555 and cannot validate cross-sample borrowing."
        ],
    },
}


def _as_counts_matrix(X):
    if sp.issparse(X):
        return X.astype(np.int32)
    return np.asarray(X, dtype=np.int32)


def _looks_like_counts(X) -> bool:
    vals = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if vals.size == 0:
        return False
    sample = vals[: min(vals.size, 1_000_000)]
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def normalize_dataset(key: str, output: Path | None = None) -> Path:
    meta = DATASET_META[key]
    path = meta["path"]
    out_path = output or path
    if not path.exists():
        raise FileNotFoundError(path)

    adata = ad.read_h5ad(path)
    if "counts" not in adata.layers:
        if not _looks_like_counts(adata.X):
            raise ValueError(
                f"{path} does not have layers['counts'] and .X does not look like raw counts."
            )
        adata.layers["counts"] = _as_counts_matrix(adata.X).copy()
    else:
        adata.layers["counts"] = _as_counts_matrix(adata.layers["counts"]).copy()

    # Keep .X as counts for benchmark fixtures whose contract says .X is raw counts.
    adata.X = adata.layers["counts"].copy()

    sample_source = meta["sample_source"]
    if "sample" not in adata.obs:
        if sample_source in adata.obs:
            adata.obs["sample"] = adata.obs[sample_source].astype(str)
        elif "orig.ident" in adata.obs:
            adata.obs["sample"] = adata.obs["orig.ident"].astype(str)
        else:
            adata.obs["sample"] = meta["dataset_key"]
    if "sampleID" not in adata.obs:
        adata.obs["sampleID"] = adata.obs["sample"].astype(str)

    if "condition" not in adata.obs:
        if "group" in adata.obs:
            adata.obs["condition"] = adata.obs["group"].astype(str)
        else:
            adata.obs["condition"] = meta["condition"]
    if "cell_type" not in adata.obs:
        adata.obs["cell_type"] = ""
    if "cell_subtype" not in adata.obs:
        adata.obs["cell_subtype"] = ""
    adata.obs["dataset"] = meta["dataset_key"]

    if "gene_symbol" not in adata.var:
        adata.var["gene_symbol"] = adata.var_names.astype(str)
    adata.var_names_make_unique()

    adata.uns.setdefault("sclucid", {})["dataset"] = {
        "dataset_key": meta["dataset_key"],
        "title": meta["title"],
        "geo_accession": meta["geo_accession"],
        "pmid": meta["pmid"],
        "species": meta["species"],
        "tissue": meta["tissue"],
        "citation": meta["citation"],
        "processing": {
            "tool": "scripts/normalize_legacy_dataset_contract.py",
            "raw_counts_layer": "counts",
            "normalized_layer": None,
        },
        "limitations": meta.get("limitations", []),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        choices=[*DATASET_META.keys(), "all"],
        default="all",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory. If omitted, files are updated in place.",
    )
    args = parser.parse_args()

    keys = DATASET_META.keys() if args.dataset == "all" else [args.dataset]
    for key in keys:
        output = None
        if args.output_dir is not None:
            output = args.output_dir / DATASET_META[key]["path"].name
        out_path = normalize_dataset(key, output=output)
        print(f"{key}: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
