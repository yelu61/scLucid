#!/usr/bin/env python3
"""Export a Bioconductor scDblFinder reference for the Kang doublet benchmark."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


def _counts_matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _stratified_subset(
    adata: ad.AnnData,
    labels: pd.Series,
    max_cells: int | None,
    seed: int,
) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep: list[int] = []
    label_counts = labels.astype(str).value_counts()
    for label, count in label_counts.items():
        idx = np.flatnonzero(labels.astype(str).to_numpy() == label)
        n_label = max(1, int(round(max_cells * count / adata.n_obs)))
        n_label = min(n_label, len(idx))
        keep.extend(rng.choice(idx, size=n_label, replace=False).tolist())
    if len(keep) > max_cells:
        keep = rng.choice(np.array(keep), size=max_cells, replace=False).tolist()
    keep = sorted(set(int(i) for i in keep))
    return adata[keep].copy()


def _write_r_script(path: Path) -> None:
    path.write_text(r"""
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

args <- commandArgs(trailingOnly=TRUE)
matrix_path <- args[[1]]
genes_path <- args[[2]]
cells_path <- args[[3]]
metadata_path <- args[[4]]
output_path <- args[[5]]
expected_rate <- as.numeric(args[[6]])
seed <- as.integer(args[[7]])

set.seed(seed)
counts <- readMM(matrix_path)
genes <- readLines(genes_path)
cells <- readLines(cells_path)
metadata <- read.csv(metadata_path, row.names=1, check.names=FALSE)
rownames(counts) <- genes
colnames(counts) <- cells

sce <- SingleCellExperiment(list(counts=counts), colData=metadata)
sce <- scDblFinder(sce, samples="sample", dbr=expected_rate, verbose=FALSE)

out <- data.frame(
  cell=colnames(sce),
  score=colData(sce)$scDblFinder.score,
  class=as.character(colData(sce)$scDblFinder.class),
  predicted=as.character(colData(sce)$scDblFinder.class) == "doublet",
  stringsAsFactors=FALSE,
  check.names=FALSE
)
write.csv(out, output_path, row.names=FALSE)
""".lstrip())


def export_reference(
    output: Path,
    max_cells: int | None,
    seed: int,
    work_dir: Path,
) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    spec = next(dataset for dataset in DATASETS if dataset.key == "kang2018.pbmc")
    adata = ad.read_h5ad(spec.path)
    labels = adata.obs["demuxlet_multiplets"].astype(str)
    adata = _stratified_subset(adata, labels, max_cells=max_cells, seed=seed)

    usable = adata.obs["demuxlet_multiplets"].astype(str).isin(["singlet", "doublet"])
    expected_rate = float(adata.obs.loc[usable, "doublet_ground_truth"].astype(bool).mean())

    X = _counts_matrix(adata)
    X = X.T.tocoo() if sp.issparse(X) else sp.coo_matrix(np.asarray(X).T)
    matrix_path = work_dir / "kang_counts_genes_by_cells.mtx"
    genes_path = work_dir / "kang_genes.txt"
    cells_path = work_dir / "kang_cells.txt"
    metadata_path = work_dir / "kang_metadata.csv"
    script_path = work_dir / "run_scdblfinder_reference.R"

    scipy.io.mmwrite(matrix_path, X)
    genes_path.write_text("\n".join(map(str, adata.var_names)))
    cells_path.write_text("\n".join(map(str, adata.obs_names)))
    pd.DataFrame(
        {
            "sample": adata.obs["sample"].astype(str).to_numpy(),
            "demuxlet_multiplets": adata.obs["demuxlet_multiplets"].astype(str).to_numpy(),
            "doublet_ground_truth": adata.obs["doublet_ground_truth"].astype(bool).to_numpy(),
        },
        index=adata.obs_names,
    ).to_csv(metadata_path)
    _write_r_script(script_path)

    subprocess.run(
        [
            "Rscript",
            str(script_path),
            str(matrix_path),
            str(genes_path),
            str(cells_path),
            str(metadata_path),
            str(output),
            str(expected_rate),
            str(seed),
        ],
        check=True,
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation_outputs/work/qc_doublet_evidence/r_scdblfinder_reference.csv"),
    )
    parser.add_argument("--max-cells", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("validation_outputs/work/qc_doublet_evidence/r_reference_inputs"),
    )
    args = parser.parse_args()
    path = export_reference(
        args.output, max_cells=args.max_cells, seed=args.seed, work_dir=args.work_dir
    )
    print(f"r_scdblfinder_reference: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
