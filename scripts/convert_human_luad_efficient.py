#!/usr/bin/env python3
"""Convert human LUAD GSE131907 to h5ad using memory-efficient approach."""

import gzip
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import sparse

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"


def convert():
    """Convert using pandas read_csv with chunks."""
    log.info("Converting human LUAD GSE131907 (memory-efficient)...")

    data_dir = DATA_DIR / "human_LUAD_GSE131907"
    count_file = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    meta_file = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"

    # Read metadata first (smaller)
    log.info("  Reading metadata...")
    meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)
    log.info(f"    Metadata: {meta_df.shape}")

    # Read count matrix in chunks and transpose on the fly
    log.info("  Reading count matrix (this will take several minutes)...")

    # Read the full matrix - this requires substantial memory
    # For 208k cells x ~20k genes, we need ~8GB RAM for dense
    # Use sparse representation instead

    with gzip.open(count_file, "rt") as f:
        header = f.readline().strip().split("\t")
        cell_names = header[1:]
        n_cells = len(cell_names)
        log.info(f"    Cells: {n_cells}")

    # Read gene by gene and build sparse matrix
    from scipy.sparse import lil_matrix

    # First pass: determine dimensions
    log.info("  Counting genes...")
    with gzip.open(count_file, "rt") as f:
        next(f)  # Skip header
        gene_names = []
        for line in f:
            gene = line.split("\t")[0]
            gene_names.append(gene)
        n_genes = len(gene_names)

    log.info(f"    Genes: {n_genes}")
    log.info(f"    Matrix size: {n_cells} x {n_genes}")

    # Build sparse matrix row by row (gene by gene)
    log.info("  Building sparse matrix...")
    data = []
    row_idx = []
    col_idx = []

    with gzip.open(count_file, "rt") as f:
        next(f)  # Skip header
        for i, line in enumerate(f):
            if i % 1000 == 0:
                log.info(f"    Processed {i}/{n_genes} genes...")
            parts = line.strip().split("\t")
            values = [int(x) for x in parts[1:]]

            # Only store non-zero values
            for j, val in enumerate(values):
                if val > 0:
                    data.append(val)
                    row_idx.append(j)  # cell index
                    col_idx.append(i)  # gene index

    log.info(f"  Creating sparse matrix ({len(data)} non-zero values)...")
    X = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(n_cells, n_genes))

    log.info("  Creating AnnData...")
    adata = sc.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = cell_names

    # Add metadata
    log.info("  Adding metadata...")
    for col in meta_df.columns:
        adata.obs[col] = meta_df.loc[adata.obs_names, col]

    col_mapping = {"Sample": "sampleID", "Cell_type": "cell_type", "Cluster": "cluster"}
    for old_col, new_col in col_mapping.items():
        if old_col in adata.obs.columns:
            adata.obs[new_col] = adata.obs[old_col]

    adata.obs["species"] = "human"
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final: {adata.shape}")
    log.info(f"  Samples: {adata.obs['sampleID'].nunique()}")

    output_path = DATA_DIR / "human_lung_cancer.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved: {output_path}")

    return output_path


if __name__ == "__main__":
    convert()
