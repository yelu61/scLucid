#!/usr/bin/env python3
"""Convert human LUAD GSE131907 to h5ad with chunked reading for memory efficiency."""

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"


def convert_chunked(chunk_size=10000):
    """Convert human LUAD GSE131907 with chunked reading."""
    log.info("Converting human LUAD GSE131907 (chunked mode)...")

    data_dir = DATA_DIR / "human_LUAD_GSE131907"
    count_file = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    meta_file = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"

    # First pass: count lines and get dimensions
    log.info("  Analyzing file...")
    with gzip.open(count_file, "rt") as f:
        header = f.readline().strip().split("\t")
        cell_names = header[1:]
        n_cells = len(cell_names)
        n_genes = sum(1 for _ in f)

    log.info(f"    Cells: {n_cells}, Genes: {n_genes}")

    # Read metadata first
    log.info("  Reading metadata...")
    meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)
    log.info(f"    Metadata shape: {meta_df.shape}")

    # Read gene names
    with gzip.open(count_file, "rt") as f:
        f.readline()  # Skip header
        gene_names = []
        for line in f:
            gene = line.split("\t")[0]
            gene_names.append(gene)

    # Initialize matrix in chunks
    log.info("  Reading count matrix in chunks...")
    chunks = []

    with gzip.open(count_file, "rt") as f:
        f.readline()  # Skip header

        for chunk_idx in range(0, n_genes, chunk_size):
            chunk_data = []
            lines_read = 0

            for line in f:
                parts = line.strip().split("\t")
                values = [int(x) for x in parts[1:]]
                chunk_data.append(values)
                lines_read += 1

                if lines_read >= chunk_size:
                    break

            if chunk_data:
                chunk_array = np.array(chunk_data)
                chunks.append(chunk_array)
                log.info(f"    Read chunk {chunk_idx // chunk_size + 1}: {chunk_array.shape}")

    # Combine chunks
    log.info("  Combining chunks...")
    full_matrix = np.vstack(chunks)

    # Create AnnData
    log.info("  Creating AnnData...")
    adata = sc.AnnData(X=sparse.csr_matrix(full_matrix.T))
    adata.var_names = gene_names
    adata.obs_names = cell_names

    # Add metadata
    log.info("  Adding metadata...")
    common_cells = adata.obs_names.intersection(meta_df.index)
    log.info(f"    Common cells: {len(common_cells)} / {len(adata.obs_names)}")

    for col in meta_df.columns:
        adata.obs[col] = meta_df.loc[adata.obs_names, col]

    # Standardize
    col_mapping = {
        "Sample": "sampleID",
        "Cell_type": "cell_type",
        "Cluster": "cluster",
    }
    for old_col, new_col in col_mapping.items():
        if old_col in adata.obs.columns:
            adata.obs[new_col] = adata.obs[old_col]

    adata.obs["species"] = "human"
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Sample count: {adata.obs['sampleID'].nunique()}")
    log.info(f"  Cell type count: {adata.obs['cell_type'].nunique()}")

    output_path = DATA_DIR / "human_lung_cancer.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


if __name__ == "__main__":
    convert_chunked(chunk_size=5000)
