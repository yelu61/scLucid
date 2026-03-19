#!/usr/bin/env python3
"""Convert human LUAD GSE131907 to h5ad with full metadata."""

import gzip
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"


def convert():
    """Convert human LUAD GSE131907 with cell annotation."""
    log.info("Converting human LUAD GSE131907...")

    data_dir = DATA_DIR / "human_LUAD_GSE131907"

    # Read count matrix
    count_file = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    log.info(f"  Reading count matrix from {count_file.name}...")
    log.info("  This is a large file, please wait...")

    # Read with gzip
    with gzip.open(count_file, "rt") as f:
        # Read header
        header = f.readline().strip().split("\t")
        cell_names = header[1:]  # First column is gene names
        log.info(f"    Found {len(cell_names)} cells")

    # Read full matrix
    count_df = pd.read_csv(count_file, sep="\t", index_col=0)
    log.info(f"    Count matrix shape: {count_df.shape}")

    # Create AnnData
    adata = sc.AnnData(X=count_df.values.T)
    adata.var_names = count_df.index.tolist()
    adata.obs_names = count_df.columns.tolist()

    # Read metadata
    meta_file = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    log.info(f"  Reading metadata from {meta_file.name}...")

    meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)
    log.info(f"    Metadata shape: {meta_df.shape}")
    log.info(f"    Metadata columns: {list(meta_df.columns)}")

    # Merge metadata
    common_cells = adata.obs_names.intersection(meta_df.index)
    log.info(f"    Common cells: {len(common_cells)} / {len(adata.obs_names)}")

    if len(common_cells) > 0:
        # Add all metadata columns
        for col in meta_df.columns:
            adata.obs[col] = meta_df.loc[adata.obs_names, col]

        log.info(f"    Added metadata: {list(meta_df.columns)}")

    # Standardize column names
    col_mapping = {
        "Sample": "sampleID",
        "Cell_type": "cell_type",
        "Cluster": "cluster",
    }
    for old_col, new_col in col_mapping.items():
        if old_col in adata.obs.columns:
            adata.obs[new_col] = adata.obs[old_col]

    adata.obs["species"] = "human"

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Obs columns: {list(adata.obs.columns)}")
    log.info(f"  Sample distribution:")
    if "sampleID" in adata.obs:
        print(adata.obs["sampleID"].value_counts().head(10))
    if "cell_type" in adata.obs:
        log.info(f"  Cell types:")
        print(adata.obs["cell_type"].value_counts().head(10))

    output_path = DATA_DIR / "human_lung_cancer.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


if __name__ == "__main__":
    convert()
