#!/usr/bin/env python3
"""
Convert raw test datasets to standardized h5ad format.

This script converts raw data (10x mtx, CSV, TXT) to h5ad format for faster loading.
"""

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io
from scipy.sparse import csr_matrix

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"


def convert_mouse_melanoma():
    """Convert mouse melanoma GSE119352 to h5ad."""
    log.info("Converting mouse melanoma GSE119352...")

    raw_dir = DATA_DIR / "mouse_melanoma_GSE119352" / "GSE119352_RAW"

    # Sample files
    samples = {
        "Control": "GSM3371684_Control",
        "aPD1": "GSM3371685_aPD1",
        "aCTLA4": "GSM3371686_aCTLA4",
        "aPD1-aCTLA4": "GSM3371687_aPD1-aCTLA4",
    }

    adata_list = []
    for sample_name, prefix in samples.items():
        log.info(f"  Loading {sample_name}...")

        # Read 10x files
        matrix = io.mmread(raw_dir / f"{prefix}_matrix.mtx.gz")
        barcodes = pd.read_csv(
            raw_dir / f"{prefix}_barcodes.tsv.gz",
            header=None, sep="\t"
        )[0].values
        features = pd.read_csv(
            raw_dir / f"{prefix}_genes.tsv.gz",
            header=None, sep="\t"
        )

        # Create AnnData
        adata = sc.AnnData(X=matrix.T.tocsr())
        adata.obs_names = [f"{b}-{sample_name}" for b in barcodes]
        adata.var_names = features[1].values  # gene symbols

        # Add metadata
        adata.obs["sampleID"] = sample_name
        adata.obs["species"] = "mouse"
        adata.obs["treatment"] = sample_name

        adata_list.append(adata)

    # Merge all samples
    adata = sc.concat(adata_list, join="outer")
    adata.var_names_make_unique()

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Samples: {adata.obs['sampleID'].value_counts().to_dict()}")

    # Save
    output_path = DATA_DIR / "mouse_melanoma.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


def convert_human_luad():
    """Convert human LUAD GSE131907 to h5ad."""
    log.info("Converting human LUAD GSE131907...")

    data_dir = DATA_DIR / "human_LUAD_GSE131907"

    # Read count matrix
    log.info("  Reading count matrix (this may take a while)...")
    count_file = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"

    with gzip.open(count_file, "rt") as f:
        # Read header to get cell names
        first_line = f.readline().strip()
        cell_names = first_line.split("\t")[1:]  # Skip first column (gene names)

    # Read full matrix
    count_df = pd.read_csv(count_file, sep="\t", index_col=0)

    log.info(f"  Count matrix shape: {count_df.shape}")

    # Create AnnData
    adata = sc.AnnData(X=count_df.values.T)
    adata.var_names = count_df.index.tolist()
    adata.obs_names = count_df.columns.tolist()

    # Read metadata
    meta_file = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)

    log.info(f"  Metadata shape: {meta_df.shape}")

    # Merge metadata
    common_cells = adata.obs_names.intersection(meta_df.index)
    adata = adata[common_cells].copy()

    for col in meta_df.columns:
        adata.obs[col] = meta_df.loc[common_cells, col]

    # Standardize column names
    if "Sample" in adata.obs:
        adata.obs["sampleID"] = adata.obs["Sample"]
    if "Cell_type" in adata.obs:
        adata.obs["cell_type"] = adata.obs["Cell_type"]

    adata.obs["species"] = "human"

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Cell types: {adata.obs.get('cell_type', 'N/A').value_counts().head().to_dict()}")

    # Save
    output_path = DATA_DIR / "human_lung_cancer.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


def create_pbmc4k_from_3k():
    """Create a 4k cell version from pbmc3k by adding synthetic variation."""
    log.info("Creating pbmc4k from pbmc3k...")

    pbmc3k_path = DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad"
    adata = sc.read_h5ad(pbmc3k_path)

    # Add sample info
    np.random.seed(42)
    samples = np.random.choice(["sample1", "sample2", "sample3", "sample4"], size=adata.n_obs)
    adata.obs["sampleID"] = samples
    adata.obs["species"] = "human"

    # Add basic QC metrics
    adata.obs["n_genes_by_counts"] = (adata.X > 0).sum(axis=1).A1
    adata.obs["total_counts"] = adata.X.sum(axis=1).A1

    # Store counts
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Shape: {adata.shape}")
    log.info(f"  Samples: {adata.obs['sampleID'].value_counts().to_dict()}")

    # Save
    output_path = DATA_DIR / "pbmc3k.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


def main():
    """Convert all test datasets."""
    log.info("=" * 60)
    log.info("Converting test datasets to h5ad format")
    log.info("=" * 60)

    converted = []

    # Convert pbmc3k with metadata
    try:
        path = create_pbmc4k_from_3k()
        converted.append(path)
    except Exception as e:
        log.error(f"Failed to convert pbmc3k: {e}")

    # Convert mouse melanoma
    try:
        path = convert_mouse_melanoma()
        converted.append(path)
    except Exception as e:
        log.error(f"Failed to convert mouse melanoma: {e}")

    # Convert human LUAD
    try:
        path = convert_human_luad()
        converted.append(path)
    except Exception as e:
        log.error(f"Failed to convert human LUAD: {e}")

    log.info("=" * 60)
    log.info("Conversion complete!")
    log.info(f"Converted files: {len(converted)}")
    for f in converted:
        size_mb = f.stat().st_size / (1024 * 1024)
        log.info(f"  {f.name}: {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
