"""
Unified data loading utilities for testing.

This module provides convenient loading of test datasets with support for:
- 10x Genomics format (mtx, gz)
- GEO processed format (CSV, TSV)
- Custom formats

Available datasets (all in h5ad format for fast loading):
- pbmc3k: PBMC 3k (human, normal, 2,700 cells, 4 samples) - FAST (~0.1s)
- mouse_melanoma: Mouse melanoma GSE119352 (mouse, BRAF-mutant, 14,618 cells,
  4 treatment groups, with tSNE and cluster annotations) - FAST (~0.2s)
- human_lung_cancer: Human LUAD GSE131907 (human, lung cancer, 208,506 cells,
  58 samples, full cell type annotations) - MEDIUM (~10s due to large size)

Usage:
    >>> from tests.fixtures import load_test_data
    >>> adata = load_test_data("pbmc")
    >>> adata = load_test_data("nsclc", subsample=1000)
"""

import gzip
import logging
import tarfile
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"


def load_test_data(
    dataset: Literal[
        "pbmc3k",
        "mouse_melanoma",
        "human_lung_cancer",
    ] = "pbmc3k",
    subsample: Optional[int] = None,
    sample_key: str = "sampleID",
) -> sc.AnnData:
    """
    Load a test dataset for testing and examples.

    Args:
        dataset: Which dataset to load
        subsample: If provided, randomly subsample to this many cells
        sample_key: Key name for sample IDs in adata.obs

    Returns:
        AnnData object with the test data

    Examples:
        >>> adata = load_test_data("pbmc3k")  # Recommended
        >>> adata = load_test_data("human_cancer", subsample=1000)
        >>> adata = load_test_data("mouse_cancer")
    """
    # Data paths configuration - all h5ad format for fast loading
    data_paths = {
        # Standardized h5ad datasets (fast loading)
        "pbmc3k": DATA_DIR / "pbmc3k.h5ad",
        "mouse_melanoma": DATA_DIR / "mouse_melanoma.h5ad",
        "human_lung_cancer": DATA_DIR / "human_lung_cancer.h5ad",
    }

    if dataset not in data_paths:
        raise ValueError(f"Unknown dataset: {dataset}. " f"Choose from {list(data_paths.keys())}")

    data_path = data_paths[dataset]

    if not data_path.exists():
        available = [d for d in data_paths if data_paths[d].exists()]
        raise FileNotFoundError(
            f"Data path not found: {data_path}\n"
            f"Available datasets: {available}\n"
            f"Run 'python scripts/convert_test_data_to_h5ad.py' to generate h5ad files."
        )

    log.info(f"Loading {dataset} dataset from {data_path}")

    # All datasets are now in h5ad format for fast loading
    adata = _load_h5ad_format(data_path, sample_key)

    # Store raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Subsample if requested
    if subsample and adata.n_obs > subsample:
        indices = np.random.choice(adata.n_obs, subsample, replace=False)
        adata = adata[indices].copy()
        log.info(f"Subsampled to {adata.n_obs} cells")

    log.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def _load_h5ad_format(data_path: Path, sample_key: str) -> sc.AnnData:
    """Load h5ad format data (AnnData HDF5)."""
    adata = sc.read_h5ad(data_path)

    # Ensure sample_key exists
    if sample_key not in adata.obs:
        # Try to guess sample_key from common names
        for possible_key in ["sample", "batch", "donor", "patient", "mouse"]:
            if possible_key in adata.obs.columns:
                adata.obs[sample_key] = adata.obs[possible_key]
                break
        else:
            # No sample info found, create default
            adata.obs[sample_key] = "sample1"

    # Ensure species annotation exists
    if "species" not in adata.obs:
        # Try to infer from gene names
        if adata.var_names[0][0].isupper():  # Human genes start with uppercase
            adata.obs["species"] = "human"
        else:  # Mouse genes start with uppercase but some are mixed
            adata.obs["species"] = "mouse"

    return adata


def _load_10x_format(data_path: Path, sample_key: str) -> sc.AnnData:
    """Load 10x Genomics format data."""
    adata = sc.read_10x_mtx(data_path, var_names="gene_symbols", cache=True)

    # Try to extract sample IDs from barcodes if multiplexed
    barcodes = adata.obs_names.tolist()

    # Check for multiplexing pattern (e.g., #1-AAACCTG...)
    if any("#" in barcode for barcode in barcodes[:10]):
        # Extract sample from barcode
        samples = [
            barcode.split("#")[1].split("-")[0] if "#" in barcode else "sample1"
            for barcode in barcodes
        ]
        adata.obs[sample_key] = pd.Categorical(samples)
    else:
        adata.obs[sample_key] = "sample1"

    # Add species annotation
    adata.obs["species"] = "human"

    return adata


def _load_geo_format(
    data_path: Path,
    sample_key: str,
    dataset: str,
) -> sc.AnnData:
    """Load GEO format data (CSV/TSV count matrix)."""
    # Look for count file
    count_files = (
        list(data_path.glob("*.csv"))
        + list(data_path.glob("*.csv.gz"))
        + list(data_path.glob("*.txt"))
        + list(data_path.glob("*.txt.gz"))
    )

    if not count_files:
        raise FileNotFoundError(f"No count file found in {data_path}")

    # Use first count file
    count_file = count_files[0]
    log.info(f"Reading count file: {count_file.name}")

    # Determine separator and compression
    if count_file.suffix == ".gz":
        open_func = gzip.open
        mode = "rt"
        base_name = count_file.stem
    else:
        open_func = open
        mode = "r"
        base_name = count_file.name

    # Detect separator
    with open_func(count_file, mode) as f:
        first_line = f.readline()
        sep = "\t" if "\t" in first_line else ","

    # Read count matrix
    with open_func(count_file, mode) as f:
        count_df = pd.read_csv(f, sep=sep, index_col=0)

    # Create AnnData
    adata = sc.AnnData(X=count_df.values.T)
    adata.var_names = count_df.index.tolist()
    adata.obs_names = count_df.columns.tolist()

    # Try to load metadata
    metadata_files = [f for f in data_path.glob("*metadata*.csv") if f != count_file] + [
        f for f in data_path.glob("*meta*.csv") if f != count_file
    ]

    if metadata_files:
        meta_file = metadata_files[0]
        log.info(f"Reading metadata: {meta_file.name}")
        try:
            meta_df = pd.read_csv(meta_file, index_col=0)
            # Merge metadata
            common_cells = adata.obs_names.intersection(meta_df.index)
            adata = adata[common_cells].copy()

            for col in meta_df.columns:
                if col not in [sample_key, "species"]:
                    adata.obs[col] = meta_df.loc[common_cells, col]

            # Set sample key if available
            if sample_key in meta_df.columns:
                adata.obs[sample_key] = meta_df.loc[common_cells, sample_key]
        except Exception as e:
            log.warning(f"Could not load metadata: {e}")

    # Set sample key if not in metadata
    if sample_key not in adata.obs:
        adata.obs[sample_key] = "sample1"

    # Set species
    if dataset.startswith("mouse"):
        adata.obs["species"] = "mouse"
    else:
        adata.obs["species"] = "human"

    return adata


def _load_crc_format(data_path: Path, sample_key: str) -> sc.AnnData:
    """Load CRC dataset (legacy format)."""
    count_file = data_path / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"

    if not count_file.exists():
        raise FileNotFoundError(f"CRC count file not found: {count_file}")

    with gzip.open(count_file, "rt") as f:
        count_df = pd.read_csv(f, sep="\t", index_col=0, nrows=1000)

    adata = sc.AnnData(X=count_df.values.T)
    adata.var_names = count_df.index.tolist()
    adata.obs_names = count_df.columns.tolist()
    adata.obs[sample_key] = "CRC_sample"
    adata.obs["species"] = "human"

    return adata


def _load_gc_format(data_path: Path, sample_key: str) -> sc.AnnData:
    """Load GC dataset (legacy tar format)."""
    tar_file = data_path / "GSE134520_RAW.tar"

    if not tar_file.exists():
        raise FileNotFoundError(f"GC tar file not found: {tar_file}")

    with tarfile.open(tar_file, "r") as tar:
        first_member = tar.getmembers()[0]
        f = tar.extractfile(first_member)
        with gzip.open(f, "rt") as gz:
            count_df = pd.read_csv(gz, sep="\t", index_col=0, nrows=500)

    adata = sc.AnnData(X=count_df.values.T)
    adata.var_names = count_df.index.tolist()
    adata.obs_names = count_df.columns.tolist()
    adata.obs[sample_key] = "GC_sample"
    adata.obs["species"] = "human"

    return adata


def get_test_config(dataset: str = "pbmc3k") -> dict:
    """
    Get default test configuration for a dataset.

    Args:
        dataset: Dataset name

    Returns:
        Dictionary with default QC parameters
    """
    configs = {
        # New v2 datasets
        "pbmc3k": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 200,
            "max_mt_percent": 20,
            "min_counts": 1000,
        },
        "human_cancer": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 100,
            "max_mt_percent": 25,
            "min_counts": 500,
        },
        "mouse_cancer": {
            "sample_key": "sampleID",
            "species": "mouse",
            "min_genes": 100,
            "max_mt_percent": 20,
            "min_counts": 500,
        },
        # Older datasets
        "pbmc": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 200,
            "max_mt_percent": 20,
            "min_counts": 1000,
        },
        "pbmc_small": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 200,
            "max_mt_percent": 20,
            "min_counts": 1000,
        },
        "nsclc": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 100,
            "max_mt_percent": 25,
            "min_counts": 500,
        },
        "mouse_melanoma": {
            "sample_key": "sampleID",
            "species": "mouse",
            "min_genes": 100,
            "max_mt_percent": 20,
            "min_counts": 500,
        },
        "mouse_lung": {
            "sample_key": "sampleID",
            "species": "mouse",
            "min_genes": 100,
            "max_mt_percent": 20,
            "min_counts": 500,
        },
        "crc": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 100,
            "max_mt_percent": 25,
            "min_counts": 500,
        },
        "gc": {
            "sample_key": "sampleID",
            "species": "human",
            "min_genes": 100,
            "max_mt_percent": 25,
            "min_counts": 500,
        },
    }
    return configs.get(dataset, configs["pbmc3k"])


def list_available_datasets() -> pd.DataFrame:
    """
    List all available datasets with information.

    Returns:
        DataFrame with dataset information
    """
    datasets_info = []

    dataset_descriptions = {
        # New v2 datasets (recommended)
        "pbmc3k": ("PBMC 3k", "Human", "Normal", 1, "Scanpy built-in, 2,700 cells"),
        "human_cancer": ("Human Breast Cancer", "Human", "Tumor", 4, "Multi-sample, TME"),
        "mouse_cancer": ("Mouse Melanoma", "Mouse", "Tumor", 3, "Multi-sample, immunotherapy"),
        # Older datasets
        "pbmc": ("PBMC 4-Donor", "Human", "Normal", 4, "10x Genomics, 20k cells"),
        "pbmc_small": ("PBMC Small", "Human", "Normal", 1, "Quick tests"),
        "nsclc": ("NSCLC GSE119911", "Human", "Lung Cancer", 20, "GEO, multi-patient"),
        "mouse_melanoma": ("Mouse Melanoma GSE279468", "Mouse", "Melanoma", "NA", "BRAF-mutant"),
        "mouse_lung": ("Mouse Lung GSE222901", "Mouse", "Lung Cancer", "NA", "Gprc5a-/-"),
        "crc": ("CRC Legacy", "Human", "Colorectal", "NA", "Old format"),
        "gc": ("GC Legacy", "Human", "Gastric", "NA", "Old format"),
    }

    for name, (desc, species, tissue, samples, notes) in dataset_descriptions.items():
        data_path = (
            DATA_DIR
            / {
                # New v2 datasets
                "pbmc3k": "pbmc3k/pbmc3k_raw.h5ad",
                "human_cancer": "human_breast_cancer_broad/simulated_breast_cancer.h5ad",
                "mouse_cancer": "mouse_melanoma_model/simulated_mouse_melanoma.h5ad",
                # Older datasets
                "pbmc": "pbmc_4_donor/filtered_feature_bc_matrix",
                "pbmc_small": "10x_pbmc_small/filtered_feature_bc_matrix",
                "nsclc": "nsclc_gse119911",
                "mouse_melanoma": "mouse_melanoma_gse279468",
                "mouse_lung": "mouse_lung_gse222901",
                "crc": "CRC",
                "gc": "GC",
            }[name]
        )

        available = data_path.exists()

        datasets_info.append(
            {
                "name": name,
                "description": desc,
                "species": species,
                "tissue_type": tissue,
                "n_samples": samples,
                "notes": notes,
                "available": "✓" if available else "✗",
            }
        )

    return pd.DataFrame(datasets_info)
