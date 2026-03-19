#!/usr/bin/env python3
"""Quick conversion of mouse melanoma to h5ad."""

import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import io

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"


def convert():
    """Convert mouse melanoma GSE119352 to h5ad."""
    log.info("Converting mouse melanoma GSE119352...")

    raw_dir = DATA_DIR / "mouse_melanoma_GSE119352" / "GSE119352_RAW"

    samples = {
        "Control": "GSM3371684_Control",
        "aPD1": "GSM3371685_aPD1",
        "aCTLA4": "GSM3371686_aCTLA4",
        "aPD1-aCTLA4": "GSM3371687_aPD1-aCTLA4",
    }

    adata_list = []
    for sample_name, prefix in samples.items():
        log.info(f"  Loading {sample_name}...")

        matrix = io.mmread(raw_dir / f"{prefix}_matrix.mtx.gz")
        barcodes = pd.read_csv(
            raw_dir / f"{prefix}_barcodes.tsv.gz",
            header=None, sep="\t"
        )[0].values
        features = pd.read_csv(
            raw_dir / f"{prefix}_genes.tsv.gz",
            header=None, sep="\t"
        )

        adata = sc.AnnData(X=matrix.T.tocsr())
        adata.obs_names = [f"{b}-{sample_name}" for b in barcodes]
        adata.var_names = features[1].values

        adata.obs["sampleID"] = sample_name
        adata.obs["species"] = "mouse"
        adata.obs["treatment"] = sample_name

        adata.var_names_make_unique()
        adata_list.append(adata)

    adata = sc.concat(adata_list, join="outer", index_unique="-")
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    log.info(f"  Final shape: {adata.shape}")

    output_path = DATA_DIR / "mouse_melanoma.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


if __name__ == "__main__":
    convert()
