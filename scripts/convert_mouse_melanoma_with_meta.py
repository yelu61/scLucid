#!/usr/bin/env python3
"""Convert mouse melanoma with full metadata."""

import gzip
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import io

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"


def convert():
    """Convert mouse melanoma GSE119352 with full metadata."""
    log.info("Converting mouse melanoma GSE119352 with metadata...")

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
    adata.layers["counts"] = adata.X.copy()

    # Add metadata
    log.info("  Adding metadata...")
    meta_file = DATA_DIR / "mouse_melanoma_GSE119352" / "GSE119352_scRNAseq_CD45_meta_data.tsv.gz"

    if meta_file.exists():
        meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)
        log.info(f"    Metadata shape: {meta_df.shape}")
        log.info(f"    Metadata columns: {list(meta_df.columns)}")

        # Create matching index (adata has format: barcode-sample-index)
        # Need to match: AAACATACAGCGTT-1-Control vs AAACATACAGCGTT-1-Control-0
        meta_df.index = meta_df.index + "-" + meta_df["Sample"]

        # Create mapping from adata index to meta index
        # Remove the trailing -N from adata obs_names for matching
        adata_obs_base = adata.obs_names.str.replace(r'-\\d+$', '', regex=True)

        # Find common cells
        common_mask = adata_obs_base.isin(meta_df.index)
        common_cells = adata.obs_names[common_mask]
        log.info(f"    Common cells: {len(common_cells)} / {len(adata.obs_names)}")
        log.info(f"    Common cells: {len(common_cells)}")

        if len(common_cells) > 0:
            # Add metadata columns
            for col in meta_df.columns:
                if col not in adata.obs:
                    adata.obs[col] = pd.NA

            # Map metadata using base index
            base_to_meta = {idx: idx for idx in meta_df.index}
            for cell in common_cells:
                base_idx = cell.rsplit('-', 1)[0]  # Remove trailing index
                if base_idx in meta_df.index:
                    for col in meta_df.columns:
                        adata.obs.at[cell, col] = meta_df.at[base_idx, col]

            log.info(f"    Added metadata: {list(meta_df.columns)}")

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Obs columns: {list(adata.obs.columns)}")

    output_path = DATA_DIR / "mouse_melanoma.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


if __name__ == "__main__":
    convert()
