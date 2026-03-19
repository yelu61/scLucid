#!/usr/bin/env python3
"""Convert mouse melanoma with full metadata (v2)."""

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
        # Create unique obs_names without sc.concat suffix
        adata.obs_names = [f"{b}-{sample_name}" for b in barcodes]
        adata.var_names = features[1].values

        adata.obs["sampleID"] = sample_name
        adata.obs["species"] = "mouse"
        adata.obs["treatment"] = sample_name

        adata.var_names_make_unique()
        adata_list.append(adata)

    # Concatenate without modifying indices
    adata = sc.concat(adata_list, join="outer")
    adata.layers["counts"] = adata.X.copy()

    # Add metadata
    log.info("  Adding metadata...")
    meta_file = DATA_DIR / "mouse_melanoma_GSE119352" / "GSE119352_scRNAseq_CD45_meta_data.tsv.gz"

    if meta_file.exists():
        meta_df = pd.read_csv(meta_file, sep="\t", index_col=0)
        log.info(f"    Metadata shape: {meta_df.shape}")

        # Create index matching adata format: barcode-sample
        meta_df.index = meta_df.index + "-" + meta_df["Sample"]

        # Find common cells
        common_cells = adata.obs_names.intersection(meta_df.index)
        log.info(f"    Common cells: {len(common_cells)} / {len(adata.obs_names)}")

        if len(common_cells) > 0:
            # Add metadata columns (only for cells with metadata)
            for col in meta_df.columns:
                if col in ['tSNE_1', 'tSNE_2']:
                    adata.obs[col] = float('nan')
                    adata.obs.loc[common_cells, col] = meta_df.loc[common_cells, col].astype(float)
                elif col == 'GraphCluster':
                    adata.obs[col] = -1  # Default for missing
                    adata.obs.loc[common_cells, col] = meta_df.loc[common_cells, col].astype(int)
                else:
                    adata.obs[col] = pd.NA
                    adata.obs.loc[common_cells, col] = meta_df.loc[common_cells, col]

            log.info(f"    Added metadata: {list(meta_df.columns)}")

    log.info(f"  Final shape: {adata.shape}")
    log.info(f"  Obs columns: {list(adata.obs.columns)}")
    log.info(f"  Sample distribution:\n{adata.obs['sampleID'].value_counts()}")

    output_path = DATA_DIR / "mouse_melanoma.h5ad"
    adata.write_h5ad(output_path)
    log.info(f"  Saved to {output_path}")

    return output_path


if __name__ == "__main__":
    convert()
