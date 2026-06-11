"""Data I/O example for scLucid.

This script demonstrates loading 10x-style data and stamping project metadata.
It is intentionally separate from ``01_workflow``: data loading is a setup step,
not the canonical QC/preprocess/analysis workflow.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

EXAMPLE_DATA_DIR = Path("data")
OUTPUT_DIR = Path("data/processed")


def load_single_10x_example() -> None:
    """Example for one Cell Ranger directory or .h5 file."""

    input_path = EXAMPLE_DATA_DIR / "single_sample" / "filtered_feature_bc_matrix"
    output_path = OUTPUT_DIR / "single_sample_raw.h5ad"

    adata = scl.ut.read_10x(
        input_path,
        sample_id="sample_1",
        species="human",
        tissue="blood",
        tissue_type="normal_tissue",
    )
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Single-sample AnnData saved to: {output_path}")


def load_multi_sample_10x_example() -> None:
    """Example for multiple 10x samples under one base directory."""

    samples = ["sample_1", "sample_2"]
    metadata_dicts = {
        "group": {
            "sample_1": "control",
            "sample_2": "treated",
        },
        "batch": {
            "sample_1": "batch_1",
            "sample_2": "batch_1",
        },
    }
    output_path = OUTPUT_DIR / "multi_sample_raw.h5ad"

    adata = scl.ut.read_10x(
        samples=samples,
        base_dir=EXAMPLE_DATA_DIR / "tenx_runs",
        metadata_dicts=metadata_dicts,
        species="human",
        tissue="blood",
        tissue_type="normal_tissue",
        output_file=output_path,
    )
    print(f"Multi-sample AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Saved to: {output_path}")


def prepare_existing_h5ad_example() -> None:
    """Example for an existing h5ad that needs scLucid handoff fields."""

    input_path = EXAMPLE_DATA_DIR / "pbmc3k.h5ad"
    output_path = OUTPUT_DIR / "pbmc3k_prepared.h5ad"

    adata = sc.read_h5ad(input_path)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "sampleID" not in adata.obs:
        adata.obs["sampleID"] = "pbmc3k"

    adata.uns.setdefault("sclucid", {}).setdefault("analysis_context", {}).update(
        {
            "species": "human",
            "tissue": "blood",
            "tissue_type": "normal_tissue",
            "dataset_type": "pbmc_or_blood",
        }
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Prepared h5ad saved to: {output_path}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Data I/O examples")
    print("- load_single_10x_example(): one Cell Ranger directory or .h5")
    print("- load_multi_sample_10x_example(): multiple sample folders with metadata")
    print("- prepare_existing_h5ad_example(): add counts/sample/context fields")
    print("\nRunning the existing-h5ad example by default.")
    prepare_existing_h5ad_example()


if __name__ == "__main__":
    main()
