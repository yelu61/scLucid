#!/usr/bin/env python3
"""Prepare downloaded public QC truth without redefining its evidence class."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

DEFAULT_ROOT = Path("data/external/qc_truth")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_manifests(
    dataset_dir: Path,
    *,
    dataset_id: str,
    accessions: list[str],
    source_urls: list[str],
    available_fields: list[str],
    artifacts: dict[str, Path],
) -> None:
    checksum_lines = []
    file_records = []
    for checksum_key, path in artifacts.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        checksum = _sha256(path)
        checksum_lines.append(f"{checksum}  {checksum_key}")
        file_records.append(
            {"path": str(path), "sha256": checksum, "size_bytes": path.stat().st_size}
        )
    (dataset_dir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n")
    (dataset_dir / "metadata_manifest.json").write_text(
        json.dumps({"available_fields": sorted(available_fields)}, indent=2) + "\n"
    )
    (dataset_dir / "provenance.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset_id,
                "accessions": accessions,
                "source_urls": source_urls,
                "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "files": file_records,
            },
            indent=2,
        )
        + "\n"
    )


def prepare_kang(root: Path, kang_h5ad: Path) -> dict[str, Any]:
    dataset_dir = root / "kang2018_pbmc"
    batch1 = pd.read_csv(
        dataset_dir / "labels/GSE96583_batch1.total.tsne.df.tsv.gz",
        sep="\t",
        index_col=0,
    )
    batch2 = pd.read_csv(
        dataset_dir / "labels/GSE96583_batch2.total.tsne.df.tsv.gz",
        sep="\t",
        index_col=0,
    )
    first = pd.DataFrame(
        {
            "barcode": batch1.index.astype(str),
            "capture": batch1["batch"].astype(str).to_numpy(),
            "condition": "ctrl",
            "donor": batch1["ind"].astype("string").to_numpy(),
            "cell_type": batch1["cell.type"].astype("string").to_numpy(),
            "demuxlet_class": batch1["multiplets"].astype(str).to_numpy(),
        }
    )
    second = pd.DataFrame(
        {
            "barcode": batch2.index.astype(str),
            "capture": np.where(batch2["stim"].astype(str).eq("ctrl"), "2.1", "2.2"),
            "condition": batch2["stim"].astype(str).to_numpy(),
            "donor": batch2["ind"].astype("string").to_numpy(),
            "cell_type": batch2["cell"].astype("string").to_numpy(),
            "demuxlet_class": batch2["multiplets"].astype(str).to_numpy(),
        }
    )
    labels = pd.concat([first, second], ignore_index=True)
    labels["truth_status"] = labels["demuxlet_class"].map(
        {"singlet": "KEEP", "doublet": "REMOVE", "ambs": "UNCERTAIN"}
    )
    labels["doublet_ground_truth"] = labels["demuxlet_class"].map(
        {"singlet": False, "doublet": True}
    )
    if labels.duplicated(["capture", "barcode"]).any():
        raise ValueError("Kang labels contain duplicate capture/barcode keys")
    output = dataset_dir / "labels/demuxlet_labels.tsv"
    labels.to_csv(output, sep="\t", index=False)

    local = ad.read_h5ad(kang_h5ad, backed="r")
    local_counts = local.obs["demuxlet_multiplets"].astype(str).value_counts().to_dict()
    local.file.close()
    _write_manifests(
        dataset_dir,
        dataset_id="kang2018_pbmc",
        accessions=["GSE96583"],
        source_urls=["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583"],
        available_fields=[
            "donor",
            "condition",
            "capture",
            "demuxlet_multiplets",
            "cell_type",
            "paired_design",
        ],
        artifacts={
            "labels/demuxlet_labels.tsv": output,
            f"registry_local/{kang_h5ad.name}": kang_h5ad,
        },
    )
    return {
        "dataset_id": "kang2018_pbmc",
        "labels": int(len(labels)),
        "label_counts": labels["demuxlet_class"].value_counts().to_dict(),
        "local_h5ad_counts": local_counts,
    }


def prepare_emtab2600(root: Path) -> dict[str, Any]:
    dataset_dir = root / "emtab2600_microscopy_quality"
    workbook = dataset_dir / "labels/Ilicic2016_TableS1_quality_annotations.xlsx"
    raw_labels = pd.read_excel(workbook, sheet_name="trainingsset mES")
    parsed = (
        raw_labels["cells"]
        .astype(str)
        .str.extract(r"trainings_set_mES_(2i|a2i|lif)_(\d+)_(\d+)\.counts")
    )
    if parsed.isna().any().any():
        raise ValueError("not every microscopy label matched the registered capture pattern")
    capture_labels = pd.DataFrame(
        {
            "capture_site": raw_labels["cells"]
            .astype(str)
            .str.replace(r"\.counts$", "", regex=True),
            "condition": parsed[0],
            "replicate": parsed[1].astype(int),
            "site": parsed[2].astype(int),
            "microscopy_quality_code": raw_labels["Extra"].astype(str),
        }
    )
    capture_labels["library"] = (
        capture_labels["condition"] + "_s" + capture_labels["replicate"].astype(str)
    )
    capture_labels["truth_cell_class"] = capture_labels["microscopy_quality_code"].map(
        {"G": "INTACT", "D": "DAMAGED", "E": "EMPTY", "M": "UNCERTAIN"}
    )
    if capture_labels["truth_cell_class"].isna().any():
        raise ValueError("unsupported microscopy quality code")
    capture_label_output = dataset_dir / "labels/microscopy_capture_labels.tsv"
    capture_labels.to_csv(capture_label_output, sep="\t", index=False)

    manifest = pd.read_csv(dataset_dir / "source/ena_file_manifest.tsv", sep="\t")
    manifest["capture_site"] = (
        "trainings_set_mES_" + manifest["sample_alias"].astype(str).str.extract(r":([^:]+)$")[0]
    )
    mapping = manifest.merge(
        capture_labels[["capture_site", "library", "microscopy_quality_code", "truth_cell_class"]],
        on="capture_site",
        how="left",
    )
    mapping["cell_barcode_mapping"] = mapping["run_accession"]
    mapping_output = dataset_dir / "labels/capture_site_barcode_mapping.tsv"
    mapping.to_csv(mapping_output, sep="\t", index=False)

    matrix_path = dataset_dir / "source/mESC960_input.matrix"
    counts = pd.read_csv(matrix_path, sep="\t", index_col=0)
    zeta_labels_path = dataset_dir / "source/mESC960_cellType_anno"
    zeta_labels = pd.read_csv(zeta_labels_path, sep="\t").set_index("Cell")["type"]
    if not zeta_labels.index.isin(counts.index).all():
        raise ValueError("ZetaSuite labels contain run accessions absent from the matrix")
    obs = pd.DataFrame(index=counts.index.astype(str))
    obs.index.name = "run_accession"
    obs["microscopy_quality_label"] = zeta_labels.reindex(obs.index).fillna("UNCERTAIN")
    obs["intact_or_damaged_or_empty"] = (
        obs["microscopy_quality_label"]
        .map({"high-quality": "INTACT", "Broken": "DAMAGED", "Empty": "EMPTY"})
        .fillna("UNCERTAIN")
    )
    obs["truth_status"] = (
        obs["intact_or_damaged_or_empty"]
        .map({"INTACT": "KEEP", "DAMAGED": "REMOVE", "EMPTY": "REMOVE"})
        .fillna("UNCERTAIN")
    )
    manifest_obs = manifest.set_index("run_accession").reindex(obs.index)
    obs["capture_site"] = manifest_obs["capture_site"].fillna("UNKNOWN").astype(str)
    obs["library"] = (
        manifest_obs["sample_alias"]
        .fillna("UNKNOWN")
        .astype(str)
        .str.extract(r":([^:]+?)(?:_\d+)?$")[0]
        .fillna("UNKNOWN")
    )
    obs["cell_barcode_mapping"] = obs.index
    prepared_counts = sparse.csr_matrix(counts.to_numpy(dtype=np.int32))
    var = pd.DataFrame(index=pd.Index(counts.columns.astype(str), name="ensembl_gene_id"))
    adata = ad.AnnData(X=prepared_counts, obs=obs, var=var)
    adata.layers["counts"] = prepared_counts.copy()
    adata.uns["truth_contract"] = {
        "type": "microscopy_labelled_legacy_capture_truth",
        "primary_source": "Ilicic et al. 2016, E-MTAB-2600/PRJEB6455",
        "processed_matrix_source": "ZetaSuite Zenodo 6395174",
        "uncertain_excluded_from_primary_binary_endpoint": True,
        "limitations": [
            "legacy Fluidigm C1/Smart-seq technology",
            "processed reconstruction covers 168 runs and 147 explicit labels",
            "microscopy phenotype is not molecular ground truth",
        ],
    }
    prepared = dataset_dir / "prepared/emtab2600_qc_truth.h5ad"
    adata.write_h5ad(prepared, compression="gzip")
    run_label_output = dataset_dir / "labels/microscopy_quality_labels.tsv"
    obs.reset_index().to_csv(run_label_output, sep="\t", index=False)
    _write_manifests(
        dataset_dir,
        dataset_id="emtab2600_microscopy_quality",
        accessions=["E-MTAB-2600", "PRJEB6455", "10.5281/zenodo.6395174"],
        source_urls=[
            "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600",
            "https://www.ebi.ac.uk/ena/browser/view/PRJEB6455",
            "https://doi.org/10.5281/zenodo.6395174",
            "https://doi.org/10.1186/s13059-016-0888-1",
        ],
        available_fields=[
            "capture_site",
            "microscopy_quality_label",
            "intact_or_damaged_or_empty",
            "library",
            "cell_barcode_mapping",
        ],
        artifacts={
            "source/ena_file_manifest.tsv": dataset_dir / "source/ena_file_manifest.tsv",
            "source/mESC960_input.matrix": matrix_path,
            "source/mESC960_cellType_anno": zeta_labels_path,
            "labels/Ilicic2016_TableS1_quality_annotations.xlsx": workbook,
            "labels/microscopy_capture_labels.tsv": capture_label_output,
            "labels/microscopy_quality_labels.tsv": run_label_output,
            "labels/capture_site_barcode_mapping.tsv": mapping_output,
            "prepared/emtab2600_qc_truth.h5ad": prepared,
        },
    )
    return {
        "dataset_id": "emtab2600_microscopy_quality",
        "source_capture_labels": int(len(capture_labels)),
        "ena_runs": int(len(mapping)),
        "ena_runs_with_source_capture_label": int(mapping["truth_cell_class"].notna().sum()),
        "prepared_runs": int(adata.n_obs),
        "prepared_label_counts": obs["intact_or_damaged_or_empty"].value_counts().to_dict(),
    }


def _read_lines(path: Path) -> np.ndarray:
    return np.asarray(path.read_text().splitlines(), dtype=object)


def prepare_hgmm(root: Path) -> dict[str, Any]:
    dataset_dir = root / "tenx_hgmm_6k"
    raw_dir = dataset_dir / "source/raw_gene_bc_matrices"
    human_barcodes = _read_lines(raw_dir / "hg19/barcodes.tsv")
    mouse_barcodes = _read_lines(raw_dir / "mm10/barcodes.tsv")
    if not np.array_equal(human_barcodes, mouse_barcodes):
        raise ValueError("HGMM raw human and mouse barcode universes differ")
    human = mmread(raw_dir / "hg19/matrix.mtx").tocsc()
    mouse = mmread(raw_dir / "mm10/matrix.mtx").tocsc()
    human_umi = np.asarray(human.sum(axis=0)).ravel().astype(np.int64)
    mouse_umi = np.asarray(mouse.sum(axis=0)).ravel().astype(np.int64)
    total_umi = human_umi + mouse_umi
    human_fraction = np.divide(
        human_umi,
        total_umi,
        out=np.zeros_like(human_umi, dtype=float),
        where=total_umi > 0,
    )
    filtered_dir = dataset_dir / "source/filtered_gene_bc_matrices"
    filtered_human = set(_read_lines(filtered_dir / "hg19/barcodes.tsv"))
    filtered_mouse = set(_read_lines(filtered_dir / "mm10/barcodes.tsv"))
    vendor_called = np.fromiter(
        (barcode in filtered_human or barcode in filtered_mouse for barcode in human_barcodes),
        dtype=bool,
        count=len(human_barcodes),
    )
    species_class = np.full(len(human_barcodes), "VENDOR_UNCALLED", dtype=object)
    called_nonzero = vendor_called & (total_umi > 0)
    species_class[called_nonzero & (human_fraction >= 0.9)] = "HUMAN_SINGLET"
    species_class[called_nonzero & (human_fraction <= 0.1)] = "MOUSE_SINGLET"
    species_class[called_nonzero & (human_fraction > 0.1) & (human_fraction < 0.9)] = (
        "CROSS_SPECIES_DOUBLET"
    )
    table = pd.DataFrame(
        {
            "barcode": human_barcodes,
            "library": "hgmm_6k",
            "human_umi": human_umi,
            "mouse_umi": mouse_umi,
            "total_umi": total_umi,
            "human_fraction": human_fraction,
            "vendor_filtered_cell_call": vendor_called,
            "species_class": species_class,
        }
    )
    table_output = dataset_dir / "labels/species_umi_counts.tsv"
    table.to_csv(table_output, sep="\t", index=False)

    aggregate_counts = sparse.csr_matrix(np.column_stack([human_umi, mouse_umi]))
    obs = table.set_index("barcode")
    var = pd.DataFrame(index=pd.Index(["human_umi", "mouse_umi"], name="feature"))
    adata = ad.AnnData(X=aggregate_counts, obs=obs, var=var)
    adata.layers["counts"] = aggregate_counts.copy()
    adata.uns["truth_contract"] = {
        "type": "species_aggregate_truth_object",
        "gene_level_expression_available": False,
        "vendor_call_is_comparator_not_independent_truth": True,
        "species_rule": "human_fraction >=0.9, <=0.1, otherwise cross-species doublet",
    }
    prepared = dataset_dir / "prepared/tenx_hgmm_qc_truth.h5ad"
    adata.write_h5ad(prepared, compression="gzip")
    _write_manifests(
        dataset_dir,
        dataset_id="tenx_hgmm_6k",
        accessions=["6k_HGMM_3p_v2"],
        source_urls=[
            "https://www.10xgenomics.com/datasets/6-k-1-1-mixture-of-fresh-frozen-human-hek-293-t-and-mouse-nih-3-t-3-cells-2-standard-1-2-0"
        ],
        available_fields=[
            "barcode",
            "human_umi",
            "mouse_umi",
            "filtered_cell_call",
            "raw_droplet_matrix",
            "library",
        ],
        artifacts={
            "source/hgmm_6k_raw_gene_bc_matrices.tar.gz": dataset_dir
            / "source/hgmm_6k_raw_gene_bc_matrices.tar.gz",
            "source/hgmm_6k_filtered_gene_bc_matrices.tar.gz": dataset_dir
            / "source/hgmm_6k_filtered_gene_bc_matrices.tar.gz",
            "labels/species_umi_counts.tsv": table_output,
            "prepared/tenx_hgmm_qc_truth.h5ad": prepared,
        },
    )
    return {
        "dataset_id": "tenx_hgmm_6k",
        "droplets": int(len(table)),
        "vendor_called": int(vendor_called.sum()),
        "species_counts": table.loc[vendor_called, "species_class"].value_counts().to_dict(),
    }


def prepare_cell_hashing(root: Path) -> dict[str, Any]:
    dataset_dir = root / "cell_hashing_gse108313"
    assignments = pd.read_csv(dataset_dir / "labels/hto_assignments.tsv", sep="\t")
    assignments["capture"] = "MixCellLines_GSM3501447"
    assignments["sample_hash"] = assignments["hash_id"]
    assignments["cell_line_or_sample_identity"] = pd.NA
    singlet = assignments["hto_classification_global"].eq("Singlet")
    assignments.loc[singlet, "cell_line_or_sample_identity"] = (
        assignments.loc[singlet, "hto_max_id"].astype(str).str.split("-").str[0]
    )
    rna_path = dataset_dir / "source/GSM3501446_MixCellLines-RNA.umi.txt.gz"
    rna_barcodes = pd.read_csv(rna_path, sep="\t", nrows=0, index_col=0).columns.astype(str)
    hto_path = dataset_dir / "source/GSM3501447_MixCellLines-HTO-counts.csv.gz"
    hto = pd.read_csv(hto_path, index_col=0)
    mapping = pd.DataFrame({"barcode": rna_barcodes})
    mapping["rna_matrix_present"] = True
    mapping["hto_matrix_present"] = mapping["barcode"].isin(hto.index)
    mapping = mapping.merge(assignments, on="barcode", how="left")
    mapping["hto_classification_global"] = mapping["hto_classification_global"].fillna(
        "UNKNOWN_NO_HTO"
    )
    mapping_output = dataset_dir / "labels/rna_hto_barcode_mapping.tsv"
    mapping.to_csv(mapping_output, sep="\t", index=False)

    joint = assignments["barcode"].astype(str)
    hto_counts = hto.loc[joint, hto.columns[:12]].to_numpy(dtype=np.int32)
    obs = assignments.set_index("barcode")
    var = pd.DataFrame(index=pd.Index(hto.columns[:12].astype(str), name="hto"))
    adata = ad.AnnData(X=sparse.csr_matrix(hto_counts), obs=obs, var=var)
    adata.layers["hto_counts"] = adata.X.copy()
    adata.uns["truth_contract"] = {
        "type": "orthogonal_hto_demultiplexing",
        "method": "Seurat HTODemux positive.quantile=0.99 seed=42",
        "label_class": "method_derived_orthogonal_evidence_not_gold_truth",
        "rna_expression_stored": False,
        "supports": "doublet and sample-hash calibration on overlapping RNA/HTO barcodes",
        "does_not_support": "complete unfiltered RNA droplet cell-calling performance",
    }
    prepared = dataset_dir / "prepared/cell_hashing_qc_truth.h5ad"
    adata.write_h5ad(prepared, compression="gzip")
    _write_manifests(
        dataset_dir,
        dataset_id="cell_hashing_gse108313",
        accessions=["GSE108313"],
        source_urls=["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313"],
        available_fields=[
            "capture",
            "hto_classification",
            "sample_hash",
            "rna_barcode",
            "cell_line_or_sample_identity",
        ],
        artifacts={
            "source/GSM3501446_MixCellLines-RNA.umi.txt.gz": rna_path,
            "source/GSM3501447_MixCellLines-HTO-counts.csv.gz": hto_path,
            "labels/hto_assignments.tsv": dataset_dir / "labels/hto_assignments.tsv",
            "labels/rna_hto_barcode_mapping.tsv": mapping_output,
            "prepared/cell_hashing_qc_truth.h5ad": prepared,
        },
    )
    return {
        "dataset_id": "cell_hashing_gse108313",
        "rna_barcodes": int(len(rna_barcodes)),
        "rna_hto_overlap": int(len(assignments)),
        "hto_classification_counts": assignments["hto_classification_global"]
        .value_counts()
        .to_dict(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=["kang", "emtab2600", "hgmm", "cell_hashing"])
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--kang-h5ad", type=Path, default=Path("data/kang2018.pbmc.h5ad"))
    args = parser.parse_args()
    if args.dataset == "kang":
        result = prepare_kang(args.root, args.kang_h5ad)
    elif args.dataset == "emtab2600":
        result = prepare_emtab2600(args.root)
    elif args.dataset == "hgmm":
        result = prepare_hgmm(args.root)
    else:
        result = prepare_cell_hashing(args.root)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
