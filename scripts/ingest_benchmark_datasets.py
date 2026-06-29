#!/usr/bin/env python3
"""Ingest downloaded GEO datasets into scLucid's benchmark h5ad format.

Supported datasets
------------------
- nsclc  : GSE127465 (Zilionis et al., 2019, Immunity)
- pancreas : GSE84133 human donors (Baron et al., 2016, Cell Systems)
- crc    : GSE132465 (Lee et al., 2020, Nature Genetics)

Output files are written to ``--output-dir`` and follow the naming convention
``<first_author><year>.<tissue>.h5ad``.

All outputs use a unified field contract:
- ``.X`` and ``.layers['counts']`` : raw counts (sparse int32)
- ``.obs['sample']``              : original sample/library id
- ``.obs['patient']`` / ``.obs['donor']`` : patient or donor id
- ``.obs['condition']``           : tumor / blood / normal / tumor-adjacent etc.
- ``.obs['cell_type']``           : author-provided major cell type
- ``.obs['cell_subtype']``        : author-provided subtype (if available)
- ``.obs['dataset']``             : scLucid dataset key
- ``.uns['sclucid']['dataset']``  : structured provenance metadata

Examples
--------
    python scripts/ingest_benchmark_datasets.py --dataset nsclc --output-dir data/
    python scripts/ingest_benchmark_datasets.py --dataset all --output-dir data/
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Dataset provenance metadata
# ---------------------------------------------------------------------------

DATASET_META: dict[str, dict[str, Any]] = {
    "nsclc": {
        "key": "zilionis2019.nsclc",
        "output_name": "zilionis2019.nsclc.h5ad",
        "title": "Zilionis et al. 2019 NSCLC",
        "geo": "GSE127465",
        "pmid": "30979687",
        "citation": (
            "Zilionis R, Engblom C, Pfirschke C, et al. "
            "Single-cell transcriptomics of human and mouse lung cancers "
            "reveals conserved myeloid populations across individuals and species. "
            "Immunity. 2019;50(5):1317-1334.e10."
        ),
        "species": "human",
        "tissue": "lung",
        "condition_values": {"t": "tumor", "b": "blood"},
        "sample_pattern": r"human_(p\d+)([tb]\d+)_raw_counts\.tsv\.gz$",
        "raw_dir": "data/GSE127465_RAW",
        "metadata_file": "data/GSE127465_human_cell_metadata_54773x25.tsv.gz",
        "gene_file": "data/GSE127465_gene_names_human_41861.tsv.gz",
    },
    "pancreas": {
        "key": "baron2016.pancreas",
        "output_name": "baron2016.pancreas.h5ad",
        "title": "Baron et al. 2016 Human Pancreas",
        "geo": "GSE84133",
        "pmid": "27667667",
        "citation": (
            "Baron M, Veres A, Wolock SL, et al. "
            "A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas "
            "Reveals Inter- and Intra-cell Population Structure. "
            "Cell Syst. 2016;3(4):346-360.e4."
        ),
        "species": "human",
        "tissue": "pancreas",
        "condition_values": None,
        "sample_pattern": r"(GSM2230\d+)_human(\d+)_umifm_counts\.csv\.gz$",
        "raw_dir": "data/GSE84133_RAW",
        "metadata_file": None,
        "gene_file": None,
    },
    "crc": {
        "key": "lee2020.crc",
        "output_name": "lee2020.crc.h5ad",
        "title": "Lee et al. 2020 Colorectal Cancer",
        "geo": "GSE132465",
        "pmid": "32451460",
        "citation": (
            "Lee HO, Hong Y, Etlioglu HE, et al. "
            "Lineage-dependent gene expression programs influence the immune landscape "
            "of colorectal cancer. Nat Genet. 2020;52(6):594-603."
        ),
        "species": "human",
        "tissue": "colorectum",
        "condition_values": None,
        "sample_pattern": None,
        "raw_dir": None,
        "metadata_file": "data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
        "gene_file": None,
        "matrix_file": "data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
    },
    "kang2018": {
        "key": "kang2018.pbmc",
        "output_name": "kang2018.pbmc.h5ad",
        "title": "Kang et al. 2018 IFN-stimulated PBMC",
        "geo": "GSE96583",
        "pmid": "29490908",
        "citation": (
            "Kang HM, Subramaniam M, Targ S, et al. "
            "Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. "
            "Nat Biotechnol. 2018;36(1):89-94."
        ),
        "species": "human",
        "tissue": "PBMC",
        "condition_values": None,
        "raw_dir": "data/GSE96583_RAW",
        "metadata_files": {
            "batch1": "data/GSE96583_batch1.total.tsne.df.tsv.gz",
            "batch2": "data/GSE96583_batch2.total.tsne.df.tsv.gz",
        },
        "gene_files": {
            "batch1": "data/GSE96583_batch1.genes.tsv.gz",
            "batch2": "data/GSE96583_batch2.genes.tsv.gz",
        },
        "matrix_files": {
            "A": "data/GSE96583_RAW/GSM2560245_A.mat.gz",
            "B": "data/GSE96583_RAW/GSM2560246_B.mat.gz",
            "C": "data/GSE96583_RAW/GSM2560247_C.mat.gz",
            "2.1": "data/GSE96583_RAW/GSM2560248_2.1.mtx.gz",
            "2.2": "data/GSE96583_RAW/GSM2560249_2.2.mtx.gz",
        },
        "barcode_files": {
            "A": "data/GSE96583_RAW/GSM2560245_barcodes.tsv.gz",
            "B": "data/GSE96583_RAW/GSM2560246_barcodes.tsv.gz",
            "C": "data/GSE96583_RAW/GSM2560247_barcodes.tsv.gz",
            "2.1": "data/GSE96583_RAW/GSM2560248_barcodes.tsv.gz",
            "2.2": "data/GSE96583_RAW/GSM2560249_barcodes.tsv.gz",
        },
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _repo_root() -> Path:
    """Return repository root (where this script lives in scripts/)."""
    return Path(__file__).resolve().parent.parent


def _read_tsv_gz(path: Path, **kwargs: Any) -> pd.DataFrame:
    """Read a gzipped TSV into a DataFrame."""
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", **kwargs)


def _read_csv_gz(path: Path, **kwargs: Any) -> pd.DataFrame:
    """Read a gzipped CSV into a DataFrame."""
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, **kwargs)


def _build_dataset_uns(meta: dict[str, Any]) -> dict[str, Any]:
    """Build standardized ``uns['sclucid']['dataset']`` provenance block."""
    return {
        "dataset_key": meta["key"],
        "title": meta["title"],
        "geo_accession": meta["geo"],
        "pmid": meta.get("pmid"),
        "species": meta["species"],
        "tissue": meta["tissue"],
        "citation": meta["citation"],
        "processing": {
            "tool": "scripts/ingest_benchmark_datasets.py",
            "raw_counts_layer": "counts",
            "normalized_layer": None,
        },
    }


def _make_sparse_anndata(
    counts: sp.csr_matrix,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    meta: dict[str, Any],
) -> ad.AnnData:
    """Create an AnnData object with unified field contract."""
    adata = ad.AnnData(X=counts, obs=obs, var=var, dtype=np.int32)
    adata.layers["counts"] = adata.X.copy()
    adata.obs["dataset"] = meta["key"]
    adata.var_names_make_unique()
    adata.uns.setdefault("sclucid", {})["dataset"] = _build_dataset_uns(meta)
    return adata


# ---------------------------------------------------------------------------
# NSCLC (GSE127465)
# ---------------------------------------------------------------------------


def _ingest_nsclc(root: Path, output_dir: Path, max_cells: int | None = None) -> Path:
    meta = DATASET_META["nsclc"]
    logger.info("[%s] Starting ingestion", meta["key"])

    raw_dir = root / meta["raw_dir"]
    metadata_path = root / meta["metadata_file"]
    gene_path = root / meta["gene_file"]

    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw directory not found: {raw_dir}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    if not gene_path.exists():
        raise FileNotFoundError(f"Gene file not found: {gene_path}")

    # Gene list
    gene_names = pd.read_csv(gene_path, header=None, names=["gene_symbol"])[
        "gene_symbol"
    ].tolist()

    # Metadata
    cell_meta = _read_tsv_gz(metadata_path, low_memory=False)
    # Build cell id used in count files: <library>_<barcode>
    cell_meta["cell_id"] = cell_meta["Library"] + "_" + cell_meta["Barcode"]
    cell_meta = cell_meta.set_index("cell_id")
    # Subsample if requested
    if max_cells is not None and len(cell_meta) > max_cells:
        rng = np.random.default_rng(42)
        keep_idx = rng.choice(cell_meta.index, size=max_cells, replace=False)
        cell_meta = cell_meta.loc[keep_idx]
        logger.info("[%s] Subsampled metadata to %d cells", meta["key"], len(cell_meta))

    keep_cells = set(cell_meta.index)

    # Collect count files
    files = sorted(raw_dir.glob("GSM*_human_*_raw_counts.tsv.gz"))
    logger.info("[%s] Found %d sample count files", meta["key"], len(files))

    row_indices: list[np.ndarray] = []
    col_indices: list[np.ndarray] = []
    values: list[np.ndarray] = []
    obs_records: list[dict[str, Any]] = []
    n_genes = len(gene_names)

    for fpath in files:
        m = re.search(meta["sample_pattern"], fpath.name)
        if not m:
            continue
        patient = m.group(1)
        library = patient + m.group(2)
        tissue_code = m.group(2)[0]  # 't' or 'b'
        condition = meta["condition_values"][tissue_code]

        logger.info("[%s] Reading %s", meta["key"], fpath.name)
        df = _read_tsv_gz(fpath, index_col="barcode")
        # genes are columns; df shape = (cells, genes)
        sample_genes = df.columns.tolist()
        if sample_genes != gene_names:
            # Reorder to global gene list; missing genes filled with 0
            df = df.reindex(columns=gene_names, fill_value=0)

        # Filter to cells in metadata subsample
        df.index = library + "_" + df.index.astype(str)
        df = df.loc[df.index.isin(keep_cells)]
        if df.empty:
            continue

        # Build sparse COO for this sample
        mat = df.to_numpy(dtype=np.int32)
        sample_obs_idx_start = len(obs_records)
        for local_i, cell_id in enumerate(df.index):
            row = mat[local_i]
            nz = row.nonzero()[0]
            if nz.size:
                row_indices.append(nz)
                col_indices.append(np.full_like(nz, sample_obs_idx_start + local_i))
                values.append(row[nz])
            meta_row = cell_meta.loc[cell_id]
            obs_records.append(
                {
                    "sample": library,
                    "patient": patient,
                    "condition": condition,
                    "cell_type": str(meta_row.get("Major cell type", "")),
                    "cell_subtype": str(meta_row.get("Minor subset", "")),
                    "barcoding_emulsion": str(meta_row.get("Barcoding emulsion", "")),
                    "total_counts": float(meta_row.get("Total counts", np.nan)),
                    "pct_counts_mt": float(
                        meta_row.get("Percent counts from mitochondrial genes", np.nan)
                    ),
                }
            )

    if not obs_records:
        raise RuntimeError("No cells matched between metadata and count files")

    all_rows = np.concatenate(row_indices)
    all_cols = np.concatenate(col_indices)
    all_vals = np.concatenate(values)
    counts = sp.coo_matrix(
        (all_vals, (all_rows, all_cols)), shape=(n_genes, len(obs_records))
    ).tocsr()
    counts = counts.transpose().tocsr()

    obs = pd.DataFrame(obs_records)
    obs.index = obs.index.astype(str)
    var = pd.DataFrame({"gene_symbol": gene_names})
    var.index = gene_names

    adata = _make_sparse_anndata(counts, obs, var, meta)
    out_path = output_dir / meta["output_name"]
    adata.write_h5ad(out_path, compression="gzip")
    logger.info("[%s] Wrote %s (%d cells x %d genes)", meta["key"], out_path, adata.n_obs, adata.n_vars)
    return out_path


# ---------------------------------------------------------------------------
# Pancreas (GSE84133 human)
# ---------------------------------------------------------------------------


def _ingest_pancreas(root: Path, output_dir: Path, max_cells: int | None = None) -> Path:
    meta = DATASET_META["pancreas"]
    logger.info("[%s] Starting ingestion", meta["key"])

    raw_dir = root / meta["raw_dir"]
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw directory not found: {raw_dir}")

    files = sorted(raw_dir.glob("GSM*_human*_umifm_counts.csv.gz"))
    if not files:
        raise FileNotFoundError(f"No human pancreas CSV files found in {raw_dir}")
    logger.info("[%s] Found %d human donor files", meta["key"], len(files))

    # Determine global gene set
    gene_names: list[str] | None = None
    for fpath in files:
        with gzip.open(fpath, "rt") as fh:
            header = fh.readline().strip().split(",")
        genes = header[3:]  # skip Unnamed: 0, barcode, assigned_cluster
        if gene_names is None:
            gene_names = genes
        elif genes != gene_names:
            # Intersect to be safe; in practice they should be identical
            gene_names = [g for g in gene_names if g in genes]

    if gene_names is None:
        raise RuntimeError("Could not determine gene list")
    gene_names = list(dict.fromkeys(gene_names))  # preserve order, dedupe
    n_genes = len(gene_names)
    logger.info("[%s] Global gene list: %d genes", meta["key"], n_genes)

    row_indices: list[np.ndarray] = []
    col_indices: list[np.ndarray] = []
    values: list[np.ndarray] = []
    obs_records: list[dict[str, Any]] = []

    for fpath in files:
        m = re.search(meta["sample_pattern"], fpath.name)
        if not m:
            continue
        gsm = m.group(1)
        donor = f"human{m.group(2)}"
        sample = f"{donor}_{gsm}"

        logger.info("[%s] Reading %s", meta["key"], fpath.name)
        df = _read_csv_gz(fpath)
        df = df.rename(columns={"Unnamed: 0": "orig_index"})
        df = df.set_index("barcode")

        # Extract cell type from assigned_cluster
        cell_types = df["assigned_cluster"].astype(str).tolist()
        df = df.drop(columns=["orig_index", "assigned_cluster"])
        # Reorder to global gene list
        df = df.reindex(columns=gene_names, fill_value=0)

        if max_cells is not None:
            n_keep = max(1, int(max_cells * len(df) / 10000))  # rough per-donor cap
            if n_keep < len(df):
                rng = np.random.default_rng(42)
                keep = rng.choice(df.index, size=n_keep, replace=False)
                df = df.loc[keep]
                cell_types = [cell_types[i] for i in df.index.get_indexer_for(keep)]

        mat = df.to_numpy(dtype=np.int32)
        sample_obs_idx_start = len(obs_records)
        for local_i, cell_id in enumerate(df.index):
            row = mat[local_i]
            nz = row.nonzero()[0]
            if nz.size:
                row_indices.append(nz)
                col_indices.append(np.full_like(nz, sample_obs_idx_start + local_i))
                values.append(row[nz])
            obs_records.append(
                {
                    "sample": sample,
                    "donor": donor,
                    "condition": "normal",
                    "cell_type": cell_types[local_i],
                    "cell_subtype": "",
                    "gsm": gsm,
                }
            )

    all_rows = np.concatenate(row_indices)
    all_cols = np.concatenate(col_indices)
    all_vals = np.concatenate(values)
    counts = sp.coo_matrix(
        (all_vals, (all_rows, all_cols)), shape=(n_genes, len(obs_records))
    ).tocsr()
    counts = counts.transpose().tocsr()

    obs = pd.DataFrame(obs_records)
    obs.index = obs.index.astype(str)
    var = pd.DataFrame({"gene_symbol": gene_names})
    var.index = gene_names

    adata = _make_sparse_anndata(counts, obs, var, meta)
    out_path = output_dir / meta["output_name"]
    adata.write_h5ad(out_path, compression="gzip")
    logger.info("[%s] Wrote %s (%d cells x %d genes)", meta["key"], out_path, adata.n_obs, adata.n_vars)
    return out_path


# ---------------------------------------------------------------------------
# Kang 2018 (GSE96583)
# ---------------------------------------------------------------------------


def _ingest_kang2018(root: Path, output_dir: Path, max_cells: int | None = None) -> Path:
    meta = DATASET_META["kang2018"]
    logger.info("[%s] Starting ingestion", meta["key"])

    matrix_files = {k: root / v for k, v in meta["matrix_files"].items()}
    barcode_files = {k: root / v for k, v in meta["barcode_files"].items()}
    gene_files = {k: root / v for k, v in meta["gene_files"].items()}
    metadata_files = {k: root / v for k, v in meta["metadata_files"].items()}

    for d in (matrix_files, barcode_files, gene_files, metadata_files):
        for k, p in d.items():
            if not p.exists():
                raise FileNotFoundError(f"Required file not found: {p}")

    # Load metadata and normalize columns
    meta1 = _read_tsv_gz(metadata_files["batch1"], index_col=0)
    meta2 = _read_tsv_gz(metadata_files["batch2"], index_col=0)
    meta1 = meta1.rename(columns={"cell.type": "cell_type"})
    meta2 = meta2.rename(columns={"cell": "cell_type"})
    meta1["batch_group"] = "batch1"
    meta2["batch_group"] = "batch2"
    meta1["stim"] = "ctrl"
    meta2["stim"] = meta2["stim"].astype(str)
    meta_all = pd.concat([meta1, meta2], axis=0)
    # Barcodes can collide across 10x runs; make cell ids unique per sample
    # by storing original barcode separately and prefixing the index with sample.
    meta_all["barcode"] = meta_all.index
    # The matrix files A/B/C belong to batch1; 2.1/2.2 to batch2. We will later
    # match metadata rows by the original barcode, but the final index must be unique.
    meta_all = meta_all.reset_index(drop=True)

    # Subsample if requested
    if max_cells is not None and len(meta_all) > max_cells:
        rng = np.random.default_rng(42)
        keep_idx = rng.choice(meta_all.index, size=max_cells, replace=False)
        meta_all = meta_all.loc[keep_idx].copy()
        logger.info("[%s] Subsampled metadata to %d cells", meta["key"], len(meta_all))

    # Build a lookup from (sample, barcode) to metadata row index for cells that
    # belong to each sample. Metadata rows do not carry sample id, so we assume
    # each metadata row can appear in at most one sample's barcode list.
    keep_barcodes_by_sample: dict[str, set[str]] = {sid: set() for sid in matrix_files}

    # Gene list — batch2 has a superset of batch1 genes. Use the union to avoid
    # losing batch2 counts, padding batch1 samples with zeros for missing genes.
    gene_df1 = _read_tsv_gz(gene_files["batch1"], header=None, names=["gene_id", "gene_symbol"])
    gene_df2 = _read_tsv_gz(gene_files["batch2"], header=None, names=["gene_id", "gene_symbol"])
    gene_df = pd.concat([gene_df1, gene_df2], ignore_index=True).drop_duplicates(subset=["gene_id"])
    gene_names = gene_df["gene_symbol"].tolist()
    gene_id_to_idx = {gid: i for i, gid in enumerate(gene_df["gene_id"])}
    n_genes = len(gene_names)
    logger.info("[%s] Merged gene list: %d genes", meta["key"], n_genes)

    row_indices: list[np.ndarray] = []
    col_indices: list[np.ndarray] = []
    values: list[np.ndarray] = []
    obs_records: list[dict[str, Any]] = []
    matched_barcodes: set[str] = set()

    for sample_id, matrix_path in matrix_files.items():
        logger.info("[%s] Reading %s", meta["key"], matrix_path.name)
        # .mat.gz and .mtx.gz are both MatrixMarket
        counts = mmread(str(matrix_path)).astype(np.int32)
        if sp.issparse(counts):
            counts = counts.toarray()
        # Original matrix is genes x cells; transpose to cells x genes
        counts = counts.T

        with gzip.open(barcode_files[sample_id], "rt") as fh:
            barcodes = [line.strip() for line in fh]

        if counts.shape[0] != len(barcodes):
            raise ValueError(
                f"Barcode count mismatch for {sample_id}: "
                f"{counts.shape[0]} cells vs {len(barcodes)} barcodes"
            )
        if counts.shape[1] != n_genes:
            # Remap each sample's genes to the merged gene index, padding missing with 0
            sample_genes = _read_tsv_gz(
                gene_files["batch1"] if sample_id in ("A", "B", "C") else gene_files["batch2"],
                header=None,
                names=["gene_id", "gene_symbol"],
            )
            sample_idx_map = np.array(
                [gene_id_to_idx.get(gid, -1) for gid in sample_genes["gene_id"]]
            )
            remapped = np.zeros((counts.shape[0], n_genes), dtype=np.int32)
            valid = sample_idx_map >= 0
            remapped[:, sample_idx_map[valid]] = counts[:, valid]
            counts = remapped
            logger.info("[%s] Remapped %s from %d to %d genes", meta["key"], sample_id, len(sample_genes), n_genes)

        # Resolve metadata rows for this sample by matching barcodes.
        # Because metadata rows do not carry sample id, a barcode appearing in
        # multiple samples would be ambiguous. We assign it to the first sample
        # where it appears and remove it from later samples' lookup.
        barcode_to_meta_idx = {}
        unmatched = []
        for barcode in barcodes:
            if barcode in matched_barcodes:
                continue
            matches = meta_all.index[meta_all["barcode"] == barcode]
            if len(matches) == 1:
                barcode_to_meta_idx[barcode] = matches[0]
                matched_barcodes.add(barcode)
            elif len(matches) > 1:
                # Ambiguous: pick the first unmatched metadata row
                barcode_to_meta_idx[barcode] = matches[0]
                matched_barcodes.add(barcode)
                meta_all = meta_all.drop(matches[1:])
            else:
                unmatched.append(barcode)
        if unmatched:
            logger.info("[%s] %s: %d barcodes not found in metadata", meta["key"], sample_id, len(unmatched))

        for local_i, barcode in enumerate(barcodes):
            if barcode not in barcode_to_meta_idx:
                continue
            meta_idx = barcode_to_meta_idx[barcode]
            row = counts[local_i]
            nz = row.nonzero()[0]
            if nz.size:
                row_indices.append(nz)
                col_indices.append(np.full_like(nz, len(obs_records)))
                values.append(row[nz])
            meta_row = meta_all.loc[meta_idx]
            demuxlet_call = str(meta_row.get("multiplets", ""))
            obs_records.append(
                {
                    "sample": sample_id,
                    "barcode": barcode,
                    "donor": str(meta_row.get("ind", "")),
                    "condition": str(meta_row.get("stim", "")),
                    "cell_type": str(meta_row.get("cell_type", "")),
                    "cell_subtype": str(meta_row.get("cluster", "")),
                    "demuxlet_multiplets": demuxlet_call,
                    "doublet_ground_truth": demuxlet_call == "doublet",
                    "batch_group": str(meta_row.get("batch_group", "")),
                }
            )
        logger.info("[%s] %s: matched %d / %d barcodes", meta["key"], sample_id, len(barcode_to_meta_idx), len(barcodes))

    if not obs_records:
        raise RuntimeError("No cells matched between metadata and count matrices")

    all_rows = np.concatenate(row_indices)
    all_cols = np.concatenate(col_indices)
    all_vals = np.concatenate(values)
    counts_sp = sp.coo_matrix(
        (all_vals, (all_cols, all_rows)), shape=(len(obs_records), n_genes)
    ).tocsr()

    obs = pd.DataFrame(obs_records)
    obs.index = obs.index.astype(str)
    var = pd.DataFrame({"gene_symbol": gene_names})
    var.index = gene_names

    adata = _make_sparse_anndata(counts_sp, obs, var, meta)
    out_path = output_dir / meta["output_name"]
    adata.write_h5ad(out_path, compression="gzip")
    logger.info("[%s] Wrote %s (%d cells x %d genes)", meta["key"], out_path, adata.n_obs, adata.n_vars)
    return out_path


# ---------------------------------------------------------------------------
# CRC (GSE132465)
# ---------------------------------------------------------------------------


def _ingest_crc(root: Path, output_dir: Path, max_cells: int | None = None) -> Path:
    meta = DATASET_META["crc"]
    logger.info("[%s] Starting ingestion", meta["key"])

    matrix_path = root / meta["matrix_file"]
    metadata_path = root / meta["metadata_file"]
    if not matrix_path.exists():
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Read annotation
    cell_meta = _read_tsv_gz(metadata_path)
    cell_meta = cell_meta.rename(columns={"Index": "cell_id"}).set_index("cell_id")
    if max_cells is not None and len(cell_meta) > max_cells:
        rng = np.random.default_rng(42)
        keep_idx = rng.choice(cell_meta.index, size=max_cells, replace=False)
        cell_meta = cell_meta.loc[keep_idx]
        logger.info("[%s] Subsampled metadata to %d cells", meta["key"], len(cell_meta))

    keep_cells = set(cell_meta.index)
    n_cells = len(keep_cells)

    # Read matrix in gene chunks and build sparse incrementally.
    # Matrix is genes x cells; we transpose to cells x genes.
    logger.info("[%s] Reading raw UMI matrix in sparse chunks", meta["key"])
    chunk_size = 1000
    row_indices: list[np.ndarray] = []
    col_indices: list[np.ndarray] = []
    values: list[np.ndarray] = []
    gene_names: list[str] = []

    with gzip.open(matrix_path, "rt") as fh:
        header = fh.readline().strip().split("\t")
        all_cell_ids = header[1:]  # first column is gene id
        # Map global cell id to position in kept subset
        keep_mask = np.array([cid in keep_cells for cid in all_cell_ids])
        kept_positions = np.where(keep_mask)[0]
        kept_id_to_local = {cid: i for i, cid in enumerate(all_cell_ids) if cid in keep_cells}

        gene_idx = 0
        buffer: list[tuple[str, list[int], list[int], list[int]]] = []
        for line in fh:
            parts = line.strip().split("\t")
            gene = parts[0]
            vals = np.array(parts[1:], dtype=np.int32)
            kept_vals = vals[keep_mask]
            nz = kept_vals.nonzero()[0]
            if nz.size:
                row_indices.append(nz)  # local cell index within kept subset
                col_indices.append(np.full_like(nz, gene_idx))
                values.append(kept_vals[nz])
            gene_names.append(gene)
            gene_idx += 1
            if gene_idx % chunk_size == 0:
                logger.info("[%s] Processed %d genes", meta["key"], gene_idx)

    all_rows = np.concatenate(row_indices)
    all_cols = np.concatenate(col_indices)
    all_vals = np.concatenate(values)
    counts = sp.coo_matrix(
        (all_vals, (all_rows, all_cols)), shape=(n_cells, len(gene_names))
    ).tocsr()

    # Reorder rows of counts to match metadata order
    metadata_order = cell_meta.index.tolist()
    position_in_counts = np.array([kept_id_to_local[cid] for cid in metadata_order])
    counts = counts[position_in_counts, :]

    obs = cell_meta.reset_index()
    obs = obs.rename(
        columns={
            "Patient": "patient",
            "Class": "condition",
            "Sample": "sample",
            "Cell_type": "cell_type",
            "Cell_subtype": "cell_subtype",
        }
    )
    obs["condition"] = obs["condition"].str.lower()
    obs.index = obs.index.astype(str)

    var = pd.DataFrame({"gene_symbol": gene_names})
    var.index = gene_names

    adata = _make_sparse_anndata(counts, obs, var, meta)
    out_path = output_dir / meta["output_name"]
    adata.write_h5ad(out_path, compression="gzip")
    logger.info("[%s] Wrote %s (%d cells x %d genes)", meta["key"], out_path, adata.n_obs, adata.n_vars)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        choices=["nsclc", "pancreas", "crc", "kang2018", "all"],
        required=True,
        help="Dataset to ingest.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data"),
        help="Directory for output h5ad files (default: data/).",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Subsample each dataset to at most N cells (for quick testing).",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Repository root (default: parent of scripts/).",
    )
    args = parser.parse_args(argv)

    root = args.repo_root or _repo_root()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = ["nsclc", "pancreas", "crc", "kang2018"] if args.dataset == "all" else [args.dataset]
    failures = []
    for ds in datasets:
        try:
            if ds == "nsclc":
                _ingest_nsclc(root, output_dir, args.max_cells)
            elif ds == "pancreas":
                _ingest_pancreas(root, output_dir, args.max_cells)
            elif ds == "crc":
                _ingest_crc(root, output_dir, args.max_cells)
            elif ds == "kang2018":
                _ingest_kang2018(root, output_dir, args.max_cells)
        except Exception as exc:
            logger.exception("[%s] Ingestion failed: %s", ds, exc)
            failures.append((ds, str(exc)))

    if failures:
        logger.error("Failures: %s", failures)
        return 1
    logger.info("All requested datasets ingested successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
