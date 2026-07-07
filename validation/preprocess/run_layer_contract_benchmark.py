#!/usr/bin/env python3
"""Build real-data preprocessing layer-contract evidence tables.

This runner is intentionally lightweight: it does not execute full
preprocessing. Instead, it records what each real h5ad currently promises at
the preprocessing boundary, and which evidence a full Phase 3 benchmark should
produce next.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


MARKER_PANELS: dict[str, tuple[str, ...]] = {
    "immune_t": ("CD3D", "CD3E", "TRAC", "CD4", "CD8A", "NKG7"),
    "myeloid": ("LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "LST1"),
    "b_plasma": ("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
    "epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
    "stromal": ("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"),
    "endothelial": ("PECAM1", "VWF", "KDR", "ENG"),
    "proliferation": ("MKI67", "TOP2A", "STMN1", "UBE2C"),
    "hypoxia_stress": ("VEGFA", "CA9", "DDIT3", "HSPA1A", "JUN"),
}


EXPECTED_TRANSITIONS: tuple[dict[str, str], ...] = (
    {
        "step": "input",
        "input_layer": "layers['counts']",
        "output_slot": "adata.X",
        "adata_X_semantics_before": "raw_counts",
        "adata_X_semantics_after": "raw_counts",
        "raw_semantics": "unset_or_source_defined",
        "evidence_key": "input_counts_contract",
    },
    {
        "step": "normalization",
        "input_layer": "layers['counts']",
        "output_slot": "layers['normalized']; adata.X",
        "adata_X_semantics_before": "raw_counts",
        "adata_X_semantics_after": "normalized_counts",
        "raw_semantics": "unset_or_normalized_full_gene_snapshot",
        "evidence_key": "normalization_contract",
    },
    {
        "step": "hvg_selection",
        "input_layer": "layers['normalized']",
        "output_slot": "var['highly_variable']; uns review summary",
        "adata_X_semantics_before": "normalized_counts",
        "adata_X_semantics_after": "normalized_counts",
        "raw_semantics": "normalized_full_gene_snapshot_if_set_raw",
        "evidence_key": "hvg_marker_preservation",
    },
    {
        "step": "scaling_pca",
        "input_layer": "layers['normalized'] or adata.X",
        "output_slot": "layers['scaled']; obsm['X_pca']",
        "adata_X_semantics_before": "normalized_or_hvg_subset",
        "adata_X_semantics_after": "scaled_or_normalized_depending_config",
        "raw_semantics": "normalized_full_gene_snapshot_if_set_raw",
        "evidence_key": "pca_parameter_evidence",
    },
    {
        "step": "batch_correction_optional",
        "input_layer": "obsm['X_pca']",
        "output_slot": "obsm['X_pca_harmony'] or method-specific key",
        "adata_X_semantics_before": "unchanged",
        "adata_X_semantics_after": "unchanged",
        "raw_semantics": "unchanged",
        "evidence_key": "batch_correction_diagnostic",
    },
    {
        "step": "neighbors_umap",
        "input_layer": "obsm['X_pca'] or corrected representation",
        "output_slot": "uns['neighbors']; obsp connectivities/distances; obsm['X_umap']",
        "adata_X_semantics_before": "unchanged",
        "adata_X_semantics_after": "unchanged",
        "raw_semantics": "unchanged",
        "evidence_key": "graph_embedding_handoff",
    },
)


def _value_counts(adata: ad.AnnData, column: str, limit: int = 20) -> dict[str, int]:
    if column not in adata.obs:
        return {}
    counts = adata.obs[column].astype(str).value_counts(dropna=False).head(limit)
    return {str(k): int(v) for k, v in counts.items()}


def _dataset_metadata(adata: ad.AnnData) -> dict[str, Any]:
    return dict(adata.uns.get("sclucid", {}).get("dataset", {}))


def _present_keys(adata: ad.AnnData, keys: tuple[str, ...], slot: str) -> str:
    source = getattr(adata, slot)
    return ";".join(key for key in keys if key in source)


def _primary_batch_key(adata: ad.AnnData) -> str:
    for key in ("donor", "patient", "sample", "sampleID", "condition"):
        if key in adata.obs and adata.obs[key].astype(str).nunique(dropna=False) > 1:
            return key
    return ""


def _biology_key(adata: ad.AnnData) -> str:
    for key in ("cell_type", "condition", "cell_subtype"):
        if key in adata.obs and adata.obs[key].astype(str).nunique(dropna=False) > 1:
            return key
    return ""


def _input_row(spec, adata: ad.AnnData) -> dict[str, Any]:
    dataset_meta = _dataset_metadata(adata)
    required_missing = [key for key in spec.required_obs if key not in adata.obs]
    return {
        "dataset": spec.key,
        "path": str(spec.path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "x_current_semantics": "raw_counts_expected_by_contract",
        "has_counts_layer": bool("counts" in adata.layers),
        "layers_present": ";".join(map(str, adata.layers.keys())),
        "raw_present": bool(adata.raw is not None),
        "obsm_present": ";".join(map(str, adata.obsm.keys())),
        "uns_dataset_key": dataset_meta.get("dataset_key", ""),
        "geo_accession": dataset_meta.get("geo_accession", ""),
        "pmid": dataset_meta.get("pmid", ""),
        "required_obs_missing": ";".join(required_missing),
        "obs_contract_keys": _present_keys(
            adata,
            (
                "sample",
                "sampleID",
                "patient",
                "donor",
                "condition",
                "cell_type",
                "cell_subtype",
                "dataset",
                "demuxlet_multiplets",
                "doublet_ground_truth",
            ),
            "obs",
        ),
        "sample_distribution": json.dumps(_value_counts(adata, "sample"), ensure_ascii=False),
        "condition_distribution": json.dumps(
            _value_counts(adata, "condition"), ensure_ascii=False
        ),
        "cell_type_distribution_top": json.dumps(
            _value_counts(adata, "cell_type", limit=12), ensure_ascii=False
        ),
        "review_required": bool(required_missing or "counts" not in adata.layers),
        "risk_note": (
            "Missing preprocessing input contract fields."
            if required_missing or "counts" not in adata.layers
            else ""
        ),
    }


def _transition_rows(spec, adata: ad.AnnData) -> list[dict[str, Any]]:
    has_counts = "counts" in adata.layers
    rows: list[dict[str, Any]] = []
    for transition in EXPECTED_TRANSITIONS:
        step = transition["step"]
        review_required = False
        risk_note = ""
        if step == "input" and not has_counts:
            review_required = True
            risk_note = "layers['counts'] missing; preprocessing cannot make count-layer claims."
        if step == "batch_correction_optional" and not _primary_batch_key(adata):
            risk_note = "No multi-level batch key detected; only no-correction baseline is relevant."
        rows.append(
            {
                "dataset": spec.key,
                "workflow": "run_preprocessing",
                "step": step,
                "input_layer": transition["input_layer"],
                "output_slot": transition["output_slot"],
                "adata_X_semantics_before": transition["adata_X_semantics_before"],
                "adata_X_semantics_after": transition["adata_X_semantics_after"],
                "raw_semantics": transition["raw_semantics"],
                "parameters": json.dumps(
                    {
                        "primary_batch_key": _primary_batch_key(adata),
                        "biology_key_to_protect": _biology_key(adata),
                    }
                ),
                "evidence_key": transition["evidence_key"],
                "review_required": review_required,
                "risk_note": risk_note,
            }
        )
    return rows


def _marker_rows(spec, adata: ad.AnnData) -> list[dict[str, Any]]:
    var_names = pd.Index(adata.var_names.astype(str))
    rows: list[dict[str, Any]] = []
    for panel, genes in MARKER_PANELS.items():
        present = [gene for gene in genes if gene in var_names]
        rows.append(
            {
                "dataset": spec.key,
                "marker_panel": panel,
                "genes_expected": len(genes),
                "genes_present": len(present),
                "coverage": len(present) / max(len(genes), 1),
                "present_genes": ";".join(present),
                "preprocess_roles": ";".join(spec.preprocess_roles),
                "review_required": bool(
                    (
                        "tumor_marker_preservation" in spec.preprocess_roles
                        or "marker_preservation" in spec.preprocess_roles
                    )
                    and panel in {"epithelial", "stromal", "hypoxia_stress"}
                    and not present
                ),
            }
        )
    return rows


def _batch_rows(spec, adata: ad.AnnData) -> list[dict[str, Any]]:
    keys = [key for key in ("sample", "sampleID", "patient", "donor", "condition") if key in adata.obs]
    rows: list[dict[str, Any]] = []
    for key in keys:
        n_groups = int(adata.obs[key].astype(str).nunique(dropna=False))
        rows.append(
            {
                "dataset": spec.key,
                "candidate_batch_key": key,
                "n_groups": n_groups,
                "group_sizes_top": json.dumps(_value_counts(adata, key), ensure_ascii=False),
                "biology_key_to_protect": _biology_key(adata),
                "recommended_use": (
                    "batch_diagnostic_candidate"
                    if n_groups > 1 and key != _biology_key(adata)
                    else "biology_or_single_group_key"
                ),
                "phase3_comparisons": "no_correction;diagnostic_only;harmony_opt_in;scanpy_standard",
            }
        )
    return rows


def run(output_dir: Path, datasets: set[str] | None = None) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    input_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    marker_rows: list[dict[str, Any]] = []
    batch_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.preprocess_roles or not spec.path.exists():
            continue
        adata = ad.read_h5ad(spec.path, backed="r")
        input_rows.append(_input_row(spec, adata))
        transition_rows.extend(_transition_rows(spec, adata))
        marker_rows.extend(_marker_rows(spec, adata))
        batch_rows.extend(_batch_rows(spec, adata))
        adata.file.close()

    paths = {
        "input_contract": output_dir / "preprocess_input_contract.tsv",
        "layer_contract": output_dir / "layer_contract_report.tsv",
        "marker_coverage": output_dir / "marker_panel_coverage.tsv",
        "batch_inputs": output_dir / "batch_diagnostic_inputs.tsv",
    }
    pd.DataFrame(input_rows).to_csv(paths["input_contract"], sep="\t", index=False)
    pd.DataFrame(transition_rows).to_csv(paths["layer_contract"], sep="\t", index=False)
    pd.DataFrame(marker_rows).to_csv(paths["marker_coverage"], sep="\t", index=False)
    pd.DataFrame(batch_rows).to_csv(paths["batch_inputs"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/preprocess_layer_contract"),
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        help="Optional dataset keys to include. Defaults to all preprocess-relevant datasets.",
    )
    args = parser.parse_args()
    paths = run(args.output_dir, datasets=set(args.datasets) if args.datasets else None)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
