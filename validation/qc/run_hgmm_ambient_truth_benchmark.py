#!/usr/bin/env python3
"""Evaluate scLucid ambient correction against 10x HGMM species truth."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc.ambient import correct_ambient_rna_linear

DEFAULT_INPUT = Path("data/external/qc_truth/tenx_hgmm_6k/prepared/tenx_hgmm_qc_truth.h5ad")
DEFAULT_RAW_DIR = Path("data/external/qc_truth/tenx_hgmm_6k/source/raw_gene_bc_matrices")
DEFAULT_REGISTRY = Path("validation/dataset_evidence_registry.json")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _thresholds(registry_path: Path) -> dict[str, float]:
    registry = json.loads(registry_path.read_text())
    endpoint = registry["endpoint_definitions"]["qc_ambient_correction"]
    return dict(endpoint["acceptance"])


def _species_metrics(
    before: np.ndarray,
    after: np.ndarray,
    classes: np.ndarray,
) -> dict[str, float | int]:
    human = classes == "HUMAN_SINGLET"
    mouse = classes == "MOUSE_SINGLET"
    if not human.any() or not mouse.any():
        raise ValueError("Both human and mouse singlet truth are required.")
    native_before = np.where(human, before[:, 0], before[:, 1])
    native_after = np.where(human, after[:, 0], after[:, 1])
    off_species_before = np.where(human, before[:, 1], before[:, 0])
    off_species_after = np.where(human, after[:, 1], after[:, 0])
    accurate_before = np.where(human, before[:, 0] > before[:, 1], before[:, 1] > before[:, 0])
    accurate_after = np.where(human, after[:, 0] > after[:, 1], after[:, 1] > after[:, 0])
    return {
        "n_species_singlets": int(len(classes)),
        "n_human_singlets": int(human.sum()),
        "n_mouse_singlets": int(mouse.sum()),
        "off_species_umi_before": float(off_species_before.sum()),
        "off_species_umi_after": float(off_species_after.sum()),
        "contamination_reduction": float(
            1.0 - off_species_after.sum() / max(off_species_before.sum(), 1.0)
        ),
        "native_marker_loss": float(1.0 - native_after.sum() / max(native_before.sum(), 1.0)),
        "identity_accuracy_before": float(np.mean(accurate_before)),
        "identity_accuracy_after": float(np.mean(accurate_after)),
        "identity_accuracy_loss": float(np.mean(accurate_before) - np.mean(accurate_after)),
    }


def _classify(
    metrics: Mapping[str, float | int],
    thresholds: Mapping[str, float],
) -> str:
    passed = (
        float(metrics["contamination_reduction"])
        >= float(thresholds["contamination_reduction_min"])
        and float(metrics["native_marker_loss"]) <= float(thresholds["native_marker_loss_max"])
        and float(metrics["identity_accuracy_loss"])
        <= float(thresholds["identity_accuracy_loss_max"])
    )
    return "PASS" if passed else "FAIL"


def _read_features(path: Path, species: str) -> list[str]:
    return [
        f"{species}::{index}::{line.split(chr(9))[1]}"
        for index, line in enumerate(path.read_text().splitlines())
    ]


def _load_gene_level_truth(
    aggregate_path: Path,
    raw_dir: Path,
    max_background: int,
) -> tuple[ad.AnnData, int, dict[str, str]]:
    aggregate = ad.read_h5ad(aggregate_path)
    obs = aggregate.obs.copy()
    singlet = obs["species_class"].isin(["HUMAN_SINGLET", "MOUSE_SINGLET"])
    background = obs.loc[
        (~obs["vendor_filtered_cell_call"].astype(bool)) & (obs["total_umi"] > 0)
    ].nsmallest(max_background, "total_umi")
    selected = np.sort(
        np.concatenate(
            [np.flatnonzero(singlet.to_numpy()), obs.index.get_indexer(background.index)]
        )
    )

    human_path = raw_dir / "hg19/matrix.mtx"
    mouse_path = raw_dir / "mm10/matrix.mtx"
    human = mmread(human_path).tocsc()[:, selected]
    mouse = mmread(mouse_path).tocsc()[:, selected]
    counts = sparse.hstack([human.T, mouse.T], format="csr")
    selected_obs = obs.iloc[selected].copy()
    selected_obs["background_for_correction"] = ~selected_obs["vendor_filtered_cell_call"].astype(
        bool
    )
    var_names = _read_features(raw_dir / "hg19/genes.tsv", "human")
    var_names.extend(_read_features(raw_dir / "mm10/genes.tsv", "mouse"))
    work = ad.AnnData(
        X=counts,
        obs=selected_obs,
        var=pd.DataFrame(index=pd.Index(var_names, name="species_gene")),
    )
    work.layers["counts"] = counts.copy()
    return (
        work,
        int(human.shape[0]),
        {
            "human_matrix": str(human_path),
            "mouse_matrix": str(mouse_path),
        },
    )


def run(
    input_path: Path,
    output_dir: Path,
    registry_path: Path = DEFAULT_REGISTRY,
    raw_dir: Path = DEFAULT_RAW_DIR,
    max_background: int = 10_000,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = _thresholds(registry_path)
    adata, n_human_features, raw_paths = _load_gene_level_truth(
        input_path,
        raw_dir,
        max_background,
    )
    required_obs = {"vendor_filtered_cell_call", "species_class"}
    missing = sorted(required_obs - set(adata.obs.columns))
    if missing:
        raise ValueError(f"HGMM truth object is missing obs columns: {missing}")
    counts = adata.layers.get("counts", adata.X)
    singlet_mask = adata.obs["species_class"].isin(["HUMAN_SINGLET", "MOUSE_SINGLET"]).to_numpy()
    before_matrix = counts[singlet_mask]
    before = np.column_stack(
        [
            np.asarray(before_matrix[:, :n_human_features].sum(axis=1)).ravel(),
            np.asarray(before_matrix[:, n_human_features:].sum(axis=1)).ravel(),
        ]
    )
    classes = adata.obs.loc[singlet_mask, "species_class"].astype(str).to_numpy()

    started = time.perf_counter()
    correction = correct_ambient_rna_linear(
        adata,
        layer="counts" if "counts" in adata.layers else None,
        output_layer="ambient_corrected_counts",
        empty_droplet_key="background_for_correction",
        record=False,
    )
    runtime_seconds = float(time.perf_counter() - started)
    if not correction.get("corrected"):
        status = "BLOCKED"
        metrics: dict[str, float | int] = {}
        blockers = [str(correction.get("reason", "ambient correction did not run"))]
    else:
        after_matrix = adata.layers["ambient_corrected_counts"][singlet_mask]
        after = np.column_stack(
            [
                np.asarray(after_matrix[:, :n_human_features].sum(axis=1)).ravel(),
                np.asarray(after_matrix[:, n_human_features:].sum(axis=1)).ravel(),
            ]
        )
        metrics = _species_metrics(before, after, classes)
        status = _classify(metrics, thresholds)
        blockers = []

    report = {
        "schema_version": "sclucid_hgmm_ambient_truth_benchmark_v1",
        "status": status,
        "dataset_id": "tenx_hgmm_6k",
        "endpoint_id": "qc_ambient_correction",
        "source": {
            "aggregate_truth": {"path": str(input_path), "sha256": _sha256(input_path)},
            "gene_level_raw_matrices": raw_paths,
        },
        "experimental_unit": "hgmm_6k_library",
        "execution": {
            "method": "scLucid.correct_ambient_rna_linear",
            "runtime_seconds": runtime_seconds,
            "correction_summary": correction,
        },
        "thresholds": thresholds,
        "metrics": metrics,
        "blockers": blockers,
        "limitations": [
            "Vendor filtered calls define the background pool used for correction and are not independent cell-calling truth.",
            "Species labels support cross-species contamination and identity checks, not fragile primary-tissue marker preservation.",
        ],
        "claim_boundary": (
            "PASS would support the registered cross-species ambient endpoint only for this HGMM library."
            if status == "PASS"
            else "scLucid ambient-correction performance is not established on the HGMM truth object."
        ),
    }
    artifact = output_dir / "hgmm_ambient_truth_benchmark.json"
    artifact.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--max-background", type=int, default=10_000)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/qc_hgmm"),
    )
    args = parser.parse_args()
    report = run(
        args.input,
        args.output_dir,
        args.registry,
        args.raw_dir,
        args.max_background,
    )
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
