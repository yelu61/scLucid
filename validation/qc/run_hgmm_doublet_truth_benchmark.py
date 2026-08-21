#!/usr/bin/env python3
"""Evaluate the scLucid doublet default against 10x HGMM species truth."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc.config import DoubletConfig
from validation.qc.run_doublet_evidence_benchmark import (
    _auc,
    _auprc,
    _metrics,
    _run_sclucid_doublet_method,
)

DEFAULT_INPUT = Path("data/external/qc_truth/tenx_hgmm_6k/prepared/tenx_hgmm_qc_truth.h5ad")
DEFAULT_RAW_DIR = Path("data/external/qc_truth/tenx_hgmm_6k/source/raw_gene_bc_matrices")
DEFAULT_REGISTRY = Path("validation/dataset_evidence_registry.json")
METHODS = ["scrublet", "scanpy_scrublet", "scdblfinder_python_pyscdblfinder"]


def _threshold(registry_path: Path) -> float:
    registry = json.loads(registry_path.read_text())
    acceptance = registry["endpoint_definitions"]["qc_doublet_calibration"]["acceptance"]
    return float(acceptance["auprc_regret_vs_best_registered_method_max"])


def _classify(selected_auprc: float, best_auprc: float, max_regret: float) -> str:
    return "PASS" if best_auprc - selected_auprc <= max_regret else "FAIL"


def _read_features(path: Path, species: str) -> list[str]:
    return [
        f"{species}::{index}::{line.split(chr(9))[1]}"
        for index, line in enumerate(path.read_text().splitlines())
    ]


def _load_called_hgmm(aggregate_path: Path, raw_dir: Path) -> ad.AnnData:
    aggregate = ad.read_h5ad(aggregate_path)
    called = aggregate.obs["vendor_filtered_cell_call"].astype(bool).to_numpy()
    selected = np.flatnonzero(called)
    human = mmread(raw_dir / "hg19/matrix.mtx").tocsc()[:, selected]
    mouse = mmread(raw_dir / "mm10/matrix.mtx").tocsc()[:, selected]
    counts = sparse.hstack([human.T, mouse.T], format="csr")
    var_names = _read_features(raw_dir / "hg19/genes.tsv", "human")
    var_names.extend(_read_features(raw_dir / "mm10/genes.tsv", "mouse"))
    obs = aggregate.obs.iloc[selected].copy()
    obs["sample"] = obs["library"].astype(str)
    obs["doublet_ground_truth"] = obs["species_class"].eq("CROSS_SPECIES_DOUBLET")
    work = ad.AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(var_names, name="species_gene")),
    )
    work.layers["counts"] = counts.copy()
    return work


def run(
    input_path: Path,
    raw_dir: Path,
    output_dir: Path,
    registry_path: Path,
    seed: int,
    methods: list[str],
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _load_called_hgmm(input_path, raw_dir)
    truth = adata.obs["doublet_ground_truth"].astype(bool)
    expected_rate = float(truth.mean())
    rows: list[dict[str, Any]] = []
    for method in methods:
        result = _run_sclucid_doublet_method(
            adata,
            method=method,
            expected_rate=expected_rate,
            sample_key="sample",
            seed=seed,
            use_heuristics=False,
        )
        predicted = result["predicted"].astype(bool)
        score = result["score"]
        rows.append(
            {
                "method": method,
                "status": result["status"],
                "auprc": _auprc(truth, score),
                "auc": _auc(truth, score),
                "runtime_seconds": float(result["runtime_seconds"]),
                "error": result["error"],
                **_metrics(truth, predicted),
            }
        )

    evaluable = [row for row in rows if row["status"] == "ok" and row["auprc"] is not None]
    selected_method = DoubletConfig().method
    selected = next((row for row in evaluable if row["method"] == selected_method), None)
    max_regret = _threshold(registry_path)
    blockers: list[str] = []
    if not evaluable:
        status = "BLOCKED"
        best: Mapping[str, Any] | None = None
        blockers.append("No registered doublet method produced an evaluable score.")
    elif selected is None:
        status = "BLOCKED"
        best = max(evaluable, key=lambda row: float(row["auprc"]))
        blockers.append(f"Configured default method was not evaluable: {selected_method}")
    else:
        best = max(evaluable, key=lambda row: float(row["auprc"]))
        status = _classify(float(selected["auprc"]), float(best["auprc"]), max_regret)

    report = {
        "schema_version": "sclucid_hgmm_doublet_truth_benchmark_v1",
        "status": status,
        "dataset_id": "tenx_hgmm_6k",
        "endpoint_id": "qc_doublet_calibration",
        "experimental_unit": "hgmm_6k_library",
        "truth_scope": "cross-species doublets only",
        "n_cells": int(adata.n_obs),
        "n_cross_species_doublets": int(truth.sum()),
        "truth_prevalence": expected_rate,
        "selected_method": selected_method,
        "best_method": None if best is None else best["method"],
        "selected_auprc": None if selected is None else selected["auprc"],
        "best_auprc": None if best is None else best["auprc"],
        "auprc_regret": (
            None if selected is None or best is None else float(best["auprc"] - selected["auprc"])
        ),
        "auprc_regret_max": max_regret,
        "automatic_removal_performed": False,
        "methods": rows,
        "blockers": blockers,
        "limitations": [
            "Cross-species truth does not identify same-species or homotypic doublets.",
            "Scores are ranking outputs, so probability calibration is not claimed.",
            "A PASS would apply only to this controlled HGMM library.",
        ],
        "claim_boundary": (
            "The configured scLucid doublet default met the registered AUPRC-regret threshold on cross-species HGMM truth."
            if status == "PASS"
            else "The configured scLucid doublet default is not scientifically accepted on cross-species HGMM truth."
        ),
    }
    (output_dir / "hgmm_doublet_truth_benchmark.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    pd.DataFrame(rows).to_csv(output_dir / "hgmm_doublet_method_metrics.tsv", sep="\t", index=False)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/current/qc_hgmm")
    )
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--methods", nargs="*", default=METHODS)
    args = parser.parse_args()
    report = run(
        args.input,
        args.raw_dir,
        args.output_dir,
        args.registry,
        args.seed,
        args.methods,
    )
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
