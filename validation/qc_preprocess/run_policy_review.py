#!/usr/bin/env python3
"""Run the read-only evidence-calibrated policy review on registered datasets."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scLucid as scl
from validation.dataset_registry import DATASETS


def _context_for(key: str, adata) -> scl.ProjectContext:
    sample_key = "sample" if "sample" in adata.obs else None
    dataset_type = "pbmc_or_blood" if "pbmc" in key else "tumor_tissue"
    return scl.ProjectContext(
        dataset_type=dataset_type,
        sample_key=sample_key,
        is_multi_sample=bool(sample_key and adata.obs[sample_key].nunique() > 1),
        input_provenance="filtered_counts",
    )


def run(output_dir: Path, dataset_keys: set[str]) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    review_rows: list[dict] = []
    sample_rows: list[dict] = []
    candidate_rows: list[dict] = []
    for spec in DATASETS:
        if spec.key not in dataset_keys or not spec.path.exists():
            continue
        adata = ad.read_h5ad(spec.path)
        card = scl.recommend_qc_policy(adata, context=_context_for(spec.key, adata))
        quick_map = card.details.get("quick_map", {})
        review_rows.append(
            {
                "dataset": spec.key,
                "status": card.status,
                "n_cells": adata.n_obs,
                "remove_cells": card.affected["remove_cells"],
                "review_cells": card.affected["review_cells"],
                "blocked_samples": json.dumps(card.affected["blocked_samples"]),
                "quick_map_status": quick_map.get("status"),
                "quick_map_review_clusters": json.dumps(quick_map.get("suspicious_clusters", [])),
                "next_action": card.next_action,
                "superiority_claim": "NOT_EVALUATED",
            }
        )
        for row in card.details["sample_decisions"]:
            sample_rows.append({"dataset": spec.key, **row})
        for row in card.candidates:
            candidate_rows.append({"dataset": spec.key, **row})

    paths = {
        "reviews": output_dir / "policy_review.tsv",
        "samples": output_dir / "sample_gate.tsv",
        "candidates": output_dir / "candidate_comparison.tsv",
    }
    pd.DataFrame(review_rows).to_csv(paths["reviews"], sep="\t", index=False)
    pd.DataFrame(sample_rows).to_csv(paths["samples"], sep="\t", index=False)
    pd.DataFrame(candidate_rows).to_csv(paths["candidates"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["lin2020.pdac", "pbmc3k", "schlesinger2020.pdac"],
    )
    args = parser.parse_args()
    for label, path in run(args.output_dir, set(args.datasets)).items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
