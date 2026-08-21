#!/usr/bin/env python3
"""Build a prediction-blinded expert annotation pack for locked QC validation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS
from validation.qc_preprocess.truth_pack import (
    TruthDatasetSpec,
    build_truth_pack,
    load_external_specs,
)


def _registered_specs(keys: set[str]) -> list[TruthDatasetSpec]:
    specs: list[TruthDatasetSpec] = []
    for item in DATASETS:
        if item.key not in keys:
            continue
        specs.append(
            TruthDatasetSpec(
                key=item.key,
                path=str(item.path),
                tissue=item.tissue,
                dataset_type=("pbmc_or_blood" if "pbmc" in item.key else "tumor_tissue"),
                sample_key="sample",
                source_role=item.modality_role,
            )
        )
    missing = keys - {item.key for item in specs}
    if missing:
        raise ValueError(f"Unknown registered datasets: {sorted(missing)}")
    return specs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "lin2020.pdac",
            "pbmc3k",
            "schlesinger2020.pdac",
            "kang2018.pbmc",
            "public_mixology",
        ],
    )
    parser.add_argument("--external-spec", type=Path)
    parser.add_argument("--seed", type=int, default=20260819)
    parser.add_argument("--primary-per-library", type=int, default=60)
    parser.add_argument("--challenge-per-axis", type=int, default=15)
    args = parser.parse_args()

    specs = _registered_specs(set(args.datasets))
    if args.external_spec:
        specs.extend(load_external_specs(args.external_spec))
    manifest = build_truth_pack(
        specs,
        args.output_dir,
        seed=args.seed,
        primary_per_library=args.primary_per_library,
        challenge_per_axis=args.challenge_per_axis,
    )
    summary = {
        "status": manifest["status"],
        "output_dir": str(args.output_dir),
        "n_datasets": len(manifest["datasets"]),
        "reviewer_dir": str(args.output_dir / "reviewer"),
        "sealed_dir": str(args.output_dir / "sealed"),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
