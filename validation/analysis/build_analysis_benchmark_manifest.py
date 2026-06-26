#!/usr/bin/env python3
"""Build a manifest of analysis validation runners.

This is the entry point for Phase 3 analysis benchmarks. It runs the
annotation, pseudobulk-DE, and proportion consistency benchmarks and writes a
unified manifest.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.analysis.run_annotation_accuracy_benchmark import (
    run as run_annotation_accuracy,
)
from validation.analysis.run_proportion_consistency_benchmark import (
    run as run_proportion_consistency,
)
from validation.analysis.run_pseudobulk_de_type1_error_benchmark import (
    run as run_pseudobulk_de_type1,
)


def run(output_dir: Path, datasets: set[str] | None, max_cells: int | None, seed: int) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    annotation_paths = run_annotation_accuracy(
        output_dir / "annotation_accuracy",
        datasets=datasets,
        max_cells=max_cells,
        seed=seed,
    )
    de_paths = run_pseudobulk_de_type1(output_dir / "pseudobulk_de_type1", seed=seed)
    proportion_paths = run_proportion_consistency(
        output_dir / "proportion_consistency", seed=seed
    )

    manifest_rows: list[dict[str, Any]] = []
    for name, path in annotation_paths.items():
        manifest_rows.append({"benchmark": "annotation_accuracy", "output_name": name, "path": str(path)})
    for name, path in de_paths.items():
        manifest_rows.append({"benchmark": "pseudobulk_de_type1", "output_name": name, "path": str(path)})
    for name, path in proportion_paths.items():
        manifest_rows.append(
            {"benchmark": "proportion_consistency", "output_name": name, "path": str(path)}
        )

    manifest_path = output_dir / "analysis_benchmark_manifest.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    return manifest_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/analysis_benchmarks")
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys for annotation benchmark.")
    parser.add_argument("--max-cells", type=int, default=3000, help="Pilot subset size.")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    manifest = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=args.max_cells,
        seed=args.seed,
    )
    print(f"manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
