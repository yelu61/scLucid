#!/usr/bin/env python3
"""Import a completed blinded review workbook into frozen QC label TSV files."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.qc_preprocess.truth_pack import (
    LABEL_COLUMNS,
    sha256_file,
    validate_frozen_label_frame,
)


def import_review_workbook(
    workbook_path: Path,
    pack_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Validate workbook label sheets and write immutable label-only tables."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty frozen-label directory: {output_dir}")

    reviewer_dir = pack_dir / "reviewer"
    sheet_specs = {
        "sample": ("Sample Review", reviewer_dir / "sample_evidence.tsv"),
        "cell": ("Cell Review", reviewer_dir / "cell_evidence.tsv"),
    }
    validated: dict[str, pd.DataFrame] = {}
    issues: list[str] = []
    for label, (sheet_name, evidence_path) in sheet_specs.items():
        workbook_frame = pd.read_excel(workbook_path, sheet_name=sheet_name, dtype=str).fillna("")
        missing = [column for column in LABEL_COLUMNS if column not in workbook_frame]
        if missing:
            issues.append(f"{sheet_name}: missing columns {missing}")
            continue
        labels = workbook_frame[LABEL_COLUMNS].copy()
        evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
        checked, table_issues = validate_frozen_label_frame(evidence, labels)
        issues.extend(f"{sheet_name}: {issue}" for issue in table_issues)
        validated[label] = checked

    if issues:
        return {
            "status": "BLOCKED",
            "issues": issues,
            "claim_boundary": "Labels were not frozen and no acceptance result may be calculated.",
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for label, frame in validated.items():
        output = output_dir / f"{label}_labels.tsv"
        frame[LABEL_COLUMNS].to_csv(output, sep="\t", index=False)
        outputs[label] = str(output)
    manifest = {
        "schema_version": "sclucid_frozen_qc_labels_v1",
        "status": "FROZEN",
        "source_workbook": str(workbook_path.resolve()),
        "source_workbook_sha256": sha256_file(workbook_path),
        "pack_dir": str(pack_dir.resolve()),
        "labels": {
            label: {"path": path, "sha256": sha256_file(Path(path))}
            for label, path in outputs.items()
        },
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--pack-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = import_review_workbook(args.workbook, args.pack_dir, args.output_dir)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "FROZEN" else 2


if __name__ == "__main__":
    raise SystemExit(main())
