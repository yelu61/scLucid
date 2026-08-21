#!/usr/bin/env python3
"""Run current read-only QC review against historical real-project membership.

Historical Step1 membership is an operational reference, not blinded biological
truth.  This runner therefore never emits scientific PASS and never modifies an
AnnData object.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scLucid as scl

PROJECTS = (
    {
        "project": "202604JJH",
        "dataset": "202604JJH",
        "raw": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202604JJH/1-DATA/Step0-combined_raw_data.h5ad"
        ),
        "historical": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202604JJH/1-DATA/Step1-sce_cleaned.h5ad"
        ),
    },
    {
        "project": "202507LPJ",
        "dataset": "202507LPJ",
        "raw": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202507LPJ/data/raw/Step0-combined_raw_data.h5ad"
        ),
        "historical": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202507LPJ/data/processed/Step1-sce_cleaned.h5ad"
        ),
    },
    {
        "project": "202603AK112",
        "dataset": "202603AK112_CT26",
        "raw": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202603AK112/data/raw/Step0-combined_raw_data-CT26.h5ad"
        ),
        "historical": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202603AK112/data/processed/Step1-sce_cleaned-CT26.h5ad"
        ),
    },
    {
        "project": "202603AK112",
        "dataset": "202603AK112_SCC7",
        "raw": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202603AK112/data/raw/Step0-combined_raw_data-SCC7.h5ad"
        ),
        "historical": Path(
            "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/"
            "Ongoing/202603AK112/data/processed/Step1-sce_cleaned-SCC7.h5ad"
        ),
    },
)


def compare_membership(
    raw_names: set[str],
    historical_kept: set[str],
    current_remove: set[str],
    current_review: set[str],
) -> dict[str, int | float]:
    """Compare a current three-way review with historical retained membership."""
    if not historical_kept <= raw_names:
        raise ValueError("Historical Step1 cells must be an exact subset of raw cells.")
    if not current_remove <= raw_names or not current_review <= raw_names:
        raise ValueError("Current QC decisions contain cells absent from the raw object.")
    if current_remove & current_review:
        raise ValueError("A cell cannot be both REMOVE and REVIEW.")

    historical_removed = raw_names - historical_kept
    current_keep = raw_names - current_remove - current_review
    n_raw = len(raw_names)
    n_historical_removed = len(historical_removed)
    union_flagged = current_remove | current_review
    return {
        "n_raw": n_raw,
        "n_historical_kept": len(historical_kept),
        "n_historical_removed": n_historical_removed,
        "historical_removed_fraction": n_historical_removed / n_raw if n_raw else 0.0,
        "n_current_remove": len(current_remove),
        "n_current_review": len(current_review),
        "n_current_keep": len(current_keep),
        "current_flagged_fraction": len(union_flagged) / n_raw if n_raw else 0.0,
        "historical_removed_current_remove": len(historical_removed & current_remove),
        "historical_removed_current_review": len(historical_removed & current_review),
        "historical_removed_current_keep": len(historical_removed & current_keep),
        "historical_removed_flagged_fraction": (
            len(historical_removed & union_flagged) / n_historical_removed
            if n_historical_removed
            else 0.0
        ),
        "historical_kept_current_remove": len(historical_kept & current_remove),
        "historical_kept_current_review": len(historical_kept & current_review),
    }


def _species(adata) -> str:
    names = pd.Index(adata.var_names.astype(str))
    return (
        "mouse"
        if int(names.str.match(r"^mt-").sum()) > int(names.str.match(r"^MT-").sum())
        else "human"
    )


def _context(adata) -> scl.ProjectContext:
    return scl.ProjectContext(
        dataset_type="tumor_tissue",
        species=_species(adata),
        assay="scrna",
        input_provenance="filtered_counts",
        sample_key="sampleID" if "sampleID" in adata.obs else None,
        condition_key="group" if "group" in adata.obs else None,
        batch_key="batch" if "batch" in adata.obs else None,
        is_multi_sample="sampleID" in adata.obs and adata.obs["sampleID"].nunique() > 1,
    )


def run(output_dir: Path, projects: tuple[dict[str, Any], ...] = PROJECTS) -> dict[str, Any]:
    """Execute all available registered project reviews and write the audit."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    blockers: list[str] = []

    for spec in projects:
        missing = [str(path) for path in (spec["raw"], spec["historical"]) if not path.exists()]
        if missing:
            blockers.append(f"{spec['dataset']}: missing {', '.join(missing)}")
            continue

        started = time.perf_counter()
        adata = ad.read_h5ad(spec["raw"])
        historical = ad.read_h5ad(spec["historical"], backed="r")
        try:
            card = scl.recommend_qc_policy(adata, context=_context(adata))
            comparison = compare_membership(
                set(adata.obs_names.astype(str)),
                set(historical.obs_names.astype(str)),
                set(card.policy.remove_obs_names),
                set(card.policy.review_obs_names),
            )
        finally:
            historical.file.close()

        row = {
            "project": spec["project"],
            "dataset": spec["dataset"],
            "review_status": card.status,
            "profile": card.policy.profile,
            "runtime_seconds": round(time.perf_counter() - started, 3),
            "blocked_samples": len(card.affected.get("blocked_samples", [])),
            "review_samples": len(card.affected.get("review_samples", [])),
            "missing_evidence": card.missing_evidence,
            "next_action": card.next_action,
            **comparison,
        }
        rows.append(row)
        for sample in card.policy.sample_decisions:
            sample_rows.append({"project": spec["project"], "dataset": spec["dataset"], **sample})

    report = {
        "schema_version": "sclucid_real_project_qc_concordance_v1",
        "status": "BLOCKED" if blockers else "REVIEW",
        "evidence_class": "HISTORICAL_OPERATIONAL_REFERENCE",
        "datasets_run": len(rows),
        "datasets_expected": len(projects),
        "blockers": blockers,
        "results": rows,
        "claim_boundary": {
            "supports": [
                "The current read-only QC reviewer executes on the registered real-project objects.",
                "Current KEEP/REVIEW/REMOVE decisions can be audited against historical Step1 membership.",
            ],
            "does_not_support": [
                "Historical membership is not blinded biological truth.",
                "Concordance cannot establish superiority or prove that either decision is correct.",
                "This audit does not satisfy the real-project UX gate or any scientific endpoint.",
            ],
        },
        "next_action": (
            "Resolve missing project objects before review."
            if blockers
            else "Use blinded KEEP/REMOVE/UNCERTAIN labels to adjudicate the observed decision conflicts."
        ),
    }
    (output_dir / "real_project_qc_concordance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    flat_rows = [
        {key: json.dumps(value) if isinstance(value, list) else value for key, value in row.items()}
        for row in rows
    ]
    pd.DataFrame(flat_rows).to_csv(
        output_dir / "real_project_qc_concordance.tsv", sep="\t", index=False
    )
    pd.DataFrame(sample_rows).to_csv(
        output_dir / "real_project_qc_sample_review.tsv", sep="\t", index=False
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/qc_real_projects"),
    )
    args = parser.parse_args()
    report = run(args.output_dir)
    print(
        json.dumps({"status": report["status"], "datasets_run": report["datasets_run"]}, indent=2)
    )
    return 0 if report["datasets_run"] == report["datasets_expected"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
