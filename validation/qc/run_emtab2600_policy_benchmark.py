#!/usr/bin/env python3
"""Run the public scLucid QC policy against E-MTAB-2600 microscopy truth."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc import recommend_qc_policy
from scLucid.utils.context import ProjectContext
from validation.qc.evaluate_cell_calling_damage_truth import (
    evaluate_truth_table,
    write_report,
)


def _decisions(obs_names: pd.Index, flagged: list[str]) -> np.ndarray:
    return np.where(obs_names.astype(str).isin(flagged), "REMOVE", "KEEP")


def run(input_path: Path, output_dir: Path) -> dict[str, object]:
    adata = ad.read_h5ad(input_path)
    context = ProjectContext(
        dataset_type="cell_line",
        species="mouse",
        assay="full_length_scrna",
        input_provenance="processed_object",
        cell_calling_source="legacy_c1_capture",
        sample_key="library",
        is_multi_sample=True,
        study_objective="external microscopy QC validation",
    )
    card = recommend_qc_policy(adata, context)
    policy = card.policy
    candidates = {row["name"]: row for row in policy.candidate_policies}
    required = {"expert_global", "per_sample_mad"}
    if not required <= set(candidates):
        raise ValueError(f"missing registered QC baselines: {sorted(required - set(candidates))}")

    frame = pd.DataFrame(index=adata.obs_names.copy())
    frame["library"] = adata.obs["library"].astype(str).to_numpy()
    frame["truth_cell_class"] = adata.obs["intact_or_damaged_or_empty"].astype(str).to_numpy()
    frame["predicted_qc_decision"] = "KEEP"
    frame.loc[frame.index.isin(policy.review_obs_names), "predicted_qc_decision"] = "REVIEW"
    frame.loc[frame.index.isin(policy.remove_obs_names), "predicted_qc_decision"] = "REMOVE"
    frame["baseline_expert_global_qc_decision"] = _decisions(
        frame.index, candidates["expert_global"]["flagged_obs_names"]
    )
    frame["baseline_per_sample_mad_qc_decision"] = _decisions(
        frame.index, candidates["per_sample_mad"]["flagged_obs_names"]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = output_dir / "emtab2600_policy_predictions.tsv"
    frame.reset_index(names="run_accession").to_csv(prediction_path, sep="\t", index=False)
    report = evaluate_truth_table(frame)
    json_path, markdown_path = write_report(report, output_dir)
    card_path = output_dir / "emtab2600_decision_card.json"
    card_path.write_text(json.dumps(card.to_dict(), indent=2) + "\n")
    return {
        "status": report["status"],
        "endpoint_status": report["endpoint_status"],
        "decision_card_status": card.status,
        "prediction_table": str(prediction_path),
        "report_json": str(json_path),
        "report_markdown": str(markdown_path),
        "decision_card": str(card_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "data/external/qc_truth/emtab2600_microscopy_quality/prepared/emtab2600_qc_truth.h5ad"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/qc_emtab2600"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.input, args.output_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
