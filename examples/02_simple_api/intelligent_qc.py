"""Evidence-calibrated QC review without automatic filtering.

This example compares project contexts through the current DecisionCard/QCPolicy
API. Candidate impacts are review evidence, not a scalar quality score or proof
that one strategy is universally superior.
"""

from __future__ import annotations

import scanpy as sc

import scLucid as scl


def prepare_pbmc_demo():
    adata = sc.datasets.pbmc3k()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"
    return adata


def print_review(label: str, card: scl.DecisionCard) -> None:
    """Print the first-screen decision and counterfactual candidate impacts."""
    print(f"\n{label}: {card.status}")
    print("-" * (len(label) + len(card.status) + 2))
    print(f"Reason: {card.reason}")
    print(f"Affected: {card.affected}")
    print(f"Next action: {card.next_action}")
    print("Candidate impacts:")
    for candidate in card.comparison:
        print(
            f"  - {candidate['name']}: "
            f"flagged={candidate['candidate_flagged_cells']}, "
            f"additional_vs_selector={candidate['additional_vs_selector']}"
        )
    if card.missing_evidence:
        print("Missing evidence:")
        for item in card.missing_evidence:
            print(f"  - {item}")


def review_contexts(adata):
    """Return PBMC and tumor-context reviews without mutating ``adata``."""
    normal = scl.recommend_qc_policy(
        adata,
        scl.ProjectContext(
            dataset_type="pbmc_or_blood",
            species="human",
            sample_key="sampleID",
            input_provenance="filtered_counts",
        ),
    )
    tumor = scl.recommend_qc_policy(
        adata,
        scl.ProjectContext(
            dataset_type="tumor_tissue",
            species="human",
            sample_key="sampleID",
            input_provenance="filtered_counts",
        ),
    )
    return normal, tumor


def main(adata=None):
    adata = prepare_pbmc_demo() if adata is None else adata
    normal_review, tumor_review = review_contexts(adata)
    print_review("PBMC review", normal_review)
    print_review("Tumor-context sensitivity review", tumor_review)
    print("\nNo cells were removed. Apply a reviewed QCPolicy explicitly if appropriate.")
    return normal_review, tumor_review


if __name__ == "__main__":
    main()
