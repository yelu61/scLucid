"""QC candidate-impact review on one dataset.

This example reports counterfactual removal differences from the current
DecisionCard. It does not treat retention rate or a heuristic score as evidence
that a strategy is scientifically superior; superiority requires registered
external truth and the validation runners under ``validation/qc``.
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


def review_candidates(adata, *, dataset_type: str = "pbmc_or_blood"):
    """Build a read-only policy comparison for one declared context."""
    return scl.recommend_qc_policy(
        adata,
        scl.ProjectContext(
            dataset_type=dataset_type,
            species="human",
            sample_key="sampleID",
            input_provenance="filtered_counts",
        ),
    )


def print_comparison(label: str, card: scl.DecisionCard) -> None:
    print(f"\n{label}: {card.status}")
    print("=" * 78)
    print(
        f"{'Candidate':<34} {'Status':<25} {'Flagged':>8} {'Extra':>8}"
    )
    print("-" * 78)
    statuses = {item["name"]: item["status"] for item in card.candidates}
    for item in card.comparison:
        print(
            f"{item['name']:<34} {statuses[item['name']]:<25} "
            f"{item['candidate_flagged_cells']:>8} "
            f"{item['additional_vs_selector']:>8}"
        )
    print("-" * 78)
    print(f"Selected-policy impact: {card.affected}")
    print(f"Next action: {card.next_action}")
    print("Claim boundary: candidate disagreement is descriptive, not truth-based accuracy.")


def main(adata=None):
    adata = prepare_pbmc_demo() if adata is None else adata
    print(f"Dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    pbmc_review = review_candidates(adata, dataset_type="pbmc_or_blood")
    tumor_review = review_candidates(adata, dataset_type="tumor_tissue")
    print_comparison("PBMC candidate-impact review", pbmc_review)
    print_comparison("Tumor-context sensitivity review", tumor_review)
    return pbmc_review, tumor_review


if __name__ == "__main__":
    main()
