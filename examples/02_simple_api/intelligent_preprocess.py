"""Evidence-calibrated preprocessing review and optional explicit execution."""

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


def print_review(card: scl.DecisionCard) -> None:
    print(f"\nPreprocessing review: {card.status}")
    print("-" * 45)
    print(f"Recommended: {card.recommended}")
    print(f"Reason: {card.reason}")
    print(f"Affected: {card.affected}")
    print("Candidates:")
    for candidate in card.candidates:
        selected = " selected" if candidate.get("selected") else ""
        print(f"  - {candidate.get('name', candidate.get('decision'))}: {candidate['status']}{selected}")
    print(f"Next action: {card.next_action}")


def review_policy(adata):
    """Return the read-only default exploration policy."""
    return scl.recommend_preprocess_policy(
        adata,
        scl.ProjectContext(
            dataset_type="pbmc_or_blood",
            species="human",
            sample_key="sampleID",
            input_provenance="filtered_counts",
        ),
        consumer="exploration",
    )


def main(adata=None, *, apply: bool = False):
    adata = prepare_pbmc_demo() if adata is None else adata
    review = review_policy(adata)
    print_review(review)
    if not apply:
        print("\nReview only: the input AnnData was not modified.")
        return review
    if review.status == "BLOCKED":
        raise RuntimeError(review.next_action)
    evidence = scl.apply_preprocess_policy(adata, review.policy)
    print(f"\nApplied policy {review.policy.policy_id}; RunEvidence={evidence.run_id}")
    print("Formal count models remain bound to layers['counts'].")
    return evidence


if __name__ == "__main__":
    main()
