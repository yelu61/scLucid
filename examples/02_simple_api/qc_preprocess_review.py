"""Canonical four-action QC and preprocessing example."""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

DATA_PATH = Path("data/pbmc3k.h5ad")
OUTPUT_DIR = Path("results/examples/qc_preprocess_review")


def show_card(card: scl.DecisionCard) -> None:
    """Print the compact first-screen contract."""
    print(f"\n{card.stage.upper()}: {card.status}")
    print(f"Reason: {card.reason}")
    print(f"Affected: {card.affected}")
    print(f"Next action: {card.next_action}")


def main(
    data_path: Path = DATA_PATH,
    output_dir: Path = OUTPUT_DIR,
) -> Path:
    """Run the canonical review/apply path and return the written artifact."""
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(data_path)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    context = scl.ProjectContext(
        dataset_type="pbmc_or_blood",
        species="human",
        sample_key="sample",
        is_multi_sample=True,
        input_provenance="filtered_counts",
    )

    qc_review = scl.recommend_qc_policy(adata, context)
    show_card(qc_review)
    if qc_review.status == "BLOCKED":
        raise RuntimeError(qc_review.next_action)
    qc_result = scl.apply_qc_policy(adata, qc_review.policy)

    pp_review = scl.recommend_preprocess_policy(
        qc_result.adata,
        context,
        consumer="exploration",
    )
    show_card(pp_review)
    if pp_review.status == "BLOCKED":
        raise RuntimeError(pp_review.next_action)
    pp_result = scl.apply_preprocess_policy(qc_result.adata, pp_review.policy)

    output = output_dir / "pbmc3k_qc_preprocess_result.h5ad"
    pp_result.adata.write_h5ad(output)
    print(f"\nSaved fingerprinted RunEvidence result to: {output}")
    return output


if __name__ == "__main__":
    main()
