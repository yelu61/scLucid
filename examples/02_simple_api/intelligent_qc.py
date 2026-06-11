"""Intelligent QC example.

Demonstrates data-driven QC recommendations for normal and tumor-aware contexts.
Use this as a concept example; project notebooks should still inspect thresholds
and filtering effects before deleting cells.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

OUTPUT_DIR = Path("results/examples/intelligent_qc")


def prepare_pbmc_demo():
    adata = sc.datasets.pbmc3k()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"
    return adata


def print_recommendation(label: str, recommendation) -> None:
    print(f"\n{label}")
    print("-" * len(label))
    print(f"Strategy: {recommendation.overall_strategy}")
    print(f"Confidence: {recommendation.overall_confidence:.2f}")
    print(f"Data quality score: {recommendation.data_quality_score:.1f}/100")
    print(
        "min_genes: "
        f"{recommendation.min_genes.threshold} "
        f"[{recommendation.min_genes.ci_lower}, {recommendation.min_genes.ci_upper}]"
    )
    print(
        "max_mt_percent: "
        f"{recommendation.max_mt_percent.threshold:.1f} "
        f"[{recommendation.max_mt_percent.ci_lower:.1f}, "
        f"{recommendation.max_mt_percent.ci_upper:.1f}]"
    )
    if recommendation.concerns:
        print("Concerns:")
        for concern in recommendation.concerns:
            print(f"  - {concern}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = prepare_pbmc_demo()

    normal_rec = scl.qc.recommend_intelligent_qc(
        adata,
        tissue_type="pbmc_or_blood",
        strategy="auto",
        plot=False,
        save_dir=OUTPUT_DIR / "normal",
    )
    print_recommendation("Normal/PBMC recommendation", normal_rec)

    tumor_like = adata.copy()
    tumor_like.obs["pct_counts_mt"] = tumor_like.obs["pct_counts_mt"] * 1.5 + 5
    tumor_rec = scl.qc.recommend_intelligent_qc(
        tumor_like,
        tissue_type="tumor",
        strategy="tumor_aware",
        plot=False,
        save_dir=OUTPUT_DIR / "tumor_aware",
    )
    print_recommendation("Tumor-aware recommendation", tumor_rec)

    print(f"\nRecommendation sidecars saved under: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
