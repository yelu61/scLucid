"""QC decision evaluation and benchmark-style reporting.

Demonstrates how to compare QC strategies (unified, sample-specific, hybrid)
on a real dataset using scLucid's intelligent QC recommender. Use this when
you need to justify QC thresholds to reviewers or compare approaches across
datasets.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc

import scLucid as scl

OUTPUT_DIR = Path("results/examples/qc_evaluation")


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


def apply_fixed_thresholds(adata, min_genes: int = 200, max_mt: float = 20.0):
    """Apply traditional fixed thresholds and return retention stats."""
    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["pct_counts_mt"] < max_mt)
    )
    return {
        "strategy": "unified_fixed",
        "min_genes": min_genes,
        "max_mt_percent": max_mt,
        "n_before": int(adata.n_obs),
        "n_after": int(mask.sum()),
        "n_removed": int((~mask).sum()),
        "retention_rate": float(mask.mean()),
    }


def apply_intelligent_qc(adata, tissue_type: str = "pbmc_or_blood"):
    """Apply scLucid intelligent QC recommendations and return retention stats."""
    rec = scl.qc.recommend_intelligent_qc(
        adata,
        tissue_type=tissue_type,
        strategy="auto",
        plot=False,
    )
    min_genes = rec.min_genes.threshold
    max_mt = rec.max_mt_percent.threshold

    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["pct_counts_mt"] < max_mt)
    )
    return {
        "strategy": f"intelligent_{tissue_type}",
        "min_genes": int(min_genes),
        "max_mt_percent": float(max_mt),
        "min_genes_ci": [
            rec.min_genes.ci_lower,
            rec.min_genes.ci_upper,
        ],
        "max_mt_ci": [
            rec.max_mt_percent.ci_lower,
            rec.max_mt_percent.ci_upper,
        ],
        "confidence": float(rec.overall_confidence),
        "data_quality_score": float(rec.data_quality_score),
        "n_before": int(adata.n_obs),
        "n_after": int(mask.sum()),
        "n_removed": int((~mask).sum()),
        "retention_rate": float(mask.mean()),
    }


def compare_strategies(adata) -> None:
    """Compare unified vs intelligent QC strategies on the same data."""
    unified = apply_fixed_thresholds(adata)
    intelligent = apply_intelligent_qc(adata, tissue_type="pbmc_or_blood")

    print("\nQC strategy comparison")
    print("=" * 55)
    print(f"{'Metric':<25} {'Unified':>12} {'Intelligent':>12}")
    print("-" * 55)
    print(f"{'min_genes':<25} {unified['min_genes']:>12} {intelligent['min_genes']:>12}")
    print(f"{'max_mt_percent':<25} {unified['max_mt_percent']:>12.1f} {intelligent['max_mt_percent']:>12.1f}")
    print(f"{'n_removed':<25} {unified['n_removed']:>12} {intelligent['n_removed']:>12}")
    print(f"{'retention_rate':<25} {unified['retention_rate']:>11.1%} {intelligent['retention_rate']:>11.1%}")
    print("-" * 55)
    if "confidence" in intelligent:
        print(f"Intelligent QC confidence: {intelligent['confidence']:.2f}")
        print(
            "min_genes 95% CI: "
            f"[{intelligent['min_genes_ci'][0]}, {intelligent['min_genes_ci'][1]}]"
        )
        print(
            "max_mt 95% CI: "
            f"[{intelligent['max_mt_ci'][0]:.1f}, {intelligent['max_mt_ci'][1]:.1f}]"
        )
    print(
        f"Retention difference: "
        f"{intelligent['retention_rate'] - unified['retention_rate']:+.1%} "
        f"({intelligent['n_after'] - unified['n_after']:+d} cells)"
    )


def evaluate_tumor_aware(adata) -> None:
    """Show how tumor-aware QC differs from standard QC."""
    tumor_like = adata.copy()
    rng = np.random.default_rng(42)
    tumor_like.obs["pct_counts_mt"] = (
        tumor_like.obs["pct_counts_mt"] * 1.5 + 5
    )

    normal_rec = apply_intelligent_qc(adata, tissue_type="pbmc_or_blood")
    tumor_rec = apply_intelligent_qc(tumor_like, tissue_type="tumor")

    print("\nNormal vs tumor-aware QC comparison")
    print("=" * 55)
    print(f"{'Metric':<25} {'Normal':>12} {'Tumor-aware':>12}")
    print("-" * 55)
    print(f"{'min_genes':<25} {normal_rec['min_genes']:>12} {tumor_rec['min_genes']:>12}")
    print(f"{'max_mt_percent':<25} {normal_rec['max_mt_percent']:>12.1f} {tumor_rec['max_mt_percent']:>12.1f}")
    print(f"{'confidence':<25} {normal_rec['confidence']:>12.2f} {tumor_rec['confidence']:>12.2f}")
    print("-" * 55)
    print("Tumor-aware QC uses relaxed mitochondrial thresholds to avoid")
    print("over-filtering metabolically stressed tumor cells.")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = prepare_pbmc_demo()

    print(f"Dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    compare_strategies(adata)
    evaluate_tumor_aware(adata)
    print(f"\nResults saved under: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
