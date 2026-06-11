"""Intelligent preprocessing example.

Shows how to request data-driven preprocessing recommendations, inspect them,
and optionally apply them through the standard preprocessing workflow.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

OUTPUT_DIR = Path("results/examples/intelligent_preprocess")


def prepare_pbmc_demo():
    adata = sc.datasets.pbmc3k()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"
    return adata


def print_strategy(strategy) -> None:
    print("\nIntelligent preprocessing recommendation")
    print("-" * 45)
    print(f"Cells: {strategy.data_profile.n_cells:,}")
    print(f"Genes: {strategy.data_profile.n_genes:,}")
    print(f"Strategy type: {strategy.data_profile.strategy_type}")
    print(f"HVGs: {strategy.hvg.n_top_genes} (confidence={strategy.hvg.confidence:.2f})")
    print(f"PCs: {strategy.pca.n_pcs} (confidence={strategy.pca.confidence:.2f})")
    print(
        "Neighbors: "
        f"n_neighbors={strategy.neighbors.n_neighbors}, "
        f"n_pcs={strategy.neighbors.n_pcs}"
    )
    print(
        "Resolution: "
        f"{strategy.resolution.resolution} "
        f"(~{strategy.resolution.n_clusters} clusters)"
    )
    if strategy.batch_correction:
        print(f"Batch correction needed: {strategy.batch_correction.needs_correction}")
    if strategy.concerns:
        print("Concerns:")
        for concern in strategy.concerns:
            print(f"  - {concern}")


def recommendation_only() -> None:
    adata = prepare_pbmc_demo()
    _, strategy = scl.pp.run_intelligent_preprocessing(
        adata,
        batch_key="sampleID",
        apply_recommendations=False,
        save_dir=str(OUTPUT_DIR / "recommendation_only"),
        fast_mode=True,
    )
    print_strategy(strategy)
    config = strategy.to_config()
    print(f"\nConfig preview: HVGs={config.hvg.n_top_genes}, n_pcs={config.graph.n_pcs}")


def custom_recommender_config() -> None:
    config = scl.pp.IntelligentPreprocessConfig(
        variance_explained_threshold=0.90,
        min_hvg_genes=1000,
        max_hvg_genes=5000,
        pca_method="cumulative_variance",
        pca_variance_threshold=0.95,
        resolution_search_space=[0.4, 0.8, 1.2],
        n_bootstrap=10,
    )
    recommender = scl.pp.IntelligentPreprocessRecommender(config=config)
    print("\nCustom recommender initialized")
    print(f"HVG range: {config.min_hvg_genes}-{config.max_hvg_genes}")
    print(f"Resolution search: {config.resolution_search_space}")
    return recommender


def optional_apply_recommendations() -> None:
    adata = prepare_pbmc_demo()
    processed, strategy = scl.pp.run_intelligent_preprocessing(
        adata,
        batch_key="sampleID",
        apply_recommendations=True,
        save_dir=str(OUTPUT_DIR / "applied"),
        fast_mode=True,
    )
    print_strategy(strategy)
    print(f"Processed shape: {processed.n_obs:,} cells x {processed.n_vars:,} genes")
    print("Review summary stored in adata.uns['sclucid']['preprocess']")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    recommendation_only()
    custom_recommender_config()

    # Uncomment for a complete preprocessing run using recommended parameters.
    # optional_apply_recommendations()

    print(f"\nOutputs saved under: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
