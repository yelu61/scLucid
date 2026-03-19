"""
Intelligent Preprocessing Example - scLucid Smart Parameter Selection

This example demonstrates the intelligent preprocessing capabilities of scLucid:
- Automatic HVG number selection based on variance explanation
- Automatic PCA dimension selection using elbow method
- Automatic n_neighbors/n_pcs optimization using silhouette score
- Automatic clustering resolution selection based on stability
- Automatic batch effect detection and correction recommendation
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import scanpy as sc
from scLucid.preprocess import (
    recommend_intelligent_preprocessing,
    run_intelligent_preprocessing,
    IntelligentPreprocessRecommender,
    PreprocessingWorkflowConfig,
)


def example_basic_usage():
    """Example 1: Basic intelligent preprocessing"""
    print("=" * 70)
    print("Example 1: Basic Intelligent Preprocessing")
    print("=" * 70)

    # Load example data
    adata = sc.datasets.pbmc3k()

    # Run intelligent analysis
    print("\n1. Running intelligent preprocessing analysis...")
    strategy = recommend_intelligent_preprocessing(
        adata,
        plot=True,
        save_dir=Path("./preprocess_output"),
    )

    # View recommendations
    print("\n2. Recommendations:")
    print(f"   Data Profile: {strategy.data_profile.n_cells} cells, "
          f"{strategy.data_profile.strategy_type} strategy")
    print(f"   HVGs: {strategy.hvg.n_top_genes} "
          f"(explains {strategy.hvg.variance_explained:.1%} variance)")
    print(f"   PCA: {strategy.pca.n_pcs} components "
          f"({strategy.pca.variance_explained:.1%} variance)")
    print(f"   Neighbors: n_neighbors={strategy.neighbors.n_neighbors}, "
          f"n_pcs={strategy.neighbors.n_pcs}")
    print(f"   Resolution: {strategy.resolution.resolution} "
          f"(~{strategy.resolution.n_clusters} clusters)")

    # Apply recommendations
    print("\n3. Applying recommendations...")
    config = strategy.to_config()
    print(f"   Config created with HVGs={config.hvg.n_top_genes}, "
          f"n_pcs={config.graph.n_pcs}")

    return strategy


def example_with_batches():
    """Example 2: Preprocessing with batch effect assessment"""
    print("\n" + "=" * 70)
    print("Example 2: Preprocessing with Batch Effect Assessment")
    print("=" * 70)

    # Load data
    adata = sc.datasets.pbmc3k()

    # Simulate batch information
    adata.obs["batch"] = ["batch1" if i < 1500 else "batch2" for i in range(adata.n_obs)]

    # Run intelligent analysis with batch key
    print("\n1. Analyzing with batch information...")
    strategy = recommend_intelligent_preprocessing(
        adata,
        batch_key="batch",
        plot=True,
        save_dir=Path("./preprocess_output_batches"),
    )

    # Check batch correction recommendation
    if strategy.batch_correction:
        print(f"\n2. Batch Effect Assessment:")
        print(f"   Needs correction: {strategy.batch_correction.needs_correction}")
        print(f"   Severity: {strategy.batch_correction.severity_score:.2f}")
        if strategy.batch_correction.needs_correction:
            print(f"   Recommended method: {strategy.batch_correction.recommended_method}")
            print(f"   Alternatives: {strategy.batch_correction.alternative_methods}")

    return strategy


def example_one_step():
    """Example 3: One-step intelligent preprocessing"""
    print("\n" + "=" * 70)
    print("Example 3: One-Step Intelligent Preprocessing")
    print("=" * 70)

    # Load data
    adata = sc.datasets.pbmc3k()

    # Run everything in one step
    print("\n1. Running complete intelligent preprocessing...")
    adata_processed, strategy = run_intelligent_preprocessing(
        adata,
        apply_recommendations=True,
        save_dir="./preprocess_output_onestep",
    )

    print(f"\n2. Preprocessing complete!")
    print(f"   Input: {adata.n_vars} genes")
    print(f"   Output: {adata_processed.n_vars} HVGs")
    print(f"   Clusters: {adata_processed.obs.get('leiden', 'N/A').nunique() if 'leiden' in adata_processed.obs else 'Run clustering separately'}")

    return adata_processed, strategy


def example_custom_config():
    """Example 4: Custom configuration for intelligent preprocessing"""
    print("\n" + "=" * 70)
    print("Example 4: Custom Configuration")
    print("=" * 70)

    # Create custom config
    config = IntelligentPreprocessConfig(
        variance_explained_threshold=0.90,  # Require 90% variance
        min_hvg_genes=1000,
        max_hvg_genes=5000,
        pca_method="cumulative_variance",
        pca_variance_threshold=0.98,
        resolution_search_space=[0.4, 0.8, 1.2, 1.6, 2.0],
        n_bootstrap=50,
    )

    print("\n1. Custom configuration:")
    print(f"   Variance threshold: {config.variance_explained_threshold:.0%}")
    print(f"   HVG range: {config.min_hvg_genes}-{config.max_hvg_genes}")
    print(f"   PCA method: {config.pca_method}")
    print(f"   Resolution search: {config.resolution_search_space}")

    # Use with recommender
    recommender = IntelligentPreprocessRecommender(config=config)
    print(f"\n2. Recommender initialized with custom config")

    return recommender


def example_inspect_before_apply():
    """Example 5: Inspect recommendations before applying"""
    print("\n" + "=" * 70)
    print("Example 5: Inspect Before Apply")
    print("=" * 70)

    adata = sc.datasets.pbmc3k()

    # Get strategy only (don't apply)
    print("\n1. Getting recommendations...")
    _, strategy = run_intelligent_preprocessing(
        adata,
        apply_recommendations=False,
        save_dir="./preprocess_output_inspect",
    )

    # Inspect in detail
    print(f"\n2. Detailed recommendations:")
    print(f"   Overall confidence: {strategy.overall_confidence:.2f}")

    print(f"\n   HVG recommendation:")
    print(f"     - n_genes: {strategy.hvg.n_top_genes}")
    print(f"     - confidence: {strategy.hvg.confidence:.2f}")
    print(f"     - 95% CI: [{strategy.hvg.ci_lower}, {strategy.hvg.ci_upper}]")
    print(f"     - method: {strategy.hvg.method}")

    print(f"\n   Neighbors recommendation:")
    print(f"     - n_neighbors: {strategy.neighbors.n_neighbors}")
    print(f"     - n_pcs: {strategy.neighbors.n_pcs}")
    print(f"     - silhouette: {strategy.neighbors.silhouette_score:.3f}")

    if strategy.concerns:
        print(f"\n3. Concerns:")
        for concern in strategy.concerns:
            print(f"   - {concern}")

    if strategy.recommendations:
        print(f"\n4. Recommendations:")
        for rec in strategy.recommendations:
            print(f"   - {rec}")

    return strategy


def main():
    """Run all examples"""
    print("\n")
    print("*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + " " * 15 + "Intelligent Preprocessing Examples" + " " * 20 + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)
    print()

    # Check if scanpy is available
    try:
        import scanpy
        print("✓ scanpy is available")
    except ImportError:
        print("✗ scanpy not available")
        return 1

    # Create output directories
    Path("./preprocess_output").mkdir(exist_ok=True)

    # Run examples
    try:
        example_basic_usage()
        example_with_batches()
        example_custom_config()
        example_inspect_before_apply()
        # example_one_step()  # Skip full execution in demo

        print("\n" + "=" * 70)
        print("✓ All examples completed successfully!")
        print("=" * 70)
        print("\nCheck ./preprocess_output/ for diagnostic plots")

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
