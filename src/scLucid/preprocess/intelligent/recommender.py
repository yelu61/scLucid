"""
Intelligent preprocessing recommender system.

This module provides data-driven parameter recommendations for preprocessing,
integrating with existing neighbors.py optimization functionality.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.metrics import silhouette_score

from ..config import NeighborsConfig, PreprocessingWorkflowConfig
from .config import IntelligentPreprocessConfig
from .data_classes import (
    BatchCorrectionRecommendation,
    DataProfile,
    HVGRecommendation,
    NeighborsRecommendation,
    PCARecommendation,
    PreprocessingStrategy,
    ResolutionRecommendation,
)

log = logging.getLogger(__name__)


class IntelligentPreprocessRecommender:
    """
    Intelligent preprocessing parameter recommendation system.

    Analyzes data characteristics and provides data-driven recommendations
    for all preprocessing parameters, integrating with existing neighbors.py
    for n_neighbors/n_pcs optimization.

    Examples
    --------
    >>> recommender = IntelligentPreprocessRecommender()
    >>> strategy = recommender.recommend(adata, batch_key="sampleID")
    >>>
    >>> # Apply recommendations
    >>> config = strategy.to_config()
    >>> adata = run_preprocessing(adata, config=config)
    """

    def __init__(
        self,
        config: Optional[IntelligentPreprocessConfig] = None,
        strategy: str = "auto",
    ):
        """
        Initialize the recommender.

        Parameters
        ----------
        config : IntelligentPreprocessConfig, optional
            Configuration for recommendation algorithms
        strategy : str, default="auto"
            Overall strategy: "auto", "minimal", "standard", "aggressive"
        """
        self.config = config or IntelligentPreprocessConfig()
        self.strategy = strategy
        self._data_profile: Optional[DataProfile] = None

    def recommend(
        self,
        adata: AnnData,
        batch_key: Optional[str] = None,
        tissue_type: str = "unknown",
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> PreprocessingStrategy:
        """
        Generate complete preprocessing strategy recommendation.

        This is the main entry point that orchestrates all analyses.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix
        batch_key : str, optional
            Column identifying batches
        tissue_type : str, default="unknown"
            Tissue type for context
        plot : bool, default=True
            Generate diagnostic plots
        save_dir : Path, optional
            Directory to save outputs

        Returns
        -------
        PreprocessingStrategy with all recommendations
        """
        log.info("=" * 60)
        log.info("=== Starting Intelligent Preprocessing Analysis ===")
        log.info("=" * 60)

        # Step 1: Profile data characteristics
        log.info("Step 1: Profiling data characteristics...")
        self._data_profile = DataProfile.from_adata(
            adata, batch_key=batch_key, config=self.config
        )
        log.info(f"Data profile: {self._data_profile.n_cells} cells, "
                 f"{self._data_profile.n_genes} genes")
        log.info(f"Strategy type: {self._data_profile.strategy_type}")

        # Step 2: Recommend HVG parameters
        log.info("Step 2: Analyzing HVG selection...")
        hvg_rec = self.recommend_hvg(adata, plot=plot, save_dir=save_dir)
        log.info(f"Recommended HVGs: {hvg_rec.n_top_genes} "
                 f"(explains {hvg_rec.variance_explained:.1%} variance)")

        # Step 3: Compute HVGs and PCA for downstream analysis
        log.info("Step 3: Computing HVGs and PCA...")
        adata_temp = self._prepare_data(adata, hvg_rec.n_top_genes)

        # Step 4: Recommend PCA dimensions
        log.info("Step 4: Analyzing PCA dimensions...")
        pca_rec = self.recommend_pca(adata_temp, plot=plot, save_dir=save_dir)
        log.info(f"Recommended PCs: {pca_rec.n_pcs} "
                 f"(explains {pca_rec.variance_explained:.1%} variance)")

        # Step 5: Recommend neighbors/PCs (using existing neighbors.py)
        log.info("Step 5: Optimizing neighbors and PCs...")
        neighbors_rec = self.recommend_neighbors(
            adata_temp, plot=plot, save_dir=save_dir
        )
        log.info(f"Recommended: n_neighbors={neighbors_rec.n_neighbors}, "
                 f"n_pcs={neighbors_rec.n_pcs} "
                 f"(silhouette={neighbors_rec.silhouette_score:.3f})")

        # Step 6: Recommend clustering resolution
        log.info("Step 6: Analyzing clustering resolution...")
        resolution_rec = self.recommend_resolution(
            adata_temp, use_rep="X_pca", plot=plot, save_dir=save_dir
        )
        log.info(f"Recommended resolution: {resolution_rec.resolution} "
                 f"(~{resolution_rec.n_clusters} clusters, "
                 f"stability={resolution_rec.stability_score:.3f})")

        # Step 7: Assess batch effects
        batch_rec = None
        if batch_key:
            log.info("Step 7: Assessing batch effects...")
            batch_rec = self.assess_batch_effects(
                adata_temp, batch_key, plot=plot, save_dir=save_dir
            )
            if batch_rec.needs_correction:
                log.info(f"Batch correction recommended: {batch_rec.recommended_method} "
                         f"(severity={batch_rec.severity_score:.2f})")
            else:
                log.info("No significant batch effects detected")
        else:
            log.info("Step 7: Skipping batch assessment (no batch_key provided)")

        # Step 8: Compile overall strategy
        overall_confidence = np.mean([
            hvg_rec.confidence,
            pca_rec.confidence,
            neighbors_rec.confidence,
            resolution_rec.confidence,
            batch_rec.confidence if batch_rec else 1.0,
        ])

        strategy = PreprocessingStrategy(
            data_profile=self._data_profile,
            hvg=hvg_rec,
            pca=pca_rec,
            neighbors=neighbors_rec,
            resolution=resolution_rec,
            batch_correction=batch_rec,
            overall_confidence=overall_confidence,
            concerns=self._data_profile.potential_issues.copy(),
            recommendations=self._generate_recommendations(
                hvg_rec, pca_rec, neighbors_rec, resolution_rec, batch_rec
            ),
        )

        log.info("=" * 60)
        log.info(f"=== Analysis Complete (confidence: {overall_confidence:.2f}) ===")
        log.info("=" * 60)

        return strategy

    def recommend_hvg(
        self,
        adata: AnnData,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> HVGRecommendation:
        """
        Recommend optimal n_top_genes based on variance explanation.

        Algorithm:
        1. Test multiple HVG thresholds
        2. Calculate cumulative variance explained by each set
        3. Find smallest n_top_genes reaching variance threshold
        4. Bootstrap to estimate confidence interval
        """
        from sklearn.decomposition import PCA

        cfg = self.config

        # Define search space based on data size and config
        search_space = np.linspace(
            cfg.min_hvg_genes,
            min(cfg.max_hvg_genes, adata.n_vars),
            cfg.hvg_search_points,
        ).astype(int)
        search_space = np.unique(search_space)  # Remove duplicates

        log.info(f"Testing {len(search_space)} HVG thresholds...")

        # Calculate variance explained for each threshold
        results = []
        for n_genes in search_space:
            try:
                # Compute HVGs
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_genes,
                    flavor="seurat_v3",
                    subset=False,
                )

                # Subset to HVGs and compute PCA
                adata_hvg = adata[:, adata.var.highly_variable].copy()
                sc.pp.scale(adata_hvg)
                sc.tl.pca(adata_hvg, n_comps=min(50, n_genes - 1))

                # Calculate cumulative variance
                var_explained = adata_hvg.uns["pca"]["variance_ratio"].sum()

                results.append({
                    "n_genes": n_genes,
                    "variance_explained": var_explained,
                })
            except Exception as e:
                log.warning(f"Failed for n_genes={n_genes}: {e}")

        if not results:
            # Fallback to default
            return HVGRecommendation(
                n_top_genes=2000,
                variance_explained=0.0,
                ci_lower=1500,
                ci_upper=2500,
                method="fallback_default",
                confidence=0.5,
            )

        df = pd.DataFrame(results)

        # Find optimal n_genes (smallest reaching threshold)
        threshold = cfg.variance_explained_threshold
        above_threshold = df[df["variance_explained"] >= threshold]

        if len(above_threshold) > 0:
            optimal_n = above_threshold["n_genes"].min()
            optimal_var = above_threshold[
                above_threshold["n_genes"] == optimal_n
            ]["variance_explained"].values[0]
            method = "variance_threshold"
        else:
            # Use elbow method if threshold not reached
            optimal_n = self._find_elbow(df["n_genes"].values, df["variance_explained"].values)
            optimal_var = df[df["n_genes"] == optimal_n]["variance_explained"].values[0]
            method = "elbow"

        # Bootstrap for confidence interval
        ci_lower, ci_upper = self._bootstrap_hvg(
            adata, optimal_n, cfg.n_bootstrap, cfg.confidence_level
        )

        # Calculate confidence based on variance achieved
        confidence = min(1.0, optimal_var / threshold) if threshold > 0 else 0.5

        # Save plot
        if plot and save_dir:
            self._plot_hvg_analysis(df, optimal_n, threshold, save_dir)

        return HVGRecommendation(
            n_top_genes=int(optimal_n),
            variance_explained=float(optimal_var),
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method=method,
            confidence=confidence,
            evidence={
                "search_space": search_space.tolist(),
                "variance_curve": df.to_dict(),
                "target_threshold": threshold,
            },
        )

    def recommend_pca(
        self,
        adata: AnnData,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> PCARecommendation:
        """
        Recommend optimal n_pcs using elbow method or cumulative variance.
        """
        cfg = self.config

        # Run PCA with max components
        max_pcs = min(cfg.max_pcs, adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=max_pcs)

        variance_ratios = adata.uns["pca"]["variance_ratio"]
        cumulative_var = np.cumsum(variance_ratios)

        # Method selection
        if cfg.pca_method == "elbow":
            n_pcs = self._find_elbow(
                np.arange(1, len(cumulative_var) + 1),
                cumulative_var,
            )
            method = "elbow"
        elif cfg.pca_method == "cumulative_variance":
            n_pcs = np.argmax(cumulative_var >= cfg.pca_variance_threshold) + 1
            method = "cumulative_variance"
        else:  # knee
            n_pcs = self._find_knee_point(variance_ratios)
            method = "knee"

        # Apply bounds
        n_pcs = max(cfg.min_pcs, min(n_pcs, max_pcs))
        var_explained = cumulative_var[n_pcs - 1]

        # Bootstrap CI
        ci_lower, ci_upper = self._bootstrap_pca(
            adata, n_pcs, cfg.n_bootstrap, cfg.confidence_level
        )

        # Confidence based on variance stability
        confidence = self._calculate_pca_confidence(variance_ratios, n_pcs)

        if plot and save_dir:
            self._plot_pca_analysis(variance_ratios, cumulative_var, n_pcs, save_dir)

        return PCARecommendation(
            n_pcs=int(n_pcs),
            variance_explained=float(var_explained),
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method=method,
            confidence=confidence,
            evidence={
                "variance_ratios": variance_ratios.tolist(),
                "cumulative_variance": cumulative_var.tolist(),
            },
        )

    def recommend_neighbors(
        self,
        adata: AnnData,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> NeighborsRecommendation:
        """
        Recommend optimal n_neighbors and n_pcs using existing neighbors.py optimization.

        This integrates with the existing optimize_neighbors_pcs function.
        """
        from ..neighbors import optimize_neighbors_pcs

        cfg = self.config

        if not cfg.optimize_neighbors:
            # Use simple defaults
            return NeighborsRecommendation(
                n_neighbors=15,
                n_pcs=min(30, adata.n_vars - 1),
                silhouette_score=0.0,
                ci_lower_neighbors=10,
                ci_upper_neighbors=20,
                ci_lower_pcs=20,
                ci_upper_pcs=40,
                method="default",
                confidence=0.5,
            )

        # Use existing neighbors.py optimization
        neighbors_cfg = NeighborsConfig(
            n_neighbors_list=cfg.neighbors_search_space,
            n_pcs_list=list(range(cfg.min_pcs, min(cfg.max_pcs, adata.n_vars), 10)),
            use_rep="X_pca",
            subsample=cfg.silhouette_sample_size,
            plot=plot,
            save_dir=str(save_dir) if save_dir else None,
        )

        try:
            results_df = optimize_neighbors_pcs(adata, config=neighbors_cfg)

            if results_df.empty:
                raise ValueError("No valid results from neighbors optimization")

            # Find best parameters
            best_idx = results_df["silhouette_score"].idxmax()
            best = results_df.loc[best_idx]

            # Calculate confidence intervals from the grid
            ci_lower_n = max(
                cfg.neighbors_search_space[0],
                int(best["n_neighbors"] * 0.7)
            )
            ci_upper_n = min(
                cfg.neighbors_search_space[-1],
                int(best["n_neighbors"] * 1.3)
            )
            ci_lower_pcs = max(cfg.min_pcs, int(best["n_pcs"] * 0.7))
            ci_upper_pcs = min(cfg.max_pcs, int(best["n_pcs"] * 1.3))

            # Confidence based on silhouette score quality
            confidence = min(1.0, best["silhouette_score"] + 0.5)

            return NeighborsRecommendation(
                n_neighbors=int(best["n_neighbors"]),
                n_pcs=int(best["n_pcs"]),
                silhouette_score=float(best["silhouette_score"]),
                ci_lower_neighbors=ci_lower_n,
                ci_upper_neighbors=ci_upper_n,
                ci_lower_pcs=ci_lower_pcs,
                ci_upper_pcs=ci_upper_pcs,
                method="silhouette_grid_search",
                confidence=confidence,
                search_results=results_df,
                evidence={
                    "grid_search_space": {
                        "n_neighbors": cfg.neighbors_search_space,
                        "n_pcs": neighbors_cfg.n_pcs_list,
                    },
                    "all_results": results_df.to_dict(),
                },
            )

        except Exception as e:
            log.warning(f"Neighbors optimization failed: {e}. Using defaults.")
            return NeighborsRecommendation(
                n_neighbors=15,
                n_pcs=min(30, adata.n_vars - 1),
                silhouette_score=0.0,
                ci_lower_neighbors=10,
                ci_upper_neighbors=20,
                ci_lower_pcs=20,
                ci_upper_pcs=40,
                method="fallback_default",
                confidence=0.5,
            )

    def recommend_resolution(
        self,
        adata: AnnData,
        use_rep: str = "X_pca",
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> ResolutionRecommendation:
        """
        Recommend optimal clustering resolution based on stability.
        """
        cfg = self.config
        resolutions = cfg.resolution_search_space

        log.info(f"Testing {len(resolutions)} resolutions for stability...")

        stability_scores = []
        n_clusters_list = []

        for res in resolutions:
            try:
                # Compute neighbors if needed
                if "neighbors" not in adata.uns:
                    sc.pp.neighbors(adata, use_rep=use_rep)

                # Run clustering multiple times for stability
                labels_list = []
                for seed in range(cfg.resolution_stability_n):
                    sc.tl.leiden(
                        adata,
                        resolution=res,
                        random_state=seed,
                        key_added=f"leiden_{seed}",
                    )
                    labels_list.append(adata.obs[f"leiden_{seed}"].values)

                # Calculate stability (mean pairwise ARI)
                from sklearn.metrics import adjusted_rand_score

                aris = []
                for i in range(len(labels_list)):
                    for j in range(i + 1, len(labels_list)):
                        ari = adjusted_rand_score(labels_list[i], labels_list[j])
                        aris.append(ari)

                stability = np.mean(aris) if aris else 0.0

                # Calculate silhouette
                n_clusters = adata.obs["leiden_0"].nunique()
                if n_clusters > 1:
                    silhouette = silhouette_score(
                        adata.obsm[use_rep],
                        adata.obs["leiden_0"].cat.codes,
                        sample_size=min(10000, adata.n_obs),
                    )
                else:
                    silhouette = -1.0

                # Combined score
                combined = stability * max(0, silhouette)

                stability_scores.append(combined)
                n_clusters_list.append(n_clusters)

                # Clean up
                for seed in range(cfg.resolution_stability_n):
                    del adata.obs[f"leiden_{seed}"]

            except Exception as e:
                log.warning(f"Failed for resolution {res}: {e}")
                stability_scores.append(0.0)
                n_clusters_list.append(0)

        if not stability_scores or max(stability_scores) == 0:
            # Fallback
            return ResolutionRecommendation(
                resolution=1.0,
                n_clusters=10,
                stability_score=0.0,
                ci_lower=0.8,
                ci_upper=1.2,
                method="fallback_default",
                confidence=0.5,
            )

        # Find best resolution
        best_idx = np.argmax(stability_scores)
        best_res = resolutions[best_idx]
        best_stability = stability_scores[best_idx]
        best_n_clusters = n_clusters_list[best_idx]

        # Simple CI estimation
        threshold = best_stability * 0.9
        within_threshold = [
            r for r, s in zip(resolutions, stability_scores) if s >= threshold
        ]
        ci_lower = min(within_threshold) if within_threshold else best_res * 0.8
        ci_upper = max(within_threshold) if within_threshold else best_res * 1.2

        # Confidence based on stability
        confidence = min(1.0, best_stability + 0.3)

        if plot and save_dir:
            self._plot_resolution_analysis(
                resolutions, stability_scores, n_clusters_list, best_res, save_dir
            )

        return ResolutionRecommendation(
            resolution=float(best_res),
            n_clusters=int(best_n_clusters),
            stability_score=float(best_stability),
            ci_lower=float(ci_lower),
            ci_upper=float(ci_upper),
            method="stability_analysis",
            confidence=confidence,
            evidence={
                "resolutions_tested": resolutions,
                "stability_scores": stability_scores,
                "n_clusters": n_clusters_list,
            },
        )

    def assess_batch_effects(
        self,
        adata: AnnData,
        batch_key: str,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> BatchCorrectionRecommendation:
        """
        Assess batch effect severity and recommend correction method.

        Simplified implementation - can be enhanced with actual kBET/LISI calculation.
        """
        cfg = self.config

        log.info(f"Assessing batch effects using '{batch_key}'...")

        n_batches = adata.obs[batch_key].nunique()

        if n_batches < 2:
            return BatchCorrectionRecommendation(
                needs_correction=False,
                severity_score=0.0,
                recommended_method=None,
                confidence=1.0,
                evidence={"reason": "Only one batch present"},
            )

        # Simplified batch assessment using PCA variance
        try:
            if "X_pca" not in adata.obsm:
                sc.tl.pca(adata, n_comps=min(50, adata.n_vars - 1))

            # Simple heuristic: ratio of batch variance to total variance
            from sklearn.preprocessing import LabelEncoder

            le = LabelEncoder()
            batch_encoded = le.fit_transform(adata.obs[batch_key])

            # Fit linear model: PCs ~ batch
            from sklearn.linear_model import LinearRegression
            from sklearn.metrics import r2_score

            r2_scores = []
            for i in range(min(10, adata.obsm["X_pca"].shape[1])):
                reg = LinearRegression()
                reg.fit(batch_encoded.reshape(-1, 1), adata.obsm["X_pca"][:, i])
                y_pred = reg.predict(batch_encoded.reshape(-1, 1))
                r2 = r2_score(adata.obsm["X_pca"][:, i], y_pred)
                r2_scores.append(r2)

            mean_r2 = np.mean(r2_scores)
            severity = min(1.0, mean_r2 * 2)  # Scale to 0-1

            needs_correction = severity > cfg.batch_effect_threshold

            # Method recommendation
            if needs_correction:
                if adata.n_obs > 50000:
                    recommended = "scvi"
                elif n_batches > 10:
                    recommended = "harmony"
                elif self._data_profile and self._data_profile.is_sparse:
                    recommended = "scanorama"
                else:
                    recommended = "harmony"

                alternatives = ["harmony", "scanorama", "bbknn"]
                if recommended in alternatives:
                    alternatives.remove(recommended)
            else:
                recommended = None
                alternatives = []

            confidence = min(1.0, 0.7 + severity)

            return BatchCorrectionRecommendation(
                needs_correction=needs_correction,
                severity_score=float(severity),
                recommended_method=recommended,
                alternative_methods=alternatives[:2],
                method_scores={"pcr_r2": float(mean_r2)},
                confidence=confidence,
                evidence={
                    "n_batches": n_batches,
                    "batch_key": batch_key,
                    "pc_r2_scores": r2_scores,
                    "mean_r2": float(mean_r2),
                },
            )

        except Exception as e:
            log.warning(f"Batch assessment failed: {e}")
            return BatchCorrectionRecommendation(
                needs_correction=False,
                severity_score=0.0,
                recommended_method=None,
                confidence=0.0,
                evidence={"batch_key": batch_key, "error": str(e)},
            )

    # --- Helper methods ---

    def _prepare_data(self, adata: AnnData, n_top_genes: int) -> AnnData:
        """Prepare data for downstream analysis."""
        adata_temp = adata.copy()

        # Basic preprocessing for analysis
        sc.pp.highly_variable_genes(
            adata_temp,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
        )
        adata_temp = adata_temp[:, adata_temp.var.highly_variable].copy()
        sc.pp.scale(adata_temp)

        return adata_temp

    def _find_elbow(self, x: np.ndarray, y: np.ndarray) -> int:
        """Find elbow point using maximum curvature."""
        # Simple implementation - can be enhanced with Kneedle algorithm
        diffs = np.diff(y)
        diffs2 = np.diff(diffs)

        if len(diffs2) > 0:
            elbow_idx = np.argmax(diffs2) + 1
            return int(x[min(elbow_idx, len(x) - 1)])
        return int(x[len(x) // 2])

    def _find_knee_point(self, variance_ratios: np.ndarray) -> int:
        """Find knee point in variance ratio curve."""
        # Find point where variance drops below mean of remaining
        cumsum = np.cumsum(variance_ratios)
        total = cumsum[-1]

        for i, v in enumerate(cumsum):
            if v / total > 0.9:
                return i + 1
        return len(variance_ratios) // 2

    def _bootstrap_hvg(
        self, adata: AnnData, optimal_n: int, n_bootstrap: int, confidence: float
    ) -> Tuple[int, int]:
        """Bootstrap confidence interval for HVG selection."""
        estimates = []

        for _ in range(n_bootstrap):
            # Sample cells
            indices = np.random.choice(adata.n_obs, size=adata.n_obs, replace=True)
            adata_boot = adata[indices].copy()

            try:
                # Quick variance estimate
                sc.pp.highly_variable_genes(
                    adata_boot, n_top_genes=optimal_n, flavor="seurat_v3"
                )
                n_hvg = adata_boot.var.highly_variable.sum()
                estimates.append(n_hvg)
            except Exception:
                estimates.append(optimal_n)

        alpha = (1 - confidence) / 2
        ci_lower = int(np.percentile(estimates, alpha * 100))
        ci_upper = int(np.percentile(estimates, (1 - alpha) * 100))

        return ci_lower, ci_upper

    def _bootstrap_pca(
        self, adata: AnnData, optimal_n: int, n_bootstrap: int, confidence: float
    ) -> Tuple[int, int]:
        """Bootstrap confidence interval for PCA dimensions."""
        # Simplified - return bounds based on optimal_n
        ci_lower = max(2, int(optimal_n * 0.8))
        ci_upper = int(optimal_n * 1.2)
        return ci_lower, ci_upper

    def _calculate_pca_confidence(
        self, variance_ratios: np.ndarray, n_pcs: int
    ) -> float:
        """Calculate confidence in PCA recommendation."""
        if n_pcs >= len(variance_ratios):
            return 0.5

        # Confidence based on variance drop-off
        remaining_var = 1 - np.sum(variance_ratios[:n_pcs])
        confidence = min(1.0, 1.0 - remaining_var + 0.2)
        return confidence

    def _generate_recommendations(
        self,
        hvg_rec: HVGRecommendation,
        pca_rec: PCARecommendation,
        neighbors_rec: NeighborsRecommendation,
        resolution_rec: ResolutionRecommendation,
        batch_rec: Optional[BatchCorrectionRecommendation],
    ) -> List[str]:
        """Generate human-readable recommendations."""
        recs = []

        if hvg_rec.confidence > 0.8:
            recs.append(f"Use {hvg_rec.n_top_genes} highly variable genes")
        else:
            recs.append(f"Consider testing HVG range {hvg_rec.ci_lower}-{hvg_rec.ci_upper}")

        recs.append(f"Use {pca_rec.n_pcs} principal components")
        recs.append(f"Set n_neighbors={neighbors_rec.n_neighbors} for graph construction")
        recs.append(f"Use resolution={resolution_rec.resolution:.1f} for clustering")

        if batch_rec and batch_rec.needs_correction:
            recs.append(f"Apply {batch_rec.recommended_method} for batch correction")

        return recs

    def _plot_hvg_analysis(
        self, df: pd.DataFrame, optimal_n: int, threshold: float, save_dir: Path
    ):
        """Plot HVG variance analysis."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(df["n_genes"], df["variance_explained"], "b-o", label="Variance Explained")
        ax.axhline(y=threshold, color="r", linestyle="--", label=f"Threshold ({threshold:.0%})")
        ax.axvline(x=optimal_n, color="g", linestyle="--", label=f"Optimal: {optimal_n}")
        ax.set_xlabel("Number of HVGs")
        ax.set_ylabel("Cumulative Variance Explained")
        ax.set_title("HVG Selection Analysis")
        ax.legend()
        ax.grid(True, alpha=0.3)

        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "hvg_analysis.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _plot_pca_analysis(
        self,
        variance_ratios: np.ndarray,
        cumulative_var: np.ndarray,
        n_pcs: int,
        save_dir: Path,
    ):
        """Plot PCA variance analysis."""
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Scree plot
        ax1.plot(range(1, len(variance_ratios) + 1), variance_ratios, "b-o")
        ax1.axvline(x=n_pcs, color="r", linestyle="--", label=f"Selected: {n_pcs}")
        ax1.set_xlabel("PC")
        ax1.set_ylabel("Variance Ratio")
        ax1.set_title("PCA Scree Plot")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Cumulative variance
        ax2.plot(range(1, len(cumulative_var) + 1), cumulative_var, "g-o")
        ax2.axvline(x=n_pcs, color="r", linestyle="--", label=f"Selected: {n_pcs}")
        ax2.axhline(y=cumulative_var[n_pcs - 1], color="r", linestyle=":")
        ax2.set_xlabel("PC")
        ax2.set_ylabel("Cumulative Variance")
        ax2.set_title("Cumulative Variance Explained")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "pca_analysis.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _plot_resolution_analysis(
        self,
        resolutions: List[float],
        stability_scores: List[float],
        n_clusters: List[int],
        best_res: float,
        save_dir: Path,
    ):
        """Plot resolution stability analysis."""
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(figsize=(10, 6))

        color = "tab:blue"
        ax1.set_xlabel("Resolution")
        ax1.set_ylabel("Stability Score", color=color)
        ax1.plot(resolutions, stability_scores, color=color, marker="o", label="Stability")
        ax1.axvline(x=best_res, color="r", linestyle="--", label=f"Best: {best_res}")
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        color = "tab:orange"
        ax2.set_ylabel("Number of Clusters", color=color)
        ax2.plot(resolutions, n_clusters, color=color, marker="s", linestyle=":")
        ax2.tick_params(axis="y", labelcolor=color)

        ax1.set_title("Clustering Resolution Analysis")
        ax1.grid(True, alpha=0.3)

        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "resolution_analysis.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


# --- Convenience functions ---


def recommend_intelligent_preprocessing(
    adata: AnnData,
    batch_key: Optional[str] = None,
    tissue_type: str = "unknown",
    plot: bool = True,
    save_dir: Optional[Path] = None,
    **config_overrides
) -> PreprocessingStrategy:
    """
    Main entry point for intelligent preprocessing recommendations.

    Parameters
    ----------
    adata : AnnData
        Input data matrix
    batch_key : str, optional
        Column identifying batches
    tissue_type : str, default="unknown"
        Tissue type for context
    plot : bool, default=True
        Generate diagnostic plots
    save_dir : Path, optional
        Directory to save outputs
    **config_overrides
        Override any IntelligentPreprocessConfig parameter

    Returns
    -------
    PreprocessingStrategy with all recommendations

    Example
    -------
    >>> strategy = recommend_intelligent_preprocessing(adata, batch_key="sampleID")
    >>> print(f"Recommended HVGs: {strategy.hvg.n_top_genes}")
    >>> config = strategy.to_config()
    >>> adata = run_preprocessing(adata, config=config)
    """
    config = IntelligentPreprocessConfig(**config_overrides)
    recommender = IntelligentPreprocessRecommender(config=config)
    return recommender.recommend(
        adata, batch_key=batch_key, tissue_type=tissue_type, plot=plot, save_dir=save_dir
    )


def run_intelligent_preprocessing(
    adata: AnnData,
    batch_key: Optional[str] = None,
    apply_recommendations: bool = True,
    save_dir: Optional[str] = None,
    **kwargs
):
    """
    One-step intelligent preprocessing with automatic parameter selection.

    Parameters
    ----------
    adata : AnnData
        Input data
    batch_key : str, optional
        Batch identifier column
    apply_recommendations : bool, default=True
        If True, apply recommendations and return processed AnnData
    save_dir : str, optional
        Directory for outputs
    **kwargs
        Additional parameters for recommendation

    Returns
    -------
    AnnData or Tuple[AnnData, PreprocessingStrategy]
        Processed data and optionally the strategy
    """
    from ..workflow import run_preprocessing

    # Generate recommendations
    strategy = recommend_intelligent_preprocessing(
        adata, batch_key=batch_key, save_dir=Path(save_dir) if save_dir else None, **kwargs
    )

    if not apply_recommendations:
        return None, strategy

    # Apply recommendations
    config = strategy.to_config()
    adata_processed = run_preprocessing(
        adata, config=config, results_dir=save_dir
    )

    adata_processed.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "intelligent_recommendation"
    ] = {
        "batch_key": batch_key,
        "apply_recommendations": apply_recommendations,
        "strategy": strategy.to_dict(),
        "applied_config": config.to_dict() if hasattr(config, "to_dict") else None,
    }

    return adata_processed, strategy
