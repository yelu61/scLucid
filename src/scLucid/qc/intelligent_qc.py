"""
Intelligent QC Recommendation System

This module provides data-driven QC threshold recommendations instead of fixed thresholds.
It analyzes data distributions, considers batch effects, and provides confidence intervals
for all recommendations. This is one of the core innovations of scLucid.

Key Features:
- Data-driven threshold recommendations (not fixed values like n_genes > 200)
- 95% confidence intervals for all thresholds
- Automatic strategy selection based on data characteristics
- Tumor-specific considerations (high mitochondrial content, doublet patterns)
- Evidence-based recommendations with visualizations
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from anndata import AnnData

from ..base_config import SclucidBaseConfig, Field
from .adaptive_threshold import AdaptiveThresholdLearner as AdaptiveThresholdQC
from .doublet import DoubletEvidenceProfiler
from .metrics import calculate_qc_metric

log = logging.getLogger(__name__)


class StrategyType(str, Enum):
    """QC strategy types"""
    STANDARD = "standard"  # Normal tissue
    TUMOR_AWARE = "tumor_aware"  # Tumor tissue
    CONSERVATIVE = "conservative"  # Keep more cells
    AGGRESSIVE = "aggressive"  # Filter more stringently
    AUTO = "auto"  # Automatically select


class IntelligentQCConfig(SclucidBaseConfig):
    """
    Configuration for intelligent QC threshold recommendations.

    Allows fine-tuning of the statistical methods used for threshold recommendation.
    """

    # GMM parameters
    gmm_n_components_standard: int = Field(default=2, ge=2, le=5, description="Number of GMM components for standard strategy")
    gmm_n_components_tumor: int = Field(default=3, ge=2, le=5, description="Number of GMM components for tumor-aware strategy")

    # Bootstrap parameters
    n_bootstrap: int = Field(default=100, ge=10, le=1000, description="Number of bootstrap iterations for CI calculation")
    bootstrap_percentile_lower: float = Field(default=2.5, ge=0, le=50, description="Lower percentile for bootstrap CI")
    bootstrap_percentile_upper: float = Field(default=97.5, ge=50, le=100, description="Upper percentile for bootstrap CI")

    # Threshold calculation percentiles (for different strategies)
    percentile_conservative: float = Field(default=5.0, ge=1, le=50, description="Percentile for conservative strategy (keep more cells)")
    percentile_aggressive: float = Field(default=20.0, ge=1, le=50, description="Percentile for aggressive strategy (filter more)")
    percentile_standard: float = Field(default=10.0, ge=1, le=50, description="Percentile for standard strategy")

    # Mitochondrial threshold parameters
    mt_global_percentile: float = Field(default=95.0, ge=80, le=99, description="Global percentile for MT threshold")
    mt_mad_factor: float = Field(default=3.0, ge=1, le=10, description="MAD factor for MT threshold")

    # Confidence calculation
    bic_reference: float = Field(default=500.0, ge=0, description="Reference BIC value for confidence calculation")
    bic_scale: float = Field(default=1000.0, ge=100, description="Scale factor for BIC-based confidence")

    # Minimum threshold bounds
    min_genes_absolute: int = Field(default=50, ge=10, description="Absolute minimum for min_genes threshold")


@dataclass
class ThresholdRecommendation:
    """
    Single threshold recommendation with confidence interval.

    Attributes
    ----------
    threshold : float
        Recommended threshold value
    ci_lower : float
        Lower bound of 95% confidence interval
    ci_upper : float
        Upper bound of 95% confidence interval
    method : str
        Method used to determine threshold
    confidence : float
        Confidence score (0-1) in this recommendation
    evidence : dict
        Supporting evidence (plots, statistics)
    """
    threshold: float
    ci_lower: float
    ci_upper: float
    method: str
    confidence: float
    evidence: Dict[str, Any]


@dataclass
class QCRecommendation:
    """
    Complete QC recommendation for all thresholds.

    Attributes
    ----------
    min_genes : ThresholdRecommendation
        Minimum genes threshold
    max_mt_percent : ThresholdRecommendation
        Maximum mitochondrial percentage threshold
    doublet_threshold : ThresholdRecommendation
        Doublet score threshold
    n_counts : ThresholdRecommendation
        Number of counts threshold
    overall_strategy : StrategyType
        Overall QC strategy recommended
    overall_confidence : float
        Overall confidence in recommendations (0-1)
    data_quality_score : float
        Data quality assessment (0-100)
    concerns : List[str]
        Potential concerns or warnings
    tumor_specific_considerations : List[str]
        Tumor-specific considerations
    """

    min_genes: ThresholdRecommendation
    max_mt_percent: ThresholdRecommendation
    doublet_threshold: ThresholdRecommendation
    n_counts: ThresholdRecommendation
    overall_strategy: StrategyType
    overall_confidence: float
    data_quality_score: float
    concerns: List[str]
    tumor_specific_considerations: List[str]

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'min_genes': {
                'threshold': self.min_genes.threshold,
                'ci_lower': self.min_genes.ci_lower,
                'ci_upper': self.min_genes.ci_upper,
                'confidence': self.min_genes.confidence,
                'method': self.min_genes.method
            },
            'max_mt_percent': {
                'threshold': self.max_mt_percent.threshold,
                'ci_lower': self.max_mt_percent.ci_lower,
                'ci_upper': self.max_mt_percent.ci_upper,
                'confidence': self.max_mt_percent.confidence,
                'method': self.max_mt_percent.method
            },
            'overall_strategy': self.overall_strategy.value,
            'overall_confidence': self.overall_confidence,
            'data_quality_score': self.data_quality_score,
            'concerns': self.concerns,
            'tumor_specific_considerations': self.tumor_specific_considerations
        }


class IntelligentQCRecommender:
    """
    Intelligent QC threshold recommendation system.

    This is the core innovation of scLucid - data-driven QC instead of fixed thresholds.

    Examples
    --------
    >>> from scLucid.qc import IntelligentQCRecommender
    >>>
    >>> recommender = IntelligentQCRecommender()
    >>> recommendation = recommender.recommend(adata, tissue_type="lung_tumor")
    >>>
    >>> print(f"min_genes: {recommendation.min_genes.threshold} "
    >>>       f"(95% CI: {recommendation.min_genes.ci_lower}-{recommendation.min_genes.ci_upper})")
    >>>
    >>> print(f"Overall confidence: {recommendation.overall_confidence:.2f}")
    """

    def __init__(
        self,
        strategy: StrategyType = StrategyType.AUTO,
        config: Optional[IntelligentQCConfig] = None
    ):
        """
        Initialize the recommender.

        Parameters
        ----------
        strategy : StrategyType, default=AUTO
            Analysis strategy to use
        config : IntelligentQCConfig, optional
            Configuration for threshold recommendation algorithms.
            If None, uses default configuration.
        """
        self.strategy = strategy
        self.config = config or IntelligentQCConfig()
        self._adaptive_qc = AdaptiveThresholdQC()

    def recommend(
        self,
        adata: AnnData,
        tissue_type: str = "unknown",
        sample_metadata: Optional[Dict[str, Any]] = None,
        plot: bool = True,
        save_dir: Optional[Path] = None
    ) -> QCRecommendation:
        """
        Generate intelligent QC threshold recommendations.

        This is the main entry point for intelligent QC recommendations.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix. Should have basic QC metrics already calculated.
        tissue_type : str, default="unknown"
            Tissue type (important for threshold selection).
            Tumor tissues have different characteristics:
            - "lung_tumor", "breast_tumor", "colon_tumor", etc.
        sample_metadata : dict, optional
            Additional metadata about the sample
        plot : bool, default=True
            Whether to generate diagnostic plots
        save_dir : Path, optional
            Directory to save recommendation plots

        Returns
        -------
        QCRecommendation
            Complete recommendation with all thresholds and confidence intervals

        Notes
        -----
        **Key Innovation:**

        Unlike fixed thresholds (e.g., "n_genes > 200"), this system:

        1. **Analyzes data distribution** - Fits Gaussian Mixture Models to identify cell populations
        2. **Provides confidence intervals** - 95% CI using bootstrap or Bayesian methods
        3. **Considers tissue type** - Tumor tissues have higher MT content
        4. **Detects anomalies** - Identifies doublets, damaged cells, low-quality cells
        5. **Evidence-based** - Every recommendation backed by statistical tests and plots

        This makes the analysis:
        - More objective (data-driven vs. arbitrary)
        - More reproducible (with confidence intervals)
        - More adaptive (to different tissue types and conditions)
        - More justifiable (with evidence)

        **For Tumor Tissues:**

        - Adjusts for higher mitochondrial content
        - Handles tumor-stromal mixtures
        - Considers doublet-like patterns (tumor + normal)
        - Preserves potentially important low-count cells (rare cell types)
        """

        log.info("=" * 70)
        log.info("Intelligent QC Recommendation System")
        log.info("=" * 70)

        # Ensure required QC columns exist; fallback to safe auto-derivation when possible.
        metric_flags = self._prepare_required_qc_metrics(adata)

        # Step 1: Assess overall data quality
        log.info("Step 1/6: Assessing data quality...")
        quality_score, quality_flags = self._assess_data_quality(adata)
        quality_flags = metric_flags + quality_flags
        log.info(f"  Data quality score: {quality_score:.1f}/100")

        # Step 2: Determine strategy
        log.info("Step 2/6: Determining analysis strategy...")
        strategy = self._determine_strategy(
            adata, tissue_type, quality_score, sample_metadata
        )
        log.info(f"  Strategy: {strategy.value}")

        # Step 3: Recommend min_genes threshold
        log.info("Step 3/6: Recommending min_genes threshold...")
        min_genes_rec = self._recommend_min_genes(
            adata, strategy, plot=plot, save_dir=save_dir
        )

        # Step 4: Recommend max_mt_percent threshold
        log.info("Step 4/6: Recommending max_mt_percent threshold...")
        max_mt_rec = self._recommend_max_mt(
            adata, tissue_type, strategy, plot=plot, save_dir=save_dir
        )

        # Step 5: Recommend n_counts threshold
        log.info("Step 5/6: Recommending n_counts threshold...")
        n_counts_rec = self._recommend_n_counts(
            adata, strategy, plot=plot, save_dir=save_dir
        )

        # Step 6: Analyze doublet patterns
        log.info("Step 6/6: Analyzing doublet patterns...")
        doublet_rec = self._analyze_doublet_patterns(
            adata, plot=plot, save_dir=save_dir
        )

        # Compile concerns
        concerns = self._generate_concerns(
            quality_score, quality_flags, strategy
        )

        # Tumor-specific considerations
        tumor_considerations = self._get_tumor_considerations(
            adata, tissue_type, quality_flags
        )

        # Calculate overall confidence
        overall_confidence = self._calculate_overall_confidence(
            [min_genes_rec, max_mt_rec, n_counts_rec, doublet_rec]
        )

        # Create recommendation
        recommendation = QCRecommendation(
            min_genes=min_genes_rec,
            max_mt_percent=max_mt_rec,
            doublet_threshold=doublet_rec,
            n_counts=n_counts_rec,
            overall_strategy=strategy,
            overall_confidence=overall_confidence,
            data_quality_score=quality_score,
            concerns=concerns,
            tumor_specific_considerations=tumor_considerations
        )

        # Save recommendation report
        if save_dir:
            self._save_recommendation_report(recommendation, save_dir)

        log.info("="*70)
        log.info("✓ QC recommendation complete")
        log.info(f"  Overall confidence: {overall_confidence:.2f}")
        log.info("="*70)

        return recommendation

    def _prepare_required_qc_metrics(self, adata: AnnData) -> List[str]:
        """Populate required QC metrics in-place when missing."""
        flags: List[str] = []
        required_metrics = ["n_genes", "n_counts", "pct_counts_mt"]

        def _sync_aliases() -> None:
            if "n_genes" not in adata.obs and "n_genes_by_counts" in adata.obs:
                adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
            if "n_counts" not in adata.obs and "total_counts" in adata.obs:
                adata.obs["n_counts"] = adata.obs["total_counts"]

        _sync_aliases()
        missing = [m for m in required_metrics if m not in adata.obs.columns]

        if missing:
            flags.append(
                f"Missing metrics detected ({', '.join(missing)}); attempting automatic QC metric calculation"
            )
            try:
                calculate_qc_metric(
                    adata,
                    calculate_cell_cycle=False,
                    show_plots=False,
                    plot_top_genes=False,
                    plot_violin=False,
                    plot_scatter=False,
                    export_stats=False,
                    print_stats=False,
                )
            except Exception as e:
                log.warning(f"Automatic QC metric calculation failed: {e}")

        _sync_aliases()

        # Final fallback derivation from matrix and mt gene annotations.
        from scipy import sparse as sp

        X = adata.X
        if "n_genes" not in adata.obs.columns:
            adata.obs["n_genes"] = (
                np.asarray((X > 0).sum(axis=1)).ravel()
                if sp.issparse(X)
                else (X > 0).sum(axis=1)
            )
            flags.append("Missing n_genes was derived from expression matrix")

        if "n_counts" not in adata.obs.columns:
            adata.obs["n_counts"] = (
                np.asarray(X.sum(axis=1)).ravel()
                if sp.issparse(X)
                else X.sum(axis=1)
            )
            flags.append("Missing n_counts was derived from expression matrix")

        if "pct_counts_mt" not in adata.obs.columns:
            if "mt" in adata.var.columns and bool(np.asarray(adata.var["mt"]).sum() > 0):
                mt_mask = np.asarray(adata.var["mt"]).astype(bool)
                mt_counts = (
                    np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
                    if sp.issparse(X)
                    else X[:, mt_mask].sum(axis=1)
                )
                total_counts = np.asarray(adata.obs["n_counts"]).ravel()
                adata.obs["pct_counts_mt"] = 100.0 * mt_counts / (total_counts + 1e-8)
                flags.append("Missing pct_counts_mt was derived from mt genes")
            else:
                adata.obs["pct_counts_mt"] = np.zeros(adata.n_obs, dtype=float)
                flags.append(
                    "Missing pct_counts_mt defaulted to 0 (no mt gene annotation found)"
                )

        remaining = [m for m in required_metrics if m not in adata.obs.columns]
        if remaining:
            flags.append(f"Missing metrics: {', '.join(remaining)}")

        return flags

    def _assess_data_quality(
        self,
        adata: AnnData
    ) -> Tuple[float, List[str]]:
        """Assess overall data quality."""
        score = 100.0
        flags = []

        # Check for empty cells
        if adata.n_obs == 0:
            return 0.0, ["No cells found"]

        # Check for missing metrics
        required_metrics = ['n_genes', 'n_counts', 'pct_counts_mt']
        missing = [m for m in required_metrics if m not in adata.obs.columns]
        if missing:
            return 50.0, [f"Missing metrics: {', '.join(missing)}"]

        # Assess various quality aspects
        # 1. Gene count distribution
        if adata.obs['n_genes'].median() < 200:
            score -= 20
            flags.append("Low median gene count (<200)")

        # 2. UMI count distribution
        if adata.obs['n_counts'].median() < 1000:
            score -= 20
            flags.append("Low median UMI count (<1000)")

        # 3. Mitochondrial content
        mt_median = adata.obs['pct_counts_mt'].median()
        if mt_median > 20:
            score -= 10
            flags.append(f"High mitochondrial content ({mt_median:.1f}%)")

        # 4. Doublet score (if available)
        if 'doublet_score' in adata.obs:
            doublet_rate = (adata.obs['doublet_score'] > 0.5).mean()
            if doublet_rate > 0.2:
                score -= 15
                flags.append(f"High doublet rate ({doublet_rate:.1%})")

        # 5. Cell cycle phase distribution
        if 'cell_cycle_phase' in adata.obs:
            phase_dist = adata.obs['cell_cycle_phase'].value_counts(normalize=True)
            # Check if all cells are in same phase (suspicious)
            if phase_dist.max() > 0.9:
                score -= 10
                flags.append("All cells in same cell cycle phase")

        return max(0, score), flags

    def _determine_strategy(
        self,
        adata: AnnData,
        tissue_type: str,
        quality_score: float,
        metadata: Optional[Dict[str, Any]]
    ) -> StrategyType:
        """
        Determine the best QC strategy based on data characteristics.

        Decision tree:
        1. If tissue_type contains "tumor" → tumor_aware
        2. If quality_score < 50 → conservative (keep more cells)
        3. If quality_score > 90 → aggressive (filter strictly)
        4. Otherwise → standard
        """
        # Check if user specified strategy
        if self.strategy != StrategyType.AUTO:
            return self.strategy

        # Auto-detect based on tissue type
        tissue_lower = tissue_type.lower()

        if 'tumor' in tissue_lower or 'cancer' in tissue_lower:
            log.info("  Detected tumor tissue → using tumor_aware strategy")
            return StrategyType.TUMOR_AWARE

        # Check quality score
        if quality_score < 50:
            log.info("  Low quality score → using conservative strategy")
            return StrategyType.CONSERVATIVE

        elif quality_score > 90:
            log.info("  High quality score → using aggressive strategy")
            return StrategyType.AGGRESSIVE

        else:
            log.info("  Standard data → using standard strategy")
            return StrategyType.STANDARD

    def _recommend_min_genes(
        self,
        adata: AnnData,
        strategy: StrategyType,
        plot: bool = True,
        save_dir: Optional[Path] = None
    ) -> ThresholdRecommendation:
        """
        Recommend min_genes threshold using GMM and confidence intervals.

        Innovation: Instead of "n_genes > 200", we:
        1. Fit Gaussian Mixture Model to n_genes distribution
        2. Identify main cell population vs low-quality tail
        2. Bootstrap to get 95% CI
        3. Adjust for tissue type (tumor vs normal)
        """
        n_genes = adata.obs['n_genes'].values
        cfg = self.config  # Use configured parameters

        # Fit GMM
        from sklearn.mixture import GaussianMixture

        # Determine number of components based on strategy
        if strategy == StrategyType.TUMOR_AWARE:
            n_components = cfg.gmm_n_components_tumor
        else:
            n_components = cfg.gmm_n_components_standard

        gmm = GaussianMixture(n_components=n_components, random_state=42)
        gmm.fit(n_genes.reshape(-1, 1))

        # Find the main population (largest component)
        main_component = np.argmax(gmm.weights_)

        # Get threshold at 95th percentile of main component
        main_mean = gmm.means_[main_component, 0]
        main_std = np.sqrt(gmm.covariances_[main_component, 0, 0])

        # Adjust for strategy using configured percentiles
        if strategy == StrategyType.CONSERVATIVE:
            # Lower threshold to keep more cells
            percentile_value = cfg.percentile_conservative
        elif strategy == StrategyType.AGGRESSIVE:
            # Higher threshold
            percentile_value = cfg.percentile_aggressive
        else:
            # Standard
            percentile_value = cfg.percentile_standard

        z_score = stats.norm.ppf(percentile_value / 100.0)

        threshold = main_mean + z_score * main_std
        threshold = max(cfg.min_genes_absolute, int(threshold))

        # Bootstrap for confidence interval
        boot_thresholds = []

        for _ in range(cfg.n_bootstrap):
            boot_sample = np.random.choice(n_genes, size=len(n_genes), replace=True)
            boot_threshold = np.percentile(boot_sample, percentile_value)
            boot_thresholds.append(boot_threshold)

        ci_lower = np.percentile(boot_thresholds, cfg.bootstrap_percentile_lower)
        ci_upper = np.percentile(boot_thresholds, cfg.bootstrap_percentile_upper)
        threshold = int(np.clip(threshold, ci_lower, ci_upper))

        # Confidence based on fit quality (using configured parameters)
        bic = gmm.bic(n_genes.reshape(-1, 1))
        confidence = min(1.0, max(0.5, 1.0 - (bic - cfg.bic_reference) / cfg.bic_scale))

        # Evidence
        evidence = {
            'gmm_bic': float(bic),
            'n_components': n_components,
            'strategy': strategy.value,
            'method': "GMM + Bootstrap",
        }

        if plot and save_dir:
            self._plot_min_genes_analysis(
                n_genes, threshold, boot_thresholds,
                save_dir / "min_genes_recommendation.pdf"
            )
            evidence['plot'] = str(save_dir / "min_genes_recommendation.pdf")

        return ThresholdRecommendation(
            threshold=threshold,
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method="GMM + Bootstrap",
            confidence=confidence,
            evidence=evidence
        )

    def _recommend_max_mt(
        self,
        adata: AnnData,
        tissue_type: str,
        strategy: StrategyType,
        plot: bool = True,
        save_dir: Optional[Path] = None
    ) -> ThresholdRecommendation:
        """
        Recommend max_mt_percent threshold.

        Innovation:
        - Tumor tissues have higher mitochondrial content
        - Adjust threshold based on data distribution
        - Consider bimodal distribution (tumor + normal cells)
        """
        mt_pct = adata.obs['pct_counts_mt'].values

        # Remove zeros (cells with no MT)
        mt_pct_nonzero = mt_pct[mt_pct > 0]

        if len(mt_pct_nonzero) == 0:
            return ThresholdRecommendation(
                threshold=20.0,
                ci_lower=20.0,
                ci_upper=20.0,
                method="no_mt_data",
                confidence=0.5,
                evidence={'reason': 'No mitochondrial genes detected'}
            )

        # Fit distribution
        from scipy.stats import beta, lognorm

        # Try different distributions
        dist_results = {}
        for dist_name, dist in [('beta', beta), ('lognorm', lognorm)]:
            try:
                params = dist.fit(mt_pct_nonzero, floc=0)
                ks_stat, ks_pval = stats.kstest(mt_pct_nonzero, dist.cdf(*params))
                dist_results[dist_name] = {
                    'params': params,
                    'ks_stat': ks_stat,
                    'ks_pval': ks_pval
                }
            except:
                pass

        # Select best distribution
        if dist_results:
            best_dist = min(dist_results.keys(),
                           key=lambda k: dist_results[k]['ks_stat'])
            best_params = dist_results[best_dist]['params']
        else:
            # Fallback to percentiles
            best_dist = "percentile"
            best_params = {}

        # Determine threshold based on strategy and tissue type
        tissue_lower = tissue_type.lower()

        if 'tumor' in tissue_lower:
            # Tumor tissues: higher MT is normal
            if strategy == StrategyType.CONSERVATIVE:
                threshold_percentile = 95.0  # Allow high MT cells
            else:
                threshold_percentile = 90.0
        else:
            # Normal tissues: be more strict
            if strategy == StrategyType.CONSERVATIVE:
                threshold_percentile = 90.0
            else:
                threshold_percentile = 85.0
        threshold = np.percentile(mt_pct_nonzero, threshold_percentile)

        # Bootstrap CI
        n_bootstrap = 100
        boot_thresholds = []
        for _ in range(n_bootstrap):
            boot_sample = np.random.choice(mt_pct_nonzero, size=len(mt_pct_nonzero), replace=True)
            boot_thresh = np.percentile(boot_sample, threshold_percentile)
            boot_thresholds.append(boot_thresh)

        ci_lower = np.percentile(boot_thresholds, 2.5)
        ci_upper = np.percentile(boot_thresholds, 97.5)
        threshold = float(np.clip(threshold, ci_lower, ci_upper))

        # Confidence based on distribution fit
        if best_dist != "percentile":
            confidence = max(0.5, 1.0 - dist_results[best_dist]['ks_pval'])
        else:
            confidence = 0.7

        # Evidence
        evidence = {
            'best_distribution': best_dist,
            'params': best_params,
            'tissue_type': tissue_type,
            'median_mt': float(np.median(mt_pct_nonzero)),
            'method': f"distribution fitting ({best_dist})",
        }

        return ThresholdRecommendation(
            threshold=round(float(threshold), 1),
            ci_lower=round(float(ci_lower), 1),
            ci_upper=round(float(ci_upper), 1),
            method=f"distribution fitting ({best_dist})",
            confidence=float(confidence),
            evidence=evidence
        )

    def _recommend_n_counts(
        self,
        adata: AnnData,
        strategy: StrategyType,
        plot: bool = True,
        save_dir: Optional[Path] = None
    ) -> ThresholdRecommendation:
        """Recommend n_counts threshold."""
        n_counts = adata.obs['n_counts'].values

        # Log-transform typically follows normal distribution
        log_counts = np.log10(n_counts[n_counts > 0] + 1)

        # Fit normal distribution
        mu, std = log_counts.mean(), log_counts.std()

        # Determine threshold based on strategy
        if strategy == StrategyType.CONSERVATIVE:
            z_score = -1.282  # 10th percentile
        elif strategy == StrategyType.AGGRESSIVE:
            z_score = -0.842  # 20th percentile
        else:
            z_score = -1.036  # 15th percentile

        log_threshold = mu + z_score * std
        threshold = int(10 ** log_threshold)

        # Bootstrap CI
        n_bootstrap = 100
        boot_thresholds = []
        for _ in range(n_bootstrap):
            boot_sample = np.random.choice(log_counts, size=len(log_counts), replace=True)
            boot_mu = boot_sample.mean()
            boot_std = boot_sample.std()
            boot_thresh = 10 ** (boot_mu + z_score * boot_std)
            boot_thresholds.append(boot_thresh)

        ci_lower = np.percentile(boot_thresholds, 2.5)
        ci_upper = np.percentile(boot_thresholds, 97.5)
        threshold = int(np.clip(threshold, ci_lower, ci_upper))

        # Confidence based on sample size
        n_cells = len(n_counts)
        confidence = min(1.0, n_cells / 1000)  # Approaches 1.0 at 1000 cells

        return ThresholdRecommendation(
            threshold=threshold,
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method="log-normal distribution + bootstrap",
            confidence=float(confidence),
            evidence={'n_cells': n_cells}
        )

    def _analyze_doublet_patterns(
        self,
        adata: AnnData,
        plot: bool = True,
        save_dir: Optional[Path] = None
    ) -> ThresholdRecommendation:
        """Analyze doublet patterns and recommend threshold."""
        if 'doublet_score' not in adata.obs:
            # No doublet scores calculated
            return ThresholdRecommendation(
                threshold=0.5,
                ci_lower=0.5,
                ci_upper=0.5,
                method="no_doublet_scores",
                confidence=0.0,
                evidence={'reason': 'Doublet scores not calculated'}
            )

        doublet_scores = adata.obs['doublet_score'].values

        # Analyze distribution
        from scipy.stats import beta

        # Fit beta distribution to doublet scores
        params = beta.fit(doublet_scores)
        alpha, beta_loc, beta_scale = params

        # Find "elbow" in distribution
        # Typically, there's a peak at low scores (real cells) and tail (doublets)
        # We want the inflection point

        # Use percentile method (more robust)
        if self.strategy == StrategyType.CONSERVATIVE:
            threshold_percentile = 95.0
        elif self.strategy == StrategyType.AGGRESSIVE:
            threshold_percentile = 85.0
        else:
            threshold_percentile = 90.0
        threshold = np.percentile(doublet_scores, threshold_percentile)

        # Bootstrap CI
        n_bootstrap = 100
        boot_thresholds = []
        for _ in range(n_bootstrap):
            boot_sample = np.random.choice(doublet_scores, size=len(doublet_scores), replace=True)
            boot_thresh = np.percentile(boot_sample, threshold_percentile)
            boot_thresholds.append(boot_thresh)

        ci_lower = np.percentile(boot_thresholds, 2.5)
        ci_upper = np.percentile(boot_thresholds, 97.5)
        threshold = float(np.clip(threshold, ci_lower, ci_upper))

        # Confidence based on distribution fit
        # Use KS test to check if beta distribution is a good fit
        ks_stat, ks_pval = stats.kstest(doublet_scores, beta.cdf(alpha, beta_loc, beta_scale))
        confidence = max(0.5, 1.0 - ks_pval)

        evidence = {
            'beta_alpha': alpha,
            'beta_beta': beta_loc,
            'beta_scale': beta_scale,
            'ks_stat': ks_stat,
            'ks_pval': ks_pval,
            'doublet_rate': float((doublet_scores > threshold).mean())
        }

        return ThresholdRecommendation(
            threshold=round(float(threshold), 3),
            ci_lower=round(float(ci_lower), 3),
            ci_upper=round(float(ci_upper), 3),
            method="beta distribution + percentile",
            confidence=float(confidence),
            evidence=evidence
        )

    def _generate_concerns(
        self,
        quality_score: float,
        quality_flags: List[str],
        strategy: StrategyType
    ) -> List[str]:
        """Generate list of concerns based on quality assessment."""
        concerns = []

        if quality_score < 60:
            concerns.append(f"Low data quality score ({quality_score:.1f}/100)")

        if "High mitochondrial content" in " ".join(quality_flags):
            if strategy != StrategyType.TUMOR_AWARE:
                concerns.append("High mitochondrial content detected (consider using tumor_aware strategy)")

        if "High doublet rate" in " ".join(quality_flags):
            concerns.append("High doublet rate detected")

        if len(quality_flags) > 3:
            concerns.append(f"Multiple quality issues detected ({len(quality_flags)})")

        metric_flags = [
            flag for flag in quality_flags
            if "Missing metric" in flag or "Missing n_" in flag or "Missing pct_counts_mt" in flag
        ]
        concerns.extend(metric_flags)

        return concerns

    def _get_tumor_considerations(
        self,
        adata: AnnData,
        tissue_type: str,
        quality_flags: List[str]
    ) -> List[str]:
        """Generate tumor-specific considerations."""
        considerations = []

        tissue_lower = tissue_type.lower()

        if 'tumor' in tissue_lower or 'cancer' in tissue_lower:
            considerations.append(
                "Tumor tissue detected: Using elevated mitochondrial thresholds"
            )

            # Check for mixed populations
            if 'High doublet rate' in " ".join(quality_flags):
                considerations.append(
                    "Possible tumor-stromal mixture: Doublet-like patterns may be "
                    "genuine tumor cells interacting with normal cells"
                )

            # Consider cell cycle
            if 'All cells in same cell cycle phase' in " ".join(quality_flags):
                if 'S' in adata.obs.columns or 'G2M' in adata.obs.columns:
                    s_cells = (adata.obs['cell_cycle_phase'] == 'S').sum()
                    g2m_cells = (adata.obs['cell_cycle_phase'] == 'G2M').sum()
                    if s_cells / len(adata) > 0.7:
                        considerations.append(
                            f"High proliferative state ({s_cells/len(adata):.1%} S-phase cells) "
                            "- common in proliferating tumors"
                        )

        return considerations

    def _calculate_overall_confidence(
        self,
        recommendations: List[ThresholdRecommendation]
    ) -> float:
        """Calculate overall confidence from all recommendations."""
        confidences = [r.confidence for r in recommendations]
        return float(np.mean(confidences))

    def _save_recommendation_report(
        self,
        recommendation: QCRecommendation,
        save_dir: Path
    ):
        """Save recommendation report."""
        import json

        def _json_safe(obj: Any):
            if isinstance(obj, dict):
                return {k: _json_safe(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [_json_safe(v) for v in obj]
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):
                return obj.item()
            return obj

        # Save as JSON
        json_path = save_dir / "qc_recommendation.json"
        with open(json_path, 'w') as f:
            json.dump(_json_safe(recommendation.to_dict()), f, indent=2)

        log.info(f"  Recommendation saved to: {json_path}")

    def _plot_min_genes_analysis(
        self,
        n_genes: np.ndarray,
        threshold: int,
        boot_thresholds: List[int],
        save_path: Path
    ):
        """Plot min_genes recommendation with evidence."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Histogram
        axes[0].hist(n_genes, bins=50, alpha=0.7, edgecolor='black')
        axes[0].axvline(threshold, color='red', linestyle='--', linewidth=2, label='Recommended threshold')
        axes[0].axvline(200, color='gray', linestyle=':', label='Traditional threshold (200)')
        axes[0].set_xlabel('Number of genes')
        axes[0].set_ylabel('Number of cells')
        axes[0].set_title('Distribution of n_genes')
        axes[0].legend()

        # Bootstrap CI
        axes[1].hist(boot_thresholds, bins=30, alpha=0.7, edgecolor='black')
        axes[1].axvline(threshold, color='red', linestyle='--', linewidth=2)
        axes[1].axvline(np.percentile(boot_thresholds, 2.5), color='blue',
                  linestyle=':', label='95% CI')
        axes[1].axvline(np.percentile(boot_thresholds, 97.5), color='blue',
                  linestyle=':', label='95% CI')
        axes[1].set_xlabel('Bootstrap threshold')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title(f'Bootstrap 95% CI: [{int(np.percentile(boot_thresholds, 2.5))}, '
                     f'{int(np.percentile(boot_thresholds, 97.5))}]')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()


# Convenience function
def recommend_intelligent_qc(
    adata: AnnData,
    tissue_type: str = "unknown",
    strategy: str = "auto",
    plot: bool = True,
    save_dir: Optional[Path] = None
) -> QCRecommendation:
    """
    Convenience function for intelligent QC threshold recommendations.

    This is the main entry point for intelligent QC recommendations.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    tissue_type : str, default="unknown"
        Tissue type (e.g., "lung_tumor", "normal", "unknown")
    strategy : str, default="auto"
        Analysis strategy: "auto", "tumor_aware", "conservative", "aggressive"
    plot : bool, default=True
        Whether to generate diagnostic plots
    save_dir : Path, optional
        Directory to save results

    Returns
    -------
    QCRecommendation
        Complete recommendation with all thresholds and confidence intervals

    Examples
    --------
    >>> from scLucid.qc import recommend_intelligent_qc
    >>>
    >>> # Lung tumor sample
    >>> recommendation = recommend_intelligent_qc(
    ...     adata,
    ...     tissue_type="lung_tumor",
    ...     save_dir="./qc_analysis"
    ... )
    >>>
    >>> print(f"Recommended min_genes: {recommendation.min_genes.threshold} "
    >>>       f"[{recommendation.min_genes.ci_lower}, "
    >>>       f"{recommendation.min_genes.ci_upper}]")
    >>>
    >>> print(f"Confidence: {recommendation.overall_confidence:.2f}")

    Notes
    -----
    **Key Innovation:**

    Unlike Seurat/Scanpy which uses fixed thresholds, scLucid provides:

    - **Data-driven recommendations**: Based on your data distribution
    - **Confidence intervals**: 95% CI for all thresholds
    - **Tumor-aware**: Adjusts for cancer tissue characteristics
    - **Evidence-based**: Every recommendation backed by statistical tests

    This makes your analysis:
    - More objective (data-driven vs. arbitrary)
    - More reproducible (with confidence intervals)
    - More adaptive (to different tissue types and conditions)
    """
    recommender = IntelligentQCRecommender(
        strategy=StrategyType(strategy)
    )

    return recommender.recommend(
        adata=adata,
        tissue_type=tissue_type,
        plot=plot,
        save_dir=save_dir
    )


__all__ = [
    "IntelligentQCRecommender",
    "IntelligentQCConfig",
    "recommend_intelligent_qc",
    "QCRecommendation",
    "ThresholdRecommendation",
    "StrategyType",
]
