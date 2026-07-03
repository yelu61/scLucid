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
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from ...base_config import Field, SclucidBaseConfig
from ...utils.context import is_tumor_context, resolve_cell_type_key
from .adaptive_threshold import (
    AdaptiveThresholdLearner as AdaptiveThresholdQC,
)
from .adaptive_threshold import (
    fit_bimodal_gmm_threshold_model,
    fit_count_mixture_threshold_model,
    fit_gmm_threshold_model,
)
from ..metrics import calculate_qc_metric

log = logging.getLogger(__name__)


def _make_json_safe(obj: Any) -> Any:
    """Recursively convert numpy/tuple types to plain Python for JSON/HDF5."""
    if isinstance(obj, dict):
        return {k: _make_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_make_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


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
    gmm_n_components_standard: int = Field(
        default=2, ge=2, le=5, description="Number of GMM components for standard strategy"
    )
    gmm_n_components_tumor: int = Field(
        default=3, ge=2, le=5, description="Number of GMM components for tumor-aware strategy"
    )

    # Count mixture parameters (for n_genes / count metrics)
    min_genes_model: Literal["auto", "zinb", "poisson_gamma", "gmm", "percentile"] = Field(
        default="auto",
        description="Distribution family for min_genes threshold. auto tries ZINB/NB and falls back to GMM/percentile.",
    )
    min_genes_percentile: float = Field(
        default=10.0, ge=1, le=50, description="Percentile for count-model min_genes threshold"
    )

    # Mitochondrial threshold parameters
    mt_model: Literal[
        "auto",
        "bimodal_gmm",
        "beta",
        "lognorm",
        "percentile",
        "sample_aware",
        "celltype_aware",
        "multicomponent",
    ] = Field(
        default="auto",
        description=(
            "Distribution family for max_mt_percent threshold. auto uses bimodal "
            "GMM for tumor tissues when supported; sample_aware / celltype_aware "
            "compute per-stratum baselines."
        ),
    )
    mt_global_percentile: float = Field(
        default=95.0, ge=80, le=99, description="Global percentile for MT threshold"
    )
    mt_mad_factor: float = Field(
        default=3.0, ge=1, le=10, description="MAD factor for MT threshold"
    )
    mt_review_mad_factor: float = Field(
        default=2.0,
        ge=1.0,
        le=5.0,
        description="MAD factor defining the MT% review band above the stratum baseline.",
    )
    mt_bimodal_min_separation: float = Field(
        default=2.0,
        ge=1.0,
        le=5.0,
        description="Minimum component separation (in std) to accept bimodal GMM for MT%.",
    )
    mt_max_components: int = Field(
        default=4,
        ge=2,
        le=6,
        description="Maximum GMM components for tumor-aware MT% modelling.",
    )
    mt_tumor_programs: Tuple[str, ...] = Field(
        default_factory=lambda: (
            "epithelial_malignant_like",
            "cell_cycle",
            "hypoxia_stress",
            "emt_stromal",
            "oxphos",
            "glycolysis",
            "mt_biogenesis",
        ),
        description="Tumor program panels used to interpret high-MT cells.",
    )

    # Stratified MT% analysis
    sample_key: str = Field(
        default="sampleID",
        description="Observation key used for sample-aware MT% baselines.",
    )

    # Bootstrap parameters
    n_bootstrap: int = Field(
        default=200, ge=50, le=2000, description="Number of bootstrap iterations for CI calculation"
    )
    bootstrap_percentile_lower: float = Field(
        default=2.5, ge=0, le=50, description="Lower percentile for bootstrap CI"
    )
    bootstrap_percentile_upper: float = Field(
        default=97.5, ge=50, le=100, description="Upper percentile for bootstrap CI"
    )

    # Data-quality scoring guardrails. These are heuristic defaults and are
    # recorded in recommendation evidence so projects can calibrate them.
    quality_min_median_genes: int = Field(default=200, ge=0)
    quality_min_median_counts: int = Field(default=1000, ge=0)
    quality_high_mt_median: float = Field(default=20.0, ge=0, le=100)
    quality_doublet_score_threshold: float = Field(default=0.5, ge=0, le=1)
    quality_high_doublet_rate: float = Field(default=0.2, ge=0, le=1)
    quality_dominant_cell_cycle_phase_fraction: float = Field(default=0.9, ge=0, le=1)
    quality_penalty_low_genes: float = Field(default=20.0, ge=0, le=100)
    quality_penalty_low_counts: float = Field(default=20.0, ge=0, le=100)
    quality_penalty_high_mt: float = Field(default=10.0, ge=0, le=100)
    quality_penalty_high_doublet: float = Field(default=15.0, ge=0, le=100)
    quality_penalty_cell_cycle_dominance: float = Field(default=10.0, ge=0, le=100)

    tumor_program_min_relative_retention: float = Field(default=0.3, ge=0, le=10)
    tumor_program_min_present_fraction: float = Field(default=0.5, ge=0, le=1)

    # Threshold calculation percentiles (for different strategies)
    percentile_conservative: float = Field(
        default=5.0,
        ge=1,
        le=50,
        description="Percentile for conservative strategy (keep more cells)",
    )
    percentile_aggressive: float = Field(
        default=20.0, ge=1, le=50, description="Percentile for aggressive strategy (filter more)"
    )
    percentile_standard: float = Field(
        default=10.0, ge=1, le=50, description="Percentile for standard strategy"
    )

    # Confidence calculation
    bic_reference: float = Field(
        default=500.0, ge=0, description="Reference BIC value for confidence calculation"
    )
    bic_scale: float = Field(
        default=1000.0, ge=100, description="Scale factor for BIC-based confidence"
    )

    # Minimum threshold bounds
    min_genes_absolute: int = Field(
        default=50, ge=10, description="Absolute minimum for min_genes threshold"
    )
    random_state: int = Field(
        default=42, description="Random seed used for bootstrap confidence intervals"
    )


@dataclass
class ThresholdRecommendation:
    """
    Single threshold recommendation with confidence interval.

    Attributes:
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

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a JSON-serializable dictionary."""
        return _make_json_safe({
            "threshold": self.threshold,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "confidence": self.confidence,
            "method": self.method,
            "evidence": self.evidence,
        })


@dataclass
class QCRecommendation:
    """
    Complete QC recommendation for all thresholds.

    Attributes:
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
            "min_genes": self.min_genes.to_dict(),
            "max_mt_percent": self.max_mt_percent.to_dict(),
            "doublet_threshold": self.doublet_threshold.to_dict(),
            "n_counts": self.n_counts.to_dict(),
            "overall_strategy": self.overall_strategy.value,
            "overall_confidence": self.overall_confidence,
            "data_quality_score": self.data_quality_score,
            "concerns": self.concerns,
            "tumor_specific_considerations": self.tumor_specific_considerations,
        }


class IntelligentQCRecommender:
    """
    Intelligent QC threshold recommendation system.

    This is the core innovation of scLucid - data-driven QC instead of fixed thresholds.

    Examples:
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
        config: Optional[IntelligentQCConfig] = None,
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
        save_dir: Optional[Path] = None,
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

        Returns:
        -------
        QCRecommendation
            Complete recommendation with all thresholds and confidence intervals

        Notes:
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

        if save_dir:
            save_dir = Path(save_dir)

        # Recommendation should not mutate the caller's AnnData; compute any
        # missing QC helpers on a working copy instead.
        required_aliases_present = all(
            col in adata.obs.columns for col in ["n_genes", "n_counts", "pct_counts_mt"]
        )
        adata_for_recommendation = adata if required_aliases_present else adata.copy()

        # Ensure required QC columns exist; fallback to safe auto-derivation when possible.
        metric_flags = self._prepare_required_qc_metrics(adata_for_recommendation)

        # Step 1: Assess overall data quality
        log.info("Step 1/6: Assessing data quality...")
        quality_score, quality_flags = self._assess_data_quality(adata_for_recommendation)
        quality_flags = metric_flags + quality_flags
        log.info(f"  Data quality score: {quality_score:.1f}/100")

        # Step 2: Determine strategy
        log.info("Step 2/6: Determining analysis strategy...")
        strategy = self._determine_strategy(
            adata_for_recommendation, tissue_type, quality_score, sample_metadata
        )
        log.info(f"  Strategy: {strategy.value}")

        # Step 3: Recommend min_genes threshold
        log.info("Step 3/6: Recommending min_genes threshold...")
        min_genes_rec = self._recommend_min_genes(
            adata_for_recommendation, strategy, plot=plot, save_dir=save_dir
        )

        # Step 4: Recommend max_mt_percent threshold
        log.info("Step 4/6: Recommending max_mt_percent threshold...")
        max_mt_rec = self._recommend_max_mt(
            adata_for_recommendation,
            tissue_type,
            strategy,
            plot=plot,
            save_dir=save_dir,
        )

        # Step 5: Recommend n_counts threshold
        log.info("Step 5/6: Recommending n_counts threshold...")
        n_counts_rec = self._recommend_n_counts(
            adata_for_recommendation, strategy, plot=plot, save_dir=save_dir
        )

        # Step 6: Analyze doublet patterns
        log.info("Step 6/6: Analyzing doublet patterns...")
        doublet_rec = self._analyze_doublet_patterns(
            adata_for_recommendation, plot=plot, save_dir=save_dir
        )

        # Compile concerns
        concerns = self._generate_concerns(quality_score, quality_flags, strategy)

        # Tumor-specific considerations
        tumor_considerations = self._get_tumor_considerations(
            adata_for_recommendation, tissue_type, quality_flags
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
            tumor_specific_considerations=tumor_considerations,
        )

        # Save recommendation report
        if save_dir:
            self._save_recommendation_report(recommendation, save_dir)

        log.info("=" * 70)
        log.info("✓ QC recommendation complete")
        log.info(f"  Overall confidence: {overall_confidence:.2f}")
        log.info("=" * 70)

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
                np.asarray((X > 0).sum(axis=1)).ravel() if sp.issparse(X) else (X > 0).sum(axis=1)
            )
            flags.append("Missing n_genes was derived from expression matrix")

        if "n_counts" not in adata.obs.columns:
            adata.obs["n_counts"] = (
                np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else X.sum(axis=1)
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
                flags.append("Missing pct_counts_mt defaulted to 0 (no mt gene annotation found)")

        remaining = [m for m in required_metrics if m not in adata.obs.columns]
        if remaining:
            flags.append(f"Missing metrics: {', '.join(remaining)}")

        return flags

    def _assess_data_quality(self, adata: AnnData) -> Tuple[float, List[str]]:
        """Assess overall data quality."""
        cfg = self.config
        score = 100.0
        flags = []

        # Check for empty cells
        if adata.n_obs == 0:
            return 0.0, ["No cells found"]

        # Check for missing metrics
        required_metrics = ["n_genes", "n_counts", "pct_counts_mt"]
        missing = [m for m in required_metrics if m not in adata.obs.columns]
        if missing:
            return 50.0, [f"Missing metrics: {', '.join(missing)}"]

        # Assess various quality aspects
        # 1. Gene count distribution
        if adata.obs["n_genes"].median() < cfg.quality_min_median_genes:
            score -= cfg.quality_penalty_low_genes
            flags.append(f"Low median gene count (<{cfg.quality_min_median_genes})")

        # 2. UMI count distribution
        if adata.obs["n_counts"].median() < cfg.quality_min_median_counts:
            score -= cfg.quality_penalty_low_counts
            flags.append(f"Low median UMI count (<{cfg.quality_min_median_counts})")

        # 3. Mitochondrial content
        mt_median = adata.obs["pct_counts_mt"].median()
        if mt_median > cfg.quality_high_mt_median:
            score -= cfg.quality_penalty_high_mt
            flags.append(f"High mitochondrial content ({mt_median:.1f}%)")

        # 4. Doublet score (if available)
        if "doublet_score" in adata.obs:
            doublet_rate = (adata.obs["doublet_score"] > cfg.quality_doublet_score_threshold).mean()
            if doublet_rate > cfg.quality_high_doublet_rate:
                score -= cfg.quality_penalty_high_doublet
                flags.append(f"High doublet rate ({doublet_rate:.1%})")

        # 5. Cell cycle phase distribution
        if "cell_cycle_phase" in adata.obs:
            phase_dist = adata.obs["cell_cycle_phase"].value_counts(normalize=True)
            # Check if all cells are in same phase (suspicious)
            if phase_dist.max() > cfg.quality_dominant_cell_cycle_phase_fraction:
                score -= cfg.quality_penalty_cell_cycle_dominance
                flags.append("All cells in same cell cycle phase")

        return max(0, score), flags

    def _determine_strategy(
        self,
        adata: AnnData,
        tissue_type: str,
        quality_score: float,
        metadata: Optional[Dict[str, Any]],
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
        if is_tumor_context(tissue_type):
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
        save_dir: Optional[Path] = None,
    ) -> ThresholdRecommendation:
        """
        Recommend min_genes threshold using count mixture models and confidence intervals.

        Uses Zero-Inflated Negative Binomial / Negative Binomial (Poisson-Gamma)
        when possible because n_genes is a discrete, right-skewed count metric.
        Falls back to GMM or percentile for small / pathological datasets.
        """
        n_genes = adata.obs["n_genes"].values
        cfg = self.config
        n_cells = len(n_genes)

        # Strategy-aware percentile
        if strategy == StrategyType.CONSERVATIVE:
            percentile_value = cfg.percentile_conservative
        elif strategy == StrategyType.AGGRESSIVE:
            percentile_value = cfg.percentile_aggressive
        else:
            percentile_value = cfg.percentile_standard

        # Small datasets: use percentile directly for stability
        if n_cells < 500 or cfg.min_genes_model == "percentile":
            return self._recommend_min_genes_percentile(n_genes, cfg, strategy)

        # Map config names to count-model names
        model = cfg.min_genes_model
        if model == "poisson_gamma":
            model = "nb"
        elif model == "auto":
            model = "auto"

        count_fit = fit_count_mixture_threshold_model(
            n_genes,
            direction="lower",
            percentile=percentile_value,
            model=model,
            random_state=cfg.random_state,
            fallback=True,
        )

        threshold = float(count_fit["threshold"])
        threshold = max(float(cfg.min_genes_absolute), threshold)

        # Stratified bootstrap that preserves the zero / non-zero structure
        rng = np.random.default_rng(cfg.random_state)
        n_bootstrap = max(cfg.n_bootstrap, 200)
        boot_thresholds = []
        for _ in range(n_bootstrap):
            boot_sample = rng.choice(n_genes, size=len(n_genes), replace=True)
            boot_fit = fit_count_mixture_threshold_model(
                boot_sample,
                direction="lower",
                percentile=percentile_value,
                model=count_fit.get("model", "auto") if count_fit.get("model") != "gmm" else "auto",
                random_state=cfg.random_state + _,
                fallback=True,
            )
            if boot_fit.get("is_success"):
                boot_thresholds.append(float(boot_fit["threshold"]))
            else:
                boot_thresholds.append(np.percentile(boot_sample, percentile_value))

        ci_lower = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_lower))
        ci_upper = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_upper))
        threshold = float(np.clip(threshold, ci_lower, ci_upper))

        # Confidence: based on AIC advantage over GMM fallback (if available)
        confidence = 0.7
        if not count_fit.get("fallback_used") and count_fit.get("all_aic"):
            all_aic = count_fit["all_aic"]
            best_aic = count_fit.get("aic")
            gmm_aic = all_aic.get("gmm")
            if best_aic is not None and gmm_aic is not None and gmm_aic > best_aic:
                delta_aic = min(gmm_aic - best_aic, 100.0)
                confidence = min(1.0, max(0.5, 0.5 + delta_aic / 40.0))
            elif best_aic is not None:
                confidence = min(1.0, max(0.5, 0.7 + 10.0 / max(best_aic, 1.0)))

        method_name = count_fit.get("method", "count mixture")
        evidence = {
            "count_model": count_fit.get("model"),
            "aic": count_fit.get("aic"),
            "all_aic": count_fit.get("all_aic"),
            "params": count_fit.get("params"),
            "n_components": count_fit.get("gmm_fit", {}).get("n_components")
            if count_fit.get("fallback_used")
            else None,
            "n_bootstrap": n_bootstrap,
            "strategy": strategy.value,
            "method": method_name,
            "fallback_used": count_fit.get("fallback_used", False),
            "method_limitation": (
                "Count-mixture threshold is a data-driven heuristic for discrete QC metrics; "
                "it should be reviewed against the actual n_genes distribution."
            ),
        }

        if plot and save_dir:
            save_dir.mkdir(parents=True, exist_ok=True)
            self._plot_min_genes_analysis(
                n_genes, int(threshold), boot_thresholds, save_dir / "min_genes_recommendation.pdf"
            )
            evidence["plot"] = str(save_dir / "min_genes_recommendation.pdf")

        return ThresholdRecommendation(
            threshold=int(threshold),
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method=method_name,
            confidence=float(confidence),
            evidence=evidence,
        )

    def _recommend_min_genes_percentile(
        self,
        n_genes: np.ndarray,
        cfg: IntelligentQCConfig,
        strategy: StrategyType,
    ) -> ThresholdRecommendation:
        """Fallback percentile-based recommendation for small datasets."""
        if strategy == StrategyType.CONSERVATIVE:
            pct = cfg.percentile_conservative
        elif strategy == StrategyType.AGGRESSIVE:
            pct = cfg.percentile_aggressive
        else:
            pct = cfg.percentile_standard

        threshold = int(np.percentile(n_genes, pct))
        threshold = max(cfg.min_genes_absolute, threshold)

        # Simple bootstrap CI (data-adaptive count)
        n_bootstrap = max(50, min(1000, len(n_genes) // 10))
        rng = np.random.default_rng(cfg.random_state)
        boot = [
            np.percentile(rng.choice(n_genes, size=len(n_genes), replace=True), pct)
            for _ in range(n_bootstrap)
        ]
        ci_lower = int(np.percentile(boot, cfg.bootstrap_percentile_lower))
        ci_upper = int(np.percentile(boot, cfg.bootstrap_percentile_upper))

        return ThresholdRecommendation(
            threshold=threshold,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            method="percentile (small dataset fallback)",
            confidence=0.6,  # Lower confidence for small data
            evidence={"n_cells": len(n_genes), "percentile": pct},
        )

    def _recommend_max_mt(
        self,
        adata: AnnData,
        tissue_type: str,
        strategy: StrategyType,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> ThresholdRecommendation:
        """Recommend max_mt_percent threshold with tumor metabolic awareness.

        For tumor tissues this method supports:

        - **Sample/cell-type aware baselines**: each stratum gets its own
          median + nMAD band so that metabolically stressed populations are not
          mechanically removed against a global low-MT baseline.
        - **Multi-component GMM**: up to ``cfg.mt_max_components`` components
          capture continuous metabolic gradients (OXPHOS, glycolysis, hypoxia,
          proliferation).
        - **Review band + hard threshold**: cells above the hard threshold are
          flagged for filtering; cells between the stratum baseline +
          ``mt_review_mad_factor`` MAD and the hard threshold are flagged for
          review but not automatically removed.

        The returned ``ThresholdRecommendation.evidence`` contains
        ``review_band_lower`` and ``review_required`` so callers can build a
        two-tier QC policy.
        """
        mt_pct = np.asarray(adata.obs["pct_counts_mt"].values, dtype=float)
        mt_pct_nonzero = mt_pct[mt_pct > 0]

        if len(mt_pct_nonzero) == 0:
            return ThresholdRecommendation(
                threshold=20.0,
                ci_lower=20.0,
                ci_upper=20.0,
                method="no_mt_data",
                confidence=0.5,
                evidence={"reason": "No mitochondrial genes detected"},
            )

        cfg = self.config
        is_tumor_dataset = is_tumor_context(tissue_type)

        # Strategy-aware percentile baseline
        if is_tumor_dataset or strategy == StrategyType.TUMOR_AWARE:
            threshold_percentile = 95.0 if strategy == StrategyType.CONSERVATIVE else 90.0
        else:
            threshold_percentile = 90.0 if strategy == StrategyType.CONSERVATIVE else 85.0

        threshold = float(np.percentile(mt_pct_nonzero, threshold_percentile))
        method_name = f"distribution fitting (percentile {threshold_percentile:.0f})"
        best_dist = "percentile"
        best_params: Dict[str, Any] = {}
        dist_results: Dict[str, Any] = {}
        bimodal_info: Dict[str, Any] = {}
        multicomponent_info: Dict[str, Any] = {}
        tumor_program_signal: Dict[str, Any] = {}
        review_required = False
        review_band_lower = float(np.percentile(mt_pct_nonzero, 80.0))
        stratum_baselines: Dict[str, Any] = {}

        # 1. Try beta / lognorm single-distribution fits (existing path)
        from scipy.stats import beta, lognorm

        for dist_name, dist in [("beta", beta), ("lognorm", lognorm)]:
            try:
                params = dist.fit(mt_pct_nonzero, floc=0)
                ks_stat, ks_pval = stats.kstest(mt_pct_nonzero, dist.cdf(*params))
                dist_results[dist_name] = {"params": params, "ks_stat": ks_stat, "ks_pval": ks_pval}
            except Exception:
                pass

        if dist_results:
            best_dist = min(dist_results.keys(), key=lambda k: dist_results[k]["ks_stat"])
            best_params = {"params": dist_results[best_dist]["params"]}
            method_name = f"distribution fitting ({best_dist})"

        # 2. Stratum-aware baseline when requested or in tumor-aware auto mode.
        requested_stratum_aware = cfg.mt_model in {"sample_aware", "celltype_aware"}
        auto_stratum_aware = (
            is_tumor_dataset
            and cfg.mt_model == "auto"
            and (cfg.sample_key in adata.obs.columns or resolve_cell_type_key(adata) is not None)
        )
        use_stratum_aware = requested_stratum_aware or auto_stratum_aware
        if use_stratum_aware:
            stratum_key = None
            cell_type_key = resolve_cell_type_key(adata)
            if cfg.mt_model == "celltype_aware" and cell_type_key is not None:
                stratum_key = cell_type_key
            elif cfg.mt_model == "sample_aware" and cfg.sample_key in adata.obs.columns:
                stratum_key = cfg.sample_key
            elif cell_type_key is not None:
                stratum_key = cell_type_key
            elif cfg.sample_key in adata.obs.columns:
                stratum_key = cfg.sample_key

            if stratum_key is not None:
                baselines, global_review_lower = self._compute_mt_stratum_baselines(
                    adata, mt_pct, stratum_key, cfg.mt_mad_factor, cfg.mt_review_mad_factor
                )
                stratum_baselines = baselines
                review_band_lower = global_review_lower
                # The hard threshold is the most permissive of the percentile
                # baseline and the highest stratum-specific hard threshold.
                stratum_hard_thresholds = [
                    info["hard_threshold"]
                    for info in baselines.values()
                    if np.isfinite(info["hard_threshold"])
                ]
                if stratum_hard_thresholds:
                    threshold = max(threshold, float(np.max(stratum_hard_thresholds)))
                    method_name = f"{cfg.mt_model}_mad"

        # 3. For tumor datasets, try bimodal GMM to identify a high-MT population
        use_bimodal = (is_tumor_dataset or strategy == StrategyType.TUMOR_AWARE) and cfg.mt_model in (
            "auto",
            "bimodal_gmm",
        )
        if use_bimodal and len(mt_pct_nonzero) >= 100:
            bimodal_fit = fit_bimodal_gmm_threshold_model(
                mt_pct_nonzero,
                direction="upper",
                random_state=cfg.random_state,
                min_separation=cfg.mt_bimodal_min_separation,
            )
            if bimodal_fit.get("is_bimodal"):
                crossing = float(bimodal_fit["threshold"])
                threshold = max(threshold, crossing)
                method_name = "bimodal_gmm"
                bimodal_info = {
                    "crossing": crossing,
                    "separation": bimodal_fit.get("separation"),
                    "component_means": bimodal_fit.get("component_means"),
                    "component_stds": bimodal_fit.get("component_stds"),
                }

        # 4. Multi-component GMM for tumor metabolic gradients.
        use_multicomponent = (
            is_tumor_dataset or strategy == StrategyType.TUMOR_AWARE
        ) and cfg.mt_model in ("auto", "multicomponent")
        if use_multicomponent and len(mt_pct_nonzero) >= 100:
            multi_fit = fit_gmm_threshold_model(
                mt_pct_nonzero,
                direction="upper",
                random_state=cfg.random_state,
                max_components=cfg.mt_max_components,
            )
            if multi_fit.get("n_components", 1) > 1:
                multi_threshold = float(multi_fit["threshold"])
                threshold = max(threshold, multi_threshold)
                method_name = "multicomponent_gmm"
                multicomponent_info = {
                    "n_components": multi_fit.get("n_components"),
                    "bic": multi_fit.get("bic"),
                    "null_bic": multi_fit.get("null_bic"),
                }

        # 5. Tumor-program co-expression check for high-MT cells
        if is_tumor_dataset and threshold is not None and np.isfinite(threshold):
            high_mt_mask = mt_pct > threshold
            if high_mt_mask.sum() >= 10:
                tumor_program_signal = self._assess_high_mt_tumor_programs(
                    adata, high_mt_mask, cfg.mt_tumor_programs
                )
                if tumor_program_signal.get("signal_detected"):
                    review_required = True

        # Bootstrap CI
        n_bootstrap = cfg.n_bootstrap
        boot_thresholds = []
        rng = np.random.default_rng(cfg.random_state)
        for _ in range(n_bootstrap):
            boot_sample = rng.choice(mt_pct_nonzero, size=len(mt_pct_nonzero), replace=True)
            boot_thresh = np.percentile(boot_sample, threshold_percentile)
            boot_thresholds.append(boot_thresh)

        ci_lower = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_lower))
        ci_upper = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_upper))
        threshold = float(np.clip(threshold, ci_lower, ci_upper))

        # Confidence
        if method_name in {"bimodal_gmm", "multicomponent_gmm"}:
            separation = bimodal_info.get("separation", 0.0)
            confidence = min(1.0, max(0.5, 0.6 + 0.05 * separation))
        elif best_dist != "percentile" and dist_results:
            confidence = min(1.0, max(0.5, 0.5 + 0.5 * dist_results[best_dist]["ks_pval"]))
        else:
            confidence = 0.7

        evidence: Dict[str, Any] = {
            "best_distribution": best_dist,
            "params": best_params,
            "tissue_type": tissue_type,
            "median_mt": float(np.median(mt_pct_nonzero)),
            "method": method_name,
            "threshold_percentile": threshold_percentile,
            "bimodal": bimodal_info,
            "multicomponent": multicomponent_info,
            "stratum_baselines": stratum_baselines,
            "review_band_lower": round(review_band_lower, 1),
            "tumor_program_signal": tumor_program_signal,
            "review_required": review_required,
        }

        risk_note = ""
        if review_required:
            risk_note = (
                "High-MT cells co-express tumor/stress/proliferation/metabolic programs; "
                "avoid mechanical deletion and review biologically."
            )
            evidence["risk_note"] = risk_note

        return ThresholdRecommendation(
            threshold=round(float(threshold), 1),
            ci_lower=round(float(ci_lower), 1),
            ci_upper=round(float(ci_upper), 1),
            method=method_name,
            confidence=float(confidence),
            evidence=evidence,
        )

    def _compute_mt_stratum_baselines(
        self,
        adata: AnnData,
        mt_pct: np.ndarray,
        stratum_key: str,
        hard_mad_factor: float,
        review_mad_factor: float,
    ) -> Tuple[Dict[str, Any], float]:
        """Compute per-stratum MT% baselines and a global review-band lower bound.

        For each stratum with at least 25 cells, the hard threshold is
        ``median + hard_mad_factor * MAD`` and the review-band lower bound is
        ``median + review_mad_factor * MAD``.  The function returns a dict of
        per-stratum summaries and the lowest review-band lower bound across
        strata, which is used as the global ``review_band_lower``.
        """
        strata = adata.obs.groupby(stratum_key, observed=False).indices
        baselines: Dict[str, Any] = {}
        review_lowers = []
        for name, idx in strata.items():
            values = mt_pct[idx]
            values = values[np.isfinite(values) & (values >= 0)]
            if len(values) < 25:
                continue
            med = float(np.median(values))
            mad = float(np.median(np.abs(values - med)))
            # MAD to std-like scaling
            scaled_mad = mad * 1.4826 if mad > 0 else 0.0
            hard = med + hard_mad_factor * scaled_mad
            review = med + review_mad_factor * scaled_mad
            baselines[str(name)] = {
                "n_cells": len(values),
                "median": round(med, 3),
                "mad": round(mad, 3),
                "hard_threshold": round(hard, 3),
                "review_threshold": round(review, 3),
            }
            review_lowers.append(review)

        global_review_lower = float(np.min(review_lowers)) if review_lowers else float(
            np.percentile(mt_pct[np.isfinite(mt_pct)], 80.0)
        )
        return baselines, global_review_lower

    def _assess_high_mt_tumor_programs(
        self,
        adata: AnnData,
        high_mt_mask: np.ndarray,
        program_panel_names: Tuple[str, ...],
    ) -> Dict[str, Any]:
        """Check whether high-MT cells retain meaningful biology vs. being debris.

        Instead of looking for enrichment of tumor programs in high-MT cells
        (which are often damaged and have globally lower expression), we compare
        the *relative retention* of the marker signal in the high-MT population
        versus the low-MT population.  If high-MT cells retain a non-trivial
        fraction of the marker signal, the high-MT peak should not be
        mechanically removed.

        Returns a dictionary with ``signal_detected`` (bool),
        ``enriched_programs``, and per-program statistics.
        """
        result: Dict[str, Any] = {"signal_detected": False, "enriched_programs": []}
        try:
            from validation.gene_panels import TUMOR_PROGRAM_PANELS
        except Exception:
            return result

        X = adata.X
        if hasattr(X, "toarray"):
            X_arr = X.toarray()
        else:
            X_arr = np.asarray(X)

        var_names = pd.Index(adata.var_names.astype(str))
        low_mt_mask = ~high_mt_mask
        n_high = int(high_mt_mask.sum())
        n_low = int(low_mt_mask.sum())
        if n_high == 0 or n_low == 0:
            return result

        programs_enriched = []
        program_stats = []

        for panel_name in program_panel_names:
            genes = TUMOR_PROGRAM_PANELS.get(panel_name, ())
            present = [g for g in genes if g in var_names]
            if len(present) < 2:
                continue
            idx = var_names.get_indexer(present)
            high_expr = np.asarray(X_arr[high_mt_mask, :][:, idx].mean(axis=0)).ravel()
            low_expr = np.asarray(X_arr[low_mt_mask, :][:, idx].mean(axis=0)).ravel()
            # Total signal per population, normalised by population size
            high_total = float(high_expr.sum())
            low_total = float(low_expr.sum())
            if low_total <= 0:
                continue
            # Relative retention of this program in high-MT vs low-MT cells
            relative_retention = (high_total / n_high) / (low_total / n_low)
            # Absolute presence: high-MT cells should express at least some marker genes
            present_fraction = float((high_expr > 0).sum()) / len(present)
            program_stats.append(
                {
                    "program": panel_name,
                    "high_total_signal": high_total,
                    "low_total_signal": low_total,
                    "relative_retention": relative_retention,
                    "present_fraction": present_fraction,
                }
            )
            if (
                relative_retention >= self.config.tumor_program_min_relative_retention
                and present_fraction >= self.config.tumor_program_min_present_fraction
            ):
                programs_enriched.append(panel_name)

        result["program_stats"] = program_stats
        result["thresholds"] = {
            "min_relative_retention": self.config.tumor_program_min_relative_retention,
            "min_present_fraction": self.config.tumor_program_min_present_fraction,
        }
        if programs_enriched:
            result["signal_detected"] = True
            result["enriched_programs"] = programs_enriched
        return result

    def _recommend_n_counts(
        self,
        adata: AnnData,
        strategy: StrategyType,
        plot: bool = True,
        save_dir: Optional[Path] = None,
    ) -> ThresholdRecommendation:
        """Recommend n_counts threshold with a count-distribution model."""
        n_counts = np.asarray(adata.obs["n_counts"].values, dtype=float)
        valid_counts = n_counts[np.isfinite(n_counts) & (n_counts >= 0)]
        positive_counts = valid_counts[valid_counts > 0]
        if positive_counts.size == 0:
            return ThresholdRecommendation(
                threshold=0,
                ci_lower=0,
                ci_upper=0,
                method="count_mixture_unavailable",
                confidence=0.0,
                evidence={"reason": "no_positive_counts", "n_cells": int(len(n_counts))},
            )

        if strategy == StrategyType.CONSERVATIVE:
            threshold_percentile = 10.0
        elif strategy == StrategyType.AGGRESSIVE:
            threshold_percentile = 20.0
        else:
            threshold_percentile = 15.0

        cfg = self.config
        count_fit = fit_count_mixture_threshold_model(
            valid_counts,
            direction="lower",
            percentile=threshold_percentile,
            model="auto",
            random_state=cfg.random_state,
            fallback=True,
        )
        threshold = float(count_fit["threshold"])

        # Bootstrap CI
        n_bootstrap = cfg.n_bootstrap
        boot_thresholds = []
        rng = np.random.default_rng(cfg.random_state)
        for _ in range(n_bootstrap):
            boot_sample = rng.choice(valid_counts, size=len(valid_counts), replace=True)
            boot_thresh = np.percentile(boot_sample, threshold_percentile)
            if np.isfinite(boot_thresh):
                boot_thresholds.append(float(boot_thresh))

        if boot_thresholds:
            ci_lower = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_lower))
            ci_upper = float(np.percentile(boot_thresholds, cfg.bootstrap_percentile_upper))
        else:
            ci_lower = ci_upper = threshold

        threshold = int(np.clip(threshold, ci_lower, ci_upper))

        n_cells = int(len(valid_counts))
        sample_score = 1.0 - np.exp(-n_cells / 1000.0)
        ci_width = max(float(ci_upper - ci_lower), 0.0)
        ci_precision = 1.0 / (1.0 + ci_width / max(float(threshold), 1.0))
        all_aic = count_fit.get("all_aic") or {}
        if len(all_aic) >= 2:
            ranked_aic = sorted(float(v) for v in all_aic.values() if np.isfinite(v))
            delta_aic = ranked_aic[1] - ranked_aic[0] if len(ranked_aic) >= 2 else 0.0
            model_score = float(np.clip(delta_aic / 20.0, 0.0, 1.0))
        else:
            model_score = 0.5 if count_fit.get("is_success") else 0.0
        fallback_penalty = 0.85 if count_fit.get("fallback_used") else 1.0
        confidence = fallback_penalty * (
            0.35 * sample_score + 0.35 * model_score + 0.30 * ci_precision
        )
        confidence = float(np.clip(confidence, 0.0, 1.0))

        return ThresholdRecommendation(
            threshold=threshold,
            ci_lower=int(ci_lower),
            ci_upper=int(ci_upper),
            method=f"{count_fit.get('method', 'count_mixture')} + bootstrap",
            confidence=float(confidence),
            evidence={
                "n_cells": n_cells,
                "model": count_fit.get("model"),
                "params": count_fit.get("params"),
                "aic": count_fit.get("aic"),
                "all_aic": count_fit.get("all_aic"),
                "fallback_used": count_fit.get("fallback_used"),
                "threshold_percentile": threshold_percentile,
                "bootstrap_ci_method": "empirical_resampled_percentile",
                "confidence_components": {
                    "sample_score": float(sample_score),
                    "model_score": float(model_score),
                    "ci_precision": float(ci_precision),
                },
            },
        )

    def _analyze_doublet_patterns(
        self, adata: AnnData, plot: bool = True, save_dir: Optional[Path] = None
    ) -> ThresholdRecommendation:
        """Analyze doublet patterns and recommend threshold."""
        if "doublet_score" not in adata.obs:
            # No doublet scores calculated
            return ThresholdRecommendation(
                threshold=0.5,
                ci_lower=0.5,
                ci_upper=0.5,
                method="no_doublet_scores",
                confidence=0.0,
                evidence={"reason": "Doublet scores not calculated"},
            )

        doublet_scores = np.asarray(adata.obs["doublet_score"].values, dtype=float)
        doublet_scores = doublet_scores[np.isfinite(doublet_scores)]
        if doublet_scores.size == 0:
            return ThresholdRecommendation(
                threshold=0.5,
                ci_lower=0.5,
                ci_upper=0.5,
                method="doublet_scores_unavailable",
                confidence=0.0,
                evidence={"reason": "no_finite_doublet_scores"},
            )

        # Use percentile method (more robust)
        if self.strategy == StrategyType.CONSERVATIVE:
            threshold_percentile = 95.0
        elif self.strategy == StrategyType.AGGRESSIVE:
            threshold_percentile = 85.0
        else:
            threshold_percentile = 90.0
        percentile_threshold = float(np.percentile(doublet_scores, threshold_percentile))
        threshold = percentile_threshold
        method = "doublet_score_percentile"
        gmm_fit: Dict[str, Any] = {}
        if doublet_scores.size >= 100 and float(np.std(doublet_scores)) > 0:
            try:
                gmm_fit = fit_gmm_threshold_model(
                    doublet_scores,
                    direction="upper",
                    random_state=self.config.random_state,
                    max_components=3,
                )
                if gmm_fit.get("n_components", 1) > 1 and np.isfinite(
                    gmm_fit.get("threshold", np.nan)
                ):
                    threshold = float(gmm_fit["threshold"])
                    method = "doublet_score_gmm"
            except Exception as exc:
                gmm_fit = {"is_success": False, "error": str(exc)}

        # Bootstrap CI
        cfg = self.config
        n_bootstrap = cfg.n_bootstrap
        boot_thresholds = []
        rng = np.random.default_rng(cfg.random_state)
        for _ in range(n_bootstrap):
            boot_sample = rng.choice(doublet_scores, size=len(doublet_scores), replace=True)
            if method == "doublet_score_gmm" and boot_sample.size >= 100:
                try:
                    boot_fit = fit_gmm_threshold_model(
                        boot_sample,
                        direction="upper",
                        random_state=cfg.random_state,
                        max_components=3,
                    )
                    boot_thresh = (
                        float(boot_fit["threshold"])
                        if boot_fit.get("n_components", 1) > 1
                        else float(np.percentile(boot_sample, threshold_percentile))
                    )
                except Exception:
                    boot_thresh = float(np.percentile(boot_sample, threshold_percentile))
            else:
                boot_thresh = float(np.percentile(boot_sample, threshold_percentile))
            boot_thresholds.append(boot_thresh)

        ci_lower = np.percentile(boot_thresholds, cfg.bootstrap_percentile_lower)
        ci_upper = np.percentile(boot_thresholds, cfg.bootstrap_percentile_upper)
        threshold = float(np.clip(threshold, ci_lower, ci_upper))

        ci_width = max(float(ci_upper - ci_lower), 0.0)
        ci_precision = 1.0 / (1.0 + ci_width / max(float(threshold), 1e-6))
        null_bic = gmm_fit.get("null_bic")
        bic = gmm_fit.get("bic")
        if method == "doublet_score_gmm" and null_bic is not None and bic is not None:
            mixture_score = float(np.clip((float(null_bic) - float(bic)) / 20.0, 0.0, 1.0))
        elif method == "doublet_score_gmm":
            mixture_score = 0.7
        else:
            mixture_score = 0.45
        sample_score = 1.0 - np.exp(-len(doublet_scores) / 1000.0)
        confidence = float(
            np.clip(0.35 * sample_score + 0.35 * mixture_score + 0.30 * ci_precision, 0.0, 1.0)
        )

        evidence = {
            "distribution_model": method,
            "gmm": gmm_fit,
            "percentile_threshold": percentile_threshold,
            "threshold_percentile": threshold_percentile,
            "doublet_rate": float((doublet_scores > threshold).mean()),
            "confidence_components": {
                "sample_score": float(sample_score),
                "mixture_score": float(mixture_score),
                "ci_precision": float(ci_precision),
            },
        }

        return ThresholdRecommendation(
            threshold=round(float(threshold), 3),
            ci_lower=round(float(ci_lower), 3),
            ci_upper=round(float(ci_upper), 3),
            method=f"{method} + bootstrap",
            confidence=float(confidence),
            evidence=evidence,
        )

    def _generate_concerns(
        self, quality_score: float, quality_flags: List[str], strategy: StrategyType
    ) -> List[str]:
        """Generate list of concerns based on quality assessment."""
        concerns = []

        if quality_score < 60:
            concerns.append(f"Low data quality score ({quality_score:.1f}/100)")

        if "High mitochondrial content" in " ".join(quality_flags):
            if strategy != StrategyType.TUMOR_AWARE:
                concerns.append(
                    "High mitochondrial content detected (consider using tumor_aware strategy)"
                )

        if "High doublet rate" in " ".join(quality_flags):
            concerns.append("High doublet rate detected")

        if len(quality_flags) > 3:
            concerns.append(f"Multiple quality issues detected ({len(quality_flags)})")

        metric_flags = [
            flag
            for flag in quality_flags
            if "Missing metric" in flag or "Missing n_" in flag or "Missing pct_counts_mt" in flag
        ]
        concerns.extend(metric_flags)

        return concerns

    def _get_tumor_considerations(
        self, adata: AnnData, tissue_type: str, quality_flags: List[str]
    ) -> List[str]:
        """Generate tumor-specific considerations."""
        considerations = []

        if is_tumor_context(tissue_type):
            considerations.append("Tumor tissue detected: Using elevated mitochondrial thresholds")

            # Check for mixed populations
            if "High doublet rate" in " ".join(quality_flags):
                considerations.append(
                    "Possible tumor-stromal mixture: Doublet-like patterns may be "
                    "genuine tumor cells interacting with normal cells"
                )

            # Consider cell cycle
            if "All cells in same cell cycle phase" in " ".join(quality_flags):
                if "S" in adata.obs.columns or "G2M" in adata.obs.columns:
                    s_cells = (adata.obs["cell_cycle_phase"] == "S").sum()
                    if s_cells / len(adata) > 0.7:
                        considerations.append(
                            f"High proliferative state ({s_cells/len(adata):.1%} S-phase cells) "
                            "- common in proliferating tumors"
                        )

        return considerations

    def _calculate_overall_confidence(
        self, recommendations: List[ThresholdRecommendation]
    ) -> float:
        """Calculate overall confidence from all recommendations."""
        confidences = [r.confidence for r in recommendations]
        return float(np.mean(confidences))

    def _save_recommendation_report(self, recommendation: QCRecommendation, save_dir: Path):
        """Save recommendation report."""
        import json

        save_dir.mkdir(parents=True, exist_ok=True)

        # Save as JSON
        json_path = save_dir / "qc_recommendation.json"
        with open(json_path, "w") as f:
            json.dump(_make_json_safe(recommendation.to_dict()), f, indent=2)

        log.info(f"  Recommendation saved to: {json_path}")

    def _plot_min_genes_analysis(
        self, n_genes: np.ndarray, threshold: int, boot_thresholds: List[int], save_path: Path
    ):
        """Plot min_genes recommendation with evidence."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Histogram
        axes[0].hist(n_genes, bins=50, alpha=0.7, edgecolor="black")
        axes[0].axvline(
            threshold, color="red", linestyle="--", linewidth=2, label="Recommended threshold"
        )
        axes[0].axvline(200, color="gray", linestyle=":", label="Traditional threshold (200)")
        axes[0].set_xlabel("Number of genes")
        axes[0].set_ylabel("Number of cells")
        axes[0].set_title("Distribution of n_genes")
        axes[0].legend()

        # Bootstrap CI
        axes[1].hist(boot_thresholds, bins=30, alpha=0.7, edgecolor="black")
        axes[1].axvline(threshold, color="red", linestyle="--", linewidth=2)
        axes[1].axvline(
            np.percentile(boot_thresholds, 2.5), color="blue", linestyle=":", label="95% CI"
        )
        axes[1].axvline(
            np.percentile(boot_thresholds, 97.5), color="blue", linestyle=":", label="95% CI"
        )
        axes[1].set_xlabel("Bootstrap threshold")
        axes[1].set_ylabel("Frequency")
        axes[1].set_title(
            f"Bootstrap 95% CI: [{int(np.percentile(boot_thresholds, 2.5))}, "
            f"{int(np.percentile(boot_thresholds, 97.5))}]"
        )
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()


# Convenience function
def recommend_intelligent_qc(
    adata: AnnData,
    tissue_type: str = "unknown",
    strategy: str = "auto",
    plot: bool = True,
    save_dir: Optional[Path] = None,
    **kwargs,
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

    Returns:
    -------
    QCRecommendation
        Complete recommendation with all thresholds and confidence intervals

    Examples:
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

    Notes:
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
    recommender = IntelligentQCRecommender(strategy=StrategyType(strategy))

    return recommender.recommend(adata=adata, tissue_type=tissue_type, plot=plot, save_dir=save_dir)


__all__ = [
    "IntelligentQCRecommender",
    "IntelligentQCConfig",
    "recommend_intelligent_qc",
    "QCRecommendation",
    "ThresholdRecommendation",
    "StrategyType",
]
