"""Ensemble doublet detection pipeline and evidence profiling.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from ..config import DoubletConfig, MarkerConfig
from .algorithms import _run_doubletdetection, _run_scrublet, _run_solo
from .core import (
    FINAL_PRED_COL,
    HEURISTIC_PRED_COL,
    HEURISTIC_SCORE_COL,
    LINEAGE_SCORES_KEY,
    _create_doublet_marker_config_from_manager,
    create_custom_marker_dict,
    generate_doublet_rates,
)
from .heuristic import _plot_doublet_summary, _run_heuristic

log = logging.getLogger(__name__)

def _merge_doublet_predictions(
    adata: AnnData,
    algorithm_score_col: str,
    heuristic_score_col: str,
    strategy: str = "weighted_average",
    algo_weight: float = 0.6,
    expected_rate: Optional[Union[float, Dict[str, float]]] = 0.1,
    score_threshold: Optional[float] = None,
) -> pd.Series:
    """
    Merge algorithmic and heuristic doublet scores for a final, more robust prediction.
    This function combines two continuous score series instead of binary predictions.

    Args:
        adata: AnnData object containing the scores.
        algorithm_score_col: Column name in adata.obs for the algorithm's score (e.g., 'scrublet_score').
        heuristic_score_col: Column name in adata.obs for the heuristic confidence score.
        strategy: The merge strategy ('weighted_average', 'max_score', 'heuristic_boost').
        algo_weight: The weight for the algorithm's score in 'weighted_average' strategy.

    Returns:
        A boolean pandas Series with the final merged doublet predictions.
    """
    algo_scores = adata.obs[algorithm_score_col].fillna(0)
    heur_scores = adata.obs[heuristic_score_col].fillna(0)

    final_score = pd.Series(0.0, index=adata.obs_names)

    if strategy == "weighted_average":
        # A simple weighted average. algo_weight determines the trust in the algorithm.
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)
    elif strategy == "max_score":
        # Takes the highest score from either method, useful if either method is considered reliable on its own.
        final_score = pd.DataFrame({"algo": algo_scores, "heur": heur_scores}).max(axis=1)
    elif strategy == "heuristic_boost":
        # Uses the algorithm score as a base and the heuristic score as a "booster".
        # This is useful for finding doublets missed by the algorithm but strongly suggested by heuristics.
        final_score = algo_scores + (heur_scores * 0.5)  # Boost factor can be tuned
    else:
        log.warning(
            f"Unknown enhanced merge strategy '{strategy}', falling back to 'weighted_average'."
        )
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)

    # Normalize the final combined score to a [0, 1] range for consistent thresholding.
    if final_score.max() > 0:
        final_score /= final_score.max()

    if score_threshold is not None:
        threshold = score_threshold
        log.info(
            f"Using user-provided doublet score threshold of {threshold:.3f} for merged predictions."
        )
    else:
        if expected_rate is None:
            log.warning("expected_doublet_rate is None, using a default of 0.1 for thresholding.")
            expected_rate = 0.1

        if isinstance(expected_rate, dict):  # Handle per-sample rates by taking the mean
            expected_rate = np.mean(list(expected_rate.values()))

        threshold = final_score.quantile(1 - expected_rate)
        log.info(
            f"Using a final score threshold of {threshold:.3f} based on expected doublet rate for merged predictions."
        )

    return final_score > threshold



def _export_doublet_stats(
    adata: AnnData,
    sample_key: str = "sampleID",
    save_dir: Optional[Union[str, Path]] = None,
    export_csv: bool = True,
    export_xlsx: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Export comprehensive doublet statistics per sample and globally.

    This function generates detailed statistical summaries of doublet detection
    results, including counts, percentages, and score distributions.

    Args:
        adata: AnnData object with doublet predictions
        sample_key: Key for sample identification
        save_dir: Directory to save statistics files
        export_csv: Whether to export as CSV files
        export_xlsx: Whether to export as Excel file

    Returns:
        Dictionary containing sample-wise and global statistics DataFrames
    """
    # Identify all doublet-related columns
    doublet_cols = [
        col
        for col in adata.obs.columns
        if any(keyword in col.lower() for keyword in ["doublet", "scrublet", "heuristic"])
    ]

    if not doublet_cols:
        log.warning("No doublet-related columns found in adata.obs")
        return {}

    log.info(f"Found doublet columns: {doublet_cols}")

    # Calculate per-sample statistics
    sample_stats = []
    unique_samples = adata.obs[sample_key].unique()
    if not pd.api.types.is_categorical_dtype(adata.obs[sample_key]):
        unique_samples = sorted(unique_samples)

    for sample in unique_samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs.loc[sample_mask]

        stats = {
            "sample": sample,
            "total_cells": len(sample_data),
        }

        for col in doublet_cols:
            if col in sample_data.columns:
                col_data = sample_data[col].dropna()
                if (
                    pd.api.types.is_numeric_dtype(col_data)
                    and not pd.api.types.is_bool_dtype(col_data)
                    and col_data.nunique() > 2
                ):
                    # Continuous column (scores)
                    stats[f"{col}_mean"] = col_data.mean()
                    stats[f"{col}_median"] = col_data.median()
                    stats[f"{col}_std"] = col_data.std()
                elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                    # Boolean/binary column (predictions)
                    positive_count = col_data.astype(bool).sum()
                    stats[f"{col}_count"] = positive_count
                    stats[f"{col}_percentage"] = (
                        (positive_count / len(sample_data) * 100) if len(sample_data) > 0 else 0
                    )

        sample_stats.append(stats)

    sample_df = pd.DataFrame(sample_stats).set_index("sample")

    global_stats = {"metric": "global", "total_cells": adata.n_obs}
    for col in doublet_cols:
        if col in adata.obs.columns:
            col_data = adata.obs[col].dropna()
            if (
                pd.api.types.is_numeric_dtype(col_data)
                and not pd.api.types.is_bool_dtype(col_data)
                and col_data.nunique() > 2
            ):
                global_stats[f"{col}_mean"] = col_data.mean()
                global_stats[f"{col}_median"] = col_data.median()
                global_stats[f"{col}_std"] = col_data.std()
            elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                positive_count = col_data.astype(bool).sum()
                global_stats[f"{col}_count"] = positive_count
                global_stats[f"{col}_percentage"] = (
                    (positive_count / adata.n_obs * 100) if adata.n_obs > 0 else 0
                )

    global_df = pd.DataFrame([global_stats])

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        if export_csv:
            sample_df.to_csv(save_dir / "doublet_stats_per_sample.csv")
            global_df.to_csv(save_dir / "doublet_stats_global.csv", index=False)
            log.info(f"Exported CSV files to {save_dir}")
        if export_xlsx:
            with pd.ExcelWriter(save_dir / "doublet_stats.xlsx") as writer:
                sample_df.to_excel(writer, sheet_name="per_sample")
                global_df.to_excel(writer, sheet_name="global", index=False)
            log.info(f"Exported Excel file to {save_dir / 'doublet_stats.xlsx'}")

    return {"sample": sample_df, "global": global_df}



def predict_doublets(
    adata: AnnData, config: DoubletConfig, sample_key: str = "sampleID", **kwargs
) -> AnnData:
    """
    Enhanced doublet prediction with a clear, config-driven workflow.
    This version integrates a quantitative heuristic score with the algorithmic score for improved accuracy.

    Args:
        adata: AnnData object containing single-cell expression data.
        config: A `DoubletConfig` object that controls the entire workflow.
        sample_key: Key for sample identification in adata.obs.

    Returns:
        AnnData object with doublet predictions added to .obs and .obsm.
    """
    # === 1. CONFIGURATION SETUP ===
    base_config = DoubletConfig()

    if config is not None:
        config_dict = config.to_dict()  # Pydantic's built-in serialization
        for key, value in config_dict.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)

    if kwargs:
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")

    cfg = base_config
    # Pydantic configs validate automatically
    log.info("--- Running Final Doublet Prediction Workflow ---")

    # Validate input data
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not found in adata.obs")

    samples = adata.obs[sample_key].unique()
    if len(samples) == 0:
        raise ValueError(f"No samples found for key '{sample_key}'")

    log.info(f"Starting doublet prediction for {adata.n_obs} cells across {len(samples)} samples")
    log.info(f"Configuration: method={cfg.method}, merge_strategy={cfg.merge_strategy}, ")

    # Initialize result columns
    algo_score_col = f"{cfg.method}_score"
    algo_pred_col = f"{cfg.method}_predicted"
    adata.obs[algo_score_col] = np.nan
    adata.obs[algo_pred_col] = False

    # Use a dispatcher for multi-algorithm support ---
    ALGORITHM_DISPATCHER = {
        "scrublet": _run_scrublet,
        "solo": _run_solo,
        "doubletdetection": _run_doubletdetection,  # Future-ready
    }
    if cfg.method not in ALGORITHM_DISPATCHER:
        raise ValueError(
            f"Method '{cfg.method}' is not supported. Available: {list(ALGORITHM_DISPATCHER.keys())}"
        )

    # === 2. ALGORITHMIC DETECTION (Per-Sample) ===
    if cfg.run_algorithm:
        log.info(f"Running {cfg.method} doublet detection...")

        for sample in samples:
            log.info(f"Processing sample '{sample}' with {cfg.method}...")
            sample_mask = adata.obs[sample_key] == sample
            data_view = adata[sample_mask]

            if data_view.n_obs < 50:
                log.warning(
                    f"Skipping {sample}: fewer than 50 cells (insufficient for reliable doublet detection)."
                )
                continue

            scores, predicted = ALGORITHM_DISPATCHER[cfg.method](data_view, sample, cfg)

            if scores is not None and predicted is not None:
                adata.obs.loc[sample_mask, algo_score_col] = scores
                adata.obs.loc[sample_mask, algo_pred_col] = predicted
    else:
        log.info("Skipping algorithmic detection as per configuration (run_algorithm=False).")

    # === 3. HEURISTIC DETECTION (Global) ===
    adata.obs[HEURISTIC_PRED_COL] = False
    adata.obs[HEURISTIC_SCORE_COL] = 0.0
    if cfg.use_heuristics:
        log.info("Running quantitative heuristic analysis...")
        # Call the new heuristic function and receive its multiple outputs
        heuristic_pred, lineage_scores_df, heuristic_scores = _run_heuristic(adata, cfg)

        # Store all the new results in the AnnData object
        adata.obsm["lineage_module_scores"] = lineage_scores_df  # Store detailed scores in .obsm
        adata.obs[HEURISTIC_PRED_COL] = heuristic_pred  # Store the binary call for simple stats
        adata.obs[HEURISTIC_SCORE_COL] = heuristic_scores  # Store the informative continuous score
        log.info(
            f"Heuristic analysis complete. Found {heuristic_pred.sum()} potential doublets based on score threshold."
        )

    # === 4. MERGE RESULTS ===
    log.info("Merging algorithmic and heuristic scores for final prediction...")
    merged_pred = _merge_doublet_predictions(
        adata,
        algorithm_score_col=algo_score_col,
        heuristic_score_col=HEURISTIC_SCORE_COL,
        strategy=cfg.merge_strategy,
        expected_rate=cfg.expected_doublet_rate,
        algo_weight=cfg.algorithm_weight,
        score_threshold=cfg.score_threshold,
    )
    adata.obs[FINAL_PRED_COL] = merged_pred

    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    adata.uns["sclucid"]["qc"]["doublet_params"].update(
        {
            "merge_strategy": cfg.merge_strategy,
            "algorithm_weight": cfg.algorithm_weight,
            "expected_doublet_rate": cfg.expected_doublet_rate,
            "score_threshold": cfg.score_threshold,
            "method": cfg.method,
        }
    )

    # === 5. SUMMARY STATISTICS ===
    log.info("\n" + "=" * 50)
    log.info("DOUBLET DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Algorithm results
    algo_count = adata.obs[algo_pred_col].sum()
    log.info(f"Algorithm ({cfg.method}): {algo_count} doublets ({algo_count / total_cells:.2%})")

    # Heuristic results
    if cfg.use_heuristics:
        heur_count = adata.obs[HEURISTIC_PRED_COL].sum()
        log.info(f"Heuristic: {heur_count} doublets ({heur_count / total_cells:.2%})")

        # Overlap analysis
        overlap_count = (adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL]).sum()
        log.info(f"Overlap: {overlap_count} doublets ({overlap_count / total_cells:.2%})")

    # Final merged results
    final_count = adata.obs[FINAL_PRED_COL].sum()
    log.info(f"Final merged: {final_count} doublets ({final_count / total_cells:.2%})")

    # Per-sample breakdown
    log.info("\nPer-sample statistics:")
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_total = sample_mask.sum()
        sample_doublets = adata.obs[FINAL_PRED_COL][sample_mask].sum()
        sample_rate = sample_doublets / sample_total
        log.info(f"  {sample}: {sample_doublets}/{sample_total} doublets ({sample_rate:.2%})")

    log.info("=" * 50)

    # === 6. Reporting & Visualization ===
    if cfg.plot_summary:
        save_path = Path(cfg.save_dir) if cfg.save_dir else None
        _plot_doublet_summary(
            adata=adata,
            sample_key=sample_key,
            save_dir=save_path,
            show=cfg.show_plots,
            plot_bar=cfg.plot_bar,
            plot_scatter=cfg.plot_scatter,
            plot_upset=cfg.plot_upset,
        )

    if cfg.export_stats and cfg.save_dir:
        _export_doublet_stats(adata, sample_key, Path(cfg.save_dir))

    log.info("Doublet prediction workflow completed.")

    return adata


class DoubletEvidenceProfiler:
    """
    Generate interpretable evidence profiles for doublet predictions.

    This class creates detailed reports explaining WHY each cell was
    flagged as a doublet, combining multiple lines of evidence.
    """

    def __init__(self, adata: AnnData):
        self.adata = adata
        self.evidence_table = None

    def generate_evidence_table(self) -> pd.DataFrame:
        """
        Create a comprehensive evidence table for each cell.

        Returns:
            DataFrame with one row per cell, columns for different evidence types
        """
        evidence = pd.DataFrame(index=self.adata.obs_names)

        # Evidence 1: Algorithmic score
        if "scrublet_score" in self.adata.obs:
            evidence["scrublet_score"] = self.adata.obs["scrublet_score"]
            evidence["scrublet_evidence"] = pd.cut(
                evidence["scrublet_score"],
                bins=[-np.inf, 0.2, 0.4, 0.6, np.inf],
                labels=["Weak", "Moderate", "Strong", "Very Strong"],
            )

        # Evidence 2: Lineage co-expression
        if "lineage_module_scores" in self.adata.obsm:
            lineage_scores = self.adata.obsm["lineage_module_scores"]

            # Count how many lineages are significantly expressed
            threshold = 0.5
            n_lineages = (lineage_scores > threshold).sum(axis=1)
            evidence["n_coexpressed_lineages"] = n_lineages

            # Identify the top 2 co-expressed lineages
            top_lineages = lineage_scores.apply(
                lambda row: (
                    lineage_scores.columns[np.argsort(row.values)[-2:]].tolist()
                    if row.max() > threshold
                    else []
                ),
                axis=1,
            )
            evidence["top_coexpressed_lineages"] = top_lineages.apply(
                lambda x: " + ".join(x) if len(x) >= 2 else "None"
            )

            # Strength of co-expression (product of top 2 scores)
            evidence["coexpression_strength"] = lineage_scores.apply(
                lambda row: np.prod(sorted(row.values)[-2:]) if row.max() > threshold else 0, axis=1
            )

        # Evidence 3: Gene count anomaly
        if "n_genes_by_counts" in self.adata.obs:
            # Z-score of gene counts
            gene_counts = self.adata.obs["n_genes_by_counts"]
            z_scores = (gene_counts - gene_counts.mean()) / gene_counts.std()
            evidence["gene_count_zscore"] = z_scores
            evidence["gene_count_anomaly"] = z_scores > 2  # High gene count

        # Evidence 4: Total UMI anomaly
        if "total_counts" in self.adata.obs:
            umi_counts = self.adata.obs["total_counts"]
            z_scores = (umi_counts - umi_counts.mean()) / umi_counts.std()
            evidence["umi_count_zscore"] = z_scores
            evidence["umi_count_anomaly"] = z_scores > 2

        # Evidence 5: Mitochondrial percentage (doublets often have lower MT%)
        if "pct_counts_mt" in self.adata.obs:
            mt_pct = self.adata.obs["pct_counts_mt"]
            # Doublets typically have LOWER MT% than singlets
            z_scores = (mt_pct - mt_pct.mean()) / mt_pct.std()
            evidence["mt_pct_zscore"] = z_scores
            evidence["low_mt_evidence"] = z_scores < -1  # Unusually low MT%

        # Combined evidence score (weighted combination)
        weights = {
            "scrublet_score": 0.3,
            "coexpression_strength": 0.3,
            "gene_count_zscore": 0.2,
            "umi_count_zscore": 0.1,
            "mt_pct_zscore": 0.1,  # Negative weight (lower is more suspicious)
        }

        evidence["combined_evidence_score"] = 0
        for feature, weight in weights.items():
            if feature in evidence.columns:
                # Normalize to [0, 1]
                normalized = (evidence[feature] - evidence[feature].min()) / (
                    evidence[feature].max() - evidence[feature].min() + 1e-10
                )
                if feature == "mt_pct_zscore":
                    normalized = 1 - normalized  # Invert for MT%
                evidence["combined_evidence_score"] += weight * normalized

        # Final classification with confidence
        evidence["doublet_confidence"] = pd.cut(
            evidence["combined_evidence_score"],
            bins=[0, 0.3, 0.5, 0.7, 1.0],
            labels=["Low", "Moderate", "High", "Very High"],
        )

        self.evidence_table = evidence
        return evidence

    def generate_doublet_report(self, cell_id: str, save_path: Optional[str] = None) -> str:
        """
        Generate a detailed textual report for a specific cell.

        Args:
            cell_id: Cell barcode
            save_path: Optional path to save the report

        Returns:
            Formatted report string
        """
        if self.evidence_table is None:
            self.generate_evidence_table()

        if cell_id not in self.evidence_table.index:
            raise ValueError(f"Cell {cell_id} not found")

        row = self.evidence_table.loc[cell_id]

        report = f"""
╔══════════════════════════════════════════════════════════════╗
║              DOUBLET EVIDENCE REPORT                         ║
║  Cell ID: {cell_id:<48}║
╚══════════════════════════════════════════════════════════════╝

┌──────────────────────────────────────────────────────────────┐
│ OVERALL ASSESSMENT                                           │
└──────────────────────────────────────────────────────────────┘
  Doublet Confidence: {row.get('doublet_confidence', 'N/A')}
  Combined Evidence Score: {row.get('combined_evidence_score', 0):.3f}

┌──────────────────────────────────────────────────────────────┐
│ EVIDENCE BREAKDOWN                                           │
└──────────────────────────────────────────────────────────────┘

1. ALGORITHMIC EVIDENCE
   • Scrublet Score: {row.get('scrublet_score', 0):.3f}
   • Strength: {row.get('scrublet_evidence', 'N/A')}

2. LINEAGE CO-EXPRESSION EVIDENCE
   • Number of Co-expressed Lineages: {row.get('n_coexpressed_lineages', 0)}
   • Top Co-expressed: {row.get('top_coexpressed_lineages', 'None')}
   • Co-expression Strength: {row.get('coexpression_strength', 0):.3f}

3. TRANSCRIPT COMPLEXITY EVIDENCE
   • Gene Count Z-score: {row.get('gene_count_zscore', 0):.2f}
   • Gene Count Anomaly: {'Yes' if row.get('gene_count_anomaly', False) else 'No'}
   • UMI Count Z-score: {row.get('umi_count_zscore', 0):.2f}
   • UMI Count Anomaly: {'Yes' if row.get('umi_count_anomaly', False) else 'No'}

4. QUALITY METRICS
   • MT% Z-score: {row.get('mt_pct_zscore', 0):.2f}
   • Low MT% Evidence: {'Yes' if row.get('low_mt_evidence', False) else 'No'}

┌──────────────────────────────────────────────────────────────┐
│ INTERPRETATION                                               │
└──────────────────────────────────────────────────────────────┘
"""

        # Add interpretation based on evidence
        if row.get("doublet_confidence") in ["High", "Very High"]:
            report += """
⚠️  This cell shows STRONG evidence of being a doublet:
"""
            if row.get("n_coexpressed_lineages", 0) >= 2:
                report += (
                    f"   • Co-expresses {row.get('n_coexpressed_lineages')} distinct lineages\n"
                )
                report += f"     ({row.get('top_coexpressed_lineages')})\n"

            if row.get("gene_count_anomaly", False):
                report += "   • Unusually high gene count (possible merged cells)\n"

            if row.get("scrublet_score", 0) > 0.5:
                report += "   • High algorithmic doublet score\n"

            report += "\n➤ RECOMMENDATION: Remove this cell from downstream analysis\n"

        elif row.get("doublet_confidence") == "Moderate":
            report += """
⚡ This cell shows MODERATE evidence of being a doublet:
   • Consider context-specific filtering
   • May be a transient cell state or true biological heterogeneity

➤ RECOMMENDATION: Review in biological context before filtering
"""
        else:
            report += """
✓ This cell shows LOW evidence of being a doublet:
   • Likely a true singlet

➤ RECOMMENDATION: Keep for downstream analysis
"""

        report += "\n" + "═" * 64 + "\n"

        if save_path:
            with open(save_path, "w") as f:
                f.write(report)
            log.info(f"Saved doublet report to {save_path}")

        return report

    def plot_evidence_heatmap(self, top_n: int = 100, save_path: Optional[str] = None):
        """
        Create a heatmap of evidence features for top doublets.
        """
        if self.evidence_table is None:
            self.generate_evidence_table()

        # Select top doublets by combined score
        top_doublets = self.evidence_table.nlargest(top_n, "combined_evidence_score")

        # Select numeric evidence columns
        evidence_cols = [
            "scrublet_score",
            "coexpression_strength",
            "gene_count_zscore",
            "umi_count_zscore",
            "mt_pct_zscore",
        ]
        evidence_cols = [col for col in evidence_cols if col in top_doublets.columns]

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.15)))

        # Normalize data for better visualization
        from sklearn.preprocessing import StandardScaler

        scaler = StandardScaler()
        normalized_data = scaler.fit_transform(top_doublets[evidence_cols])

        sns.heatmap(
            normalized_data,
            xticklabels=[col.replace("_", " ").title() for col in evidence_cols],
            yticklabels=False,  # Too many cells to label
            cmap="RdYlBu_r",
            center=0,
            cbar_kws={"label": "Standardized Score"},
            ax=ax,
        )

        ax.set_title(f"Evidence Heatmap for Top {top_n} Doublets")
        ax.set_xlabel("Evidence Type")
        ax.set_ylabel(f"Cells (n={top_n})")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            log.info(f"Saved evidence heatmap to {save_path}")

        return fig

    def export_evidence_summary(self, output_dir: str, top_n_reports: int = 50):
        """
        Export comprehensive evidence summaries.

        Creates:
        - evidence_table.csv: Full evidence table
        - top_doublets_reports/: Individual reports for top doublets
        - evidence_heatmap.png: Heatmap visualization
        """
        from pathlib import Path

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Export full table
        if self.evidence_table is None:
            self.generate_evidence_table()

        self.evidence_table.to_csv(output_path / "evidence_table.csv")
        log.info(f"Exported evidence table to {output_path / 'evidence_table.csv'}")

        # Generate individual reports for top doublets
        reports_dir = output_path / "top_doublets_reports"
        reports_dir.mkdir(exist_ok=True)

        top_doublets = self.evidence_table.nlargest(top_n_reports, "combined_evidence_score")

        for i, (cell_id, row) in enumerate(top_doublets.iterrows(), 1):
            report = self.generate_doublet_report(cell_id)
            report_path = reports_dir / f"rank_{i:03d}_{cell_id}.txt"
            with open(report_path, "w") as f:
                f.write(report)

        log.info(f"Generated {top_n_reports} individual reports in {reports_dir}")

        # Generate heatmap
        self.plot_evidence_heatmap(
            top_n=min(100, top_n_reports), save_path=output_path / "evidence_heatmap.png"
        )


def predict_doublets_with_profiling(
    adata: AnnData,
    config: DoubletConfig,
    sample_key: str = "sampleID",
    generate_reports: bool = True,
    top_n_reports: int = 50,
    **kwargs,
) -> AnnData:
    """
    Enhanced doublet prediction with evidence profiling.

    This wrapper adds biological interpretability to doublet predictions.
    """
    # Run standard doublet detection
    adata = predict_doublets(adata, config, sample_key, **kwargs)

    if generate_reports:
        log.info("Generating doublet evidence profiles...")

        profiler = DoubletEvidenceProfiler(adata)
        profiler.generate_evidence_table()

        # Export comprehensive reports
        if config.save_dir:
            profiler.export_evidence_summary(
                output_dir=Path(config.save_dir) / "evidence_profiles", top_n_reports=top_n_reports
            )

        # Add evidence table to AnnData
        adata.obs = adata.obs.join(
            profiler.evidence_table[["combined_evidence_score", "doublet_confidence"]]
        )

    return adata
