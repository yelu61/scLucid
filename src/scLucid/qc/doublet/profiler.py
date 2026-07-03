"""Doublet evidence profiling utilities.

This module provides interpretable, cell-level evidence tables and reports for
doublet predictions. It is kept separate from :mod:`ensemble` so that the core
merging logic does not depend on plotting/reporting utilities.
"""

from __future__ import annotations

import json
import logging
import warnings
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData

from ..config import DoubletConfig
from .core import (
    ALGORITHM_SCORE_COL,
    COMBINED_SCORE_COL,
    EXPECTED_HETEROTYPIC_RATE_COL,
    EXPECTED_HOMOTYPIC_RATE_COL,
    EXPECTED_TOTAL_RATE_COL,
    HETEROTYPIC_RISK_COL,
    HEURISTIC_SCORE_COL,
    HOMOTYPIC_RISK_COL,
)
from .ensemble import predict_doublets

log = logging.getLogger(__name__)

__all__ = ["DoubletEvidenceProfiler", "predict_doublets_with_profiling"]


# --- Helper Functions ---
def _safe_zscore(values: pd.Series) -> pd.Series:
    """Return z-scores, falling back to zeros when variance is zero."""
    std = values.std()
    if pd.isna(std) or std == 0:
        return pd.Series(0.0, index=values.index)
    return (values - values.mean()) / std


# --- Main Functions ---
class DoubletEvidenceProfiler:
    """
    Generate interpretable evidence profiles for doublet predictions.

    This class creates detailed reports explaining WHY each cell was
    flagged as a doublet, combining multiple lines of evidence.
    """

    def __init__(self, adata: AnnData):
        self.adata = adata
        self.evidence_table: Optional[pd.DataFrame] = None

    def generate_evidence_table(self) -> pd.DataFrame:
        """
        Create a comprehensive evidence table for each cell.

        Returns:
            DataFrame with one row per cell, columns for different evidence types
        """
        evidence = pd.DataFrame(index=self.adata.obs_names)

        # Evidence 1: Algorithmic score
        if ALGORITHM_SCORE_COL in self.adata.obs:
            evidence[ALGORITHM_SCORE_COL] = self.adata.obs[ALGORITHM_SCORE_COL]
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
                lambda row: np.prod(sorted(row.values)[-2:]) if row.max() > threshold else 0,
                axis=1,
            )

        if HEURISTIC_SCORE_COL in self.adata.obs:
            evidence["heuristic_evidence_score"] = self.adata.obs[HEURISTIC_SCORE_COL]
        if COMBINED_SCORE_COL in self.adata.obs:
            evidence[COMBINED_SCORE_COL] = self.adata.obs[COMBINED_SCORE_COL]
        if HETEROTYPIC_RISK_COL in self.adata.obs:
            evidence[HETEROTYPIC_RISK_COL] = self.adata.obs[HETEROTYPIC_RISK_COL]
        if HOMOTYPIC_RISK_COL in self.adata.obs:
            evidence[HOMOTYPIC_RISK_COL] = self.adata.obs[HOMOTYPIC_RISK_COL]
        if EXPECTED_TOTAL_RATE_COL in self.adata.obs:
            evidence[EXPECTED_TOTAL_RATE_COL] = self.adata.obs[EXPECTED_TOTAL_RATE_COL]
            evidence[EXPECTED_HETEROTYPIC_RATE_COL] = self.adata.obs.get(
                EXPECTED_HETEROTYPIC_RATE_COL, np.nan
            )
            evidence[EXPECTED_HOMOTYPIC_RATE_COL] = self.adata.obs.get(
                EXPECTED_HOMOTYPIC_RATE_COL, np.nan
            )
        if "external_doublet_evidence" in self.adata.obs:
            evidence["external_doublet_evidence"] = (
                self.adata.obs["external_doublet_evidence"].fillna(False).astype(bool)
            )

        # Evidence 3: Gene count anomaly
        if "n_genes_by_counts" in self.adata.obs:
            gene_counts = self.adata.obs["n_genes_by_counts"]
            z_scores = _safe_zscore(gene_counts)
            evidence["gene_count_zscore"] = z_scores
            evidence["gene_count_anomaly"] = z_scores > 2  # High gene count

        # Evidence 4: Total UMI anomaly
        if "total_counts" in self.adata.obs:
            umi_counts = self.adata.obs["total_counts"]
            z_scores = _safe_zscore(umi_counts)
            evidence["umi_count_zscore"] = z_scores
            evidence["umi_count_anomaly"] = z_scores > 2

        # Evidence 5: Mitochondrial percentage is descriptive QC context only.
        if "pct_counts_mt" in self.adata.obs:
            mt_pct = self.adata.obs["pct_counts_mt"]
            z_scores = _safe_zscore(mt_pct)
            evidence["mt_pct_zscore"] = z_scores

        # Combined evidence score (weighted combination)
        weights = {
            ALGORITHM_SCORE_COL: 0.25,
            "scrublet_score": 0.15,
            "heuristic_evidence_score": 0.20,
            HETEROTYPIC_RISK_COL: 0.20,
            HOMOTYPIC_RISK_COL: 0.15,
            "gene_count_zscore": 0.05,
        }

        evidence["combined_evidence_score"] = 0
        for feature, weight in weights.items():
            if feature in evidence.columns:
                # Normalize to [0, 1]
                normalized = (evidence[feature] - evidence[feature].min()) / (
                    evidence[feature].max() - evidence[feature].min() + 1e-10
                )
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
  Doublet Confidence: {row.get("doublet_confidence", "N/A")}
  Combined Evidence Score: {row.get("combined_evidence_score", 0):.3f}

┌──────────────────────────────────────────────────────────────┐
│ EVIDENCE BREAKDOWN                                           │
└──────────────────────────────────────────────────────────────┘

1. ALGORITHMIC EVIDENCE
   • Algorithm Score: {row.get(ALGORITHM_SCORE_COL, row.get("scrublet_score", 0)):.3f}
   • Scrublet Score: {row.get("scrublet_score", 0):.3f}
   • Strength: {row.get("scrublet_evidence", "N/A")}

2. LINEAGE CO-EXPRESSION EVIDENCE
   • Number of Co-expressed Lineages: {row.get("n_coexpressed_lineages", 0)}
   • Top Co-expressed: {row.get("top_coexpressed_lineages", "None")}
   • Co-expression Strength: {row.get("coexpression_strength", 0):.3f}
   • Heuristic Evidence Score: {row.get("heuristic_evidence_score", 0):.3f}

3. TRANSCRIPT COMPLEXITY EVIDENCE
   • Gene Count Z-score: {row.get("gene_count_zscore", 0):.2f}
   • Gene Count Anomaly: {"Yes" if row.get("gene_count_anomaly", False) else "No"}
   • UMI Count Z-score: {row.get("umi_count_zscore", 0):.2f}
   • UMI Count Anomaly: {"Yes" if row.get("umi_count_anomaly", False) else "No"}

4. RISK DECOMPOSITION
   • Heterotypic Risk: {row.get(HETEROTYPIC_RISK_COL, 0):.3f}
   • Homotypic Risk: {row.get(HOMOTYPIC_RISK_COL, 0):.3f}
   • Combined Score: {row.get(COMBINED_SCORE_COL, row.get("combined_evidence_score", 0)):.3f}
   • Scores are evidence indices, not calibrated probabilities.

5. QUALITY METRICS
   • MT% Z-score: {row.get("mt_pct_zscore", 0):.2f}
   • MT% is shown as descriptive QC context, not as doublet evidence.

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
            ALGORITHM_SCORE_COL,
            "heuristic_evidence_score",
            HETEROTYPIC_RISK_COL,
            HOMOTYPIC_RISK_COL,
            COMBINED_SCORE_COL,
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

    def export_evidence_summary(
        self,
        output_dir: str,
        top_n_reports: int = 50,
        max_table_rows: Optional[int] = 100_000,
    ):
        """
        Export comprehensive evidence summaries.

        Creates:
        - evidence_table.csv: Evidence table, capped to the highest-risk rows
          when ``max_table_rows`` is set.
        - evidence_export_summary.json: Export provenance and truncation status.
        - top_doublets_reports/: Individual reports for top doublets
        - evidence_heatmap.png: Heatmap visualization
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Export full table
        if self.evidence_table is None:
            self.generate_evidence_table()

        total_rows = int(len(self.evidence_table))
        export_table = self.evidence_table
        truncated = False
        if max_table_rows is not None and total_rows > max_table_rows:
            truncated = True
            export_table = self.evidence_table.nlargest(max_table_rows, "combined_evidence_score")
        export_table.to_csv(output_path / "evidence_table.csv")
        export_summary = {
            "schema_version": "doublet_evidence_export_summary_v1",
            "total_rows": total_rows,
            "exported_rows": int(len(export_table)),
            "truncated": truncated,
            "max_table_rows": max_table_rows,
            "sort_key": "combined_evidence_score" if truncated else None,
        }
        (output_path / "evidence_export_summary.json").write_text(
            json.dumps(export_summary, indent=2),
            encoding="utf-8",
        )
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
    Transitional convenience wrapper for doublet prediction with evidence profiling.

    Prefer calling :func:`predict_doublets` and then instantiating
    :class:`DoubletEvidenceProfiler` explicitly when building new workflows.
    """
    warnings.warn(
        "predict_doublets_with_profiling is a transitional convenience wrapper. "
        "Use predict_doublets plus DoubletEvidenceProfiler for new code.",
        FutureWarning,
        stacklevel=2,
    )
    # Run standard doublet detection
    adata = predict_doublets(adata, config, sample_key, **kwargs)

    if generate_reports:
        log.info("Generating doublet evidence profiles...")

        profiler = DoubletEvidenceProfiler(adata)
        profiler.generate_evidence_table()

        # Export comprehensive reports
        if config.save_dir:
            profiler.export_evidence_summary(
                output_dir=Path(config.save_dir) / "evidence_profiles",
                top_n_reports=top_n_reports,
            )

        # Add evidence table to AnnData
        adata.obs = adata.obs.join(
            profiler.evidence_table[["combined_evidence_score", "doublet_confidence"]]
        )

    return adata
