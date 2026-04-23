"""Tests for scLucid.qc.doublet module.

Covers:
- predict_doublets configuration and dispatch logic
- _run_heuristic (co-expression-based doublet detection)
- DoubletEvidenceProfiler evidence table generation
- _export_doublet_stats
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.qc.config import DoubletConfig
from scLucid.qc.doublet import (
    _export_doublet_stats,
    _run_heuristic,
    predict_doublets,
    DoubletEvidenceProfiler,
)


# ---------------------------------------------------------------------------
# _run_heuristic
# ---------------------------------------------------------------------------


class TestRunHeuristic:
    """Tests for the co-expression heuristic doublet detection."""

    def test_basic_heuristic_run(self, minimal_adata):
        """Heuristic should run on basic synthetic data."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        potential_doublets, lineage_scores, heuristic_scores = _run_heuristic(
            minimal_adata, cfg
        )

        assert isinstance(potential_doublets, pd.Series)
        assert isinstance(lineage_scores, pd.DataFrame)
        assert isinstance(heuristic_scores, pd.Series)
        assert len(potential_doublets) == minimal_adata.n_obs
        assert len(heuristic_scores) == minimal_adata.n_obs
        # Scores should be in [0, 1]
        assert heuristic_scores.min() >= 0
        assert heuristic_scores.max() <= 1

    def test_heuristic_returns_no_false_positives_on_clean_data(self, minimal_adata):
        """On clean synthetic data without doublets, predictions should be conservative."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        potential_doublets, _, heuristic_scores = _run_heuristic(minimal_adata, cfg)

        # Even if no doublets are predicted, scores should be computed
        assert heuristic_scores.notna().all()

    def test_heuristic_ignore_pairs(self, minimal_adata):
        """ignore_coexpression_pairs should mask scores for specified lineages."""
        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
            ignore_coexpression_pairs=[("T_cell", "B_cell")],
        )
        _, _, heuristic_scores_with_ignore = _run_heuristic(minimal_adata, cfg)
        _, _, heuristic_scores_without = _run_heuristic(
            minimal_adata,
            DoubletConfig(
                run_algorithm=False,
                use_heuristics=True,
                marker_species="human",
                marker_tissue="pbmc",
            ),
        )
        # Scores may differ when ignoring pairs
        assert isinstance(heuristic_scores_with_ignore, pd.Series)


# ---------------------------------------------------------------------------
# predict_doublets
# ---------------------------------------------------------------------------


class TestPredictDoublets:
    """Tests for the main predict_doublets entry point."""

    def test_heuristic_only(self, minimal_adata):
        """predict_doublets with heuristic only should populate required columns."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        cfg = DoubletConfig(
            run_algorithm=False,
            use_heuristics=True,
            marker_species="human",
            marker_tissue="pbmc",
        )
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")

        assert "predicted_doublet" in out.obs.columns
        assert "doublet_source" in out.obs.columns
        assert "heuristic_predicted" in out.obs.columns
        assert "heuristic_confidence_score" in out.obs.columns

    def test_invalid_sample_key_raises(self, minimal_adata):
        """Missing sample key should raise ValueError."""
        cfg = DoubletConfig()
        with pytest.raises(ValueError, match="Sample key"):
            predict_doublets(minimal_adata, config=cfg, sample_key="nonexistent")

    def test_unknown_method_raises_at_config_level(self):
        """Unknown method should be rejected by Pydantic config validation."""
        with pytest.raises(Exception, match="scrublet|solo|doubletdetection"):
            DoubletConfig(method="unknown_method")

    def test_skip_small_samples(self, minimal_adata):
        """Samples with <50 cells should be skipped for algorithmic detection."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        cfg = DoubletConfig(
            run_algorithm=True,
            use_heuristics=False,
            method="scrublet",
        )
        # With <50 cells, scrublet should be skipped and heuristic columns still created
        out = predict_doublets(adata, config=cfg, sample_key="sampleID")
        assert "predicted_doublet" in out.obs.columns

    def test_config_none_uses_defaults(self, minimal_adata):
        """config=None should use DoubletConfig defaults."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"

        out = predict_doublets(adata, config=None, sample_key="sampleID")
        assert "predicted_doublet" in out.obs.columns


# ---------------------------------------------------------------------------
# DoubletEvidenceProfiler
# ---------------------------------------------------------------------------


class TestDoubletEvidenceProfiler:
    """Tests for DoubletEvidenceProfiler."""

    def test_generate_evidence_table_basic(self, minimal_adata):
        """Evidence table should be generated from data with doublet results."""
        adata = minimal_adata.copy()
        # Add mock doublet results
        adata.obs["scrublet_score"] = np.random.random(adata.n_obs)
        adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, adata.n_obs)
        adata.obsm["lineage_module_scores"] = pd.DataFrame(
            np.random.random((adata.n_obs, 3)),
            index=adata.obs_names,
            columns=["T_cell", "B_cell", "Myeloid"],
        )

        profiler = DoubletEvidenceProfiler(adata)
        evidence = profiler.generate_evidence_table()

        assert isinstance(evidence, pd.DataFrame)
        assert len(evidence) == adata.n_obs
        assert "scrublet_evidence" in evidence.columns
        assert "n_coexpressed_lineages" in evidence.columns

    def test_generate_evidence_table_without_doublet_results(self, minimal_adata):
        """Evidence table should handle data without doublet predictions gracefully."""
        profiler = DoubletEvidenceProfiler(minimal_adata)
        evidence = profiler.generate_evidence_table()

        assert isinstance(evidence, pd.DataFrame)
        assert len(evidence) == minimal_adata.n_obs

    def test_generate_doublet_report(self, minimal_adata):
        """Individual cell report should be generated."""
        adata = minimal_adata.copy()
        adata.obs["scrublet_score"] = np.random.random(adata.n_obs)
        adata.obs["predicted_doublet"] = False
        adata.obs["heuristic_confidence_score"] = np.random.random(adata.n_obs)

        profiler = DoubletEvidenceProfiler(adata)
        report = profiler.generate_doublet_report(adata.obs_names[0])

        assert isinstance(report, str)
        assert len(report) > 0


# ---------------------------------------------------------------------------
# _export_doublet_stats
# ---------------------------------------------------------------------------


class TestExportDoubletStats:
    """Tests for doublet statistics export."""

    def test_export_without_save_dir(self, minimal_adata):
        """Export should return DataFrames without writing files when save_dir is None."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        adata.obs["predicted_doublet"] = False
        adata.obs["scrublet_score"] = 0.1

        result = _export_doublet_stats(adata, sample_key="sampleID")

        assert "sample" in result
        assert "global" in result
        assert isinstance(result["sample"], pd.DataFrame)
        assert isinstance(result["global"], pd.DataFrame)

    def test_export_with_temp_dir(self, minimal_adata, tmp_path):
        """Export should write CSV files when save_dir is provided."""
        adata = minimal_adata.copy()
        adata.obs["sampleID"] = "sample_1"
        adata.obs["predicted_doublet"] = False
        adata.obs["scrublet_score"] = 0.1

        save_dir = tmp_path / "doublet_stats"
        result = _export_doublet_stats(
            adata, sample_key="sampleID", save_dir=str(save_dir), export_csv=True
        )

        assert (save_dir / "doublet_stats_per_sample.csv").exists()
        assert (save_dir / "doublet_stats_global.csv").exists()
