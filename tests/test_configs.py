"""
Tests for Pydantic-based configuration system.

Validates all config classes have correct validation and serialization.
"""

import json
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

# Analysis Configs
from scLucid.analysis.config import (
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    EnrichmentConfig,
)
from scLucid.base_config import SclucidBaseConfig, WorkflowConfigBase

# Preprocess Configs
from scLucid.preprocess.config import (
    HVGConfig,
    NormalizationConfig,
    PreprocessingWorkflowConfig,
)

# QC Configs
from scLucid.qc.config import (
    DoubletConfig,
    QCThresholds,
    QCWorkflowConfig,
)


@pytest.mark.unit
class TestBaseConfig:
    """Test base configuration functionality."""

    def test_base_config_creation(self):
        """Test creating base config."""
        config = SclucidBaseConfig(save_dir="./test_output", verbose=True)
        assert config.save_dir == "./test_output"
        assert config.verbose is True

    def test_base_config_to_dict(self):
        """Test serialization to dict."""
        config = SclucidBaseConfig(save_dir="./test", verbose=False)
        d = config.to_dict()

        assert isinstance(d, dict)
        assert d["save_dir"] == "./test"
        assert d["verbose"] is False

    def test_base_config_to_json(self):
        """Test serialization to JSON."""
        config = SclucidBaseConfig(save_dir="./test")
        json_str = config.to_json()

        assert isinstance(json_str, str)
        data = json.loads(json_str)
        assert data["save_dir"] == "./test"

    def test_base_config_from_dict(self):
        """Test deserialization from dict."""
        data = {"save_dir": "./from_dict", "verbose": True, "plot": True}
        config = SclucidBaseConfig.from_dict(data)

        assert config.save_dir == "./from_dict"
        assert config.verbose is True

    def test_base_config_from_json(self):
        """Test deserialization from JSON."""
        json_str = '{"save_dir": "./from_json", "verbose": false}'
        config = SclucidBaseConfig.from_json(json_str)

        assert config.save_dir == "./from_json"
        assert config.verbose is False

    def test_base_config_save_load_file(self):
        """Test saving and loading from file."""
        config = SclucidBaseConfig(save_dir="./test_save")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            config.to_json_file(temp_path)
            loaded = SclucidBaseConfig.from_json_file(temp_path)

            assert loaded.save_dir == "./test_save"
        finally:
            Path(temp_path).unlink(missing_ok=True)


@pytest.mark.unit
class TestWorkflowConfigBase:
    """Test workflow base configuration."""

    def test_n_jobs_validation(self):
        """Test n_jobs validation."""
        # Valid values
        config = WorkflowConfigBase(n_jobs=-1)
        assert config.n_jobs == -1

        config = WorkflowConfigBase(n_jobs=4)
        assert config.n_jobs == 4

        # Invalid: n_jobs=0 should raise error
        with pytest.raises(ValueError):
            WorkflowConfigBase(n_jobs=0)

        # Invalid: negative other than -1
        with pytest.raises(ValueError):
            WorkflowConfigBase(n_jobs=-2)


@pytest.mark.unit
class TestQCThresholds:
    """Test QC thresholds configuration."""

    def test_valid_thresholds(self):
        """Test creating valid thresholds."""
        config = QCThresholds(min_genes=200, max_genes=5000, pc_mt=20.0)
        assert config.min_genes == 200
        assert config.max_genes == 5000
        assert config.pc_mt == 20.0

    def test_invalid_gene_thresholds(self):
        """Test that min > max raises error."""
        with pytest.raises(ValueError):
            QCThresholds(min_genes=1000, max_genes=500)

    def test_invalid_count_thresholds(self):
        """Test that min > max for counts raises error."""
        with pytest.raises(ValueError):
            QCThresholds(min_counts=5000, max_counts=1000)

    def test_pc_mt_range(self):
        """Test that pc_mt must be 0-100."""
        with pytest.raises(ValueError):
            QCThresholds(pc_mt=150)

        with pytest.raises(ValueError):
            QCThresholds(pc_mt=-5)

    def test_nmads_positive(self):
        """Test that nmads must be positive."""
        with pytest.raises(ValueError):
            QCThresholds(nmads=0)

        with pytest.raises(ValueError):
            QCThresholds(nmads=-1)


@pytest.mark.unit
class TestDoubletConfig:
    """Test doublet detection configuration."""

    def test_valid_doublet_config(self):
        """Test creating valid config."""
        config = DoubletConfig(method="scrublet", scr_n_pcs=30)
        assert config.method == "scrublet"
        assert config.scr_n_pcs == 30

    def test_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError):
            DoubletConfig(method="invalid_method")

    def test_scr_n_pcs_validation(self):
        """Test scr_n_pcs must be > 1."""
        with pytest.raises(ValueError):
            DoubletConfig(method="scrublet", scr_n_pcs=1)

    def test_dd_n_components_validation(self):
        """Test dd_n_components must be > 1."""
        with pytest.raises(ValueError):
            DoubletConfig(method="doubletdetection", dd_n_components=1)

    def test_solo_epochs_validation(self):
        """Test solo_n_epochs must be >= 100."""
        with pytest.raises(ValueError):
            DoubletConfig(method="solo", solo_n_epochs=50)

    def test_expected_doublet_rate_validation(self):
        """Test expected doublet rate must be 0-1."""
        with pytest.raises(ValueError):
            DoubletConfig(expected_doublet_rate=1.5)

        with pytest.raises(ValueError):
            DoubletConfig(expected_doublet_rate=-0.1)

        # Valid: dict of rates
        config = DoubletConfig(expected_doublet_rate={"sample1": 0.05, "sample2": 0.08})
        assert config.expected_doublet_rate["sample1"] == 0.05


@pytest.mark.unit
class TestNormalizationConfig:
    """Test normalization configuration."""

    def test_valid_normalization(self):
        """Test valid config."""
        config = NormalizationConfig(method="standard", target_sum=1e4)
        assert config.method == "standard"
        assert config.target_sum == 1e4

    def test_invalid_target_sum(self):
        """Test target_sum must be positive."""
        with pytest.raises(ValueError):
            NormalizationConfig(target_sum=0)

        with pytest.raises(ValueError):
            NormalizationConfig(target_sum=-1000)

    def test_max_fraction_range(self):
        """Test max_fraction must be 0-1."""
        with pytest.raises(ValueError):
            NormalizationConfig(max_fraction=1.5)

    def test_reserved_layer_names(self):
        """Test that reserved layer names are rejected."""
        with pytest.raises(ValueError):
            NormalizationConfig(output_layer="X")

        with pytest.raises(ValueError):
            NormalizationConfig(output_layer="raw")


@pytest.mark.unit
class TestHVGConfig:
    """Test HVG configuration."""

    def test_valid_hvg_config(self):
        """Test valid config."""
        config = HVGConfig(n_top_genes=2000, method="scanpy")
        assert config.n_top_genes == 2000

    def test_n_top_genes_range(self):
        """Test n_top_genes bounds."""
        with pytest.raises(ValueError):
            HVGConfig(n_top_genes=50)  # Below minimum

        with pytest.raises(ValueError):
            HVGConfig(n_top_genes=25000)  # Above maximum

    def test_span_validation(self):
        """Test span must be in valid range."""
        with pytest.raises(ValueError):
            HVGConfig(span=0.005)  # Too small

        with pytest.raises(ValueError):
            HVGConfig(span=1.5)  # Too large


@pytest.mark.unit
class TestClusteringConfig:
    """Test clustering configuration."""

    def test_valid_clustering(self):
        """Test valid config."""
        config = ClusteringConfig(method="leiden", resolution=1.0)
        assert config.method == "leiden"
        assert config.resolution == 1.0

    def test_invalid_method(self):
        """Test invalid method."""
        with pytest.raises(ValueError):
            ClusteringConfig(method="invalid")

    def test_resolution_positive(self):
        """Test resolution must be positive."""
        with pytest.raises(ValueError):
            ClusteringConfig(resolution=-0.5)

    def test_n_clusters_for_kmeans(self):
        """Test n_clusters must be >= 2 for kmeans."""
        with pytest.raises(ValueError):
            ClusteringConfig(method="kmeans", n_clusters=1)


@pytest.mark.unit
class TestAnnotationConfig:
    """Test annotation configuration."""

    def test_valid_annotation(self):
        """Test valid config."""
        config = AnnotationConfig(run_celltypist=False, final_method="combined")
        assert config.run_celltypist is False
        assert config.final_method == "combined"

    def test_valid_hierarchical_annotation(self):
        """Hierarchical annotation config should accept lineage/state-specific fields."""
        config = AnnotationConfig(
            final_method="hierarchical",
            target_lineage="T cells",
            lineage_marker_config="base_human",
            state_signature_names=["T_cell_activation"],
            nomenclature_style="modular",
        )
        assert config.final_method == "hierarchical"
        assert config.target_lineage == "T cells"
        assert config.nomenclature_style == "modular"

    def test_min_confidence_range(self):
        """Test min_confidence must be 0-1."""
        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=1.5)

        with pytest.raises(ValueError):
            AnnotationConfig(min_confidence=-0.1)

    def test_invalid_final_method(self):
        """Test invalid final_method."""
        with pytest.raises(ValueError):
            AnnotationConfig(final_method="invalid")


@pytest.mark.unit
class TestEnrichmentConfig:
    """Test enrichment configuration."""

    def test_valid_enrichment(self):
        """Test valid config."""
        config = EnrichmentConfig(method="ora", organism="human")
        assert config.method == "ora"
        assert config.organism == "human"

    def test_max_padj_range(self):
        """Test max_padj must be 0-1."""
        with pytest.raises(ValueError):
            EnrichmentConfig(max_padj=1.5)

    def test_gsea_permutations_minimum(self):
        """Test gsea_permutations must be >= 100."""
        with pytest.raises(ValueError):
            EnrichmentConfig(gsea_permutations=50)


@pytest.mark.unit
class TestConfigSerialization:
    """Test config serialization round-trips."""

    def test_qc_workflow_roundtrip(self):
        """Test QC workflow config round-trip."""
        original = QCWorkflowConfig(
            sample_key="sample",
            species="mouse",
            thresholds=QCThresholds(min_genes=300),
        )

        # Serialize and deserialize
        json_str = original.to_json()
        restored = QCWorkflowConfig.from_json(json_str)

        assert restored.sample_key == original.sample_key
        assert restored.species == original.species

    def test_preprocess_workflow_roundtrip(self):
        """Test preprocess workflow config round-trip."""
        original = PreprocessingWorkflowConfig(
            n_jobs=4,
            normalization=NormalizationConfig(target_sum=1e4),
        )

        json_str = original.to_json()
        restored = PreprocessingWorkflowConfig.from_json(json_str)

        assert restored.n_jobs == original.n_jobs

    def test_analysis_workflow_roundtrip(self):
        """Test analysis workflow config round-trip."""
        original = AnalysisWorkflowConfig(
            clustering=ClusteringConfig(resolution=0.8),
        )

        json_str = original.to_json()
        restored = AnalysisWorkflowConfig.from_json(json_str)

        assert restored.clustering.resolution == 0.8


@pytest.mark.unit
class TestConfigWithExtras:
    """Test that configs handle extra fields gracefully."""

    def test_extra_fields_ignored(self):
        """Test that extra fields are ignored (not rejected)."""
        # This should not raise an error due to extra="ignore"
        config = SclucidBaseConfig(
            save_dir="./test",
            unknown_field="should_be_ignored",
        )
        assert config.save_dir == "./test"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
