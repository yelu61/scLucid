"""Smoke tests for previously untested modules."""

import pytest


class TestCellCycle:
    def test_score_cell_cycle_import(self):
        from scLucid.qc.cycle import score_cell_cycle

        assert callable(score_cell_cycle)

    def test_score_cell_cycle_with_explicit_genes(self, qc_test_adata):
        from scLucid.qc.cycle import score_cell_cycle

        # Use genes from the test data's var_names that exist
        avail = set(qc_test_adata.var_names)
        s_genes = [g for g in avail if g.startswith("marker_")][:3]
        g2m_genes = [g for g in avail if g.startswith("marker_")][3:6]
        if len(s_genes) >= 3 and len(g2m_genes) >= 3:
            result = score_cell_cycle(
                qc_test_adata,
                s_genes=s_genes,
                g2m_genes=g2m_genes,
                plot=False,
            )
            assert result is not None
        else:
            pytest.skip("Not enough marker genes in test data")

    def test_available_species(self):
        from scLucid.qc.cycle import _get_available_species

        species = _get_available_species()
        assert "human" in species
        assert "mouse" in species


class TestSettings:
    def test_set_figure_params(self):
        from scLucid.settings import set_figure_params

        set_figure_params(style="default", font_style="nature", dpi=100)
        import matplotlib

        assert matplotlib.rcParams["figure.dpi"] == 100

    def test_interactive_mode(self):
        from scLucid.settings import is_interactive_mode, set_interactive_mode

        orig = is_interactive_mode()
        set_interactive_mode(False)
        assert not is_interactive_mode()
        set_interactive_mode(orig)

    def test_setup_logging(self):
        from scLucid.settings import setup_logging

        setup_logging(level="WARNING")


class TestWorkflowCheckpoint:
    def test_imports(self):
        from scLucid.workflow_checkpoint import (
            StepStatus,
            WorkflowCheckpoint,
            WorkflowState,
            WorkflowStep,
        )

        assert WorkflowStep.QC is not None
        assert StepStatus.PENDING is not None

    def test_workflow_state_init(self):
        from scLucid.workflow_checkpoint import WorkflowState

        state = WorkflowState(workflow_id="test_001")
        assert state.workflow_id == "test_001"

    def test_checkpoint_info_init(self, temp_output_dir):
        from datetime import datetime
        from pathlib import Path

        from scLucid.workflow_checkpoint import CheckpointInfo

        info = CheckpointInfo(
            step="qc",
            timestamp=datetime.now(),
            adata_path=Path(temp_output_dir) / "test.h5ad",
            config_hash="abc123",
            metadata={"min_genes": 200},
        )
        assert info.step == "qc"
        assert info.metadata["min_genes"] == 200

    def test_workflow_checkpoint_init(self, temp_output_dir):
        from scLucid.workflow_checkpoint import WorkflowCheckpoint

        ckpt = WorkflowCheckpoint(
            results_dir=temp_output_dir,
            workflow_id="test_wf",
        )
        assert str(ckpt.results_dir) == temp_output_dir
        assert ckpt.workflow_id == "test_wf"
