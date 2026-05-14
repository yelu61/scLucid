"""Tests for marker manager helpers and built-in resource composition."""

from scLucid.utils import Manager, get_marker_manager
from scLucid.utils.manager import load_gene_set_manager, load_gene_sets


def test_get_marker_manager_loads_scoped_state_resources():
    """Requested state markers should remain discoverable after resource reorganization."""
    mgr = get_marker_manager(
        species="human",
        states=["T cell exhaustion-like", "Stress-high"],
    )

    assert "T cell exhaustion-like" in mgr.CELLS
    assert "Stress-high" in mgr.CELLS
    assert mgr.CELLS["T cell exhaustion-like"].metadata["scope"] == "lineage_restricted"
    assert "T cells" in mgr.CELLS["T cell exhaustion-like"].metadata["applies_to"]
    assert mgr.CELLS["Stress-high"].metadata["scope"] == "all"


def test_marker_manager_parses_metadata_and_negative_markers():
    """Top-level metadata and negative markers should be first-class Manager data."""
    mgr = get_marker_manager(species="mouse")

    assert mgr.metadata["species"] == "mouse"
    assert "version" in mgr.metadata
    assert "metadata" not in mgr.CELLS
    assert "Epcam" in mgr.CELLS["Immune"].negative_markers
    assert "Ptprc" in mgr.CELLS["Epithelial"].negative_markers


def test_marker_manager_preserves_extra_definition_fields_as_metadata():
    """Resource fields such as cancer marker descriptions should not be dropped."""
    mgr = get_marker_manager(species="human", cancer_type="Lung Cancer")

    assert mgr.CELLS["LUAD"].metadata["description"] == "Lung adenocarcinoma"


def test_functional_signatures_load_as_marker_manager_view():
    """Functional signatures should be available through Manager instead of a separate class."""
    mgr = get_marker_manager(species="human", view="state_annotation", include_functional=True)

    assert "Cytotoxicity" in mgr.CELLS
    assert mgr.CELLS["Cytotoxicity"].metadata["kind"] == "functional_program"
    assert mgr.CELLS["Cytotoxicity"].metadata["category"] == "Immune_Function"


def test_geneset_resources_load_through_manager_helpers():
    """Legacy geneset JSON resources should still be loadable through Manager helpers."""
    genesets = load_gene_sets(species="human", name="cell_cycle")
    mgr = load_gene_set_manager(species="human", name="cancer_hallmarks", kind="geneset")

    assert "s_genes" in genesets
    assert "HALLMARK_APOPTOSIS" in mgr.CELLS


def test_new_marker_resource_scaffolds_are_manager_readable():
    """Resource scaffolds should stay compatible with the unified Manager."""
    for resource_name in ["cancer_human", "cell_state_mouse"]:
        mgr = Manager(resource_name, case_sensitive=True)
        assert len(mgr.CELLS) > 0

    cancer_mgr = Manager("cancer_human", case_sensitive=True)
    assert "Epithelial tumor identity" in cancer_mgr.CELLS
    assert cancer_mgr.CELLS["Epithelial tumor identity"].metadata["kind"] == "tumor_evidence"

    mouse_states = get_marker_manager(species="mouse", states=["Stress-high"])
    assert "Stress-high" in mouse_states.CELLS
