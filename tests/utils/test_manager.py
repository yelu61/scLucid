"""Tests for marker manager helpers and built-in resource composition."""

import re

from scLucid.utils import Manager, get_marker_manager
from scLucid.utils.manager import load_gene_set_manager, load_gene_sets


MARKER_SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


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
    mgr = get_marker_manager(
        species="human",
        cancer_type="Lung Cancer",
        view="tumor_interpretation",
    )

    assert mgr.CELLS["LUAD"].metadata["description"] == "Lung adenocarcinoma"


def test_functional_signatures_load_as_marker_manager_view():
    """Functional signatures should be available through Manager instead of a separate class."""
    mgr = get_marker_manager(species="human", view="state_annotation", include_functional=True)

    assert "Cytotoxicity" in mgr.CELLS
    assert mgr.CELLS["Cytotoxicity"].metadata["kind"] == "functional_program"
    assert mgr.CELLS["Cytotoxicity"].metadata["category"] == "immune_function"


def test_geneset_resources_load_through_manager_helpers():
    """Legacy geneset JSON resources should still be loadable through Manager helpers."""
    genesets = load_gene_sets(species="human", name="cell_cycle")
    mgr = load_gene_set_manager(species="human", name="cancer_hallmarks", kind="geneset")

    assert "s_genes" in genesets
    assert "HALLMARK_APOPTOSIS" in mgr.CELLS


def test_new_marker_resource_scaffolds_are_manager_readable():
    """Resource scaffolds should stay compatible with the unified Manager."""
    for resource_name in ["registry_human", "registry_mouse", "tumor_human", "tissue_human"]:
        mgr = Manager(resource_name, case_sensitive=True)
        assert len(mgr.CELLS) > 0

    cancer_mgr = Manager("tumor_human", case_sensitive=True)
    assert "Epithelial tumor identity" in cancer_mgr.CELLS
    assert cancer_mgr.CELLS["Epithelial tumor identity"].metadata["kind"] == "tumor_evidence"

    mouse_states = get_marker_manager(species="mouse", states=["Stress-high"])
    assert "Stress-high" in mouse_states.CELLS


def test_marker_resources_do_not_store_display_aliases_as_gene_symbols():
    """Built-in marker resources should keep display aliases out of marker lists."""
    for resource_name in [
        "registry_human",
        "registry_mouse",
        "tissue_human",
        "tumor_human",
    ]:
        mgr = Manager(resource_name, case_sensitive=True)
        for cell in mgr.CELLS.values():
            for gene in list(cell.markers) + list(cell.negative_markers):
                assert gene == gene.strip(), f"{resource_name}:{cell.name}:{gene!r}"
                assert MARKER_SYMBOL_RE.match(gene), f"{resource_name}:{cell.name}:{gene!r}"


def test_marker_resource_metadata_routes_context_away_from_global_annotation():
    """Cancer/tissue context entries should not appear as global cell-type labels."""
    cancer_mgr = get_marker_manager(
        species="human",
        cancer_type="Lung Cancer",
        view="global_annotation",
    )
    tumor_mgr = get_marker_manager(
        species="human",
        cancer_type="Lung Cancer",
        view="tumor_interpretation",
    )

    assert "Lung Cancer" not in cancer_mgr.CELLS
    assert "LUAD" not in cancer_mgr.CELLS
    assert "LUSC" not in cancer_mgr.CELLS
    assert "SCLC" not in cancer_mgr.CELLS
    assert "Lung Cancer" in tumor_mgr.CELLS
    assert "LUSC" in tumor_mgr.CELLS
    assert tumor_mgr.CELLS["Lung Cancer"].metadata["kind"] == "cancer_context"
    assert tumor_mgr.CELLS["LUSC"].metadata["kind"] == "cancer_context"


def test_refactored_marker_resource_metadata_contract():
    """Core resources should expose species/schema metadata and key routing fields."""
    human = Manager("registry_human", case_sensitive=True)
    tissue = Manager("tissue_human", root_key="Pancreas", case_sensitive=True)
    cancer = Manager("tumor_human", case_sensitive=True)

    assert human.metadata["species"] == "human"
    assert human.metadata["schema"] == "scLucid_marker_registry_v2"
    assert human.CELLS["Immune"].metadata["kind"] == "cell_type"
    assert human.CELLS["Epithelial"].metadata["use_for_malignancy_interpretation"] is True
    assert "PTPRC" in human.CELLS["Epithelial"].negative_markers

    assert tissue.metadata["species"] == "human"
    assert tissue.CELLS["Pancreas Tissue"].metadata["kind"] == "tissue_context"
    assert tissue.CELLS["Pancreas Tissue"].metadata["use_for_global_annotation"] is False
    assert tissue.CELLS["Acinar cells"].metadata["kind"] == "tissue_context"

    assert cancer.metadata["schema"] == "scLucid_marker_tumor_resource_v2"
    assert cancer.CELLS["Ovarian Cancer"].markers == [
        "PAX8",
        "WT1",
        "MUC16",
        "WFDC2",
        "EPCAM",
        "KRT7",
    ]
    assert "CA19-9" in cancer.CELLS["Pancreatic Cancer"].metadata[
        "excluded_non_gene_markers"
    ]


def test_registry_views_separate_identity_state_artifact_and_tumor_layers():
    """Unified registries should expose explicit views for each annotation layer."""
    compartment = get_marker_manager(species="human", view="compartment_annotation")
    lineage = get_marker_manager(species="human", view="lineage_annotation")
    subtype = get_marker_manager(species="human", tissue="Lung", view="subtype_annotation")
    state = get_marker_manager(species="human", view="state_annotation")
    program = get_marker_manager(species="human", view="program_scoring")
    artifact = get_marker_manager(species="human", view="artifact_annotation")

    assert {"Immune", "Stromal", "Neural"}.issubset(compartment.CELLS)
    assert "Epithelial" in lineage.CELLS
    assert "Lung Goblet cell" in subtype.CELLS
    assert "Cytotoxicity" in program.CELLS
    assert "Treg program" in program.CELLS
    assert program.CELLS["Treg program"].metadata["alias_of"] == "Treg"
    assert "Stress-high" in artifact.CELLS
    assert "Ribosomal-high" in artifact.CELLS
