"""Tests for marker manager helpers and built-in resource composition."""

import json
import re
from pathlib import Path

import tomllib

from scLucid.utils import (
    Manager,
    canonicalize_marker_label,
    get_gene_display_aliases,
    get_marker_aliases,
    get_marker_manager,
)
from scLucid.utils.manager import load_gene_set_manager, load_gene_sets

MARKER_SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
RESOURCE_DIR = Path(__file__).parents[2] / "src" / "scLucid" / "resources"
DOCS_DIR = Path(__file__).parents[2] / "docs"


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
    assert mgr.CELLS["HALLMARK_APOPTOSIS"].metadata["granularity"] == "program"
    assert mgr.CELLS["HALLMARK_APOPTOSIS"].metadata["review_status"] == "needs_review"


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
    assert human.CELLS["Immune"].metadata["granularity"] == "compartment"
    assert human.CELLS["T cells"].metadata["granularity"] == "lineage"
    assert human.CELLS["CD8+ T"].metadata["granularity"] == "subtype"
    assert human.CELLS["Epithelial"].metadata["use_for_malignancy_interpretation"] is True
    assert "PTPRC" in human.CELLS["Epithelial"].negative_markers

    assert tissue.metadata["species"] == "human"
    assert tissue.CELLS["Pancreas Tissue"].metadata["kind"] == "tissue_context"
    assert tissue.CELLS["Pancreas Tissue"].metadata["use_for_global_annotation"] is False
    assert tissue.metadata["schema"] == "scLucid_marker_tissue_resource_v2"
    assert tissue.CELLS["Acinar cells"].metadata["kind"] == "cell_type"
    assert tissue.CELLS["Acinar cells"].metadata["granularity"] == "tissue_subtype"
    assert tissue.CELLS["Acinar cells"].metadata["use_for_global_annotation"] is True

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
    state_view = get_marker_manager(species="human", view="state_annotation")
    program = get_marker_manager(species="human", view="program_scoring")
    artifact = get_marker_manager(species="human", view="artifact_annotation")

    assert {"Immune", "Stromal", "Neural"}.issubset(compartment.CELLS)
    assert "Epithelial" in lineage.CELLS
    assert "Lung Goblet cell" in subtype.CELLS
    assert "T cell exhaustion-like" in state_view.CELLS
    assert "CD8+ T | exhausted-like" in state_view.CELLS
    assert "CD8+ T | progenitor exhausted-like" in state_view.CELLS
    assert "CD8+ T | terminal exhausted-like" in state_view.CELLS
    assert "CD8+ T | exhausted-like" not in subtype.CELLS
    assert "Macro_NLRP3" not in subtype.CELLS
    assert "Macro_ISG15" not in subtype.CELLS
    assert "CD4+ Th1" in subtype.CELLS
    assert "CD4+ Th2" in subtype.CELLS
    assert "Tfr" in subtype.CELLS
    assert "AtM/TAAB B" in subtype.CELLS
    assert "FOLR2+ resident macrophage" in subtype.CELLS
    assert "PI16+ fibroblast" in subtype.CELLS
    assert "Cytotoxicity" in program.CELLS
    assert "Treg program" in program.CELLS
    assert "T cell activation" in program.CELLS
    assert "T_cell_activation" not in program.CELLS
    assert "HALLMARK_APOPTOSIS" in program.CELLS
    assert "EMT_epithelial" in program.CELLS
    assert "T cells" not in program.CELLS
    assert program.CELLS["Treg program"].metadata["alias_of"] == "Treg"
    assert "Stress-high" in artifact.CELLS
    assert "Ribosomal-high" in artifact.CELLS


def test_marker_resources_expose_manager_routing_contract():
    """Built-in marker resources should declare schemas that match manager routing."""
    resources = Path(__file__).parents[2] / "src" / "scLucid" / "resources"
    expected = {
        "marker_registry_human.toml": "scLucid_marker_registry_v2",
        "marker_registry_mouse.toml": "scLucid_marker_registry_v2",
        "marker_tissue_human.toml": "scLucid_marker_tissue_resource_v2",
        "marker_tumor_human.toml": "scLucid_marker_tumor_resource_v2",
    }

    for filename, schema in expected.items():
        data = tomllib.loads((resources / filename).read_text())
        metadata = data["metadata"]
        assert metadata["schema"] == schema
        assert metadata["species"] in {"human", "mouse"}
        assert "resource_type" in metadata
        assert "curation_status" in metadata


def test_lineage_negative_markers_are_systematic():
    """Major marker-manager lineages should include exclusion markers for conflicts."""
    human = Manager("registry_human", case_sensitive=True)
    mouse = Manager("registry_mouse", case_sensitive=True)

    assert {"MS4A1", "LYZ", "EPCAM", "PECAM1"}.issubset(
        set(human.CELLS["T cells"].negative_markers)
    )
    assert {"CD3D", "NKG7", "LYZ", "EPCAM"}.issubset(
        set(human.CELLS["B cells"].negative_markers)
    )
    assert {"PTPRC", "EPCAM", "COL1A1"}.issubset(
        set(human.CELLS["Endothelial cells"].negative_markers)
    )
    assert {"Ptprc", "Epcam", "Col1a1"}.issubset(
        set(mouse.CELLS["Endothelial"].negative_markers)
    )


def test_tissue_and_tumor_children_route_to_specific_views():
    """Manager inference should route resource children without hand-coded labels."""
    subtype = get_marker_manager(species="human", tissue="Pancreas", view="subtype_annotation")
    tumor = get_marker_manager(
        species="human",
        cancer_type="Lung Cancer",
        view="tumor_interpretation",
    )

    assert "Acinar cells" in subtype.CELLS
    assert subtype.CELLS["Acinar cells"].metadata["granularity"] == "tissue_subtype"
    assert "LUSC" in tumor.CELLS
    assert tumor.CELLS["LUSC"].metadata["granularity"] == "cancer_subtype"
    assert tumor.CELLS["LUSC"].metadata["use_for_global_annotation"] is False


def test_marker_alias_resource_keeps_nomenclature_out_of_marker_lists():
    """Alias resource should store label/gene synonyms without becoming marker evidence."""
    alias_data = tomllib.loads((RESOURCE_DIR / "marker_aliases.toml").read_text())

    assert alias_data["metadata"]["schema"] == "scLucid_marker_aliases_v1"
    assert any(
        item["canonical"] == "T cell exhaustion-like"
        and "Tex" in item["aliases"]
        for item in alias_data["label_aliases"]
    )
    assert any(
        item["symbol"] == "CD274"
        and "PD-L1" in item["display_aliases"]
        for item in alias_data["gene_aliases"]
    )
    assert canonicalize_marker_label("Tex") == "T cell exhaustion-like"
    assert canonicalize_marker_label("T_cell_activation") == "T cell activation"
    assert "dysfunctional T" in get_marker_aliases("T cell exhaustion-like")
    assert "PD-L1" in get_gene_display_aliases("CD274")
    assert "CD62L" in get_gene_display_aliases("SELL")


def test_marker_entries_expose_required_routing_metadata():
    """Manager-normalized marker entries should expose the routing contract."""
    required = {
        "kind",
        "granularity",
        "scope",
        "review_status",
        "use_for_global_annotation",
        "use_for_state_annotation",
        "use_for_malignancy_interpretation",
    }
    for resource_name in ["registry_human", "registry_mouse", "tissue_human", "tumor_human"]:
        mgr = Manager(resource_name, case_sensitive=True)
        for name, cell in mgr.CELLS.items():
            missing = required.difference(cell.metadata)
            assert not missing, f"{resource_name}:{name} missing {sorted(missing)}"


def test_manager_contract_audit_and_summary_are_machine_readable():
    """Manager should expose a compact curation audit for resource maintenance."""
    reference_data = tomllib.loads((RESOURCE_DIR / "references.toml").read_text())
    known_ids = {item["source_id"] for item in reference_data["references"]}
    mgr = Manager("registry_human", case_sensitive=True)

    assert mgr.validate_resource_contract(known_source_ids=known_ids) == []

    summary = mgr.audit_summary()
    assert summary["n_entries"] == len(mgr.CELLS)
    assert summary["by_kind"]["cell_type"] > 0
    assert summary["views"]["global_annotation"] > 0

    rows = mgr.get_marker_table()
    t_cell_row = next(row for row in rows if row["name"] == "T cells")
    assert t_cell_row["doublet_lineage"] is True
    assert t_cell_row["path"].endswith("T cells")


def test_metadata_selectors_support_include_and_exclude_routing():
    """Metadata selectors should support downstream-specific resource views."""
    mgr = Manager("registry_human", case_sensitive=True)
    immune_lineages = mgr.select_by_metadata(
        kind="cell_type",
        doublet_lineage=True,
        exclude={"review_status": ["conflict", "deprecated"]},
        include_children=False,
    )
    no_states = mgr.exclude_by_metadata(kind=["state", "functional_program"], include_children=False)

    assert "T cells" in immune_lineages.CELLS
    assert "B cells" in immune_lineages.CELLS
    assert "T cell exhaustion-like" not in immune_lineages.CELLS
    assert "T cell exhaustion-like" not in no_states.CELLS
    assert "T cells" in no_states.CELLS


def test_program_scoring_genesets_expose_useful_categories():
    """Cancer genesets should carry practical scoring categories into Manager metadata."""
    program = get_marker_manager(species="human", view="program_scoring")

    assert program.CELLS["HALLMARK_APOPTOSIS"].metadata["category"] == "cell_fate_stress"
    assert program.CELLS["HALLMARK_APOPTOSIS"].metadata["source_collection"] == "MSigDB Hallmark"
    assert "SRC0128" in program.CELLS["HALLMARK_APOPTOSIS"].metadata["source_ids"]
    assert "cell_identity_annotation" in program.CELLS["HALLMARK_APOPTOSIS"].metadata["not_for"]
    assert program.CELLS["Immune_infiltrated"].metadata["category"] == "TME"
    assert program.CELLS["Immune_infiltrated"].metadata["source_collection"] == "scLucid cancer signatures seed"
    assert program.CELLS["ITH_cell_cycle"].metadata["category"] == "ITH_Hallmarks"
    assert (
        program.CELLS["ITH_cell_cycle"].metadata["source_collection"]
        == "Gavish 2023 transcriptional ITH hallmarks"
    )
    assert program.CELLS["ITH_cell_cycle"].metadata["zotero_item_key"] == "RNXMRWIT"
    assert program.CELLS["Stromal_barrier"].metadata["kind"] == "geneset"
    assert "Immune_desert" not in program.CELLS
    assert program.CELLS["TCellSI_cytotoxicity"].metadata["category"] == "TCellSI"
    assert program.CELLS["TCellSI_cytotoxicity"].metadata["source_ids"] == ["SRC0135"]
    assert "cell_identity_annotation" in program.CELLS["TCellSI_cytotoxicity"].metadata["not_for"]
    assert program.CELLS["Tryptophan_degradation"].metadata["category"] == "Myeloid_Macrophage"
    assert program.CELLS["TREM2_macrophage_signature"].metadata["category"] == "Myeloid_Macrophage"
    assert program.CELLS["LRRC15_myofibroblast_program"].metadata["category"] == "CAF_Fibroblast"
    assert program.CELLS["PI16_fibroblast_program"].metadata["category"] == "CAF_Fibroblast"
    assert program.CELLS["ITH_glycolysis"].metadata["category"] == "ITH_Hallmarks"
    assert program.CELLS["ITH_cuproptosis_metabolism"].metadata["source_ids"] == ["SRC0016"]


def test_gavish_ith_tumor_states_route_to_tumor_interpretation():
    """Primary-tumor ITH states should support tumor review without global leakage."""
    tumor = get_marker_manager(species="human", view="tumor_interpretation")
    global_view = get_marker_manager(species="human", view="global_annotation")

    assert "Malignant ITH stemness state" in tumor.CELLS
    assert "Malignant ITH interferon/MHC-II state" in tumor.CELLS
    assert "Malignant ITH secretory/mucin state" in tumor.CELLS
    assert tumor.CELLS["Malignant ITH stemness state"].metadata["source_ids"] == ["SRC0016"]
    assert "Malignant ITH stemness state" not in global_view.CELLS


def test_state_requests_accept_marker_aliases():
    """Historical state/program names should resolve through marker_aliases."""
    mgr = get_marker_manager(species="human", states=["T_cell_activation"])

    assert "T cell activation" in mgr.CELLS
    assert "T_cell_activation" not in mgr.CELLS


def test_marker_source_ids_resolve_to_references():
    """Curated marker entries should only reference known bibliography IDs."""
    reference_data = tomllib.loads((RESOURCE_DIR / "references.toml").read_text())
    known_ids = {item["source_id"] for item in reference_data["references"]}

    for resource_name in ["registry_human", "registry_mouse", "tissue_human", "tumor_human"]:
        mgr = Manager(resource_name, case_sensitive=True)
        for name, cell in mgr.CELLS.items():
            for source_id in cell.metadata.get("source_ids", []):
                assert source_id in known_ids, f"{resource_name}:{name}:{source_id}"


def test_registry_views_cover_doublet_detection_and_plotting():
    """Doublet and plotting consumers should be routed by Manager views."""
    doublet = get_marker_manager(species="human", view="doublet_detection")
    plotting = get_marker_manager(
        species="human",
        cancer_type="Lung Cancer",
        view="plotting",
        include_functional=True,
    )

    assert {"T cells", "B cells", "NK cells", "Macrophages"}.issubset(doublet.CELLS)
    assert "Cycling" not in doublet.CELLS
    assert "Immune" in plotting.CELLS
    assert "Cytotoxicity" in plotting.CELLS
    assert "Diploid immune anchor" not in plotting.CELLS


def test_marker_curation_candidates_cover_all_batches_and_targets():
    """Each batch curation document should have at least one classified queue item."""
    candidates_path = DOCS_DIR / "marker_resources" / "marker_curation_candidates.jsonl"
    rows = [
        json.loads(line)
        for line in candidates_path.read_text().splitlines()
        if line.strip()
    ]
    required_fields = {
        "batch_id",
        "entry_name",
        "target_resource",
        "kind",
        "granularity",
        "markers",
        "negative_markers",
        "source_ids",
        "evidence_tier",
        "review_status",
        "notes",
    }
    valid_targets = {
        "marker_registry_human.toml",
        "marker_registry_mouse.toml",
        "marker_tissue_human.toml",
        "marker_tumor_human.toml",
        "marker_aliases.toml",
        "genesets_cancer_signatures.json",
        "genesets_cancer_hallmarks.json",
        "skip",
    }
    reference_data = tomllib.loads((RESOURCE_DIR / "references.toml").read_text())
    known_ids = {item["source_id"] for item in reference_data["references"]}

    assert {row["batch_id"] for row in rows} == {f"{i:02d}" for i in range(1, 13)}
    for row in rows:
        assert required_fields.issubset(row)
        for target in row["target_resource"].split(";"):
            assert target in valid_targets
        assert row["review_status"] in {"needs_review", "reviewed", "conflict", "skip"}
        for source_id in row["source_ids"]:
            assert source_id in known_ids
        if row["target_resource"] == "skip":
            assert not row["markers"]


def test_marker_curation_literature_index_covers_all_batch_sources():
    """The literature index should make all batch-md source papers auditable."""
    index_path = (
        DOCS_DIR / "marker_resources" / "marker_curation_literature_index.jsonl"
    )
    rows = [
        json.loads(line)
        for line in index_path.read_text().splitlines()
        if line.strip()
    ]
    assert {row["batch_id"] for row in rows} == {f"{i:02d}" for i in range(1, 13)}
    assert len(rows) == 141
    assert sum(bool(row["reference_source_ids"]) for row in rows) >= 137
    assert sum(not bool(row["reference_source_ids"]) for row in rows) <= 4
    assert sum(bool(row.get("zotero_item_key")) for row in rows) >= 133
    assert any(
        row["doi"] == "10.1038/s41586-023-06130-4"
        and row["reference_source_ids"] == ["SRC0016"]
        and row.get("zotero_item_key") == "RNXMRWIT"
        for row in rows
    )
    assert any(
        row["doi"] == "10.1038/s41577-025-01238-2"
        and row["reference_source_ids"] == ["SRC0132"]
        for row in rows
    )
