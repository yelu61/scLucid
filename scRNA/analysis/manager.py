"""
Cell type marker manager for single-cell RNA-seq analysis.

This module provides a robust system for managing and querying hierarchical
cell type markers, enabling flexible and reproducible cell type annotation.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional

import tomllib
from anndata import AnnData
from rich.console import Console
from rich.panel import Panel
from rich.style import Style
from rich.text import Text
from rich.tree import Tree


def _get_marker_path(name_or_path: str) -> Path:
    """
    Finds a marker file, searching for built-in resources first, then local paths.
    """
    full_resource_name = f"marker_{name_or_path}.toml"
    try:
        resource_path = resources.files("scRNA").joinpath(
            f"resources/{full_resource_name}"
        )
        if resource_path.is_file():
            return resource_path
    except (ModuleNotFoundError, FileNotFoundError):
        pass

    path = Path(name_or_path)
    if path.is_file():
        return path

    raise FileNotFoundError(
        f"Marker configuration '{name_or_path}' could not be found as a local file "
        f"or as a built-in resource (expected filename: '{full_resource_name}')."
    )


@dataclass
class CellType:
    """Data class representing a cell type with its markers and metadata."""

    name: str
    color: Optional[str]
    markers: List[str]
    level: Literal["major", "minor"]
    parent: Optional[CellType] = None
    minor: List[CellType] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Post-initialization processing."""
        if self.color == "":
            self.color = None

    # ADD THIS METHOD INSIDE THE CLASS
    def to_dict(self) -> dict:
        """Converts this CellType object back to a dictionary."""
        d = {
            "name": self.name,
            "color": self.color,
            "markers": self.markers,
        }
        # Recursively convert children to dictionaries as well
        if self.minor:
            d["minor"] = [m.to_dict() for m in self.minor]
        return d


class Manager:
    """Manages hierarchical cell type markers for analysis."""

    def __init__(self, config: str, root_key: Optional[str] = None) -> None:
        """Initializes the marker manager from a configuration."""
        config_file_path = _get_marker_path(config)
        with open(config_file_path, "rb") as f:
            data = tomllib.load(f)

        data_to_parse = data
        if root_key:
            if root_key not in data:
                raise KeyError(f"root_key '{root_key}' not found in {config_file_path}")
            data_to_parse = {root_key: data[root_key]}

        self.CELLS: Dict[str, CellType] = {}
        self.CLUSTERS: Dict[str, List[CellType]] = {}
        self._parse_level(data_to_parse)

    def _parse_level(self, level_data: dict, parent_obj: Optional[CellType] = None):
        """Recursively parses a level of the marker hierarchy."""
        for major_name, definitions in level_data.items():
            if parent_obj is None:
                self.CLUSTERS[major_name] = []

            for cell_def in definitions:
                cell_name = cell_def["name"]

                # Create or update the cell type object
                if cell_name in self.CELLS:
                    cell_obj = self.CELLS[cell_name]
                    cell_obj.color = (
                        cell_def.get("color")
                        if cell_def.get("color")
                        else cell_obj.color
                    )
                    cell_obj.markers.extend(
                        m
                        for m in cell_def.get("markers", [])
                        if m not in cell_obj.markers
                    )
                else:
                    cell_obj = CellType(
                        name=cell_name,
                        color=cell_def.get("color"),
                        markers=cell_def.get("markers", []),
                        level="minor" if parent_obj else "major",
                        parent=parent_obj,
                    )
                    self.CELLS[cell_name] = cell_obj

                # Establish parent-child relationship
                if parent_obj:
                    parent_obj.minor.append(cell_obj)
                else:
                    self.CLUSTERS[major_name].append(cell_obj)

                # Recurse if there are further subtypes
                if "minor" in cell_def:
                    self._parse_level(
                        {cell_name: cell_def["minor"]}, parent_obj=cell_obj
                    )

    def __getitem__(self, key: str) -> CellType:
        """Allows dictionary-style access to cell types."""
        if key not in self.CELLS:
            raise KeyError(f"Cell type '{key}' not found in manager.")
        return self.CELLS[key]

    def intersect_with(self, adata: AnnData) -> None:
        """Filters all marker lists to include only genes present in the AnnData object."""
        genes_in_data = set(adata.var_names)
        for cell in self.CELLS.values():
            cell.markers = [m for m in cell.markers if m in genes_in_data]

    def query(
        self, info: Literal["color", "markers"], key: str | Sequence[str]
    ) -> dict:
        """Queries for information for one or more cell types."""
        keys_to_query = [key] if isinstance(key, str) else key
        missing = [k for k in keys_to_query if k not in self.CELLS]
        if missing:
            raise KeyError(f"Cell types not found: {', '.join(missing)}")

        return {k: getattr(self.CELLS[k], info) for k in keys_to_query}

    def show_tree(self) -> None:
        """Displays the entire cell type hierarchy as a rich tree."""
        console = Console()
        tree = Tree("🔖 [bold cyan]Cell Type Hierarchy[/bold cyan]", guide_style="cyan")

        def add_to_tree(parent_node: Tree, cell_obj: CellType):
            """Helper function to recursively build the rich tree."""
            panel_content = Text(", ".join(cell_obj.markers), justify="left")
            node_panel = Panel(
                panel_content,
                title=f"[bold]{cell_obj.name}[/bold]",
                title_align="left",
                border_style=Style(color=cell_obj.color or "white"),
                width=60,
            )
            child_node = parent_node.add(node_panel)
            for minor_obj in cell_obj.minor:
                add_to_tree(child_node, minor_obj)

        for major_name, cell_objects in self.CLUSTERS.items():
            major_node = tree.add(f"[bold magenta]{major_name}[/bold magenta]")
            for cell_obj in cell_objects:
                add_to_tree(major_node, cell_obj)

        console.print(tree)

    def merge_from(self, other_manager: "Manager") -> None:
        """Merges all definitions from another Manager instance into this one."""
        for major_name, cell_list in other_manager.CLUSTERS.items():
            self._parse_level({major_name: [cell.to_dict() for cell in cell_list]})


def get_marker_manager(
    species: str,
    tissue: Optional[str] = None,
    states: Optional[List[str]] = None,
) -> Manager:
    """
    Factory function to build a Manager by combining base, tissue, and state markers.

    Args:
        species: The species ('human' or 'mouse').
        tissue: The tissue name, corresponding to a top-level key in the tissue-specific file.
        states: A list of cell states, corresponding to keys in the cell-state file.

    Returns:
        A fully configured and combined Manager instance.
    """
    print(
        f"Building manager for species: {species}, tissue: {tissue}, states: {states}"
    )

    mgr = Manager(f"base_{species}")

    if tissue:
        mgr_tissue = Manager(f"tissue_specific_{species}", root_key=tissue)
        mgr.merge_from(mgr_tissue)
        print(f"-> Merged tissue-specific markers for '{tissue}'.")

    if states:
        # Load the entire cell state file once
        mgr_states_all = Manager(f"cell_state_{species}")

        # Now, selectively merge only the requested states
        for state_name in states:
            if state_name in mgr_states_all.CELLS:
                # Find the parent category ("Cell States")
                parent_category = "Cell States"

                # Re-create a temporary manager with just this one state to merge
                state_obj = mgr_states_all.CELLS[state_name]
                temp_mgr = Manager(f"cell_state_{species}", root_key=parent_category)
                temp_mgr.CELLS = {state_name: state_obj}
                temp_mgr.CLUSTERS = {parent_category: [state_obj]}

                mgr.merge_from(temp_mgr)
                print(f"-> Merged cell state markers for '{state_name}'.")
            else:
                print(
                    f"Warning: State '{state_name}' not found in {f'cell_state_{species}.toml'}"
                )

    print("✅ Manager built successfully.")
    return mgr
