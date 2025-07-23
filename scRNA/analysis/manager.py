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
from rich.table import Table
from rich.text import Text
from rich.tree import Tree


def _get_marker_path(name_or_path: str) -> Path:
    """
    Finds a marker file, searching for built-in resources first, then local paths.

    Args:
        name_or_path: Name of a built-in set (e.g., 'base_human') or a direct file path.

    Returns:
        A Path object to a valid TOML file.

    Raises:
        FileNotFoundError: If the configuration cannot be located.
    """
    # First, try to resolve it as a built-in resource.
    full_resource_name = f"marker_{name_or_path}.toml" # <-- This is the key change
    try:
        # Note: The file name in the resources directory should not have "marker_" prefix.
        # It should be exactly 'base_human.toml', 'tissue_specific_human.toml', etc.
        resource_path = resources.files("scRNA").joinpath(f"resources/{full_resource_name}")
        if resource_path.is_file():
            return resource_path
    except (ModuleNotFoundError, FileNotFoundError):
        # This is not an error, just means it's not a built-in resource.
        pass

    # If not found as a resource, try as a direct file path.
    path = Path(name_or_path)
    if path.is_file():
        return path

    raise FileNotFoundError(
        f"Marker configuration '{name_or_path}' could not be found as a local file "
        f"or as a built-in resource (expected filename: '{full_resource_name}')."
    )


@dataclass
class CellType:
    """
    Data class representing a cell type with its markers and metadata.

    Attributes:
        name: The unique name of the cell type.
        color: A hex color code for visualization.
        markers: A list of marker genes.
        level: The hierarchical level ('major' or 'minor').
        minor: A list of names of direct subtypes.
    """
    name: str
    color: Optional[str]
    markers: List[str]
    level: Literal["major", "minor"]
    minor: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Post-initialization processing."""
        if self.color == "":
            self.color = None


class Manager:
    """
    Manages hierarchical cell type markers for analysis.

    This class loads marker definitions from TOML files, supports unlimited
    nesting of cell types, and provides methods for querying, visualization,
    and dynamic modification.
    """

    def __init__(self, config: str, root_key: Optional[str] = None) -> None:
        """
        Initializes the marker manager from a configuration.

        Args:
            config: Path to a TOML file or the name of a built-in marker set.
            root_key: If specified, only loads definitions from this top-level
                      category within the TOML file.
        """
        config_file_path = _get_marker_path(config)
        with open(config_file_path, "rb") as f:
            data = tomllib.load(f)

        data_to_parse = data
        if root_key:
            if root_key not in data:
                raise KeyError(f"root_key '{root_key}' not found in {config_file_path}")
            data_to_parse = {root_key: data[root_key]}

        self.CELLS: Dict[str, CellType] = {}
        self.CLUSTERS: Dict[str, List[str]] = {}
        self._parse_level(data_to_parse)

    def _parse_level(self, level_data: dict, parent_obj: Optional[CellType] = None):
        """
        Recursively parses a level of the marker hierarchy.

        Args:
            level_data: A dictionary representing the current level from the TOML data.
            parent_obj: The parent CellType object for the current level, if any.
        """
        for major_name, definitions in level_data.items():
            if parent_obj is None:  # This is a top-level category
                self.CLUSTERS[major_name] = []

            for cell_def in definitions:
                cell_name = cell_def["name"]
                if cell_name in self.CELLS:
                    # If cell type already exists, update it instead of overwriting
                    # This can happen when merging files with overlapping but refined definitions
                    cell_obj = self.CELLS[cell_name]
                    cell_obj.color = cell_def.get("color", cell_obj.color)
                    cell_obj.markers.extend(m for m in cell_def.get("markers", []) if m not in cell_obj.markers)
                else:
                    cell_obj = CellType(
                        name=cell_name,
                        color=cell_def.get("color"),
                        markers=cell_def.get("markers", []),
                        level="minor" if parent_obj else "major",
                    )
                    self.CELLS[cell_name] = cell_obj

                if parent_obj:
                    parent_obj.minor.append(cell_name)
                else:
                    self.CLUSTERS[major_name].append(cell_name)

                if "minor" in cell_def:
                    self._parse_level({cell_name: cell_def["minor"]}, parent_obj=cell_obj)

    def __getitem__(self, key: str) -> CellType:
        """Allows dictionary-style access to cell types."""
        if key not in self.CELLS:
            raise KeyError(f"Cell type '{key}' not found in manager.")
        return self.CELLS[key]

    def intersect_with(self, adata: AnnData) -> None:
        """
        Filters all marker lists to include only genes present in the AnnData object.
        """
        genes_in_data = set(adata.var_names)
        for cell in self.CELLS.values():
            cell.markers = [m for m in cell.markers if m in genes_in_data]

    def query(self, info: Literal["color", "markers"], key: str | Sequence[str]) -> dict:
        """
        Queries for information (colors or markers) for one or more cell types.
        """
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
            for minor_name in cell_obj.minor:
                add_to_tree(child_node, self.CELLS[minor_name])

        for major_name, cell_names in self.CLUSTERS.items():
            major_node = tree.add(f"[bold magenta]{major_name}[/bold magenta]")
            for cell_name in cell_names:
                add_to_tree(major_node, self.CELLS[cell_name])
        
        console.print(tree)

    def add_cell_type(self, **kwargs) -> None:
        """A simple wrapper to add a single cell type definition."""
        self._parse_level({"adhoc": [kwargs]})

    def merge_from(self, other_manager: Manager) -> None:
        """
        Merges all definitions from another Manager instance into this one.

        Args:
            other_manager: Another Manager instance to merge from.
        """
        for cell_name, cell_obj in other_manager.CELLS.items():
            parent_name = None
            for major, minors in other_manager.CLUSTERS.items():
                if cell_name in minors:
                    parent_name = major
                    break
            self.add_cell_type(
                name=cell_obj.name,
                markers=cell_obj.markers,
                color=cell_obj.color,
                level=cell_obj.level,
                parent=parent_name
            )


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
    print(f"Building manager for species: {species}, tissue: {tissue}, states: {states}")

    # 1. Load the base manager for the specified species
    mgr = Manager(f"base_{species}")

    # 2. Optionally, load and merge tissue-specific markers
    if tissue:
        mgr_tissue = Manager(f"tissue_specific_{species}", root_key=tissue)
        mgr.merge_from(mgr_tissue)
        print(f"-> Merged tissue-specific markers for '{tissue}'.")

    # 3. Optionally, load and merge cell state markers
    if states:
        for state in states:
            mgr_state = Manager(f"cell_state_{species}", root_key=state)
            mgr.merge_from(mgr_state)
            print(f"-> Merged cell state markers for '{state}'.")

    print("✅ Manager built successfully.")
    return mgr