"""
Cell type marker manager for single-cell RNA-seq analysis.

This module provides functionality for managing and querying cell type markers,
enabling robust cell type annotation and visualization.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Literal

import tomllib
from anndata import AnnData
from rich.console import Console, Group
from rich.panel import Panel
from rich.style import Style
from rich.table import Table
from rich.text import Text
from rich.tree import Tree


def _get_marker_path(name_or_path: str) -> Path:
    """
    Parse built-in marker name or external file path.

    Args:
        name_or_path: If 'human' or 'mouse', etc., looks for built-in resources.
                      Otherwise, treated as a file path.
    Returns:
        A Path object pointing to a TOML file.
    """
    # Try treating it as a file path
    path = Path(name_or_path)
    if path.is_file():
        return path

    # If not a file, try as a built-in resource
    # .joinpath() is a safe way to build paths
    # .files() returns a Traversable object pointing to the package resources
    try:
        resource_path = resources.files("scRNA").joinpath(
            f"resources/manager_{name_or_path}.toml"
        )
        if resource_path.is_file():
            return resource_path
    except (ModuleNotFoundError, FileNotFoundError):
        # Fails if the scRNA package or resource file doesn't exist
        pass

    raise FileNotFoundError(
        f"Marker configuration '{name_or_path}' could not be found as a local file "
        "or as a built-in resource."
    )


@dataclass
class CellType:
    """
    Data class representing a cell type with its markers and metadata.

    Attributes:
        name: Name of the cell type
        color: Hex color code for visualization
        markers: List of marker genes
        level: Classification level ("major" or "minor")
        minor: List of minor cell types (for major cell types)
    """

    name: str
    color: str | None
    markers: list[str]
    level: Literal["major", "minor"]
    minor: list[str] | None = None

    def __post_init__(self) -> None:
        """Post-initialization processing"""
        if self.color == "":
            self.color = None


class Manager:
    """
    Manager for cell type markers and their hierarchical relationships.

    Provides functionality for loading, querying, and visualizing cell type marker genes.
    """

    def __init__(self, config: str | Path) -> None:
        """
        Initialize the marker manager from a TOML configuration file.

        Args:
            config_file: Path to the TOML configuration file containing cell type definitions
        """
        config_file_path = _get_marker_path(config)

        with open(config_file_path, "rb") as f:
            data = tomllib.load(f)

        self.CELLS = {}
        self.CLUSTERS = {}

        for major_cells in data.values():
            for major_cell in major_cells:
                minor_cell_name = None
                if minor_cells := major_cell.get("minor"):
                    minor_cell_name = [minor_cell["name"] for minor_cell in minor_cells]
                    self.CELLS.update(
                        {
                            minor_cell["name"]: CellType(
                                name=minor_cell["name"],
                                color=minor_cell["color"],
                                markers=minor_cell["markers"],
                                level="minor",
                            )
                            for minor_cell in minor_cells
                        }
                    )
                    self.CLUSTERS[major_cell["name"]] = minor_cell_name
                else:
                    self.CLUSTERS[major_cell["name"]] = []
                self.CELLS[major_cell["name"]] = CellType(
                    name=major_cell["name"],
                    color=major_cell["color"],
                    markers=major_cell["markers"],
                    level="major",
                    minor=minor_cell_name,
                )

    def __getitem__(self, key: str) -> CellType:
        """
        Get a cell type by name.

        Args:
            key: Name of the cell type

        Returns:
            CellType object

        Raises:
            KeyError: If the cell type is not found
        """
        if key not in self.CELLS:
            raise KeyError(f"{key} is not in CELLS")
        else:
            return self.CELLS[key]

    def intersect_with(self, adata: AnnData) -> None:
        """
        Intersect marker genes with genes present in the AnnData object.

        This updates the marker lists to only include genes that are present in the dataset.

        Args:
            adata: AnnData object containing gene expression data
        """
        genes = adata.var_names
        for cell_name, cell_type in self.CELLS.items():
            self.CELLS[cell_name].markers = list(
                set(cell_type.markers).intersection(genes)
            )

    def query(
        self,
        info: Literal["color", "markers"] = "color",
        key: str | Sequence[str] | None = None,
        cluster: str | None = None,
        include_major_type: bool = True,
    ) -> dict[str, str | list[str]]:
        """
        Query information about cell types.

        Args:
            info: Type of information to query ("color" or "markers")
            key: Cell type name(s) to query
            cluster: Major cell type to query (returns info for all subtypes)
            include_major_type: Whether to include the major type when querying a cluster

        Returns:
            Dictionary mapping cell types to requested information

        Raises:
            KeyError: If the specified cell type or cluster is not found
            ValueError: If neither key nor cluster is specified
            TypeError: If key is not a string or sequence of strings
        """
        if cluster is not None:
            if cluster not in self.CLUSTERS.keys():
                raise KeyError(f"{cluster} is not in CLUSTERS")
            clusters_ = (
                [cluster] + self.CLUSTERS[cluster]
                if include_major_type
                else self.CLUSTERS[cluster]
            )
            return {
                cell_type: getattr(self.CELLS[cell_type], info)
                for cell_type in clusters_
            }
        elif key is not None:
            if isinstance(key, str):
                if key not in self.CELLS:
                    raise KeyError(f"{key} is not in CELLS")
                else:
                    return {key: getattr(self.CELLS[key], info)}
            elif isinstance(key, Sequence):
                if missing_key := [i for i in key if i not in self.CELLS]:
                    raise KeyError(f"{missing_key} is not in CELLS")
                return {
                    cell_type: getattr(self.CELLS[cell_type], info) for cell_type in key
                }
            else:
                raise TypeError(f"index must be str or Sequence[str], not {type(key)}")
        else:
            raise ValueError("either cluster or key must be specified")

    def show_clusters(self) -> None:
        """
        Display a table of cell type clusters with their colors.

        Shows major cell types and their corresponding minor (sub) cell types.
        """
        table = Table(title="Cluster")
        table.add_column("Major type", justify="center")
        table.add_column("Minor type", justify="center", max_width=80)
        for major_name, minors_names in self.CLUSTERS.items():
            major_cell = self.CELLS[major_name]
            minor_cells = [self.CELLS[minor_name] for minor_name in minors_names]
            minor_text = Text()
            for minor_cell in minor_cells:
                minor_text.append(
                    f"{minor_cell.name}, ",
                    style=Style(color=minor_cell.color, bold=True),
                )
            table.add_row(
                Text(
                    major_cell.name,
                    style=Style(bgcolor=major_cell.color, bold=True, color="black"),
                ),
                minor_text,
            )
            table.add_section()

        console = Console()
        console.print(table)

    def show_markers(self) -> None:
        """
        Display a table of cell types and their marker genes.

        Shows marker genes for both major cell types and their subtypes.
        """
        table = Table(title="Marker")
        table.add_column("Cell type", justify="left")
        table.add_column("Markers", justify="center", max_width=80)
        for major_name, minor_names in self.CLUSTERS.items():
            major_cell = self.CELLS[major_name]
            minor_cells = [self.CELLS[minor_name] for minor_name in minor_names]
            table.add_row(
                Text(
                    major_cell.name,
                    style=Style(bold=True, bgcolor=major_cell.color, color="black"),
                    justify="center",
                ),
                ", ".join(major_cell.markers),
            )
            table.add_section()
            for minor_cell in minor_cells:
                table.add_row(
                    Text(
                        minor_cell.name,
                        style=Style(bold=True, bgcolor=minor_cell.color, color="black"),
                        justify="center",
                    ),
                    ", ".join(minor_cell.markers),
                )
                table.add_section()
        console = Console()
        console.print(table)

    def show_tree(self) -> None:
        """
        Display a hierarchical tree of cell types and their marker genes.

        Visualizes the hierarchical structure of cell types and their markers.
        """
        tree = Tree("Cluster")
        for major_name, minor_names in self.CLUSTERS.items():
            major_cell = self.CELLS[major_name]
            major = tree.add(
                Group(
                    Text(
                        major_name,
                        style=Style(bold=True, bgcolor=major_cell.color, color="black"),
                    ),
                    Panel(
                        Text(",".join(major_cell.markers), justify="center"), width=40
                    ),
                )
            )

            for minor_name in minor_names:
                minor_cell = self.CELLS[minor_name]
                major.add(
                    Group(
                        Text(
                            minor_name,
                            style=Style(
                                bold=True, bgcolor=minor_cell.color, color="black"
                            ),
                        ),
                        Panel(
                            Text(",".join(minor_cell.markers), justify="center"),
                            width=40,
                        ),
                    )
                )

        console = Console()
        console.print(tree)

    def add_cell_type(
        self,
        name: str,
        markers: list[str],
        color: str = None,
        level: Literal["major", "minor"] = "major",
        parent: str = None,
    ) -> None:
        """
        Add a new cell type to the manager.

        Args:
            name: Name of the cell type
            markers: List of marker genes
            color: Hex color code for visualization
            level: Classification level ("major" or "minor")
            parent: Parent cell type for minor cell types

        Raises:
            ValueError: If adding a minor cell type without a parent,
                       or if the parent does not exist
        """
        if level == "minor" and parent is None:
            raise ValueError("Minor cell types must have a parent specified")

        if level == "minor" and parent not in self.CELLS:
            raise ValueError(f"Parent cell type '{parent}' does not exist")

        # Create the cell type
        cell_type = CellType(
            name=name,
            color=color,
            markers=markers,
            level=level,
            minor=[] if level == "major" else None,
        )

        # Add to CELLS dictionary
        self.CELLS[name] = cell_type

        # Update CLUSTERS dictionary for major cell types
        if level == "major":
            self.CLUSTERS[name] = []

        # Update parent's minor list and CLUSTERS for minor cell types
        if level == "minor":
            if self.CELLS[parent].minor is None:
                self.CELLS[parent].minor = []
            self.CELLS[parent].minor.append(name)
            self.CLUSTERS[parent].append(name)
