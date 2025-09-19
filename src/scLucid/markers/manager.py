"""
Cell type marker manager for single-cell RNA-seq analysis.

This module provides a robust system for managing and querying hierarchical
cell type markers, enabling flexible and reproducible cell type annotation.
The marker system supports species-specific, tissue-specific, and cell state-specific
markers, allowing for customizable annotation strategies.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union

import tomllib
from anndata import AnnData
from rich.console import Console
from rich.panel import Panel
from rich.style import Style
from rich.text import Text
from rich.tree import Tree

# Configure logging
log = logging.getLogger(__name__)

# Define allowed marker file formats
MARKER_FORMATS = ["toml", "json"]

# Define known species for validation
KNOWN_SPECIES = [
    "human",
    "mouse",
]


def _get_marker_path(name_or_path: str) -> Path:
    """
    Finds a marker file, searching for built-in resources first, then local paths.

    This function first checks if the file exists as a built-in resource within
    the package, then checks if it's a valid path to a local file.

    Args:
        name_or_path: Either a marker name (e.g., "base_human") or a path to a marker file

    Returns:
        Path object pointing to the marker file

    Raises:
        FileNotFoundError: If the marker file cannot be found
    """
    # Try first with .toml extension
    full_resource_name = f"marker_{name_or_path}.toml"
    log.debug(f"Looking for built-in resource: {full_resource_name}")

    try:
        resource_path = resources.files("scLucid").joinpath(
            f"resources/{full_resource_name}"
        )
        if resource_path.is_file():
            log.debug(f"Found built-in resource: {resource_path}")
            return resource_path
    except (ModuleNotFoundError, FileNotFoundError):
        log.debug(f"Built-in resource not found: {full_resource_name}")

    # Next, try with .json extension
    full_resource_name = f"marker_{name_or_path}.json"
    log.debug(f"Looking for built-in resource: {full_resource_name}")

    try:
        resource_path = resources.files("scLucid").joinpath(
            f"resources/{full_resource_name}"
        )
        if resource_path.is_file():
            log.debug(f"Found built-in resource: {resource_path}")
            return resource_path
    except (ModuleNotFoundError, FileNotFoundError):
        log.debug(f"Built-in resource not found: {full_resource_name}")

    # Try as a direct path
    path = Path(name_or_path)
    if path.is_file():
        log.debug(f"Found local file: {path}")
        return path

    # Try adding extensions if no extension is present
    if "." not in name_or_path:
        for ext in MARKER_FORMATS:
            test_path = Path(f"{name_or_path}.{ext}")
            if test_path.is_file():
                log.debug(f"Found local file with added extension: {test_path}")
                return test_path

    # If we get here, we couldn't find the file
    error_msg = (
        f"Marker configuration '{name_or_path}' could not be found as a local file "
        f"or as a built-in resource (expected filename: 'marker_{name_or_path}.toml' or "
        f"'marker_{name_or_path}.json')."
    )
    log.error(error_msg)
    raise FileNotFoundError(error_msg)


def _load_marker_file(file_path: Path) -> dict:
    """
    Loads a marker file in either TOML or JSON format.

    Args:
        file_path: Path to the marker file

    Returns:
        Dictionary containing the marker data

    Raises:
        ValueError: If the file format is not supported
    """
    suffix = file_path.suffix.lower()

    if suffix == ".toml":
        log.debug(f"Loading TOML file: {file_path}")
        with open(file_path, "rb") as f:
            return tomllib.load(f)
    elif suffix == ".json":
        log.debug(f"Loading JSON file: {file_path}")
        with open(file_path, "r", encoding="utf-8") as f:
            return json.load(f)
    else:
        error_msg = f"Unsupported marker file format: {suffix}. Supported formats: {MARKER_FORMATS}"
        log.error(error_msg)
        raise ValueError(error_msg)


@dataclass
class CellType:
    """
    Data class representing a cell type with its markers and metadata.

    Attributes:
        name: The name of the cell type
        color: Optional color for visualization (hex code or named color)
        markers: List of marker gene names associated with this cell type
        level: Whether this is a major or minor cell type in the hierarchy
        parent: Reference to the parent cell type (if any)
        minor: List of child cell types (subtypes)
        metadata: Additional metadata about this cell type
    """

    name: str
    color: Optional[str]
    markers: List[str]
    level: Literal["major", "minor"]
    parent: Optional[CellType] = None
    minor: List[CellType] = field(default_factory=list)
    metadata: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Post-initialization processing."""
        # Handle empty color strings
        if self.color == "":
            self.color = None

        # Ensure markers list is unique
        if self.markers:
            self.markers = list(dict.fromkeys(self.markers))

    def to_dict(self) -> dict:
        """
        Converts this CellType object back to a dictionary.

        Returns:
            Dictionary representation of the cell type
        """
        d = {
            "name": self.name,
            "color": self.color,
            "markers": self.markers,
        }

        # Add metadata if present
        if self.metadata:
            d["metadata"] = self.metadata

        # Recursively convert children to dictionaries as well
        if self.minor:
            d["minor"] = [m.to_dict() for m in self.minor]

        return d

    def add_markers(self, new_markers: List[str]) -> None:
        """
        Adds new markers to this cell type, ensuring uniqueness.

        Args:
            new_markers: List of marker gene names to add
        """
        if not new_markers:
            return

        # Add only markers that aren't already in the list
        existing = set(self.markers)
        self.markers.extend([m for m in new_markers if m not in existing])

    def get_all_markers(self, include_subtypes: bool = False) -> List[str]:
        """
        Gets all markers for this cell type, optionally including subtypes.

        Args:
            include_subtypes: Whether to include markers from all subtypes

        Returns:
            List of unique marker genes
        """
        markers = set(self.markers)

        if include_subtypes and self.minor:
            for subtype in self.minor:
                markers.update(subtype.get_all_markers(include_subtypes=True))

        return list(markers)

    def get_hierarchy_path(self) -> str:
        """
        Gets the full hierarchical path of this cell type.

        Returns:
            String representation of the path (e.g., "Immune cells/T cells/CD8+ T cells")
        """
        parts = []
        current = self

        # Walk up the hierarchy
        while current:
            parts.append(current.name)
            current = current.parent

        # Reverse to get top-down path and join with slashes
        return "/".join(reversed(parts))


class Manager:
    """
    Manages hierarchical cell type markers for single-cell RNA-seq analysis.

    This class provides functionality to load, manage, and query cell type marker
    genes from configuration files. It supports hierarchical organization of cell types
    and enables flexible annotation strategies.
    """

    def __init__(
        self, config: str, root_key: Optional[str] = None, case_sensitive: bool = False
    ) -> None:
        """
        Initializes the marker manager from a configuration file.

        Args:
            config: Path to marker configuration file or name of built-in marker set
            root_key: If provided, only load this top-level key from the configuration
            case_sensitive: Whether gene names should be case-sensitive

        Raises:
            FileNotFoundError: If the marker file cannot be found
            KeyError: If root_key is provided but not found in the configuration
        """
        log.info(f"Initializing marker manager from '{config}'")

        self.CELLS: Dict[str, CellType] = {}
        self.CLUSTERS: Dict[str, List[CellType]] = {}
        self.case_sensitive = case_sensitive
        self._source_file = config

        try:
            config_file_path = _get_marker_path(config)
            data = _load_marker_file(config_file_path)

            data_to_parse = data
            if root_key:
                if root_key not in data:
                    error_msg = f"root_key '{root_key}' not found in {config_file_path}"
                    log.error(error_msg)
                    raise KeyError(error_msg)
                log.info(f"Loading only '{root_key}' from configuration")
                data_to_parse = {root_key: data[root_key]}

            self._parse_level(data_to_parse)

            log.info(
                f"Loaded {len(self.CELLS)} cell types in {len(self.CLUSTERS)} major categories"
            )

        except Exception as e:
            log.error(f"Error initializing marker manager: {str(e)}")
            raise

    def _parse_level(
        self, level_data: dict, parent_obj: Optional[CellType] = None
    ) -> None:
        """
        Recursively parses a level of the marker hierarchy.

        Args:
            level_data: Dictionary containing cell type definitions
            parent_obj: Parent CellType object (if any)
        """
        for major_name, definitions in level_data.items():
            log.debug(f"Parsing level: {major_name}")

            if parent_obj is None:
                self.CLUSTERS[major_name] = []

            for cell_def in definitions:
                # Validate required fields
                if "name" not in cell_def:
                    log.warning(
                        f"Skipping cell type definition without 'name' field in {major_name}"
                    )
                    continue

                cell_name = cell_def["name"]
                log.debug(f"Processing cell type: {cell_name}")

                # Create or update the cell type object
                if cell_name in self.CELLS:
                    # Update existing cell type
                    cell_obj = self.CELLS[cell_name]

                    # Update color if provided
                    if "color" in cell_def and cell_def["color"]:
                        cell_obj.color = cell_def["color"]

                    # Add new markers
                    new_markers = cell_def.get("markers", [])
                    if not self.case_sensitive:
                        new_markers = [m.upper() for m in new_markers]

                    cell_obj.add_markers(new_markers)

                    # Update metadata if provided
                    if "metadata" in cell_def:
                        cell_obj.metadata.update(cell_def["metadata"])

                    log.debug(f"Updated existing cell type: {cell_name}")

                else:
                    # Create new cell type
                    markers = cell_def.get("markers", [])
                    if not self.case_sensitive:
                        markers = [m.upper() for m in markers]

                    cell_obj = CellType(
                        name=cell_name,
                        color=cell_def.get("color"),
                        markers=markers,
                        level="minor" if parent_obj else "major",
                        parent=parent_obj,
                        metadata=cell_def.get("metadata", {}),
                    )
                    self.CELLS[cell_name] = cell_obj
                    log.debug(
                        f"Created new cell type: {cell_name} with {len(markers)} markers"
                    )

                # Establish parent-child relationship
                if parent_obj:
                    # Check if this cell is already a child (avoid duplicates)
                    if cell_obj not in parent_obj.minor:
                        parent_obj.minor.append(cell_obj)
                else:
                    # Add to top-level cluster if not already there
                    if cell_obj not in self.CLUSTERS[major_name]:
                        self.CLUSTERS[major_name].append(cell_obj)

                # Recurse if there are further subtypes
                if "minor" in cell_def:
                    self._parse_level(
                        {cell_name: cell_def["minor"]}, parent_obj=cell_obj
                    )

    def __getitem__(self, key: str) -> CellType:
        """
        Allows dictionary-style access to cell types.

        Args:
            key: Name of the cell type to retrieve

        Returns:
            CellType object

        Raises:
            KeyError: If the cell type is not found
        """
        if key not in self.CELLS:
            raise KeyError(f"Cell type '{key}' not found in manager.")
        return self.CELLS[key]

    def __contains__(self, key: str) -> bool:
        """
        Checks if a cell type is in the manager.

        Args:
            key: Name of the cell type to check

        Returns:
            True if the cell type exists, False otherwise
        """
        return key in self.CELLS

    def __len__(self) -> int:
        """
        Gets the number of cell types in the manager.

        Returns:
            Number of cell types
        """
        return len(self.CELLS)

    def intersect_with(self, adata: AnnData) -> Tuple[int, int]:
        """
        Filters all marker lists to include only genes present in the AnnData object.

        Args:
            adata: AnnData object containing gene expression data

        Returns:
            Tuple of (total markers before filtering, total markers after filtering)
        """
        log.info(
            f"Intersecting markers with AnnData object containing {adata.n_vars} genes"
        )

        genes_in_data = {g.upper() for g in adata.var_names}

        total_before = 0
        total_after = 0

        for cell_type, cell in self.CELLS.items():
            n_before = len(cell.markers)
            total_before += n_before

            cell.markers = [m for m in cell.markers if m.upper() in genes_in_data]
            
            n_after = len(cell.markers)
            total_after += n_after

            if n_before > 0 and n_after == 0:
                log.warning(
                    f"Cell type '{cell_type}' has no markers left after intersection"
                )
            elif n_before > n_after:
                log.debug(
                    f"Cell type '{cell_type}': {n_after}/{n_before} markers remain after intersection"
                )

        retention = total_after / total_before * 100 if total_before > 0 else 0
        log.info(
            f"Marker intersection: {total_after}/{total_before} markers retained ({retention:.1f}%)"
        )

        return total_before, total_after

    def query(
        self,
        info: Literal["color", "markers", "metadata"],
        key: Union[str, Sequence[str]],
    ) -> dict:
        """
        Queries for information for one or more cell types.

        Args:
            info: Type of information to retrieve ('color', 'markers', or 'metadata')
            key: Cell type name(s) to query

        Returns:
            Dictionary mapping cell type names to requested information

        Raises:
            KeyError: If any requested cell type is not found
            ValueError: If info type is not recognized
        """
        if info not in ["color", "markers", "metadata"]:
            raise ValueError(
                f"Invalid info type: {info}. Must be 'color', 'markers', or 'metadata'"
            )

        keys_to_query = [key] if isinstance(key, str) else key

        missing = [k for k in keys_to_query if k not in self.CELLS]
        if missing:
            raise KeyError(f"Cell types not found: {', '.join(missing)}")

        return {k: getattr(self.CELLS[k], info) for k in keys_to_query}

    def show_tree(self, max_markers: int = 5) -> None:
        """
        Displays the entire cell type hierarchy as a rich tree.

        Args:
            max_markers: Maximum number of markers to display per cell type
        """
        console = Console()
        tree = Tree("🔖 [bold cyan]Cell Type Hierarchy[/bold cyan]", guide_style="cyan")

        def add_to_tree(parent_node: Tree, cell_obj: CellType):
            """Helper function to recursively build the rich tree."""
            # Format markers (limit to max_markers)
            if cell_obj.markers:
                marker_text = ", ".join(cell_obj.markers[:max_markers])
                if len(cell_obj.markers) > max_markers:
                    marker_text += f"... (+{len(cell_obj.markers) - max_markers} more)"
            else:
                marker_text = "[italic]No markers[/italic]"

            # Create panel content with markers
            panel_content = Text(marker_text, justify="left")

            # Add metadata if available
            if cell_obj.metadata:
                meta_text = ", ".join(f"{k}: {v}" for k, v in cell_obj.metadata.items())
                panel_content.append("\n[dim]" + meta_text + "[/dim]")

            # Create panel with cell type information
            node_panel = Panel(
                panel_content,
                title=f"[bold]{cell_obj.name}[/bold] ({len(cell_obj.markers)} markers)",
                title_align="left",
                border_style=Style(color=cell_obj.color or "white"),
                width=70,
            )

            # Add to tree and process children
            child_node = parent_node.add(node_panel)
            for minor_obj in cell_obj.minor:
                add_to_tree(child_node, minor_obj)

        # Process each major category
        for major_name, cell_objects in self.CLUSTERS.items():
            major_node = tree.add(f"[bold magenta]{major_name}[/bold magenta]")
            for cell_obj in cell_objects:
                add_to_tree(major_node, cell_obj)

        console.print(tree)

        # Print summary
        total_markers = sum(len(cell.markers) for cell in self.CELLS.values())
        console.print(
            f"[bold]Summary:[/bold] {len(self.CELLS)} cell types with {total_markers} markers "
            f"across {len(self.CLUSTERS)} major categories."
        )

    def merge_from(self, other_manager: "Manager") -> None:
        """
        Merges all definitions from another Manager instance into this one.

        Args:
            other_manager: Another Manager instance to merge from
        """
        log.info(
            f"Merging markers from another manager with {len(other_manager.CELLS)} cell types"
        )

        for major_name, cell_list in other_manager.CLUSTERS.items():
            self._parse_level({major_name: [cell.to_dict() for cell in cell_list]})

        log.info(
            f"After merging: {len(self.CELLS)} cell types in {len(self.CLUSTERS)} categories"
        )

    def save(self, output_path: str, format: Literal["toml", "json"] = "toml") -> None:
        """
        Saves the current marker configuration to a file.

        Args:
            output_path: Path where the configuration will be saved
            format: File format ('toml' or 'json')

        Raises:
            ValueError: If the format is not supported
        """
        if format not in MARKER_FORMATS:
            raise ValueError(
                f"Unsupported format: {format}. Supported formats: {MARKER_FORMATS}"
            )

        # Convert to dictionary format
        output_data = {}
        for major_name, cell_list in self.CLUSTERS.items():
            output_data[major_name] = [cell.to_dict() for cell in cell_list]

        # Save to file
        output_file = Path(output_path)
        if not output_file.parent.exists():
            output_file.parent.mkdir(parents=True)

        if format == "toml":
            # Since Python's standard library doesn't include a TOML writer,
            # we use JSON as a fallback and inform the user
            log.warning(
                "TOML writing is not supported in the standard library. Using JSON format instead."
            )
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(output_data, f, indent=2)
        elif format == "json":
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(output_data, f, indent=2)

        log.info(f"Saved marker configuration to {output_file}")

    def get_markers_by_level(
        self, level: Literal["major", "minor", "all"] = "all"
    ) -> Dict[str, List[str]]:
        """
        Gets all markers organized by cell type level.

        Args:
            level: Which level of cell types to include ('major', 'minor', or 'all')

        Returns:
            Dictionary mapping cell type names to marker lists
        """
        markers = {}

        for name, cell in self.CELLS.items():
            if level == "all" or cell.level == level:
                markers[name] = cell.markers

        return markers

    def get_all_markers(self) -> List[str]:
        """
        Gets a deduplicated list of all marker genes across all cell types.

        Returns:
            List of unique marker genes
        """
        all_markers = set()
        for cell in self.CELLS.values():
            all_markers.update(cell.markers)

        return list(all_markers)

    def filter_markers(self, min_genes_per_type: int = 3) -> Set[str]:
        """
        Filters out cell types with too few markers and returns the removed cell types.

        Args:
            min_genes_per_type: Minimum number of marker genes required for a cell type

        Returns:
            Set of cell type names that were removed due to insufficient markers
        """
        removed = set()

        for name, cell in list(self.CELLS.items()):
            if len(cell.markers) < min_genes_per_type:
                log.warning(
                    f"Removing cell type '{name}' with only {len(cell.markers)} markers"
                )

                # Remove from parent's minor list if applicable
                if cell.parent:
                    if cell in cell.parent.minor:
                        cell.parent.minor.remove(cell)

                # Remove from clusters list
                for cluster_name, cell_list in self.CLUSTERS.items():
                    if cell in cell_list:
                        cell_list.remove(cell)

                # Remove from CELLS dictionary
                del self.CELLS[name]
                removed.add(name)

        log.info(
            f"Filtered out {len(removed)} cell types with fewer than {min_genes_per_type} markers"
        )
        return removed

    def add_cell_type(
        self,
        name: str,
        markers: List[str],
        parent: Optional[str] = None,
        color: Optional[str] = None,
        metadata: Optional[Dict[str, str]] = None,
    ) -> CellType:
        """
        Adds a new cell type to the manager.

        Args:
            name: Name of the cell type
            markers: List of marker genes
            parent: Name of the parent cell type (if any)
            color: Color for visualization
            metadata: Additional metadata

        Returns:
            The newly created CellType object

        Raises:
            KeyError: If the parent cell type is specified but not found
        """
        # Check if cell type already exists
        if name in self.CELLS:
            log.warning(
                f"Cell type '{name}' already exists, updating instead of creating new"
            )
            cell_obj = self.CELLS[name]
            cell_obj.add_markers(markers)
            if color:
                cell_obj.color = color
            if metadata:
                cell_obj.metadata.update(metadata or {})
            return cell_obj

        # Process markers for case sensitivity
        if not self.case_sensitive:
            markers = [m.upper() for m in markers]

        # Create new cell type
        parent_obj = None
        if parent:
            if parent not in self.CELLS:
                raise KeyError(f"Parent cell type '{parent}' not found")
            parent_obj = self.CELLS[parent]

        # Determine level based on parent
        level = "minor" if parent_obj else "major"

        # Create cell type
        cell_obj = CellType(
            name=name,
            color=color,
            markers=markers,
            level=level,
            parent=parent_obj,
            metadata=metadata or {},
        )

        # Add to manager
        self.CELLS[name] = cell_obj

        # Add to parent if applicable
        if parent_obj:
            parent_obj.minor.append(cell_obj)
        else:
            # If no parent, add to "User-Added" cluster
            cluster_name = "User-Added"
            if cluster_name not in self.CLUSTERS:
                self.CLUSTERS[cluster_name] = []
            self.CLUSTERS[cluster_name].append(cell_obj)

        log.info(f"Added new cell type: {name} with {len(markers)} markers")
        return cell_obj

    def get_doublet_lineage_markers(self) -> Dict[str, List[str]]:
        """
        Gets a clean set of mutually exclusive lineages for doublet detection.

        This method scans the entire hierarchy and selects only the cell types
        explicitly tagged with `metadata = { doublet_lineage = true }`.

        Returns:
            Dictionary mapping cell type names to their marker lists.
        """
        lineage_markers = {}
        for name, cell in self.CELLS.items():
            if cell.metadata.get("doublet_lineage") is True:
                if cell.markers:  # Only include if it has markers
                    lineage_markers[name] = cell.markers
        
        # Fix for mouse genes: convert uppercase genes to title case if species is mouse
        if hasattr(self, '_source_file') and 'mouse' in self._source_file.lower():
            for lineage, genes in lineage_markers.items():
                lineage_markers[lineage] = [g.title() if g.isupper() else g for g in genes]
        
        log.info(f"Extracted {len(lineage_markers)} dedicated lineages for doublet detection.")
        return lineage_markers


def get_marker_manager(
    species: str,
    tissue: Optional[str] = None,
    states: Optional[List[str]] = None,
    case_sensitive: bool = True,
) -> Manager:
    """
    Factory function to build a Manager by combining base, tissue, and state markers.

    This function creates a comprehensive marker manager by layering tissue-specific
    and cell state-specific markers on top of the base markers for a species.

    Args:
        species: The species ('human', 'mouse', etc.)
        tissue: The tissue name, corresponding to a top-level key in the tissue-specific file
        states: A list of cell states, corresponding to keys in the cell-state file
        case_sensitive: Whether gene names should be case-sensitive

    Returns:
        A fully configured and combined Manager instance

    Raises:
        ValueError: If the species is not recognized
    """
    # Validate species
    species = species.lower()
    if species not in KNOWN_SPECIES:
        log.warning(f"Species '{species}' not in known list: {KNOWN_SPECIES}")

    log.info(
        f"Building marker manager for species: {species}, tissue: {tissue}, states: {states}"
    )

    try:
        # Load base markers for the species
        mgr = Manager(f"base_{species}", case_sensitive=case_sensitive)
        log.info(f"Loaded base markers for {species}")

        # Add tissue-specific markers if specified
        if tissue:
            try:
                mgr_tissue = Manager(
                    f"tissue_specific_{species}",
                    root_key=tissue,
                    case_sensitive=case_sensitive,
                )
                mgr.merge_from(mgr_tissue)
                log.info(f"Merged tissue-specific markers for '{tissue}'")
            except (FileNotFoundError, KeyError) as e:
                log.warning(
                    f"Could not load tissue-specific markers for '{tissue}': {str(e)}"
                )

        # Add cell state markers if specified
        if states:
            try:
                # Load the entire cell state file once
                mgr_states_all = Manager(
                    f"cell_state_{species}", case_sensitive=case_sensitive
                )

                # Find the parent category for states
                parent_category = "Cell States"

                # Now, selectively merge only the requested states
                merged_count = 0
                for state_name in states:
                    if state_name in mgr_states_all.CELLS:
                        # Create a temporary manager with just this one state
                        state_obj = mgr_states_all.CELLS[state_name]
                        temp_mgr = Manager(
                            f"cell_state_{species}",
                            root_key=parent_category,
                            case_sensitive=case_sensitive,
                        )
                        temp_mgr.CELLS = {state_name: state_obj}
                        temp_mgr.CLUSTERS = {parent_category: [state_obj]}

                        # Merge this state into the main manager
                        mgr.merge_from(temp_mgr)
                        log.info(f"Merged cell state markers for '{state_name}'")
                        merged_count += 1
                    else:
                        log.warning(
                            f"State '{state_name}' not found in cell state file"
                        )

                if merged_count == 0 and states:
                    log.warning(f"None of the requested states {states} were found")

            except FileNotFoundError as e:
                log.warning(f"Could not load cell state markers: {str(e)}")

        log.info(f"Manager built successfully with {len(mgr.CELLS)} cell types")
        return mgr

    except Exception as e:
        log.error(f"Error building marker manager: {str(e)}")
        raise
