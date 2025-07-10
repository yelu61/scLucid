from collections.abc import Sequence
from dataclasses import dataclass
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


@dataclass
class CellType:
    name: str
    color: str | None
    markers: list[str]
    level: Literal["major", "minor"]
    minor: list[str] | None = None

    def __post_init__(self) -> None:
        if self.color == "":
            self.color = None


class Manager:
    def __init__(self, config_file: str | Path) -> None:
        with open(config_file, "rb") as f:
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
        if key not in self.CELLS:
            raise KeyError(f"{key} is not in CELLS")
        else:
            return self.CELLS[key]
        
    def intersect_with(self, adata: AnnData):
        genes = adata.var_names
        for cell_name, cell_type in self.CELLS.items():
            self.CELLS[cell_name].markers = list(set(cell_type.markers).intersection(genes))

    def query(
        self,
        info: Literal["color", "markers"] = "color",
        key: str | Sequence[str] | None = None,
        cluster: str | None = None,
        include_major_type: bool = True,
    ) -> dict[str, str | list[str]]:
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

    def show_tree(self):
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
