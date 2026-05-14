"""Marker gene heatmap — per-cluster expression of canonical markers.

Generates a 4-cell-type synthetic dataset with planted marker signatures,
computes the per-cluster mean expression, and renders a heatmap with cell
types on the y-axis and marker genes on the x-axis. Z-score normalisation
across cell types is applied so the visual emphasises specificity rather
than absolute magnitude.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.plotting.theme import apply_theme

OUTPUT_DIR = Path("results/publication_figures")


def _make_marker_data(seed: int = 0):
    """Synthesise expression with planted markers per cell type."""
    rng = np.random.default_rng(seed)
    cell_types = ["T cell", "B cell", "Monocyte", "NK cell"]
    n_cells_per_type = 150

    # Canonical marker names per cell type (4 markers each = 16 displayed)
    markers = {
        "T cell": ["CD3D", "CD3E", "CD8A", "CD4"],
        "B cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "Monocyte": ["CD14", "LYZ", "S100A8", "S100A9"],
        "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    }
    all_markers = [m for marker_list in markers.values() for m in marker_list]

    expression = []
    labels = []
    rows_per_type = []
    for ct in cell_types:
        base = rng.poisson(1, size=(n_cells_per_type, len(all_markers))).astype(float)
        # Boost expression on this cell type's markers
        for marker in markers[ct]:
            idx = all_markers.index(marker)
            base[:, idx] += rng.poisson(10, size=n_cells_per_type)
        expression.append(base)
        labels.extend([ct] * n_cells_per_type)
        rows_per_type.append(n_cells_per_type)

    X = np.vstack(expression)
    var = pd.DataFrame(index=all_markers)
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(labels, categories=cell_types)},
        index=[f"cell_{i:04d}" for i in range(X.shape[0])],
    )
    return AnnData(X=X, obs=obs, var=var), markers


def main() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata, markers = _make_marker_data()

    # Per-cell-type log-mean expression
    cell_types = list(markers.keys())
    matrix = np.zeros((len(cell_types), adata.n_vars))
    for i, ct in enumerate(cell_types):
        mask = (adata.obs["cell_type"] == ct).to_numpy()
        matrix[i] = np.log1p(adata.X[mask].mean(axis=0))

    # Z-score per gene across cell types so the heatmap shows specificity
    matrix = (matrix - matrix.mean(axis=0, keepdims=True)) / (
        matrix.std(axis=0, keepdims=True) + 1e-8
    )

    apply_theme("nature")

    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=9)
    ax.set_xticks(range(adata.n_vars))
    ax.set_xticklabels(adata.var_names, rotation=45, ha="right", fontsize=8)

    # Group separators between cell-type marker blocks
    for boundary in range(4, adata.n_vars, 4):
        ax.axvline(boundary - 0.5, color="black", linewidth=0.4)

    ax.set_title("Cell-type marker expression (z-score)")
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("z-score")
    cbar.ax.tick_params(labelsize=8)

    out = OUTPUT_DIR / "fig02_marker_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
