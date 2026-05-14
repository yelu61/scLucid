"""UMAP scatter colored by cell type — the Figure-1 staple.

Generates a synthetic clustered AnnData, runs PCA + UMAP via scanpy, then
renders a UMAP scatter with cells colored by cell type. The resulting PDF
uses embedded TrueType fonts so every text label can be edited in
Illustrator.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from scLucid.plotting.theme import apply_theme

OUTPUT_DIR = Path("results/publication_figures")


def _make_synthetic_data(n_cells_per_type: int = 200, n_genes: int = 200, seed: int = 0):
    """Synthesise a 4-cell-type dataset with clearly separable signatures."""
    rng = np.random.default_rng(seed)
    cell_types = ["T cell", "B cell", "Monocyte", "NK cell"]

    blocks = []
    labels = []
    for idx, ct in enumerate(cell_types):
        base = rng.poisson(2, size=(n_cells_per_type, n_genes)).astype(float)
        # Plant a 30-gene signature for each cell type
        signature_start = idx * 30
        signature_end = signature_start + 30
        base[:, signature_start:signature_end] += rng.poisson(8, size=(n_cells_per_type, 30))
        blocks.append(base)
        labels.extend([ct] * n_cells_per_type)

    X = np.vstack(blocks)
    var = pd.DataFrame(index=[f"gene_{i:03d}" for i in range(n_genes)])
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(labels, categories=cell_types)},
        index=[f"cell_{i:04d}" for i in range(X.shape[0])],
    )
    return AnnData(X=X, obs=obs, var=var)


def main() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = _make_synthetic_data()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.4, random_state=0)

    apply_theme("nature")
    palette = {
        "T cell": "#E64B35",
        "B cell": "#4DBBD5",
        "Monocyte": "#00A087",
        "NK cell": "#3C5488",
    }

    fig, ax = plt.subplots(figsize=(5, 4.5))
    coords = adata.obsm["X_umap"]
    for cell_type, color in palette.items():
        mask = (adata.obs["cell_type"] == cell_type).to_numpy()
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=8,
            color=color,
            alpha=0.85,
            edgecolors="none",
            label=cell_type,
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Cell type composition")
    ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=8,
        markerscale=1.5,
    )
    ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    out = OUTPUT_DIR / "fig01_umap_annotation.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
