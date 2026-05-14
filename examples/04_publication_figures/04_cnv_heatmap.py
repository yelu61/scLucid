"""CNV profile heatmap — chromosome-ordered copy-number signal per cell group.

Synthesises malignant + normal cells with planted CNV gains/losses on a few
chromosomes, then renders a heatmap with cells on the y-axis (sorted by
group) and chromosomal bin position on the x-axis.

This is the Figure 2 staple of tumor single-cell papers: visually shows
which cells carry tumor-specific copy-number alterations vs which match the
diploid reference.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scLucid.plotting.theme import apply_theme

OUTPUT_DIR = Path("results/publication_figures")


def _make_cnv_profile(seed: int = 0):
    """Synthesise a (n_cells, n_bins) CNV matrix with two cell groups."""
    rng = np.random.default_rng(seed)
    n_bins_per_chrom = [50] * 22  # 22 autosomes
    n_bins = sum(n_bins_per_chrom)
    n_malignant = 200
    n_normal = 200

    # Normal cells: diploid baseline with small noise centered at 0
    normal = rng.normal(0.0, 0.1, size=(n_normal, n_bins))

    # Malignant cells: planted gains on chr7, losses on chr10 & chr17
    malignant = rng.normal(0.0, 0.1, size=(n_malignant, n_bins))
    chrom_starts = np.concatenate(([0], np.cumsum(n_bins_per_chrom)))
    malignant[:, chrom_starts[6] : chrom_starts[7]] += rng.normal(0.4, 0.05, (n_malignant, 50))
    malignant[:, chrom_starts[9] : chrom_starts[10]] -= rng.normal(0.5, 0.05, (n_malignant, 50))
    malignant[:, chrom_starts[16] : chrom_starts[17]] -= rng.normal(0.3, 0.05, (n_malignant, 50))

    cnv = np.vstack([normal, malignant])
    labels = ["Normal"] * n_normal + ["Malignant"] * n_malignant

    obs = pd.DataFrame(
        {"cell_group": pd.Categorical(labels, categories=["Normal", "Malignant"])},
        index=[f"cell_{i:04d}" for i in range(len(labels))],
    )
    return cnv, obs, chrom_starts


def main() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cnv, obs, chrom_starts = _make_cnv_profile()
    apply_theme("nature")

    fig = plt.figure(figsize=(7, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.04, 1], wspace=0.04)
    ax_groups = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    # Group color bar on the left
    group_palette = {"Normal": "#1f78b4", "Malignant": "#d9352a"}
    group_codes = np.array(
        [list(group_palette).index(g) for g in obs["cell_group"]]
    ).reshape(-1, 1)
    ax_groups.imshow(
        group_codes,
        aspect="auto",
        cmap=plt.matplotlib.colors.ListedColormap(list(group_palette.values())),
        interpolation="nearest",
    )
    ax_groups.set_yticks([100, 300])
    ax_groups.set_yticklabels(["Normal", "Malignant"], rotation=90, va="center", fontsize=8)
    ax_groups.set_xticks([])
    for spine in ax_groups.spines.values():
        spine.set_visible(False)

    # Main heatmap
    im = ax_heat.imshow(cnv, aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6)

    # Chromosome boundaries
    for boundary in chrom_starts[1:-1]:
        ax_heat.axvline(boundary - 0.5, color="black", linewidth=0.3)
    chrom_centers = (chrom_starts[:-1] + chrom_starts[1:]) / 2
    ax_heat.set_xticks(chrom_centers)
    ax_heat.set_xticklabels(
        [str(i + 1) for i in range(22)], fontsize=7
    )
    ax_heat.set_yticks([])
    ax_heat.set_xlabel("Chromosome")
    ax_heat.set_ylabel("Cells")
    ax_heat.set_title("CNV profile by cell group")

    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.6, pad=0.02)
    cbar.set_label("Mean modified expression")
    cbar.ax.tick_params(labelsize=8)

    out = OUTPUT_DIR / "fig04_cnv_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
