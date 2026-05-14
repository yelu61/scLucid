"""Volcano plot — differential expression between two conditions.

Synthesises bulk-aggregated single-cell expression for two conditions
("treated" vs "control"), computes per-gene log2 fold-change and
significance, and renders a volcano plot with top hits labeled. This is the
standard differential-expression figure for tumor-vs-normal or
responder-vs-non-responder comparisons.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from scipy import stats

from scLucid.plotting.theme import apply_theme

OUTPUT_DIR = Path("results/publication_figures")


def _simulate_de(seed: int = 0):
    """Simulate per-gene mean expression in two conditions."""
    rng = np.random.default_rng(seed)
    n_genes = 800

    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    # Most genes are unchanged
    mean_control = rng.lognormal(2.0, 0.6, n_genes)
    mean_treated = mean_control * rng.lognormal(0.0, 0.15, n_genes)

    # Plant 40 strongly up-regulated, 40 strongly down-regulated genes
    up_idx = rng.choice(n_genes, 40, replace=False)
    down_pool = np.setdiff1d(np.arange(n_genes), up_idx)
    down_idx = rng.choice(down_pool, 40, replace=False)
    mean_treated[up_idx] *= rng.lognormal(1.5, 0.2, len(up_idx))
    mean_treated[down_idx] *= rng.lognormal(-1.5, 0.2, len(down_idx))

    # Compute log2 fold change
    log2fc = np.log2((mean_treated + 1) / (mean_control + 1))

    # P-value from a fake replicated mean using a small noise model
    n_replicates = 6
    control_reps = rng.lognormal(np.log(mean_control + 1)[:, None], 0.18, (n_genes, n_replicates))
    treated_reps = rng.lognormal(np.log(mean_treated + 1)[:, None], 0.18, (n_genes, n_replicates))
    _, p_values = stats.ttest_ind(np.log(treated_reps + 1), np.log(control_reps + 1), axis=1)
    p_values = np.clip(p_values, 1e-300, 1.0)

    return gene_names, log2fc, p_values, up_idx, down_idx


def main() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    genes, log2fc, p_values, planted_up, planted_down = _simulate_de()
    neg_log10_p = -np.log10(p_values)

    apply_theme("nature")

    fig, ax = plt.subplots(figsize=(5.5, 5))

    # Significance + magnitude thresholds
    fc_thresh = 1.0
    p_thresh = 0.05

    significant_up = (log2fc > fc_thresh) & (p_values < p_thresh)
    significant_down = (log2fc < -fc_thresh) & (p_values < p_thresh)
    non_significant = ~(significant_up | significant_down)

    ax.scatter(
        log2fc[non_significant],
        neg_log10_p[non_significant],
        s=8,
        color="#cccccc",
        alpha=0.5,
        edgecolors="none",
    )
    ax.scatter(
        log2fc[significant_up],
        neg_log10_p[significant_up],
        s=10,
        color="#d9352a",
        alpha=0.85,
        edgecolors="none",
        label=f"Up (n={int(significant_up.sum())})",
    )
    ax.scatter(
        log2fc[significant_down],
        neg_log10_p[significant_down],
        s=10,
        color="#1f78b4",
        alpha=0.85,
        edgecolors="none",
        label=f"Down (n={int(significant_down.sum())})",
    )

    # Reference lines
    ax.axvline(fc_thresh, color="black", linestyle="--", linewidth=0.5)
    ax.axvline(-fc_thresh, color="black", linestyle="--", linewidth=0.5)
    ax.axhline(-np.log10(p_thresh), color="black", linestyle="--", linewidth=0.5)

    # Label the top 5 up- and down-regulated genes by p-value × |fc|
    sig_indices = np.where(significant_up | significant_down)[0]
    if sig_indices.size:
        priority = neg_log10_p[sig_indices] * np.abs(log2fc[sig_indices])
        top = sig_indices[np.argsort(priority)[-10:]]
        texts = [
            ax.text(log2fc[i], neg_log10_p[i], genes[i], fontsize=7) for i in top
        ]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.3))

    ax.set_xlabel("log$_2$ fold change (treated / control)")
    ax.set_ylabel("$-\\log_{10}$ p-value")
    ax.set_title("Differential expression: treated vs control")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    out = OUTPUT_DIR / "fig03_volcano_de.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
