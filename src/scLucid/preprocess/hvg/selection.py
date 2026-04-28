"""HVG set selection, suggestion, and comparison.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

try:
    from matplotlib_venn import venn2, venn3
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

from ...utils import use_layer_as_X
from ..config import HVGConfig
from .core import find_hvgs

log = logging.getLogger(__name__)

def suggest_hvg_choice(adata: AnnData, hvg_keys: List[str], mode: str) -> None:
    """
    Provides data-driven guidance for HVG set selection by analyzing overlap.
    """
    if len(hvg_keys) < 2:
        log.info("Guidance is most useful when comparing 2 or more HVG sets.")
        return

    sets = [set(adata.var_names[adata.var[key]]) for key in hvg_keys]

    # --- Quantitative Analysis ---
    intersection_set = set.intersection(*sets)
    union_set = set.union(*sets)
    jaccard_index = len(intersection_set) / len(union_set) if len(union_set) > 0 else 0

    # --- Build Report ---
    msg = ["=" * 50, "==== HVG Selection Guidance ====", "=" * 50]
    msg.append(f"Comparing {len(hvg_keys)} HVG sets: {', '.join(hvg_keys)}")
    for i, key in enumerate(hvg_keys):
        msg.append(f"- Set '{key}': {len(sets[i])} genes")

    msg.append("\n--- Overlap Analysis ---")
    msg.append(f"- Intersection (genes in all sets): {len(intersection_set)} genes")
    msg.append(f"- Union (genes in any set): {len(union_set)} genes")
    msg.append(f"- Jaccard Similarity Index: {jaccard_index:.3f} (Intersection / Union)")

    msg.append(f"\n--- Recommendation for your chosen mode ('{mode}') ---")

    # --- Generate Tailored Advice ---
    if jaccard_index > 0.7:
        msg.append("Data-driven verdict: **High Overlap**.")
        msg.append("Both methods identify a very similar core set of variable genes.")
        if mode == "intersection":
            msg.append(
                "Your choice of 'intersection' is a safe and robust strategy. You will get a high-confidence set of HVGs."
            )
        elif mode == "union":
            msg.append(
                "Your choice of 'union' is also reasonable. It will add a few extra genes without a high risk of introducing noise."
            )

    elif 0.4 <= jaccard_index <= 0.7:
        msg.append("Data-driven verdict: **Moderate Overlap**.")
        msg.append("The methods agree on a core set of genes but also identify unique ones.")
        if mode == "intersection":
            msg.append(
                "Your choice of 'intersection' is the most conservative and reproducible option. You will get a high-confidence set but may miss some subtle biological signals."
            )
        elif mode == "union":
            msg.append(
                "Your choice of 'union' is more inclusive and better for discovery, but may introduce noise. Be sure to check for over-clustering in downstream analysis."
            )

    else:  # Low overlap
        msg.append("Data-driven verdict: **Low Overlap**.")
        msg.append(
            "⚠️ **Warning:** The selected methods are identifying very different sets of genes. This could be due to strong batch effects or fundamental differences in the algorithms."
        )
        if mode == "intersection":
            msg.append(
                f"Your choice of 'intersection' will result in a very small set of {len(intersection_set)} genes. This may not be enough for stable downstream analysis. Please verify."
            )
        elif mode == "union":
            msg.append(
                "Your choice of 'union' will combine two very different gene lists, which could be risky. It is highly recommended to first visualize the UMAPs from each HVG set individually to understand why they differ so much."
            )
        msg.append(
            "\n**Suggestion:** The 'custom' method is often more robust for multi-sample datasets than the standard 'scanpy' method. Consider trusting the 'custom' set or investigating potential batch effects further."
        )

    print("\n".join(msg))


def select_hvg_sets(
    adata: AnnData,
    hvg_keys: Union[str, List[str]],
    mode: Literal["direct", "intersection", "union", "difference"] = "direct",
    subset: bool = True,
    keep_raw: bool = True,
    copy: bool = False,
    output_key: str = "highly_variable_selected",
    plot_venn: bool = True,
    show_stats: bool = True,
    show_suggestion: bool = True,
    save_dir: Optional[str] = None,
    **kwargs,
) -> AnnData:
    """
    Select HVG genes using one or more masks, with set operations, summary and visualization.
    """
    if isinstance(hvg_keys, str):
        hvg_keys = [hvg_keys]
    for k in hvg_keys:
        if k not in adata.var:
            raise KeyError(f"HVG key '{k}' not found in adata.var.")

    # --- Suggestion for HVG set choice ---
    if show_suggestion:
        # Call the new, data-driven guidance function
        suggest_hvg_choice(adata, hvg_keys, mode)

    hvg_sets = [set(adata.var_names[adata.var[k]]) for k in hvg_keys]
    set_names = hvg_keys

    # --- Combine sets ---
    if mode == "direct":
        combined_set = hvg_sets[0]
    elif mode == "intersection":
        combined_set = set.intersection(*hvg_sets)
    elif mode == "union":
        combined_set = set.union(*hvg_sets)
    elif mode == "difference":
        if len(hvg_sets) < 2:
            log.warning("Difference mode needs at least 2 masks. Falling back to direct.")
            combined_set = hvg_sets[0]
        else:
            combined_set = hvg_sets[0].copy()
            for s in hvg_sets[1:]:
                combined_set -= s
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    # --- Output stats ---
    stats_msg = []
    if show_stats:
        stats_msg.append("==== HVG Set Statistics ====")
        for i, (name, s) in enumerate(zip(set_names, hvg_sets)):
            stats_msg.append(f"Set {i + 1} [{name}]: {len(s)} genes")
        if len(hvg_sets) == 2:
            intersect = hvg_sets[0] & hvg_sets[1]
            only0 = hvg_sets[0] - hvg_sets[1]
            only1 = hvg_sets[1] - hvg_sets[0]
            union = hvg_sets[0] | hvg_sets[1]
            stats_msg.append(f"Intersection: {len(intersect)}")
            stats_msg.append(f"Only {set_names[0]}: {len(only0)}")
            stats_msg.append(f"Only {set_names[1]}: {len(only1)}")
            stats_msg.append(f"Union: {len(union)}")
        elif len(hvg_sets) == 3:
            intersect = set.intersection(*hvg_sets)
            union = set.union(*hvg_sets)
            stats_msg.append(f"Intersection (all): {len(intersect)}")
            stats_msg.append(f"Union: {len(union)}")
        stats_msg.append(f"Selected set [{mode}]: {len(combined_set)} genes")
        stats_msg = "\n".join(stats_msg)
        print(stats_msg)
        if save_dir:
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            with open(f"{save_dir}/hvg_set_stats.txt", "w") as f:
                f.write(stats_msg)

    # --- Plot Venn diagram if needed ---
    if plot_venn and HAS_VENN and (2 <= len(hvg_sets) <= 3):
        plt.figure(figsize=(6, 5))
        if len(hvg_sets) == 2:
            venn2(subsets=hvg_sets, set_labels=set_names)
        elif len(hvg_sets) == 3:
            venn3(subsets=hvg_sets, set_labels=set_names)
        plt.title(f"HVG Sets Venn Diagram ({mode})")
        if save_dir:
            save_path = Path(save_dir)
            # save_path.mkdir(parents=True, exist_ok=True)
            plt.savefig(f"{save_dir}/hvg_venn_{mode}.png", dpi=150, bbox_inches="tight")
        plt.show()
    elif plot_venn and not HAS_VENN and (2 <= len(hvg_sets) <= 3):
        log.warning("matplotlib_venn is not installed. Skipping Venn plot.")

    mask_combined = adata.var_names.isin(list(combined_set))
    adata.var[output_key] = mask_combined

    log.info(f"Created final HVG mask in '.var['{output_key}']' with {mask_combined.sum()} genes.")

    if subset:
        if keep_raw and adata.raw is None:
            adata.raw = adata.copy()

        if copy:
            adata_subset = adata[:, mask_combined].copy()
            log.info(
                f"Created a new subsetted AnnData object with {mask_combined.sum()} final HVGs."
            )
            return adata_subset
        else:
            adata._inplace_subset_var(mask_combined)
            log.info(f"Subsetted AnnData object in-place to {mask_combined.sum()} final HVGs.")
            return adata

    return adata


