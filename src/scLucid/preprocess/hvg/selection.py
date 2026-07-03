"""HVG set selection, suggestion, and comparison.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

import matplotlib.pyplot as plt
from anndata import AnnData

try:
    from matplotlib_venn import venn2, venn3

    HAS_VENN = True
except ImportError:
    HAS_VENN = False

from scLucid.plotting.plotting_utils import _is_interactive_backend

log = logging.getLogger(__name__)

__all__ = [
    "suggest_hvg_choice",
    "select_hvg_sets",
    "select_and_audit_hvgs",
]


def _is_protected_hvg_role(role: object) -> bool:
    role_text = str(role or "").lower()
    return any(
        token in role_text for token in ("protected", "marker", "program", "biology", "curated")
    )


def suggest_hvg_choice(
    adata: AnnData,
    hvg_keys: List[str],
    mode: str = "auto",
    set_roles: Optional[Dict[str, str]] = None,
) -> Dict[str, object]:
    """Return structured guidance for choosing an HVG set operation.

    ``set_roles`` lets callers distinguish ordinary algorithmic HVG masks from
    curated/protected biology masks. Small protected marker sets naturally have
    low Jaccard overlap with a 2k variance HVG set, so overlap alone would
    recommend a destructive intersection.
    """
    if len(hvg_keys) < 2:
        return {
            "requested_mode": mode,
            "recommended_mode": mode if mode != "auto" else "direct",
            "jaccard_index": 1.0,
            "overlap_level": "single",
            "risk": "low",
            "n_genes_per_set": {hvg_keys[0]: int(adata.var[hvg_keys[0]].sum())} if hvg_keys else {},
            "messages": ["HVG guidance is most useful when comparing two or more masks."],
        }

    sets = [set(adata.var_names[adata.var[key]]) for key in hvg_keys]
    intersection_set = set.intersection(*sets)
    union_set = set.union(*sets)
    jaccard_index = len(intersection_set) / len(union_set) if union_set else 0.0

    set_roles = set_roles or {}
    protected_keys = [key for key in hvg_keys if _is_protected_hvg_role(set_roles.get(key))]
    has_protected_set = bool(protected_keys)

    if has_protected_set and len(hvg_keys) >= 2:
        overlap_level = "protected_biology"
        recommended_mode = "union"
        risk = "review"
    elif jaccard_index > 0.7:
        overlap_level = "high"
        recommended_mode = "union"
        risk = "low"
    elif jaccard_index >= 0.4:
        overlap_level = "moderate"
        recommended_mode = "intersection"
        risk = "moderate"
    else:
        overlap_level = "low"
        recommended_mode = "intersection"
        risk = "high"

    effective_mode = recommended_mode if mode == "auto" else mode
    messages = [
        "=" * 50,
        "==== HVG Selection Guidance ====",
        "=" * 50,
        f"Comparing {len(hvg_keys)} HVG sets: {', '.join(hvg_keys)}",
    ]
    for key, genes in zip(hvg_keys, sets):
        messages.append(f"- Set '{key}': {len(genes)} genes")
    messages.extend(
        [
            "",
            "--- Overlap Analysis ---",
            f"- Intersection (genes in all sets): {len(intersection_set)} genes",
            f"- Union (genes in any set): {len(union_set)} genes",
            f"- Jaccard Similarity Index: {jaccard_index:.3f} (Intersection / Union)",
            "",
            f"Recommended mode: {recommended_mode} (overlap={overlap_level}, risk={risk})",
            f"Requested/effective mode: {mode} -> {effective_mode}",
        ]
    )
    if has_protected_set:
        messages.append(
            "Protected biology HVG set detected; auto mode favors union so marker/program genes are not lost."
        )
    if risk == "high":
        messages.append(
            "Warning: HVG methods disagree strongly; inspect batch effects and consider conservative intersection."
        )

    return {
        "requested_mode": mode,
        "recommended_mode": recommended_mode,
        "effective_mode": effective_mode,
        "jaccard_index": float(jaccard_index),
        "overlap_level": overlap_level,
        "risk": risk,
        "set_roles": {key: set_roles.get(key, "") for key in hvg_keys},
        "protected_hvg_keys": protected_keys,
        "n_genes_per_set": {key: len(genes) for key, genes in zip(hvg_keys, sets)},
        "n_intersection": len(intersection_set),
        "n_union": len(union_set),
        "messages": messages,
    }


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
    set_roles: Optional[Dict[str, str]] = None,
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
        suggestion = suggest_hvg_choice(adata, hvg_keys, mode, set_roles=set_roles)
        print("\n".join(suggestion.get("messages", [])))

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
        if _is_interactive_backend():
            plt.show()
        else:
            plt.close()
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


def select_and_audit_hvgs(
    adata: AnnData,
    *,
    hvg_keys: Union[str, List[str]],
    mode: Literal["auto", "direct", "intersection", "union", "difference"] = "union",
    subset: bool = True,
    keep_raw: bool = False,
    output_key: str = "highly_variable_selected",
    evaluate_stability: bool = False,
    stability_key: Optional[str] = None,
    stability_kwargs: Optional[Dict] = None,
    save_dir: Optional[str] = None,
    set_roles: Optional[Dict[str, str]] = None,
    **kwargs,
) -> tuple[AnnData, Dict[str, object]]:
    """Select final HVGs and store a compact audit summary.

    This wrapper intentionally assumes HVG masks already exist in ``adata.var``.
    Use ``find_hvgs`` for method-specific detection, then this function for the
    set-operation decision and optional stability check.
    """
    if isinstance(hvg_keys, str):
        hvg_key_list = [hvg_keys]
    else:
        hvg_key_list = list(hvg_keys)

    suggestion = (
        suggest_hvg_choice(adata, hvg_key_list, mode, set_roles=set_roles)
        if len(hvg_key_list) >= 2
        else None
    )
    effective_mode = (
        suggestion.get("recommended_mode", "direct")
        if mode == "auto" and isinstance(suggestion, dict)
        else mode
    )

    audit: Dict[str, object] = {
        "hvg_keys": hvg_key_list,
        "mode": mode,
        "effective_mode": effective_mode,
        "subset": bool(subset),
        "output_key": output_key,
    }
    if set_roles:
        audit["set_roles"] = {key: set_roles.get(key, "") for key in hvg_key_list}

    if suggestion is not None:
        audit["suggestion"] = suggestion
        if suggestion.get("messages"):
            print("\n".join(suggestion["messages"]))
        if suggestion.get("risk") == "high":
            audit.setdefault("warnings", []).append(
                "HVG methods have low overlap; selected conservative intersection in auto mode."
                if mode == "auto"
                else "HVG methods have low overlap; inspect selected mode carefully."
            )

    stability_result = None
    if evaluate_stability:
        try:
            from .stability import evaluate_hvg_stability

            key_for_stability = stability_key or hvg_key_list[0]
            stability_result = evaluate_hvg_stability(
                adata,
                hvg_key=key_for_stability,
                **(stability_kwargs or {}),
            )
            audit["stability_key"] = key_for_stability
            audit["stability_available"] = True
            stability_summary = (
                adata.uns.get("sclucid", {}).get("preprocess", {}).get("hvg_stability", {})
            )
            audit["stability_summary"] = stability_summary
            if stability_summary.get("overall_score", 1.0) < 0.5:
                audit.setdefault("warnings", []).append(
                    "HVG stability is low; consider increasing n_top_genes or comparing methods."
                )
        except Exception as exc:
            log.warning("HVG stability evaluation failed: %s", exc)
            audit["stability_available"] = False
            audit["stability_error"] = str(exc)

    result = select_hvg_sets(
        adata,
        hvg_keys=hvg_key_list,
        mode=effective_mode,
        subset=subset,
        keep_raw=keep_raw,
        output_key=output_key,
        save_dir=save_dir,
        set_roles=set_roles,
        show_suggestion=kwargs.pop("show_suggestion", False),
        **kwargs,
    )

    selected_count = (
        int(result.var[output_key].sum()) if output_key in result.var else result.n_vars
    )
    audit["n_selected"] = selected_count
    if stability_result is not None:
        audit["stability_result_type"] = type(stability_result).__name__

    result.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["hvg_selection_audit"] = audit
    return result, audit
