"""
Enhanced HVG selection and diagnostics for single-cell RNA-seq data.

Provides config-driven, reproducible, and fully traceable workflows
for highly variable gene (HVG) selection, with batch/sample awareness,
gene-type exclusion, automatic reporting, and large data support.
"""

import logging
from dataclasses import asdict, replace
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

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

from ..utils.utils import use_layer_as_X
from .config import HVGConfig

log = logging.getLogger(__name__)

__all__ = [
    "find_hvgs",
    "suggest_hvg_choice",
    "select_hvg_sets",
    "evaluate_hvg_stability",
    "plot_hvg_metrics",
]


# --- Helper Functions ---#
def _diagnose_input_for_hvg(X, max_n=10000):
    """Compute and log basic statistics for the input matrix (dense or sparse)."""
    n_cells, n_genes = X.shape
    arr = X
    if hasattr(arr, "toarray"):
        arr = arr.toarray()
    if arr.shape[0] > max_n:
        arr = arr[:max_n]
    if arr.shape[1] > max_n:
        arr = arr[:, :max_n]
    mean = np.mean(arr)
    std = np.std(arr)
    min_val = np.min(arr)
    max_val = np.max(arr)
    zero_frac = np.mean(arr == 0)
    log.info(
        f"[HVG Input] shape={n_cells}x{n_genes}, mean={mean:.2f}, std={std:.2f}, min={min_val:.2f}, max={max_val:.2f}, zero_frac={zero_frac:.2%}"
    )
    return dict(
        mean=mean,
        std=std,
        min=min_val,
        max=max_val,
        zero_frac=zero_frac,
        n_cells=n_cells,
        n_genes=n_genes,
    )


def _write_hvg_report(
    report_path: Path,
    stats: dict,
    n_hvg: int,
    config: HVGConfig,
    gene_type_counts: dict,
) -> None:
    """Write a simple markdown HVG report."""
    with open(report_path, "w") as f:
        f.write("# HVG Selection Report\n\n")
        f.write(f"**Method:** {config.method}\n\n")
        f.write(
            f"**Input shape:** {stats.get('n_cells')} cells × {stats.get('n_genes')} genes\n\n"
        )
        f.write("## Input Statistics\n")
        for k, v in stats.items():
            f.write(f"- {k}: {v:.3g}\n")
        f.write(f"\n## HVG count: {n_hvg}\n")
        f.write("\n## Parameters\n")
        for k, v in config.__dict__.items():
            f.write(f"- {k}: {v}\n")
        f.write("\n## Excluded gene types\n")
        for k, v in gene_type_counts.items():
            f.write(f"- {k}: {v}\n")


def _detect_gene_types(var_names: pd.Index) -> Dict[str, np.ndarray]:
    """Detect common gene types (mito/ribo/hb/heatshock/IEG) by naming."""
    gene_names = var_names.astype(str)
    gene_types = {
        "mitochondrial": gene_names.str.match(r"^(MT-|mt-|MT\.|mt\.)")
        | gene_names.str.match(r"^(MTRNR|MTATP|MTND|MTCO|MTCYB)"),
        "ribosomal": gene_names.str.match(r"^(RP[SL]|Rp[sl])"),
        "hemoglobin": gene_names.str.match(r"^HB[^P]")
        | gene_names.str.contains(r"^hemoglobin", case=False),
        "heat_shock": gene_names.str.match(r"^HSP")
        | gene_names.str.contains(r"^heat shock", case=False),
        "immediate_early": np.isin(
            gene_names, ["FOS", "JUN", "JUNB", "EGR1", "NR4A1", "ZFP36"]
        ),
    }
    return gene_types


def _exclude_genes(
    adata: AnnData,
    hvg_mask: np.ndarray,
    exclude_types: List[str] = ["mitochondrial", "ribosomal"],
    gene_types: Optional[Dict[str, np.ndarray]] = None,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Exclude specific gene types from HVG selection."""
    if gene_types is None:
        gene_types = _detect_gene_types(adata.var_names)
    excluded_counts = {}
    updated_mask = hvg_mask.copy()
    for gene_type in exclude_types:
        if gene_type in gene_types:
            type_mask = gene_types[gene_type]
            excluded_counts[gene_type] = int((hvg_mask & type_mask).sum())
            updated_mask = updated_mask & ~type_mask
            log.info(f"[HVG] Excluded {excluded_counts[gene_type]} {gene_type} genes")
        else:
            log.warning(f"[HVG] Gene type '{gene_type}' not recognized")
    return updated_mask, excluded_counts


def _get_sample_specific_genes(
    adata: AnnData,
    sample_key: str,
    n_specific_genes: int,
    layer: Optional[str] = None,
    method: str = "t-test",
) -> List[str]:
    """
    Identifies top genes specific to each sample using differential expression analysis.
    """
    n_samples = adata.obs[sample_key].nunique()
    if n_samples <= 1:
        log.info(
            f"[HVG] Only one sample group found for '{sample_key}'. Skipping sample-specific gene identification."
        )
        return []
    log.info(f"[HVG] Identifying sample-specific genes across {n_samples} groups...")

    original_X = None
    temp_adata = adata
    try:
        if layer is not None:
            if layer not in temp_adata.layers:
                raise KeyError(f"Layer '{layer}' not found in AnnData object")
            original_X = temp_adata.X
            temp_adata.X = temp_adata.layers[layer]
        valid_methods = ["t-test", "wilcoxon", "logreg"]
        if method not in valid_methods:
            log.warning(
                f"Method '{method}' not in {valid_methods}. Defaulting to 't-test'."
            )
            method = "t-test"
        sc.tl.rank_genes_groups(
            temp_adata,
            groupby=sample_key,
            method=method,
            n_genes=n_specific_genes,
            pts=True,
        )
        specific_genes_df = pd.DataFrame(temp_adata.uns["rank_genes_groups"]["names"])
        specific_genes = list(np.unique(specific_genes_df.values.flatten()))
        log.info(f"[HVG] Identified {len(specific_genes)} sample-specific genes")
        return specific_genes
    except Exception as e:
        log.error(f"[HVG] Could not identify sample-specific genes: {str(e)}")
        log.exception("Detailed error:")
        return []
    finally:
        if original_X is not None:
            temp_adata.X = original_X
        if "rank_genes_groups" in temp_adata.uns:
            del temp_adata.uns["rank_genes_groups"]


def _to_savable_dict(d: dict) -> dict:
    """Recursively convert a dictionary to be h5ad-savable."""
    savable = {}
    for k, v in d.items():
        if isinstance(v, tuple):
            savable[k] = list(v)  # Convert tuples to lists
        elif isinstance(v, dict):
            savable[k] = _to_savable_dict(v)
        elif isinstance(v, (str, int, float, bool, list)) or v is None:
            savable[k] = v
        else:
            # For other complex objects, convert to string representation
            savable[k] = str(v)
    return savable


# --- Main Functions ---#
def find_hvgs(
    adata: AnnData,
    config: Optional[HVGConfig] = None,
    input_layer: str = "normalized",
    **kwargs,
) -> AnnData:
    """
    Config-driven, reproducible HVG selection with diagnostics and reporting.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = HVGConfig()
    else:
        active_config = replace(config)  # Use a copy

    # Apply overrides from kwargs
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)
        else:
            log.warning(f"Ignoring unknown parameter: '{key}'")

    # Extract parameters from the final config for use in the function
    force = kwargs.get("force", False)
    report = (
        active_config.report
        if hasattr(active_config, "report")
        else kwargs.get("report", False)
    )
    plot = kwargs.get("plot", active_config.plot)
    save_dir = active_config.save_dir
    n_top_genes = active_config.n_top_genes
    method = active_config.method
    batch_key = active_config.batch_key
    flavor = active_config.flavor
    exclude_gene_types = active_config.exclude_gene_types
    span = active_config.span

    output_key = (
        f"highly_variable_{method}_{flavor}"
        if method == "scanpy"
        else f"highly_variable_{method}"
    )

    log.info(f"[HVG] Diagnosing input data from layer '{input_layer}' ...")
    X = adata.layers.get(input_layer, adata.X)
    stats = _diagnose_input_for_hvg(X)

    if output_key in adata.var and not force:
        n_existing = adata.var[output_key].sum()
        log.info(
            f"[HVG] Annotations found in '{output_key}' with {n_existing} genes. Use force=True to overwrite."
        )
        return adata

    with use_layer_as_X(adata, input_layer):
        if method == "scanpy":
            log.info(
                f"[HVG] Selecting with scanpy ({flavor}), n_top_genes={n_top_genes}, batch_key={batch_key}"
            )
            sc.pp.highly_variable_genes(
                adata,
                flavor=flavor,
                n_top_genes=n_top_genes,
                batch_key=batch_key,
                span=span,
                inplace=True,
            )
            adata.var[output_key] = adata.var["highly_variable"].copy()
            for metric in ["means", "dispersions", "dispersions_norm"]:
                if metric in adata.var:
                    adata.var[f"{output_key}_{metric}"] = adata.var[metric].copy()
            if output_key != "highly_variable":
                del adata.var["highly_variable"]

        elif method == "custom":
            log.info(f"[HVG] Custom HVG selection (n_top_genes={n_top_genes})")
            # ... (rest of the custom method logic remains the same)
            sample_key = active_config.sample_key
            min_n_samples = active_config.min_n_samples
            n_highly_expressed_genes = active_config.n_highly_expressed_genes
            n_specific_genes = active_config.n_specific_genes
            if sample_key not in adata.obs.columns:
                raise KeyError(f"Sample key '{sample_key}' not found in adata.obs")
            samples = adata.obs[sample_key].unique()
            hvg_masks = []
            for sample in samples:
                sample_mask = adata.obs[sample_key] == sample
                if sample_mask.sum() > 10:
                    sample_adata_view = adata[sample_mask, :]
                    sample_adata = sc.AnnData(
                        X=sample_adata_view.X, var=sample_adata_view.var
                    )
                    sc.pp.highly_variable_genes(
                        sample_adata, n_top_genes=n_top_genes, inplace=True
                    )
                    hvg_masks.append(sample_adata.var["highly_variable"])
            if not hvg_masks:
                raise ValueError("[HVG] No samples had enough cells to compute HVGs")
            combined_df = pd.concat(hvg_masks, axis=1)
            sample_counts = combined_df.sum(axis=1)
            combined_hvgs = sample_counts >= min_n_samples
            adata.var[f"{output_key}_sample_count"] = sample_counts
            gene_expr = np.array(adata.X.sum(axis=0)).flatten()
            top_expr_indices = np.argsort(-gene_expr)[:n_highly_expressed_genes]
            top_expr_genes = adata.var_names[top_expr_indices]
            specific_genes = _get_sample_specific_genes(
                adata, sample_key, n_specific_genes, layer=input_layer
            )
            exclude_genes = set(top_expr_genes) | set(specific_genes)
            final_hvg_mask = combined_hvgs & ~adata.var_names.isin(exclude_genes)
            adata.var[output_key] = final_hvg_mask
            adata.var[f"{output_key}_highly_expressed"] = adata.var_names.isin(
                top_expr_genes
            )
            adata.var[f"{output_key}_sample_specific"] = adata.var_names.isin(
                specific_genes
            )

        elif method == "triku":
            try:
                import triku
            except ImportError:
                raise ImportError("Please install triku: pip install triku")
            log.info("[HVG] Running triku method")
            result = triku.tl.triku(adata, return_all=True)
            adata.var[output_key] = result["highly_variable"]
            adata.var[f"{output_key}_score"] = result["score"]
        else:
            raise ValueError(f"Unknown method '{method}'.")

    # Exclude gene types
    gene_type_counts = {}
    if exclude_gene_types:
        log.info(f"[HVG] Excluding gene types: {exclude_gene_types}")
        current_mask = adata.var[output_key]
        updated_mask, excluded_counts = _exclude_genes(
            adata, current_mask, exclude_gene_types
        )
        adata.var[output_key] = updated_mask
        gene_type_counts.update(excluded_counts)

    n_hvg = int(adata.var[output_key].sum())
    log.info(f"[HVG] Final number of highly variable genes: {n_hvg}")

    # --- Store metadata in .uns ---
    # Use a helper function to ensure the dictionary is savable
    savable_params = _to_savable_dict(asdict(active_config))
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["hvg"] = {
        "output_key": output_key,
        "method": method,
        "params": savable_params,
        "n_hvg": n_hvg,
        "input_stats": stats,
        "excluded_gene_types": gene_type_counts,
    }

    if plot:
        plot_hvg_metrics(
            adata,
            output_key,
            save_path=Path(save_dir) / "hvg_metrics.png" if save_dir else None,
        )
        plt.show()

    if report and save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        report_path = Path(save_dir) / "hvg_report.md"
        _write_hvg_report(report_path, stats, n_hvg, active_config, gene_type_counts)
        log.info(f"[HVG] Report written to {report_path}")

    return adata


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
    msg.append(
        f"- Jaccard Similarity Index: {jaccard_index:.3f} (Intersection / Union)"
    )

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
        msg.append(
            "The methods agree on a core set of genes but also identify unique ones."
        )
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
            log.warning(
                "Difference mode needs at least 2 masks. Falling back to direct."
            )
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

    log.info(
        f"Created final HVG mask in '.var['{output_key}']' with {mask_combined.sum()} genes."
    )

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
            log.info(
                f"Subsetted AnnData object in-place to {mask_combined.sum()} final HVGs."
            )
            return adata

    return adata


def evaluate_hvg_stability(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_bootstrap: int = 20,
    sample_fraction: float = 0.8,
    method: str = "scanpy",
    flavor: str = "seurat",
    n_top_genes: Optional[int] = 2000,
    layer: Optional[str] = None,
    random_state: Optional[int] = 42,
    plot: bool = True,
    save_path: Optional[str] = None,
) -> AnnData:
    """
    Evaluates the stability of HVG selection through bootstrap resampling.
    Adds stability info to .uns['sclucid']['preprocess']['hvg_stability'].
    """
    import random

    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be a positive integer")
    if not 0 < sample_fraction < 1:
        raise ValueError("sample_fraction must be between 0 and 1")
    current_hvgs = set(adata.var_names[adata.var[hvg_key]])
    log.info(f"[HVG stability] Evaluating stability of {len(current_hvgs)} HVGs")
    if random_state is not None:
        random.seed(random_state)
        np.random.seed(random_state)
    gene_selection_count = {gene: 0 for gene in adata.var_names}
    n_cells_per_bootstrap = int(adata.n_obs * sample_fraction)
    report_interval = max(1, n_bootstrap // 10)
    for i in range(n_bootstrap):
        if i % report_interval == 0:
            log.info(f"[HVG stability] Bootstrap iteration {i + 1}/{n_bootstrap}")
        cell_indices = np.random.choice(
            adata.n_obs, size=n_cells_per_bootstrap, replace=False
        )
        bootstrap_adata_view = adata[cell_indices, :]
        bootstrap_adata = sc.AnnData(
            X=bootstrap_adata_view.X, var=bootstrap_adata_view.var
        )
        find_hvgs(
            bootstrap_adata,
            HVGConfig(method=method, n_top_genes=n_top_genes, flavor=flavor),
            force=True,
            plot=False,
        )
        bootstrap_hvgs = set(
            bootstrap_adata.var_names[
                bootstrap_adata.var[f"highly_variable_{method}_{flavor}"]
            ]
            if method == "scanpy"
            else bootstrap_adata.var_names[
                bootstrap_adata.var[f"highly_variable_{method}"]
            ]
        )
        for gene in bootstrap_hvgs:
            if gene in gene_selection_count:
                gene_selection_count[gene] += 1
    selection_frequency = {
        gene: count / n_bootstrap for gene, count in gene_selection_count.items()
    }
    adata.var["hvg_selection_frequency"] = pd.Series(
        [selection_frequency.get(gene, 0) for gene in adata.var_names],
        index=adata.var_names,
    )
    stability_score = np.mean(
        [selection_frequency.get(gene, 0) for gene in current_hvgs]
    )
    top_quartile = np.quantile(
        [selection_frequency.get(gene, 0) for gene in current_hvgs], 0.75
    )
    bottom_quartile = np.quantile(
        [selection_frequency.get(gene, 0) for gene in current_hvgs], 0.25
    )
    log.info("[HVG stability] Stability metrics:")
    log.info(f"  - Overall stability score: {stability_score:.3f}")
    log.info(f"  - Top 25% of HVGs selected with frequency >= {top_quartile:.3f}")
    log.info(f"  - Bottom 25% of HVGs selected with frequency <= {bottom_quartile:.3f}")
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "hvg_stability"
    ] = {
        "overall_score": stability_score,
        "top_quartile": top_quartile,
        "bottom_quartile": bottom_quartile,
        "n_bootstrap": n_bootstrap,
        "sample_fraction": sample_fraction,
        "method": method,
    }
    if plot:
        try:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sns.histplot(
                adata.var["hvg_selection_frequency"], bins=30, kde=True, ax=axes[0]
            )
            axes[0].set_title("HVG Selection Frequency Distribution")
            axes[0].set_xlabel("Selection Frequency")
            axes[0].set_ylabel("Number of Genes")
            axes[0].axvline(
                stability_score,
                color="red",
                linestyle="--",
                label=f"Avg HVG Stability: {stability_score:.3f}",
            )
            axes[0].legend()
            if "means" in adata.var:
                x = "means"
            elif f"{hvg_key}_means" in adata.var:
                x = f"{hvg_key}_means"
            else:
                if layer is None:
                    adata.var["temp_means"] = np.array(adata.X.mean(axis=0)).flatten()
                else:
                    adata.var["temp_means"] = np.array(
                        adata.layers[layer].mean(axis=0)
                    ).flatten()
                x = "temp_means"
            scatter = axes[1].scatter(
                adata.var[x],
                adata.var["hvg_selection_frequency"],
                c=adata.var[hvg_key].astype(int),
                alpha=0.6,
                cmap="coolwarm",
                s=10,
            )
            axes[1].set_title("HVG Stability vs. Mean Expression")
            axes[1].set_xlabel("Mean Expression")
            axes[1].set_ylabel("Selection Frequency")
            cbar = plt.colorbar(scatter, ax=axes[1])
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(["Not HVG", "HVG"])
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                log.info(f"Saved stability plot to {save_path}")
            if "temp_means" in adata.var:
                del adata.var["temp_means"]
        except Exception as e:
            log.warning(f"[HVG stability] Failed to create stability plot: {str(e)}")
    return adata


def plot_hvg_metrics(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_top_genes: int = 20,
    metrics: Optional[List[str]] = None,
    show_gene_labels: bool = True,
    size_by_expr: bool = True,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Creates visualizations of HVG metrics to evaluate selection quality.
    """
    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")
    available_metrics = {}
    for column in adata.var.columns:
        if any(x in column for x in ["mean", "disp", "var", "score"]):
            if "mean" in column:
                available_metrics["mean"] = column
            elif any(x in column for x in ["disp", "var"]):
                available_metrics["dispersion"] = column
            elif "norm" in column:
                available_metrics["norm_dispersion"] = column
            elif "score" in column:
                available_metrics["score"] = column
    if metrics is not None:
        for metric in metrics:
            if metric not in adata.var.columns:
                log.warning(f"Metric '{metric}' not found in adata.var columns")
    if "norm_dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["norm_dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "score" in available_metrics:
        x = available_metrics["mean"] if "mean" in available_metrics else None
        y = available_metrics["score"]
        plot_type = "score"
    else:
        method_specific_x = f"{hvg_key}_means"
        method_specific_y = f"{hvg_key}_dispersions_norm"
        if method_specific_x in adata.var and method_specific_y in adata.var:
            x = method_specific_x
            y = method_specific_y
            plot_type = "dispersion_vs_mean"
        else:
            log.warning("[HVG plot] Could not find appropriate metrics for plotting")
            adata.var["_temp_mean"] = np.array(adata.X.mean(axis=0)).flatten()
            x = "_temp_mean"
            if "hvg_selection_frequency" in adata.var:
                y = "hvg_selection_frequency"
                plot_type = "stability"
            else:
                raise ValueError(
                    "[HVG plot] Cannot create HVG plot: no appropriate metrics found"
                )
    fig, ax = plt.subplots(figsize=(10, 8))
    if size_by_expr and "mean" in available_metrics:
        sizes = np.clip(adata.var[available_metrics["mean"]] * 20, 5, 200)
    else:
        sizes = 30
    scatter = ax.scatter(
        adata.var[x],
        adata.var[y],
        s=sizes,
        c=adata.var[hvg_key].astype(int),
        cmap="coolwarm",
        alpha=0.7,
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Not HVG", "HVG"])
    if plot_type == "dispersion_vs_mean":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Dispersion (normalized)")
        ax.set_title("Highly Variable Genes: Dispersion vs. Mean")
        ax.set_xscale("log")
    elif plot_type == "score":
        ax.set_xlabel("Mean Expression" if x else "Gene Index")
        ax.set_ylabel("HVG Score")
        ax.set_title("Highly Variable Genes: Score Distribution")
    elif plot_type == "stability":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Selection Frequency")
        ax.set_title("HVG Selection Stability")
    if show_gene_labels and n_top_genes > 0:
        hvg_mask = adata.var[hvg_key]
        if hvg_mask.sum() > 0:
            top_indices = adata.var.loc[hvg_mask, y].nlargest(n_top_genes).index
            for idx in top_indices:
                gene_name = idx
                ax.annotate(
                    gene_name,
                    (adata.var.loc[idx, x], adata.var.loc[idx, y]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    bbox=dict(
                        boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8
                    ),
                )
    n_hvgs = adata.var[hvg_key].sum()
    total_genes = len(adata.var)
    ax.set_title(
        f"{ax.get_title()}\n{n_hvgs} HVGs selected ({n_hvgs / total_genes:.1%} of {total_genes} genes)"
    )
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"[HVG plot] Saved metrics plot to {save_path}")
    if "_temp_mean" in adata.var:
        del adata.var["_temp_mean"]
    return fig
