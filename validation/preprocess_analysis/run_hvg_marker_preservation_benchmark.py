#!/usr/bin/env python3
"""Benchmark HVG marker/program preservation on real datasets."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS
from validation.gene_panels import MARKER_PANELS, TUMOR_PROGRAM_PANELS, present_genes


def _matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _log_normalized_matrix(adata: ad.AnnData, target_sum: float = 1e4):
    X = _matrix(adata).astype(float)
    totals = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    if sp.issparse(X):
        Xn = X.multiply(scale[:, None]).tocsr(copy=False)
        Xn.data = np.log1p(Xn.data)
        return Xn
    return np.log1p(np.asarray(X) * scale[:, None])


def _gene_variance(X) -> np.ndarray:
    if sp.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        return np.maximum(mean_sq - mean * mean, 0.0)
    return np.asarray(X).var(axis=0)


def _ranked_hvgs(adata: ad.AnnData, n_top_genes: int) -> tuple[list[str], pd.Series]:
    X = _log_normalized_matrix(adata)
    variance = _gene_variance(X)
    names = pd.Index(adata.var_names.astype(str))
    ranks = pd.Series(variance, index=names).rank(method="first", ascending=False)
    selected_idx = np.argsort(-variance, kind="mergesort")[: min(n_top_genes, adata.n_vars)]
    return list(names[selected_idx]), ranks


def _protected_hvgs(
    standard_hvgs: list[str],
    ranks: pd.Series,
    protected_genes: set[str],
    n_top_genes: int,
) -> list[str]:
    selected = list(dict.fromkeys(standard_hvgs))
    selected_set = set(selected)
    protected_present = [gene for gene in protected_genes if gene in ranks.index]
    missing_protected = sorted(
        [gene for gene in protected_present if gene not in selected_set],
        key=lambda gene: float(ranks.loc[gene]),
    )
    for gene in missing_protected:
        if len(selected) < n_top_genes:
            selected.append(gene)
            selected_set.add(gene)
            continue
        replace_idx = None
        for idx in range(len(selected) - 1, -1, -1):
            candidate = selected[idx]
            if candidate not in protected_genes:
                replace_idx = idx
                break
        if replace_idx is None:
            break
        selected_set.remove(selected[replace_idx])
        selected[replace_idx] = gene
        selected_set.add(gene)
    return selected


def _panel_gene_set(var_names, panels: dict[str, tuple[str, ...]]) -> set[str]:
    genes: set[str] = set()
    for panel_genes in panels.values():
        genes.update(present_genes(var_names, panel_genes))
    return genes


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def _overlap_row(
    dataset: str,
    left_name: str,
    left: set[str],
    right_name: str,
    right: set[str],
) -> dict[str, Any]:
    union = left | right
    intersection = left & right
    return {
        "dataset": dataset,
        "left_set": left_name,
        "right_set": right_name,
        "n_left": len(left),
        "n_right": len(right),
        "n_intersection": len(intersection),
        "n_union": len(union),
        "jaccard": len(intersection) / len(union) if union else 0.0,
        "left_only": len(left - right),
        "right_only": len(right - left),
    }


def _suggest_mode_from_jaccard(jaccard: float) -> dict[str, Any]:
    """Mirror current scLucid overlap-only HVG auto guidance for evidence review."""
    if jaccard > 0.7:
        return {"recommended_mode": "union", "overlap_level": "high", "risk": "low"}
    if jaccard >= 0.4:
        return {"recommended_mode": "intersection", "overlap_level": "moderate", "risk": "moderate"}
    return {"recommended_mode": "intersection", "overlap_level": "low", "risk": "high"}


def _mean_inclusion(rows: list[dict[str, Any]], panel_type: str) -> float:
    values = [row["inclusion_rate"] for row in rows if row["panel_type"] == panel_type]
    return float(np.mean(values)) if values else np.nan


def _panel_rows(
    *,
    dataset: str,
    strategy: str,
    hvg_set: set[str],
    panel_type: str,
    panels: dict[str, tuple[str, ...]],
    var_names,
    n_top_genes: int,
    selection_mode: str,
    budget_preserving: bool,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for panel, genes in panels.items():
        present = present_genes(var_names, genes)
        retained = [gene for gene in present if gene in hvg_set]
        rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "selection_mode": selection_mode,
                "panel_type": panel_type,
                "panel": panel,
                "n_top_genes": n_top_genes,
                "hvg_set_size": len(hvg_set),
                "budget_preserving": budget_preserving,
                "genes_expected": len(genes),
                "genes_present": len(present),
                "genes_retained": len(retained),
                "inclusion_rate": len(retained) / max(len(present), 1),
                "present_genes": ";".join(present),
                "retained_genes": ";".join(retained),
                "review_required": bool(present and len(retained) / max(len(present), 1) < 0.5),
            }
        )
    return rows


def run(
    output_dir: Path, datasets: set[str] | None, max_cells: int | None, seed: int, n_top_genes: int
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    marker_rows: list[dict[str, Any]] = []
    program_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []
    strategy_rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []

    protected_genes = {
        gene
        for panel in (*MARKER_PANELS.values(), *TUMOR_PROGRAM_PANELS.values())
        for gene in panel
    }

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.preprocess_roles or not spec.path.exists():
            continue
        adata = _subset(ad.read_h5ad(spec.path), max_cells=max_cells, seed=seed)
        standard_hvgs, ranks = _ranked_hvgs(adata, n_top_genes=n_top_genes)
        standard_set = set(standard_hvgs)
        custom_marker_set = _panel_gene_set(adata.var_names, MARKER_PANELS)
        custom_program_set = _panel_gene_set(adata.var_names, TUMOR_PROGRAM_PANELS)
        custom_marker_program_set = custom_marker_set | custom_program_set
        protected_hvgs = _protected_hvgs(
            standard_hvgs,
            ranks,
            protected_genes=protected_genes,
            n_top_genes=min(n_top_genes, adata.n_vars),
        )
        union_set = standard_set | custom_marker_program_set
        intersection_set = standard_set & custom_marker_program_set
        auto_jaccard = _jaccard(standard_set, custom_marker_program_set)
        auto_guidance = _suggest_mode_from_jaccard(auto_jaccard)
        auto_set = union_set if auto_guidance["recommended_mode"] == "union" else intersection_set
        semantic_auto_set = union_set

        set_library = {
            "standard_variance_hvg": standard_set,
            "custom_marker_hvg": custom_marker_set,
            "custom_program_hvg": custom_program_set,
            "custom_marker_program_hvg": custom_marker_program_set,
            "direct_standard": standard_set,
            "union_standard_custom": union_set,
            "intersection_standard_custom": intersection_set,
            "auto_recommended_standard_custom": auto_set,
            "semantic_auto_protected_union": semantic_auto_set,
            "sclucid_marker_program_retained": set(protected_hvgs),
        }

        strategy_metadata = {
            "standard_variance_hvg": {
                "selection_mode": "standard_variance",
                "budget_preserving": True,
                "note": "Variance-ranked baseline.",
            },
            "custom_marker_hvg": {
                "selection_mode": "custom_marker_only",
                "budget_preserving": len(custom_marker_set) <= n_top_genes,
                "note": "Curated lineage marker panel only; diagnostic, not a standalone HVG recommendation.",
            },
            "custom_program_hvg": {
                "selection_mode": "custom_program_only",
                "budget_preserving": len(custom_program_set) <= n_top_genes,
                "note": "Curated program panel only; diagnostic, not a standalone HVG recommendation.",
            },
            "custom_marker_program_hvg": {
                "selection_mode": "custom_marker_program_only",
                "budget_preserving": len(custom_marker_program_set) <= n_top_genes,
                "note": "Combined curated marker and program panel.",
            },
            "direct_standard": {
                "selection_mode": "direct",
                "budget_preserving": True,
                "note": "Direct use of the first standard HVG mask.",
            },
            "union_standard_custom": {
                "selection_mode": "union",
                "budget_preserving": len(union_set) <= n_top_genes,
                "note": "Maximum biology retention with extra gene budget when custom genes are outside standard HVGs.",
            },
            "intersection_standard_custom": {
                "selection_mode": "intersection",
                "budget_preserving": True,
                "note": "Conservative overlap; useful as a negative-control risk example for protected marker sets.",
            },
            "auto_recommended_standard_custom": {
                "selection_mode": f"auto->{auto_guidance['recommended_mode']}",
                "budget_preserving": auto_guidance["recommended_mode"] == "intersection"
                or len(auto_set) <= n_top_genes,
                "note": (
                    "Mirrors current overlap-only auto guidance. High-risk rows indicate where semantic-aware "
                    "protected-marker guidance should override raw overlap."
                ),
            },
            "semantic_auto_protected_union": {
                "selection_mode": "auto(protected)->union",
                "budget_preserving": len(semantic_auto_set) <= n_top_genes,
                "note": (
                    "Protected-aware auto mode: favor union when a curated marker/program set is explicitly "
                    "identified as biology that must not be lost."
                ),
            },
            "sclucid_marker_program_retained": {
                "selection_mode": "budget_preserving_retention",
                "budget_preserving": len(protected_hvgs) <= min(n_top_genes, adata.n_vars),
                "note": "scLucid protected-gene policy: rescue curated marker/program genes while preserving HVG budget.",
            },
        }

        overlap_rows.extend(
            [
                _overlap_row(
                    spec.key,
                    "standard_variance_hvg",
                    standard_set,
                    "custom_marker_hvg",
                    custom_marker_set,
                ),
                _overlap_row(
                    spec.key,
                    "standard_variance_hvg",
                    standard_set,
                    "custom_program_hvg",
                    custom_program_set,
                ),
                _overlap_row(
                    spec.key,
                    "standard_variance_hvg",
                    standard_set,
                    "custom_marker_program_hvg",
                    custom_marker_program_set,
                ),
            ]
        )

        dataset_panel_rows: list[dict[str, Any]] = []
        for strategy, hvg_set in set_library.items():
            metadata = strategy_metadata[strategy]
            marker_panel_rows = _panel_rows(
                dataset=spec.key,
                strategy=strategy,
                hvg_set=hvg_set,
                panel_type="lineage_marker",
                panels=MARKER_PANELS,
                var_names=adata.var_names,
                n_top_genes=n_top_genes,
                selection_mode=metadata["selection_mode"],
                budget_preserving=metadata["budget_preserving"],
            )
            program_panel_rows = _panel_rows(
                dataset=spec.key,
                strategy=strategy,
                hvg_set=hvg_set,
                panel_type="tumor_program",
                panels=TUMOR_PROGRAM_PANELS,
                var_names=adata.var_names,
                n_top_genes=n_top_genes,
                selection_mode=metadata["selection_mode"],
                budget_preserving=metadata["budget_preserving"],
            )
            marker_rows.extend(marker_panel_rows)
            program_rows.extend(program_panel_rows)
            dataset_panel_rows.extend(marker_panel_rows)
            dataset_panel_rows.extend(program_panel_rows)
            strategy_rows.append(
                {
                    "dataset": spec.key,
                    "strategy": strategy,
                    "selection_mode": metadata["selection_mode"],
                    "hvg_set_size": len(hvg_set),
                    "n_top_genes": n_top_genes,
                    "budget_preserving": metadata["budget_preserving"],
                    "jaccard_with_standard": _jaccard(standard_set, hvg_set),
                    "auto_overlap_level": (
                        auto_guidance["overlap_level"]
                        if strategy == "auto_recommended_standard_custom"
                        else ""
                    ),
                    "auto_risk": (
                        auto_guidance["risk"]
                        if strategy == "auto_recommended_standard_custom"
                        else ""
                    ),
                    "mean_marker_inclusion_rate": _mean_inclusion(
                        marker_panel_rows, "lineage_marker"
                    ),
                    "mean_program_inclusion_rate": _mean_inclusion(
                        program_panel_rows, "tumor_program"
                    ),
                    "review_required": bool(
                        any(
                            row["review_required"] for row in marker_panel_rows + program_panel_rows
                        )
                    ),
                    "note": metadata["note"],
                }
            )

        figure_strategies = {
            "standard_variance_hvg",
            "union_standard_custom",
            "intersection_standard_custom",
            "auto_recommended_standard_custom",
            "semantic_auto_protected_union",
            "sclucid_marker_program_retained",
        }
        for row in dataset_panel_rows:
            if row["strategy"] not in figure_strategies:
                continue
            figure_rows.append(
                {
                    "figure_panel": "3B",
                    "dataset": spec.key,
                    "panel_type": row["panel_type"],
                    "panel": row["panel"],
                    "strategy": row["strategy"],
                    "selection_mode": row["selection_mode"],
                    "inclusion_rate": row["inclusion_rate"],
                    "genes_present": row["genes_present"],
                    "genes_retained": row["genes_retained"],
                    "hvg_set_size": row["hvg_set_size"],
                    "budget_preserving": row["budget_preserving"],
                    "review_required": row["review_required"],
                    "context": json.dumps(
                        {
                            "n_top_genes": n_top_genes,
                            "max_cells": max_cells,
                            "auto_recommended_mode": auto_guidance["recommended_mode"],
                            "auto_jaccard_standard_custom": auto_jaccard,
                            "auto_risk": auto_guidance["risk"],
                        }
                    ),
                }
            )

    paths = {
        "marker_preservation": output_dir / "hvg_marker_preservation.tsv",
        "program_retention": output_dir / "program_gene_retention.tsv",
        "strategy_summary": output_dir / "hvg_strategy_summary.tsv",
        "set_overlap": output_dir / "hvg_set_overlap.tsv",
        "figure3": output_dir / "figure3_hvg_data.tsv",
    }
    pd.DataFrame(marker_rows).to_csv(paths["marker_preservation"], sep="\t", index=False)
    pd.DataFrame(program_rows).to_csv(paths["program_retention"], sep="\t", index=False)
    pd.DataFrame(strategy_rows).to_csv(paths["strategy_summary"], sep="\t", index=False)
    pd.DataFrame(overlap_rows).to_csv(paths["set_overlap"], sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(paths["figure3"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/preprocess_hvg_preservation")
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys to include.")
    parser.add_argument(
        "--max-cells", type=int, default=5000, help="Pilot subset size. Use 0 for full data."
    )
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=args.max_cells,
        seed=args.seed,
        n_top_genes=args.n_top_genes,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
