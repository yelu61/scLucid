#!/usr/bin/env python3
"""Generate deterministic QC positive/negative controls with explicit truth.

This suite is intentionally small and mechanistic. It can validate input,
execution, determinism, and expected-direction behavior. It cannot establish
external scientific superiority on real tissue.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

SCHEMA_VERSION = "sclucid_controlled_qc_truth_v1"
GENERATOR_VERSION = "1.0.0"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _gene_names(n_genes: int) -> list[str]:
    if n_genes < 80:
        raise ValueError("n_genes must be at least 80.")
    special = [
        "MT-CO1",
        "MT-CO2",
        "MT-CO3",
        "MT-ND1",
        "MT-ND2",
        "MT-ND3",
        "MT-CYB",
        "MT-ATP6",
        "A_MARKER1",
        "A_MARKER2",
        "A_MARKER3",
        "B_MARKER1",
        "B_MARKER2",
        "B_MARKER3",
        "TUMOR_MARKER1",
        "TUMOR_MARKER2",
        "TUMOR_MARKER3",
        "STRESS1",
        "STRESS2",
        "APOPTOSIS1",
    ]
    return special + [f"GENE_{idx:04d}" for idx in range(n_genes - len(special))]


def _base_rate(n_genes: int, lineage: str, assay_profile: str) -> np.ndarray:
    # Keep the intact-cell background deep enough that the default 240-gene
    # fixture genuinely clears the reviewer's protocol safety floor.  This
    # prevents the locked-HQ check from passing merely because low complexity
    # triggers only one of several catastrophic axes.
    rate = np.full(n_genes, 4.0, dtype=float)
    rate[:8] = 0.12 if assay_profile == "scrna" else 0.025
    if lineage == "lineage_A":
        rate[8:11] = 3.0
    elif lineage == "lineage_B":
        rate[11:14] = 3.0
    elif lineage == "rare_tumor":
        rate[14:17] = 4.0
        rate[8:11] = 0.8
    return rate


def _draw_native(
    rng: np.random.Generator,
    n_genes: int,
    lineage: str,
    assay_profile: str,
    *,
    depth_scale: float,
) -> np.ndarray:
    heterogeneity = rng.gamma(shape=2.0, scale=0.5, size=n_genes)
    rate = _base_rate(n_genes, lineage, assay_profile) * heterogeneity * depth_scale
    return rng.poisson(rate).astype(np.int32)


def build_controlled_truth(
    *,
    seed: int = 20260820,
    n_genes: int = 240,
) -> AnnData:
    """Return a deterministic AnnData containing known QC mechanisms."""
    rng = np.random.default_rng(seed)
    genes = _gene_names(n_genes)
    ambient_profile = np.full(n_genes, 0.02, dtype=float)
    ambient_profile[8:17] = 0.25
    ambient_profile[:8] = 0.04

    library_specs = [
        ("hq_scrna", "scrna", {"lineage_A": 100, "lineage_B": 80, "rare_tumor": 12}, 25, 60, 18),
        (
            "degraded_scrna",
            "scrna",
            {"lineage_A": 50, "lineage_B": 40, "rare_tumor": 8},
            70,
            70,
            20,
        ),
        ("hq_snrna", "snrna", {"lineage_A": 70, "lineage_B": 60, "rare_tumor": 10}, 20, 50, 15),
    ]

    native_rows: list[np.ndarray] = []
    ambient_rows: list[np.ndarray] = []
    obs_rows: list[dict[str, object]] = []
    cell_index = 0

    def add_row(
        *,
        library: str,
        assay_profile: str,
        droplet_class: str,
        lineage: str,
        native: np.ndarray,
        ambient: np.ndarray,
        low_rna: bool,
        damage_fraction: float,
        source_pair: str = "",
    ) -> None:
        nonlocal cell_index
        total = int(native.sum() + ambient.sum())
        realized_ambient = float(ambient.sum() / total) if total else 0.0
        is_rare = lineage == "rare_tumor"
        if droplet_class == "intact":
            policy_label = "KEEP"
        elif droplet_class in {"damaged", "empty"}:
            policy_label = "REMOVE"
        else:
            policy_label = "UNCERTAIN"
        obs_rows.append(
            {
                "obs_name": f"truth_{cell_index:05d}",
                "library": library,
                "sample": library,
                "assay_profile": assay_profile,
                "truth_droplet_class": droplet_class,
                "truth_policy_label": policy_label,
                "truth_lineage": lineage,
                "truth_is_rare": is_rare,
                "truth_low_rna": low_rna,
                "truth_protected_keep": bool(droplet_class == "intact" and (low_rna or is_rare)),
                "truth_is_cell_capture": droplet_class != "empty",
                "truth_is_doublet": droplet_class == "doublet",
                "truth_has_ambient": bool(ambient.sum() > 0),
                "truth_ambient_fraction": realized_ambient,
                "truth_damage_fraction": damage_fraction,
                "truth_source_pair": source_pair,
            }
        )
        native_rows.append(native)
        ambient_rows.append(ambient)
        cell_index += 1

    for library, assay_profile, intact_counts, n_damaged, n_empty, n_doublet in library_specs:
        intact_pool: list[tuple[str, np.ndarray]] = []
        for lineage, n_cells in intact_counts.items():
            for index in range(n_cells):
                low_rna = index < max(2, int(round(0.12 * n_cells)))
                depth = (0.38 if low_rna else 1.0) * rng.lognormal(mean=0.0, sigma=0.18)
                native = _draw_native(
                    rng,
                    n_genes,
                    lineage,
                    assay_profile,
                    depth_scale=depth,
                )
                ambient_depth = 0.8 if index % 9 else 4.0
                ambient = rng.poisson(ambient_profile * ambient_depth).astype(np.int32)
                add_row(
                    library=library,
                    assay_profile=assay_profile,
                    droplet_class="intact",
                    lineage=lineage,
                    native=native,
                    ambient=ambient,
                    low_rna=low_rna,
                    damage_fraction=0.0,
                )
                intact_pool.append((lineage, native.copy()))

        for index in range(n_damaged):
            lineage = ("lineage_A", "lineage_B", "rare_tumor")[index % 3]
            source = _draw_native(rng, n_genes, lineage, assay_profile, depth_scale=1.0)
            retain_probability = 0.18 if library == "degraded_scrna" else 0.3
            native = rng.binomial(source, retain_probability).astype(np.int32)
            native[:8] += rng.poisson(2.5, size=8).astype(np.int32)
            native[17:20] += rng.poisson(2.0, size=3).astype(np.int32)
            ambient = rng.poisson(ambient_profile * 3.0).astype(np.int32)
            add_row(
                library=library,
                assay_profile=assay_profile,
                droplet_class="damaged",
                lineage=lineage,
                native=native,
                ambient=ambient,
                low_rna=False,
                damage_fraction=1.0 - retain_probability,
            )

        for _ in range(n_empty):
            ambient = rng.poisson(ambient_profile * 1.5).astype(np.int32)
            add_row(
                library=library,
                assay_profile=assay_profile,
                droplet_class="empty",
                lineage="none",
                native=np.zeros(n_genes, dtype=np.int32),
                ambient=ambient,
                low_rna=False,
                damage_fraction=0.0,
            )

        for _ in range(n_doublet):
            first, second = rng.choice(len(intact_pool), size=2, replace=False)
            first_lineage, first_counts = intact_pool[int(first)]
            second_lineage, second_counts = intact_pool[int(second)]
            ambient = rng.poisson(ambient_profile * 1.5).astype(np.int32)
            add_row(
                library=library,
                assay_profile=assay_profile,
                droplet_class="doublet",
                lineage="mixed",
                native=first_counts + second_counts,
                ambient=ambient,
                low_rna=False,
                damage_fraction=0.0,
                source_pair=f"{first_lineage}+{second_lineage}",
            )

    native_matrix = sparse.csr_matrix(np.vstack(native_rows), dtype=np.int32)
    ambient_matrix = sparse.csr_matrix(np.vstack(ambient_rows), dtype=np.int32)
    counts = (native_matrix + ambient_matrix).tocsr()
    counts.eliminate_zeros()
    obs = pd.DataFrame(obs_rows).set_index("obs_name")
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    var["is_mitochondrial"] = var.index.str.startswith("MT-")
    var["truth_feature_group"] = "background"
    var.loc[var.index.str.contains("_MARKER"), "truth_feature_group"] = "lineage_marker"
    var.loc[var.index.str.startswith("MT-"), "truth_feature_group"] = "mitochondrial"
    var.loc[var.index.str.startswith(("STRESS", "APOPTOSIS")), "truth_feature_group"] = (
        "damage_program"
    )

    adata = AnnData(X=counts.copy(), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.layers["native_counts"] = native_matrix
    adata.layers["ambient_counts"] = ambient_matrix
    adata.uns["sclucid_controlled_truth"] = {
        "schema_version": SCHEMA_VERSION,
        "generator_version": GENERATOR_VERSION,
        "seed": seed,
        "truth_source": "deterministic_simulation",
        "counts_identity": "layers[counts] = layers[native_counts] + layers[ambient_counts]",
        "truth_semantics": {
            "droplet_class": "Known simulated mechanism; independent of the expected policy action.",
            "policy_label": (
                "KEEP for intact cells, REMOVE for damaged/empty droplets, and UNCERTAIN "
                "for doublets because the product contract keeps doublet evidence review-only."
            ),
            "protected_keep": "Intact low-RNA or rare-tumor cells that QC should not erase.",
        },
        "primary_claim_boundary": (
            "Engineering and mechanism control only; external scientific superiority is not supported."
        ),
    }
    return adata


def write_suite(adata: AnnData, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = output_dir / "controlled_qc_truth.h5ad"
    truth_path = output_dir / "controlled_qc_truth_cells.tsv"
    manifest_path = output_dir / "controlled_qc_truth_manifest.json"
    adata.write_h5ad(h5ad_path)
    adata.obs.to_csv(truth_path, sep="\t", index=True)
    class_counts = (
        adata.obs.groupby(["library", "truth_droplet_class"], observed=True)
        .size()
        .rename("n")
        .reset_index()
        .to_dict(orient="records")
    )
    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS",
        "evidence_scope": "engineering_fixture_generation",
        "scientific_status": "SIMULATION_PASS_NOT_EXTERNAL",
        "generator_version": GENERATOR_VERSION,
        "seed": int(adata.uns["sclucid_controlled_truth"]["seed"]),
        "n_observations": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "class_counts": class_counts,
        "artifacts": {
            "h5ad": str(h5ad_path),
            "h5ad_sha256": _sha256(h5ad_path),
            "truth_cells": str(truth_path),
            "truth_cells_sha256": _sha256(truth_path),
        },
        "supports": [
            "Input/layer identity and deterministic policy execution tests.",
            "Known-direction tests for empty, ambient, damage, doublet, and rare-cell mechanisms.",
        ],
        "does_not_support": [
            "Real-tissue QC superiority.",
            "Tumor-project generalization.",
            "A release PASS for scientific QC heads without external truth.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/qc_controlled_truth"),
    )
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--n-genes", type=int, default=240)
    args = parser.parse_args()
    adata = build_controlled_truth(seed=args.seed, n_genes=args.n_genes)
    manifest = write_suite(adata, args.output_dir)
    print(json.dumps({"status": manifest["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
