#!/usr/bin/env python3
"""Build real-data evidence for Analysis inference-design contracts.

Kang2018 PBMC provides a paired donor-level ctrl/stim design for executable
proportion and pseudobulk checks. Lin2020 PDAC intentionally provides a
single-condition counterexample: the runner records a BLOCKED condition-level
inference status rather than inventing a comparison or cell-type annotation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from scLucid.analysis.config import PseudobulkDEConfig
from scLucid.analysis.differential_expression.de_core import run_pseudobulk_de
from scLucid.analysis.proportion.config import ProportionConfig
from scLucid.analysis.proportion.pseudobulk import celltype_proportion_analysis

PBMC_FILENAME = "kang2018.pbmc.h5ad"
PDAC_FILENAME = "lin2020.pdac.h5ad"
DERIVED_SAMPLE_KEY = "__sclucid_analysis_sample"


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_dict()
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_provenance(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "sha256": _sha256(path),
    }


def _valid_label_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.strip()
    return values.notna() & ~text.isin({"", "nan", "None", "NA", "<NA>"})


def build_design_snapshot(
    adata: ad.AnnData,
    *,
    dataset: str,
    sample_key: str,
    condition_key: str,
    experimental_unit_key: str,
    cell_type_key: str,
    paired_key: str | None = None,
    batch_key: str | None = None,
) -> dict[str, Any]:
    """Return a compact, claim-oriented audit of one real study design."""
    configured = [sample_key, condition_key, experimental_unit_key, cell_type_key]
    configured.extend(key for key in (paired_key, batch_key) if key)
    missing_columns = [key for key in dict.fromkeys(configured) if key not in adata.obs]
    blockers: list[str] = []
    if missing_columns:
        blockers.append(f"missing metadata columns: {missing_columns}")
        return {
            "dataset": dataset,
            "status": "BLOCKED",
            "missing_columns": missing_columns,
            "blockers": blockers,
        }

    obs = adata.obs
    sample_condition_conflicts = int(
        (
            obs.groupby(sample_key, observed=True)[condition_key].nunique(dropna=False)
            > 1
        ).sum()
    )
    sample_unit_conflicts = int(
        (
            obs.groupby(sample_key, observed=True)[experimental_unit_key].nunique(dropna=False)
            > 1
        ).sum()
    )
    design_columns = list(
        dict.fromkeys([sample_key, condition_key, experimental_unit_key])
    )
    design_rows = obs[design_columns].drop_duplicates()
    conditions = sorted(design_rows[condition_key].dropna().astype(str).unique().tolist())
    units_per_condition = {
        str(condition): int(frame[experimental_unit_key].nunique())
        for condition, frame in design_rows.groupby(condition_key, observed=True)
    }
    cell_type_mask = _valid_label_mask(obs[cell_type_key])
    usable_cell_types = sorted(obs.loc[cell_type_mask, cell_type_key].astype(str).unique().tolist())
    complete_pairs = 0
    if paired_key:
        pair_counts = obs[[paired_key, condition_key]].drop_duplicates().groupby(
            paired_key, observed=True
        )[condition_key].nunique()
        complete_pairs = int((pair_counts == len(conditions)).sum())

    if sample_condition_conflicts:
        blockers.append(
            f"{sample_condition_conflicts} sample values map to multiple conditions"
        )
    if sample_unit_conflicts:
        blockers.append(
            f"{sample_unit_conflicts} sample values map to multiple experimental units"
        )
    if len(conditions) < 2:
        blockers.append("fewer than two observed condition levels")
    if len(usable_cell_types) < 2:
        blockers.append("fewer than two usable non-empty cell-type labels")
    if units_per_condition and min(units_per_condition.values()) < 2:
        blockers.append("fewer than two independent experimental units in a condition")

    return {
        "dataset": dataset,
        "status": "BLOCKED" if blockers else "READY",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "sample_key": sample_key,
        "condition_key": condition_key,
        "experimental_unit_key": experimental_unit_key,
        "paired_key": paired_key,
        "batch_key": batch_key,
        "cell_type_key": cell_type_key,
        "n_samples": int(obs[sample_key].nunique()),
        "conditions": conditions,
        "experimental_units_per_condition": units_per_condition,
        "n_complete_pairs": complete_pairs,
        "n_usable_cell_types": len(usable_cell_types),
        "usable_cell_types": usable_cell_types,
        "n_cells_without_usable_cell_type": int((~cell_type_mask).sum()),
        "sample_condition_conflicts": sample_condition_conflicts,
        "sample_unit_conflicts": sample_unit_conflicts,
        "blockers": blockers,
    }


def _prepare_paired_pbmc(adata: ad.AnnData) -> tuple[ad.AnnData, list[str]]:
    required = {"donor", "condition", "cell_type"}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Kang2018 PBMC is missing required metadata: {missing}")
    obs = adata.obs
    valid = (
        _valid_label_mask(obs["donor"])
        & _valid_label_mask(obs["cell_type"])
        & obs["condition"].astype(str).isin(["ctrl", "stim"])
    )
    candidate = adata[valid].copy()
    donor_condition_counts = candidate.obs[["donor", "condition"]].drop_duplicates().groupby(
        "donor", observed=True
    )["condition"].nunique()
    complete_donors = sorted(
        donor_condition_counts[donor_condition_counts == 2].index.astype(str).tolist()
    )
    donor_text = candidate.obs["donor"].astype(str)
    candidate = candidate[donor_text.isin(complete_donors)].copy()
    candidate.obs["donor"] = candidate.obs["donor"].astype(str)
    candidate.obs["condition"] = candidate.obs["condition"].astype(str)
    candidate.obs["cell_type"] = candidate.obs["cell_type"].astype(str)
    candidate.obs[DERIVED_SAMPLE_KEY] = (
        candidate.obs["donor"] + "::" + candidate.obs["condition"]
    )
    return candidate, complete_donors


def _select_genes(adata: ad.AnnData, max_genes: int) -> list[str]:
    preferred = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1"]
    selected = [gene for gene in preferred if gene in adata.var_names]
    matrix = adata.layers.get("counts", adata.X)
    totals = np.asarray(matrix.sum(axis=0)).ravel()
    for idx in np.argsort(-totals):
        gene = str(adata.var_names[int(idx)])
        if gene not in selected:
            selected.append(gene)
        if len(selected) >= max_genes:
            break
    return selected[:max_genes]


def _eligible_cell_types(
    adata: ad.AnnData,
    *,
    min_cells_per_sample: int,
    max_cell_types: int,
) -> list[str]:
    coverage = adata.obs.groupby(
        ["cell_type", DERIVED_SAMPLE_KEY], observed=True
    ).size().unstack(fill_value=0)
    eligible = coverage.index[(coverage >= min_cells_per_sample).all(axis=1)].astype(str).tolist()
    abundance = adata.obs["cell_type"].astype(str).value_counts()
    eligible.sort(key=lambda label: (-int(abundance.get(label, 0)), label))
    return eligible[:max_cell_types]


def _matrix_rows(pbmc_snapshot: dict[str, Any], pdac_snapshot: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "dataset": "kang2018.pbmc",
            "context_field": "sample_key",
            "resolved_column": DERIVED_SAMPLE_KEY,
            "proportion_field": "sample_col",
            "pseudobulk_field": "sample_col",
            "statistical_role": "one aggregate observation per donor-condition",
            "source": "derived from donor + condition; original sample is a capture/batch label",
            "status": pbmc_snapshot["status"],
        },
        {
            "dataset": "kang2018.pbmc",
            "context_field": "condition_key",
            "resolved_column": "condition",
            "proportion_field": "condition_col",
            "pseudobulk_field": "condition_key",
            "statistical_role": "directional ctrl-to-stim contrast",
            "source": "tracked real-data metadata",
            "status": pbmc_snapshot["status"],
        },
        {
            "dataset": "kang2018.pbmc",
            "context_field": "experimental_unit_key / paired_key",
            "resolved_column": "donor",
            "proportion_field": "experimental_unit_col / pairing_col",
            "pseudobulk_field": "experimental_unit_col / block_col",
            "statistical_role": "independent donor and paired repeated-measures block",
            "source": "tracked real-data metadata",
            "status": pbmc_snapshot["status"],
        },
        {
            "dataset": "kang2018.pbmc",
            "context_field": "batch_key",
            "resolved_column": "batch_group",
            "proportion_field": "batch_col candidate",
            "pseudobulk_field": "design_covariates candidate",
            "statistical_role": "review candidate, not auto-applied",
            "source": "constant/aliased after paired-donor restriction; explicit review required",
            "status": "REVIEW",
        },
        {
            "dataset": "lin2020.pdac",
            "context_field": "sample_key / experimental_unit_key",
            "resolved_column": pdac_snapshot.get("sample_key"),
            "proportion_field": "sample_col / experimental_unit_col",
            "pseudobulk_field": "sample_col / experimental_unit_col",
            "statistical_role": "independent tumor specimens",
            "source": "tracked real-data metadata",
            "status": pdac_snapshot["status"],
        },
        {
            "dataset": "lin2020.pdac",
            "context_field": "condition_key",
            "resolved_column": "condition",
            "proportion_field": "condition_col",
            "pseudobulk_field": "condition_key",
            "statistical_role": "comparison factor",
            "source": "one observed level (Primary tumor); no contrast is identifiable",
            "status": pdac_snapshot["status"],
        },
    ]


def run(
    data_dir: Path,
    output_dir: Path,
    *,
    max_genes: int = 200,
    max_cell_types: int = 4,
    min_cells_per_sample: int = 10,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    pbmc_path = (data_dir / PBMC_FILENAME).resolve()
    pdac_path = (data_dir / PDAC_FILENAME).resolve()
    for path in (pbmc_path, pdac_path):
        if not path.exists():
            raise FileNotFoundError(path)

    pbmc, complete_donors = _prepare_paired_pbmc(ad.read_h5ad(pbmc_path))
    pbmc_snapshot = build_design_snapshot(
        pbmc,
        dataset="kang2018.pbmc",
        sample_key=DERIVED_SAMPLE_KEY,
        condition_key="condition",
        experimental_unit_key="donor",
        cell_type_key="cell_type",
        paired_key="donor",
        batch_key="batch_group" if "batch_group" in pbmc.obs else None,
    )
    if pbmc_snapshot["status"] != "READY":
        raise ValueError(f"Kang2018 PBMC design is not executable: {pbmc_snapshot['blockers']}")

    proportion_config = ProportionConfig(
        celltype_col="cell_type",
        sample_col=DERIVED_SAMPLE_KEY,
        condition_col="condition",
        experimental_unit_col="donor",
        pairing_col="donor",
        test_method="clr-paired-t-test",
        auto_configure=False,
        min_samples_per_condition=2,
        plot_types=[],
        export_data=False,
    )
    prop_df, prop_stats = celltype_proportion_analysis(pbmc, proportion_config)
    prop_df.index.name = DERIVED_SAMPLE_KEY
    prop_path = output_dir / "pbmc_proportion_estimates.tsv"
    prop_stats_path = output_dir / "pbmc_proportion_statistics.tsv"
    prop_df.to_csv(prop_path, sep="\t", index=True)
    prop_stats.to_csv(prop_stats_path, sep="\t", index=False)

    genes = _select_genes(pbmc, max_genes=max_genes)
    cell_types = _eligible_cell_types(
        pbmc,
        min_cells_per_sample=min_cells_per_sample,
        max_cell_types=max_cell_types,
    )
    if not cell_types:
        raise ValueError("No PBMC cell type has adequate cells in every donor-condition sample")
    pbmc_de = pbmc[:, genes].copy()
    de_config = PseudobulkDEConfig(
        sample_col=DERIVED_SAMPLE_KEY,
        condition_key="condition",
        experimental_unit_col="donor",
        block_col="donor",
        contrasts=[("ctrl", "stim")],
        groupby="cell_type",
        group_names=cell_types,
        layer="counts" if "counts" in pbmc_de.layers else None,
        min_cells_per_sample=min_cells_per_sample,
        min_counts=0,
        min_samples_per_condition=2,
        method="linear_model_logcpm",
        key_added="real_pbmc_pseudobulk_de",
    )
    de_results = run_pseudobulk_de(pbmc_de, de_config)
    de_path = output_dir / "pbmc_pseudobulk_de.tsv"
    de_results.to_csv(de_path, sep="\t", index=False)
    de_design = pbmc_de.uns["sclucid"]["analysis"]["de"][
        "real_pbmc_pseudobulk_de_design"
    ]

    pdac = ad.read_h5ad(pdac_path, backed="r")
    pdac_sample_key = "sampleID" if "sampleID" in pdac.obs else "sample"
    pdac_snapshot = build_design_snapshot(
        pdac,
        dataset="lin2020.pdac",
        sample_key=pdac_sample_key,
        condition_key="condition",
        experimental_unit_key=pdac_sample_key,
        cell_type_key="cell_type",
    )
    if getattr(pdac, "file", None) is not None:
        pdac.file.close()

    matrix_path = output_dir / "metadata_propagation_matrix.tsv"
    pd.DataFrame(_matrix_rows(pbmc_snapshot, pdac_snapshot)).to_csv(
        matrix_path, sep="\t", index=False
    )
    design_path = output_dir / "real_data_design_audit.tsv"
    pd.DataFrame(
        [
            {
                **snapshot,
                "conditions": json.dumps(snapshot.get("conditions", [])),
                "experimental_units_per_condition": json.dumps(
                    snapshot.get("experimental_units_per_condition", {})
                ),
                "usable_cell_types": json.dumps(snapshot.get("usable_cell_types", [])),
                "blockers": json.dumps(snapshot.get("blockers", [])),
            }
            for snapshot in (pbmc_snapshot, pdac_snapshot)
        ]
    ).to_csv(design_path, sep="\t", index=False)

    proportion_design = pbmc.uns["sclucid"]["proportion"]["design"]
    pbmc_inference_ready = bool(
        not de_results.empty
        and de_results["valid_for_publication_inference"].any()
        and not prop_stats.empty
        and prop_stats["valid_for_publication_inference"].any()
    )
    pdac_blocked_safely = bool(
        pdac_snapshot["status"] == "BLOCKED"
        and "fewer than two observed condition levels" in pdac_snapshot["blockers"]
    )
    manifest = {
        "schema_version": "1.0",
        "scope": "analysis_inference_contract",
        "gate_status": "PASS" if pbmc_inference_ready and pdac_blocked_safely else "REVIEW",
        "claim_boundary": (
            "Kang2018 supports a paired donor-level executable contract check. "
            "Lin2020 PDAC has one condition and empty cell-type labels, so it supports "
            "only fail-closed design evidence, not a tumor condition-effect claim."
        ),
        "datasets": {
            "kang2018.pbmc": {
                "provenance": _file_provenance(pbmc_path),
                "design": pbmc_snapshot,
                "complete_donors": complete_donors,
                "selected_cell_types": cell_types,
                "selected_genes": genes,
                "proportion_design": proportion_design,
                "pseudobulk_design": de_design,
                "inference_ready": pbmc_inference_ready,
            },
            "lin2020.pdac": {
                "provenance": _file_provenance(pdac_path),
                "design": pdac_snapshot,
                "inference_ready": False,
                "blocked_safely": pdac_blocked_safely,
            },
        },
        "artifacts": {
            "metadata_propagation_matrix": str(matrix_path),
            "real_data_design_audit": str(design_path),
            "pbmc_proportion_estimates": str(prop_path),
            "pbmc_proportion_statistics": str(prop_stats_path),
            "pbmc_pseudobulk_de": str(de_path),
        },
    }
    manifest_path = output_dir / "analysis_inference_evidence_manifest.json"
    manifest["artifacts"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(_json_safe(manifest), indent=2), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/work/analysis_inference_contract"),
    )
    parser.add_argument("--max-genes", type=int, default=200)
    parser.add_argument("--max-cell-types", type=int, default=4)
    parser.add_argument("--min-cells-per-sample", type=int, default=10)
    args = parser.parse_args()
    manifest = run(
        args.data_dir,
        args.output_dir,
        max_genes=args.max_genes,
        max_cell_types=args.max_cell_types,
        min_cells_per_sample=args.min_cells_per_sample,
    )
    print(json.dumps(_json_safe(manifest), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
