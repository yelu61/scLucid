"""Build and validate prediction-blinded QC expert-truth packs.

The reviewer-facing files contain only input-derived QC evidence and anonymous
identifiers. scLucid policy calls, original identifiers, and source paths stay
in the sealed mapping and are used only after expert labels are frozen.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd

from scLucid.qc.policy.reviewer import _metric_frame
from scLucid.utils.context import ProjectContext

TRUTH_PACK_SCHEMA_VERSION = "sclucid_qc_truth_pack_v1"
ALLOWED_LABELS = {"KEEP", "REMOVE", "UNCERTAIN"}
LABEL_COLUMNS = [
    "case_id",
    "evidence_hash",
    "expert_label",
    "reviewer_id",
    "confidence",
    "rationale",
]
FORBIDDEN_REVIEWER_COLUMNS = {
    "dataset_key",
    "sample_key",
    "original_obs_name",
    "policy_status",
    "selector_call",
    "candidate_call",
    "qc_decision",
    "predicted_doublet",
    "doublet_ground_truth",
}


@dataclass(frozen=True)
class TruthDatasetSpec:
    """Minimum source metadata needed to build a QC truth pack."""

    key: str
    path: str
    tissue: str = "unknown"
    dataset_type: str = "unknown"
    assay: str = "scrna"
    species: str = "auto"
    sample_key: str | None = None
    condition_key: str | None = None
    input_provenance: str = "filtered_counts"
    source_role: str = "development"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stable_seed(seed: int, *tokens: str) -> int:
    payload = "\0".join([str(seed), *map(str, tokens)]).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % (2**32)


def _case_id(seed: int, dataset_key: str, obs_name: str) -> str:
    payload = f"{seed}\0{dataset_key}\0{obs_name}".encode()
    return f"C-{hashlib.sha256(payload).hexdigest()[:16].upper()}"


def _evidence_hash(row: dict[str, Any]) -> str:
    payload = json.dumps(row, sort_keys=True, ensure_ascii=True, default=str, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _infer_species(adata) -> str:
    names = pd.Index(adata.var_names.astype(str))
    n_human = int(names.str.match(r"^MT-").sum())
    n_mouse = int(names.str.match(r"^mt-").sum())
    return "mouse" if n_mouse > n_human else "human"


def _percentile(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() < 2:
        return pd.Series(np.nan, index=values.index)
    return 100.0 * numeric.rank(method="average", pct=True)


def _quantiles(values: pd.Series) -> dict[str, float | None]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return {"p05": None, "p50": None, "p95": None}
    return {
        "p05": float(numeric.quantile(0.05)),
        "p50": float(numeric.quantile(0.50)),
        "p95": float(numeric.quantile(0.95)),
    }


def _select_cases(
    metrics: pd.DataFrame,
    sample_values: pd.Series,
    *,
    dataset_key: str,
    seed: int,
    primary_per_library: int,
    challenge_per_axis: int,
) -> pd.DataFrame:
    """Select a uniform primary set plus a prediction-independent challenge set."""
    rows: list[dict[str, Any]] = []
    samples = sample_values.astype(str)
    for sample in sorted(samples.unique()):
        positions = np.flatnonzero(samples.to_numpy() == sample)
        rng = np.random.default_rng(_stable_seed(seed, dataset_key, sample))
        n_primary = min(int(primary_per_library), len(positions))
        primary = set(rng.choice(positions, size=n_primary, replace=False).tolist())

        sub = metrics.iloc[positions]
        challenge: set[int] = set()
        axes = (
            ("n_genes_by_counts", True),
            ("pct_counts_mt", False),
            ("pct_counts_in_top_20_genes", False),
        )
        for column, ascending in axes:
            finite = sub[column].dropna().sort_values(ascending=ascending)
            chosen_names = finite.index[: int(challenge_per_axis)]
            challenge.update(metrics.index.get_indexer(chosen_names).tolist())
        challenge.difference_update(primary)

        for position in sorted(primary | challenge):
            rows.append(
                {
                    "position": int(position),
                    "sampling_tier": (
                        "primary_uniform" if position in primary else "secondary_metric_challenge"
                    ),
                }
            )
    return pd.DataFrame(rows).drop_duplicates("position").sort_values("position")


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def _label_template(evidence: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "case_id": evidence["case_id"],
            "evidence_hash": evidence["evidence_hash"],
            "expert_label": "",
            "reviewer_id": "",
            "confidence": "",
            "rationale": "",
        }
    )


def _assert_no_reviewer_leakage(frames: Iterable[pd.DataFrame]) -> None:
    for frame in frames:
        overlap = FORBIDDEN_REVIEWER_COLUMNS.intersection(frame.columns)
        if overlap:
            raise RuntimeError(f"Reviewer pack leaks forbidden columns: {sorted(overlap)}")


def build_truth_pack(
    specs: Iterable[TruthDatasetSpec],
    output_dir: Path,
    *,
    seed: int = 20260819,
    primary_per_library: int = 60,
    challenge_per_axis: int = 15,
) -> dict[str, Any]:
    """Create a deterministic reviewer pack and a separate sealed mapping."""
    specs = sorted(specs, key=lambda item: item.key)
    if not specs:
        raise ValueError("At least one truth-pack dataset is required.")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty truth pack: {output_dir}")

    reviewer_dir = output_dir / "reviewer"
    sealed_dir = output_dir / "sealed"
    reviewer_dir.mkdir(parents=True, exist_ok=True)
    sealed_dir.mkdir(parents=True, exist_ok=True)

    dataset_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    sealed_datasets: list[dict[str, Any]] = []
    sealed_samples: list[dict[str, Any]] = []
    sealed_cells: list[dict[str, Any]] = []

    for dataset_number, spec in enumerate(specs, start=1):
        source_path = Path(spec.path).expanduser().resolve()
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        dataset_alias = f"D{dataset_number:03d}"
        adata = ad.read_h5ad(source_path)
        species = _infer_species(adata) if spec.species == "auto" else spec.species
        sample_key = spec.sample_key if spec.sample_key in adata.obs else None
        if sample_key:
            sample_values = adata.obs[sample_key].astype(str)
        else:
            sample_values = pd.Series("single_library", index=adata.obs_names, dtype=str)
        context = ProjectContext(
            dataset_type=spec.dataset_type,
            tissue=spec.tissue,
            assay=spec.assay,
            species=species,
            sample_key=sample_key,
            condition_key=(spec.condition_key if spec.condition_key in adata.obs else None),
            is_multi_sample=bool(sample_values.nunique() > 1),
            input_provenance=spec.input_provenance,
        )
        metrics, metric_provenance = _metric_frame(adata, context)
        selections = _select_cases(
            metrics,
            sample_values,
            dataset_key=spec.key,
            seed=seed,
            primary_per_library=primary_per_library,
            challenge_per_axis=challenge_per_axis,
        )
        sample_aliases = {
            sample: f"{dataset_alias}-L{idx:03d}"
            for idx, sample in enumerate(sorted(sample_values.unique()), start=1)
        }

        dataset_row = {
            "dataset_alias": dataset_alias,
            "tissue": spec.tissue,
            "dataset_type": spec.dataset_type,
            "assay": spec.assay,
            "species": species,
            "input_provenance": spec.input_provenance,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_libraries": int(sample_values.nunique()),
            "mt_evidence_status": (
                "AVAILABLE" if metric_provenance["mt_available"] else "NOT_EVALUABLE"
            ),
        }
        dataset_row["evidence_hash"] = _evidence_hash(dataset_row)
        dataset_rows.append(dataset_row)

        for sample in sorted(sample_values.unique()):
            mask = sample_values.to_numpy() == sample
            sub = metrics.iloc[np.flatnonzero(mask)]
            sample_case = f"S-{hashlib.sha256(f'{seed}|{spec.key}|{sample}'.encode()).hexdigest()[:16].upper()}"
            genes_q = _quantiles(sub["n_genes_by_counts"])
            counts_q = _quantiles(sub["total_counts"])
            mt_q = _quantiles(sub["pct_counts_mt"])
            top_q = _quantiles(sub["pct_counts_in_top_20_genes"])
            evidence = {
                "case_id": sample_case,
                "dataset_alias": dataset_alias,
                "library_alias": sample_aliases[sample],
                "n_cells": int(mask.sum()),
                "genes_p05": genes_q["p05"],
                "genes_p50": genes_q["p50"],
                "genes_p95": genes_q["p95"],
                "counts_p05": counts_q["p05"],
                "counts_p50": counts_q["p50"],
                "counts_p95": counts_q["p95"],
                "mt_p05": mt_q["p05"],
                "mt_p50": mt_q["p50"],
                "mt_p95": mt_q["p95"],
                "top20_p05": top_q["p05"],
                "top20_p50": top_q["p50"],
                "top20_p95": top_q["p95"],
            }
            evidence["evidence_hash"] = _evidence_hash(evidence)
            sample_rows.append(evidence)
            sealed_samples.append(
                {
                    "case_id": sample_case,
                    "dataset_key": spec.key,
                    "dataset_alias": dataset_alias,
                    "sample_key": sample_key or "__single_library__",
                    "original_sample": sample,
                    "library_alias": sample_aliases[sample],
                }
            )

        for selection in selections.itertuples(index=False):
            position = int(selection.position)
            obs_name = str(adata.obs_names[position])
            sample = str(sample_values.iloc[position])
            sub_mask = sample_values.to_numpy() == sample
            sub_metrics = metrics.iloc[np.flatnonzero(sub_mask)]
            row = metrics.iloc[position]
            case_id = _case_id(seed, spec.key, obs_name)
            evidence = {
                "case_id": case_id,
                "dataset_alias": dataset_alias,
                "library_alias": sample_aliases[sample],
                "total_counts": float(row["total_counts"]),
                "n_genes_by_counts": float(row["n_genes_by_counts"]),
                "pct_counts_mt": (
                    float(row["pct_counts_mt"]) if np.isfinite(row["pct_counts_mt"]) else None
                ),
                "pct_counts_in_top_20_genes": float(row["pct_counts_in_top_20_genes"]),
                "genes_library_percentile": float(
                    _percentile(sub_metrics["n_genes_by_counts"]).loc[obs_name]
                ),
                "mt_library_percentile": (
                    float(_percentile(sub_metrics["pct_counts_mt"]).loc[obs_name])
                    if np.isfinite(row["pct_counts_mt"])
                    else None
                ),
                "top20_library_percentile": float(
                    _percentile(sub_metrics["pct_counts_in_top_20_genes"]).loc[obs_name]
                ),
            }
            evidence["evidence_hash"] = _evidence_hash(evidence)
            cell_rows.append(evidence)
            sealed_cells.append(
                {
                    "case_id": case_id,
                    "dataset_key": spec.key,
                    "dataset_alias": dataset_alias,
                    "original_obs_name": obs_name,
                    "original_sample": sample,
                    "library_alias": sample_aliases[sample],
                    "sampling_tier": selection.sampling_tier,
                }
            )

        sealed_datasets.append(
            {
                **asdict(spec),
                "path": str(source_path),
                "dataset_alias": dataset_alias,
                "source_sha256": sha256_file(source_path),
                "source_size": int(source_path.stat().st_size),
                "resolved_species": species,
                "resolved_sample_key": sample_key,
                "metric_provenance": metric_provenance,
                "context": context.model_dump(mode="json"),
            }
        )
        del adata

    datasets = pd.DataFrame(dataset_rows)
    samples = pd.DataFrame(sample_rows).sample(frac=1, random_state=seed).reset_index(drop=True)
    cells = pd.DataFrame(cell_rows).sample(frac=1, random_state=seed).reset_index(drop=True)
    sealed_sample_frame = pd.DataFrame(sealed_samples)
    sealed_cell_frame = pd.DataFrame(sealed_cells)
    _assert_no_reviewer_leakage((datasets, samples, cells))

    reviewer_paths = {
        "dataset_context": reviewer_dir / "dataset_context.tsv",
        "sample_evidence": reviewer_dir / "sample_evidence.tsv",
        "cell_evidence": reviewer_dir / "cell_evidence.tsv",
        "sample_labels": reviewer_dir / "sample_labels.tsv",
        "cell_labels": reviewer_dir / "cell_labels.tsv",
    }
    _write_tsv(datasets, reviewer_paths["dataset_context"])
    _write_tsv(samples, reviewer_paths["sample_evidence"])
    _write_tsv(cells, reviewer_paths["cell_evidence"])
    _write_tsv(_label_template(samples), reviewer_paths["sample_labels"])
    _write_tsv(_label_template(cells), reviewer_paths["cell_labels"])
    _write_tsv(sealed_sample_frame, sealed_dir / "sample_key.tsv")
    _write_tsv(sealed_cell_frame, sealed_dir / "cell_key.tsv")

    manifest = {
        "schema_version": TRUTH_PACK_SCHEMA_VERSION,
        "status": "AWAITING_BLINDED_EXPERT_LABELS",
        "seed": seed,
        "primary_per_library": primary_per_library,
        "challenge_per_axis": challenge_per_axis,
        "primary_endpoint_sampling_tier": "primary_uniform",
        "allowed_labels": sorted(ALLOWED_LABELS),
        "datasets": sealed_datasets,
        "reviewer_files": {
            name: {
                "path": str(path.relative_to(output_dir)),
                "sha256": sha256_file(path),
            }
            for name, path in reviewer_paths.items()
        },
    }
    (sealed_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    reviewer_manifest = {
        "schema_version": TRUTH_PACK_SCHEMA_VERSION,
        "status": "AWAITING_BLINDED_EXPERT_LABELS",
        "n_datasets": len(datasets),
        "n_sample_cases": len(samples),
        "n_cell_cases": len(cells),
        "allowed_labels": sorted(ALLOWED_LABELS),
        "leakage_guard": "No scLucid prediction, candidate call, original identifier, or source path is present.",
    }
    (reviewer_dir / "manifest.json").write_text(json.dumps(reviewer_manifest, indent=2) + "\n")
    (output_dir / "README.md").write_text(
        "# scLucid blinded QC truth pack\n\n"
        "Give only the `reviewer/` directory to reviewers. Keep `sealed/` unavailable "
        "until labels are frozen. Reviewers edit only `sample_labels.tsv` and "
        "`cell_labels.tsv`; allowed labels are KEEP, REMOVE, and UNCERTAIN.\n"
    )
    return manifest


def load_external_specs(path: Path) -> list[TruthDatasetSpec]:
    payload = json.loads(path.read_text())
    rows = payload.get("datasets", payload)
    if not isinstance(rows, list):
        raise ValueError("External spec must be a list or contain a 'datasets' list.")
    return [TruthDatasetSpec(**row) for row in rows]


def validate_frozen_labels(
    evidence_path: Path,
    labels_path: Path,
) -> tuple[pd.DataFrame, list[str]]:
    """Validate case coverage, evidence hashes, labels, and reviewer metadata."""
    evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    labels = pd.read_csv(labels_path, sep="\t", dtype=str).fillna("")
    return validate_frozen_label_frame(evidence, labels)


def validate_frozen_label_frame(
    evidence: pd.DataFrame,
    labels: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Validate an in-memory expert label table against frozen evidence."""
    evidence = evidence.astype(str).fillna("")
    labels = labels.astype(str).fillna("")
    issues: list[str] = []
    missing_columns = [column for column in LABEL_COLUMNS if column not in labels]
    if missing_columns:
        return labels, [f"Missing label columns: {missing_columns}"]
    if labels["case_id"].duplicated().any():
        issues.append("Duplicate case_id values are not allowed.")
    expected = set(evidence["case_id"])
    observed = set(labels["case_id"])
    if expected != observed:
        issues.append(
            f"Label case coverage mismatch: missing={len(expected - observed)}, extra={len(observed - expected)}."
        )
    expected_hash = evidence.set_index("case_id")["evidence_hash"]
    observed_hash = labels.set_index("case_id")["evidence_hash"]
    common = expected_hash.index.intersection(observed_hash.index)
    if not expected_hash.loc[common].equals(observed_hash.loc[common]):
        issues.append("Evidence hashes changed; reviewer evidence or label linkage is not frozen.")
    normalized = labels["expert_label"].str.upper().str.strip()
    invalid = sorted(set(normalized) - ALLOWED_LABELS - {""})
    if invalid:
        issues.append(f"Unsupported expert labels: {invalid}")
    missing = int((normalized == "").sum())
    if missing:
        issues.append(f"Missing expert labels: {missing}")
    labeled = normalized != ""
    if (labels.loc[labeled, "reviewer_id"].str.strip() == "").any():
        issues.append("Every labeled case requires reviewer_id.")
    labels = labels.copy()
    labels["expert_label"] = normalized
    return labels, issues


__all__ = [
    "ALLOWED_LABELS",
    "TRUTH_PACK_SCHEMA_VERSION",
    "TruthDatasetSpec",
    "build_truth_pack",
    "load_external_specs",
    "sha256_file",
    "validate_frozen_label_frame",
    "validate_frozen_labels",
]
