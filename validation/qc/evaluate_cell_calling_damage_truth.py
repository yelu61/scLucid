#!/usr/bin/env python3
"""Evaluate independent cell-calling and damaged-cell truth tables.

The evaluator consumes externally supplied truth and predictions. It never
derives truth from scLucid predictions, and it cannot verify blinding or truth
independence from the table alone.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DEFAULT_REGISTRY = Path("validation/dataset_evidence_registry.json")
CELL_CALL_ENDPOINT = "qc_cell_calling"
DAMAGE_ENDPOINT = "qc_damage_classification"
TRUTH_CLASSES = {"INTACT", "DAMAGED", "EMPTY", "UNCERTAIN"}
QC_DECISIONS = {"KEEP", "REMOVE", "REVIEW", "UNCERTAIN"}


def _load_thresholds(registry_path: Path) -> dict[str, dict[str, Any]]:
    registry = json.loads(registry_path.read_text())
    endpoints = registry.get("endpoint_definitions", {})
    thresholds: dict[str, dict[str, Any]] = {}
    for endpoint_id in (CELL_CALL_ENDPOINT, DAMAGE_ENDPOINT):
        endpoint = endpoints.get(endpoint_id)
        if not isinstance(endpoint, dict) or not isinstance(endpoint.get("acceptance"), dict):
            raise ValueError(f"Registry has no acceptance contract for {endpoint_id!r}.")
        thresholds[endpoint_id] = dict(endpoint["acceptance"])
    return thresholds


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    raise ValueError("Input must be a .csv, .tsv, or tab-delimited .txt file.")


def _normalize_truth(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip().str.upper()
    invalid = sorted(set(normalized.dropna()) - TRUTH_CLASSES)
    if invalid:
        raise ValueError(f"Unsupported truth_cell_class values: {invalid}")
    return normalized


def _normalize_decision(series: pd.Series) -> pd.Series:
    bool_values = series.dropna().map(lambda value: isinstance(value, (bool, np.bool_)))
    if not bool_values.empty and bool_values.all():
        return series.map({True: "REMOVE", False: "KEEP"}).astype("string")
    normalized = series.astype("string").str.strip().str.upper()
    normalized = normalized.replace(
        {
            "TRUE": "REMOVE",
            "1": "REMOVE",
            "FALSE": "KEEP",
            "0": "KEEP",
        }
    )
    invalid = sorted(set(normalized.dropna()) - QC_DECISIONS)
    if invalid:
        raise ValueError(f"Unsupported QC decision values: {invalid}")
    return normalized.fillna("REVIEW")


def _normalize_call(series: pd.Series) -> pd.Series:
    truthy = {"TRUE", "1", "YES", "Y", "CALLED", "CELL"}
    falsey = {"FALSE", "0", "NO", "N", "NOT_CALLED", "EMPTY", "NONCELL"}

    def convert(value: Any) -> bool | None:
        if pd.isna(value):
            return None
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, (int, np.integer)) and value in {0, 1}:
            return bool(value)
        token = str(value).strip().upper().replace(" ", "_")
        if token in truthy:
            return True
        if token in falsey:
            return False
        raise ValueError(f"Unsupported predicted cell-call value: {value!r}")

    return series.map(convert).astype("boolean")


def _normalize_flag(series: pd.Series) -> pd.Series:
    return _normalize_call(series).fillna(False).astype(bool)


def _ratio(numerator: int, denominator: int) -> float | None:
    return float(numerator / denominator) if denominator else None


def _cell_call_metrics(
    truth: pd.Series,
    called: pd.Series,
    low_rna: pd.Series,
) -> dict[str, Any]:
    known_prediction = called.notna()
    true_cell = truth.isin(["INTACT", "DAMAGED"])
    empty = truth.eq("EMPTY")
    evaluable_truth = true_cell | empty
    low_rna_cell = true_cell & low_rna
    called_bool = called.fillna(False).astype(bool)
    n_called = int((called_bool & known_prediction & evaluable_truth).sum())
    return {
        "n_rows": int(len(truth)),
        "n_true_cells": int(true_cell.sum()),
        "n_low_rna_true_cells": int(low_rna_cell.sum()),
        "n_empty_droplets": int(empty.sum()),
        "n_called": n_called,
        "n_truth_uncertain_excluded": int((~evaluable_truth).sum()),
        "n_missing_predictions": int((~known_prediction & evaluable_truth).sum()),
        "true_cell_recall": _ratio(int((true_cell & called_bool).sum()), int(true_cell.sum())),
        "low_rna_cell_recall": _ratio(
            int((low_rna_cell & called_bool).sum()),
            int(low_rna_cell.sum()),
        ),
        "empty_droplet_fdr": _ratio(int((empty & called_bool).sum()), n_called),
    }


def _damage_metrics(truth: pd.Series, decision: pd.Series) -> dict[str, Any]:
    truth_binary = truth.isin(["INTACT", "DAMAGED"])
    decided = decision.isin(["KEEP", "REMOVE"])
    evaluated = truth_binary & decided
    intact = truth.eq("INTACT") & evaluated
    damaged = truth.eq("DAMAGED") & evaluated
    remove = decision.eq("REMOVE")
    return {
        "n_truth_binary": int(truth_binary.sum()),
        "n_evaluated": int(evaluated.sum()),
        "n_intact_evaluated": int(intact.sum()),
        "n_damaged_evaluated": int(damaged.sum()),
        "decision_coverage": _ratio(int(evaluated.sum()), int(truth_binary.sum())),
        "keep_false_removal_rate": _ratio(int((intact & remove).sum()), int(intact.sum())),
        "damaged_cell_recall": _ratio(int((damaged & remove).sum()), int(damaged.sum())),
    }


def _percentile_interval(values: Sequence[float]) -> list[float] | list[None]:
    if not values:
        return [None, None]
    return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]


def _grouped_bootstrap(
    frame: pd.DataFrame,
    groups: pd.Series,
    metric_fn: Callable[[pd.DataFrame], Mapping[str, float | None]],
    metric_names: Sequence[str],
    *,
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups = pd.Index(groups.dropna().astype(str).unique())
    if len(unique_groups) < 2:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "Grouped bootstrap requires at least two libraries.",
            "n_bootstrap_requested": int(n_bootstrap),
            "ci95": {name: [None, None] for name in metric_names},
        }
    group_values = groups.astype(str).to_numpy()
    rng = np.random.default_rng(seed)
    sampled_metrics: dict[str, list[float]] = {name: [] for name in metric_names}
    effective = 0
    for _ in range(int(n_bootstrap)):
        sampled_groups = rng.choice(
            unique_groups.to_numpy(),
            size=len(unique_groups),
            replace=True,
        )
        positions: list[int] = []
        for group in sampled_groups:
            positions.extend(np.flatnonzero(group_values == group).tolist())
        metrics = metric_fn(frame.iloc[positions])
        if not any(metrics.get(name) is not None for name in metric_names):
            continue
        effective += 1
        for name in metric_names:
            value = metrics.get(name)
            if value is not None:
                sampled_metrics[name].append(float(value))
    if not effective:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "No bootstrap replicate contained evaluable truth and predictions.",
            "n_bootstrap_requested": int(n_bootstrap),
            "ci95": {name: [None, None] for name in metric_names},
        }
    return {
        "status": "EVALUATED",
        "seed": int(seed),
        "n_bootstrap_requested": int(n_bootstrap),
        "n_bootstrap_effective": effective,
        "ci95": {name: _percentile_interval(sampled_metrics[name]) for name in metric_names},
    }


def _cell_call_endpoint(
    frame: pd.DataFrame,
    *,
    truth: pd.Series,
    groups: pd.Series,
    low_rna: pd.Series,
    called_col: str,
    thresholds: Mapping[str, Any],
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    if called_col not in frame:
        return {
            "status": "BLOCKED",
            "thresholds": dict(thresholds),
            "blockers": [f"Missing required prediction column: {called_col}"],
            "claim_boundary": "Cell-calling performance is not evaluable.",
        }
    called = _normalize_call(frame[called_col])
    aggregate = _cell_call_metrics(truth, called, low_rna)
    blockers: list[str] = []
    if aggregate["n_true_cells"] == 0:
        blockers.append("No intact or damaged true cells are present.")
    if aggregate["n_low_rna_true_cells"] == 0:
        blockers.append("No independently labelled low-RNA true cells are present.")
    if aggregate["n_empty_droplets"] == 0:
        blockers.append("No independently labelled empty droplets are present.")
    if aggregate["n_missing_predictions"]:
        blockers.append("One or more cell-call predictions are missing.")
    if aggregate["n_called"] == 0:
        blockers.append("No droplets were predicted as cells; empty-droplet FDR is undefined.")

    by_library = []
    for library, positions in groups.groupby(groups, sort=True).groups.items():
        row_metrics = _cell_call_metrics(
            truth.loc[positions],
            called.loc[positions],
            low_rna.loc[positions],
        )
        by_library.append({"library": str(library), **row_metrics})

    bootstrap_frame = pd.DataFrame(
        {"truth": truth, "called": called, "low_rna": low_rna},
        index=frame.index,
    )
    bootstrap = _grouped_bootstrap(
        bootstrap_frame,
        groups,
        lambda sampled: _cell_call_metrics(
            sampled["truth"],
            sampled["called"].astype("boolean"),
            sampled["low_rna"].astype(bool),
        ),
        ("true_cell_recall", "low_rna_cell_recall", "empty_droplet_fdr"),
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if blockers:
        status = "BLOCKED"
    else:
        passed = bool(
            aggregate["true_cell_recall"] >= thresholds["true_cell_recall_min"]
            and aggregate["low_rna_cell_recall"] >= thresholds["low_rna_cell_recall_min"]
            and aggregate["empty_droplet_fdr"] <= thresholds["empty_droplet_fdr_max"]
        )
        status = "PASS" if passed else "FAIL"
    return {
        "status": status,
        "thresholds": dict(thresholds),
        "aggregate_metrics": aggregate,
        "by_library": by_library,
        "grouped_bootstrap": bootstrap,
        "blockers": blockers,
        "claim_boundary": (
            "PASS supports the registered cell-calling thresholds only for this supplied "
            "independent truth table."
            if status == "PASS"
            else "Cell-calling acceptance is not established for this truth table."
        ),
    }


def _expected_calibration_error(truth: pd.Series, probabilities: pd.Series) -> dict[str, Any]:
    valid = truth.isin(["INTACT", "DAMAGED"]) & probabilities.notna()
    numeric = pd.to_numeric(probabilities.loc[valid], errors="coerce")
    valid_numeric = numeric.notna() & numeric.between(0.0, 1.0)
    numeric = numeric.loc[valid_numeric]
    binary_truth = truth.loc[numeric.index].eq("DAMAGED").astype(float)
    if numeric.empty:
        return {
            "status": "NOT_EVALUABLE",
            "n": 0,
            "expected_calibration_error": None,
        }
    bin_ids = np.minimum((numeric.to_numpy() * 10).astype(int), 9)
    ece = 0.0
    for bin_id in range(10):
        in_bin = bin_ids == bin_id
        if not in_bin.any():
            continue
        confidence = float(numeric.to_numpy()[in_bin].mean())
        observed = float(binary_truth.to_numpy()[in_bin].mean())
        ece += float(in_bin.mean()) * abs(confidence - observed)
    return {
        "status": "EVALUATED",
        "n": int(len(numeric)),
        "expected_calibration_error": float(ece),
    }


def _comparison_metrics(
    truth: pd.Series,
    selector: pd.Series,
    baseline: pd.Series,
) -> dict[str, Any]:
    common = truth.isin(["INTACT", "DAMAGED"])
    common &= selector.isin(["KEEP", "REMOVE"])
    common &= baseline.isin(["KEEP", "REMOVE"])
    selector_metrics = _damage_metrics(truth.loc[common], selector.loc[common])
    baseline_metrics = _damage_metrics(truth.loc[common], baseline.loc[common])
    selector_recall = selector_metrics["damaged_cell_recall"]
    baseline_recall = baseline_metrics["damaged_cell_recall"]
    gain = None
    if selector_recall is not None and baseline_recall is not None:
        gain = float(selector_recall - baseline_recall)
    return {
        "n_common_evaluated": int(common.sum()),
        "selector_metrics": selector_metrics,
        "baseline_metrics": baseline_metrics,
        "absolute_recall_gain": gain,
    }


def _damage_endpoint(
    frame: pd.DataFrame,
    *,
    truth: pd.Series,
    groups: pd.Series,
    decision_col: str,
    baseline_columns: Mapping[str, str],
    probability_col: str,
    thresholds: Mapping[str, Any],
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    blockers = [
        f"Missing required prediction column: {decision_col}"
        for _ in [0]
        if decision_col not in frame
    ]
    blockers.extend(
        f"Missing required baseline column for {name}: {column}"
        for name, column in baseline_columns.items()
        if column not in frame
    )
    if blockers:
        return {
            "status": "BLOCKED",
            "thresholds": dict(thresholds),
            "blockers": blockers,
            "claim_boundary": (
                "Damaged-cell performance and superiority are not evaluable; both exact "
                "baseline call columns are required."
            ),
        }

    decision = _normalize_decision(frame[decision_col])
    baselines = {
        name: _normalize_decision(frame[column]) for name, column in baseline_columns.items()
    }
    aggregate = _damage_metrics(truth, decision)
    if aggregate["n_intact_evaluated"] == 0:
        blockers.append("No decided, independently labelled intact cells are present.")
    if aggregate["n_damaged_evaluated"] == 0:
        blockers.append("No decided, independently labelled damaged cells are present.")

    comparisons: list[dict[str, Any]] = []
    for index, (name, baseline) in enumerate(baselines.items()):
        comparison = _comparison_metrics(truth, decision, baseline)
        if comparison["absolute_recall_gain"] is None:
            blockers.append(f"Baseline {name} has no common evaluable damaged-cell truth.")
        bootstrap_frame = pd.DataFrame(
            {"truth": truth, "selector": decision, "baseline": baseline},
            index=frame.index,
        )
        bootstrap = _grouped_bootstrap(
            bootstrap_frame,
            groups,
            lambda sampled: {
                "absolute_recall_gain": _comparison_metrics(
                    sampled["truth"],
                    sampled["selector"],
                    sampled["baseline"],
                )["absolute_recall_gain"]
            },
            ("absolute_recall_gain",),
            n_bootstrap=n_bootstrap,
            seed=seed + 101 + index,
        )
        comparisons.append(
            {
                "baseline": name,
                **comparison,
                "grouped_bootstrap": bootstrap,
            }
        )

    by_library = []
    for library, positions in groups.groupby(groups, sort=True).groups.items():
        row = {
            "library": str(library),
            "selector_metrics": _damage_metrics(truth.loc[positions], decision.loc[positions]),
            "comparisons": [],
        }
        for name, baseline in baselines.items():
            row["comparisons"].append(
                {
                    "baseline": name,
                    **_comparison_metrics(
                        truth.loc[positions],
                        decision.loc[positions],
                        baseline.loc[positions],
                    ),
                }
            )
        by_library.append(row)

    truth_binary = truth.isin(["INTACT", "DAMAGED"])
    review = decision.isin(["REVIEW", "UNCERTAIN"]) & truth_binary
    truth_uncertain = truth.eq("UNCERTAIN")
    coverage = {
        "primary_truth_rows": int(truth_binary.sum()),
        "primary_decided_rows": int((truth_binary & ~review).sum()),
        "prediction_review_or_uncertain_rows": int(review.sum()),
        "prediction_review_or_uncertain_rate": _ratio(int(review.sum()), int(truth_binary.sum())),
        "truth_uncertain_rows_excluded": int(truth_uncertain.sum()),
        "empty_rows_excluded": int(truth.eq("EMPTY").sum()),
    }
    bootstrap_frame = pd.DataFrame(
        {"truth": truth, "decision": decision},
        index=frame.index,
    )
    bootstrap = _grouped_bootstrap(
        bootstrap_frame,
        groups,
        lambda sampled: _damage_metrics(sampled["truth"], sampled["decision"]),
        ("keep_false_removal_rate", "damaged_cell_recall", "decision_coverage"),
        n_bootstrap=n_bootstrap,
        seed=seed + 100,
    )
    calibration = (
        _expected_calibration_error(truth, frame[probability_col])
        if probability_col in frame
        else {
            "status": "NOT_AVAILABLE",
            "n": 0,
            "expected_calibration_error": None,
            "reason": f"Optional probability column is absent: {probability_col}",
        }
    )
    if blockers:
        status = "BLOCKED"
    else:
        passes = bool(
            aggregate["keep_false_removal_rate"] <= thresholds["keep_false_removal_rate_max"]
            and all(
                row["absolute_recall_gain"]
                >= thresholds["absolute_recall_gain_vs_each_baseline_min"]
                for row in comparisons
            )
        )
        status = "PASS" if passes else "FAIL"
    return {
        "status": status,
        "thresholds": dict(thresholds),
        "aggregate_metrics": aggregate,
        "coverage": coverage,
        "by_library": by_library,
        "comparisons": comparisons,
        "grouped_bootstrap": bootstrap,
        "calibration": calibration,
        "blockers": blockers,
        "claim_boundary": (
            "PASS supports improved damaged-cell recall at the registered intact-cell "
            "guardrail only for this truth table and these two baselines."
            if status == "PASS"
            else "Damaged-cell superiority is not established for this truth table."
        ),
    }


def evaluate_truth_table(
    frame: pd.DataFrame,
    *,
    registry_path: Path = DEFAULT_REGISTRY,
    library_col: str = "library",
    truth_col: str = "truth_cell_class",
    cell_call_col: str = "predicted_cell_call",
    qc_decision_col: str = "predicted_qc_decision",
    low_rna_col: str = "truth_low_rna",
    probability_col: str = "predicted_damage_probability",
    baseline_global_col: str = "baseline_expert_global_qc_decision",
    baseline_mad_col: str = "baseline_per_sample_mad_qc_decision",
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Evaluate both endpoints against supplied, prediction-independent truth."""
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be positive.")
    thresholds = _load_thresholds(registry_path)
    common_missing = [column for column in (library_col, truth_col) if column not in frame]
    if common_missing:
        reason = f"Missing common truth columns: {common_missing}"
        endpoints = {
            endpoint_id: {
                "status": "BLOCKED",
                "thresholds": thresholds[endpoint_id],
                "blockers": [reason],
                "claim_boundary": "Performance is not evaluable.",
            }
            for endpoint_id in (CELL_CALL_ENDPOINT, DAMAGE_ENDPOINT)
        }
        return _assemble_report(
            frame,
            registry_path,
            endpoints,
            library_col=library_col,
            truth_col=truth_col,
            seed=seed,
            n_bootstrap=n_bootstrap,
        )

    truth = _normalize_truth(frame[truth_col])
    groups = frame[library_col].astype("string")
    if groups.isna().any() or groups.str.strip().eq("").any():
        raise ValueError("library values must be non-missing and non-empty.")
    low_rna = (
        _normalize_flag(frame[low_rna_col])
        if low_rna_col in frame
        else pd.Series(False, index=frame.index)
    )
    cell_call = _cell_call_endpoint(
        frame,
        truth=truth,
        groups=groups,
        low_rna=low_rna,
        called_col=cell_call_col,
        thresholds=thresholds[CELL_CALL_ENDPOINT],
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if low_rna_col not in frame:
        cell_call.setdefault("blockers", []).append(
            f"Missing required low-RNA truth column: {low_rna_col}"
        )
        cell_call["status"] = "BLOCKED"
        cell_call["claim_boundary"] = "Cell-calling performance is not evaluable."
    damage = _damage_endpoint(
        frame,
        truth=truth,
        groups=groups,
        decision_col=qc_decision_col,
        baseline_columns={
            "expert_global": baseline_global_col,
            "per_sample_mad": baseline_mad_col,
        },
        probability_col=probability_col,
        thresholds=thresholds[DAMAGE_ENDPOINT],
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    return _assemble_report(
        frame,
        registry_path,
        {CELL_CALL_ENDPOINT: cell_call, DAMAGE_ENDPOINT: damage},
        library_col=library_col,
        truth_col=truth_col,
        seed=seed,
        n_bootstrap=n_bootstrap,
    )


def _assemble_report(
    frame: pd.DataFrame,
    registry_path: Path,
    endpoints: Mapping[str, Mapping[str, Any]],
    *,
    library_col: str,
    truth_col: str,
    seed: int,
    n_bootstrap: int,
) -> dict[str, Any]:
    statuses = {name: result["status"] for name, result in endpoints.items()}
    if all(status == "PASS" for status in statuses.values()):
        status = "PASS"
    elif "BLOCKED" in statuses.values():
        status = "BLOCKED"
    else:
        status = "FAIL"
    libraries = (
        sorted(frame[library_col].dropna().astype(str).unique().tolist())
        if library_col in frame
        else []
    )
    return {
        "schema_version": "sclucid_cell_calling_damage_truth_evaluation_v1",
        "status": status,
        "endpoint_status": statuses,
        "input": {
            "n_rows": int(len(frame)),
            "libraries": libraries,
            "library_column": library_col,
            "truth_column": truth_col,
        },
        "threshold_source": str(registry_path),
        "bootstrap": {
            "unit": "library",
            "seed": int(seed),
            "n_bootstrap": int(n_bootstrap),
        },
        "truth_independence": {
            "status": "NOT_VERIFIED_FROM_TABLE",
            "requirement": (
                "Truth must be obtained independently of all evaluated prediction columns."
            ),
        },
        "endpoints": dict(endpoints),
        "claim_boundary": {
            "supported": (
                ["Registered endpoint thresholds were met on the supplied truth table."]
                if status == "PASS"
                else []
            ),
            "unsupported": [
                "The evaluator does not prove that truth was blinded or prediction-independent.",
                "A result on one supplied table does not establish cross-dataset or tumor-project generalization.",
                "Predictions are never used to create or revise truth labels.",
            ],
        },
    }


def write_report(report: Mapping[str, Any], output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "cell_calling_damage_truth_evaluation.json"
    markdown_path = output_dir / "cell_calling_damage_truth_evaluation.md"
    json_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = [
        "# Cell-calling and damaged-cell truth evaluation",
        "",
        f"Overall status: **{report['status']}**",
        "",
        "| Endpoint | Status |",
        "| --- | --- |",
    ]
    for endpoint, status in report["endpoint_status"].items():
        lines.append(f"| `{endpoint}` | **{status}** |")
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            *[f"- {item}" for item in report["claim_boundary"]["unsupported"]],
            "",
        ]
    )
    markdown_path.write_text("\n".join(lines))
    return json_path, markdown_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--library-col", default="library")
    parser.add_argument("--truth-col", default="truth_cell_class")
    parser.add_argument("--cell-call-col", default="predicted_cell_call")
    parser.add_argument("--qc-decision-col", default="predicted_qc_decision")
    parser.add_argument("--low-rna-col", default="truth_low_rna")
    parser.add_argument("--probability-col", default="predicted_damage_probability")
    parser.add_argument(
        "--baseline-global-col",
        default="baseline_expert_global_qc_decision",
    )
    parser.add_argument(
        "--baseline-mad-col",
        default="baseline_per_sample_mad_qc_decision",
    )
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    report = evaluate_truth_table(
        _read_table(args.input),
        registry_path=args.registry,
        library_col=args.library_col,
        truth_col=args.truth_col,
        cell_call_col=args.cell_call_col,
        qc_decision_col=args.qc_decision_col,
        low_rna_col=args.low_rna_col,
        probability_col=args.probability_col,
        baseline_global_col=args.baseline_global_col,
        baseline_mad_col=args.baseline_mad_col,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
    )
    json_path, markdown_path = write_report(report, args.output_dir)
    print(
        json.dumps(
            {"status": report["status"], "json": str(json_path), "markdown": str(markdown_path)}
        )
    )


if __name__ == "__main__":
    main()
