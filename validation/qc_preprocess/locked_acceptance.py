"""Locked scientific acceptance calculations for evidence-calibrated policies.

This module evaluates predictions against blinded external labels. It does not
create expert truth and it never converts unlabeled datasets into evidence of
superiority.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd


def _binary_metrics(truth: pd.Series, predicted_remove: pd.Series) -> dict[str, float]:
    keep = truth == "KEEP"
    remove = truth == "REMOVE"
    false_removal = float(predicted_remove[keep].mean()) if keep.any() else float("nan")
    recall = float(predicted_remove[remove].mean()) if remove.any() else float("nan")
    return {
        "n_keep": int(keep.sum()),
        "n_remove": int(remove.sum()),
        "false_removal_rate": false_removal,
        "low_quality_recall": recall,
    }


def _grouped_bootstrap_difference(
    truth: pd.Series,
    selector: pd.Series,
    baseline: pd.Series,
    groups: pd.Series,
    *,
    n_bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    unique_groups = pd.Index(groups.dropna().astype(str).unique())
    if len(unique_groups) < 2:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "Grouped bootstrap requires at least two sample/library groups.",
            "ci95": [None, None],
        }
    rng = np.random.default_rng(seed)
    differences: list[float] = []
    group_values = groups.astype(str)
    for _ in range(int(n_bootstrap)):
        sampled = rng.choice(unique_groups.to_numpy(), size=len(unique_groups), replace=True)
        positions: list[int] = []
        for group in sampled:
            positions.extend(np.flatnonzero(group_values.to_numpy() == group).tolist())
        sampled_truth = truth.iloc[positions]
        remove = sampled_truth == "REMOVE"
        if not remove.any():
            continue
        selector_recall = float(selector.iloc[positions][remove.to_numpy()].mean())
        baseline_recall = float(baseline.iloc[positions][remove.to_numpy()].mean())
        differences.append(selector_recall - baseline_recall)
    if not differences:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "No expert REMOVE labels were present in bootstrap samples.",
            "ci95": [None, None],
        }
    return {
        "status": "EVALUATED",
        "n_bootstrap_effective": len(differences),
        "ci95": [
            float(np.percentile(differences, 2.5)),
            float(np.percentile(differences, 97.5)),
        ],
    }


def evaluate_qc_policy_acceptance(
    policy: Any,
    expert_truth: pd.Series,
    sample_or_library: pd.Series,
    *,
    baseline_names: Sequence[str] = ("expert_global", "per_sample_mad"),
    max_keep_false_removal: float = 0.02,
    min_absolute_recall_gain: float = 0.05,
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Evaluate the locked QC endpoint using KEEP/REMOVE/UNCERTAIN truth."""
    labels = expert_truth.astype(str).str.upper()
    invalid = sorted(set(labels.dropna()) - {"KEEP", "REMOVE", "UNCERTAIN"})
    if invalid:
        raise ValueError(f"Unsupported expert labels: {invalid}")
    common = labels.index.intersection(sample_or_library.index)
    labels = labels.loc[common]
    groups = sample_or_library.loc[common]
    evaluated = labels.isin(["KEEP", "REMOVE"])
    labels = labels[evaluated]
    groups = groups[evaluated]
    if not (labels == "KEEP").any() or not (labels == "REMOVE").any():
        return {
            "status": "NOT_EVALUABLE",
            "reason": "Both expert KEEP and REMOVE labels are required; UNCERTAIN is excluded.",
        }

    selector_names = set(policy.remove_obs_names)
    selector = pd.Series(labels.index.isin(selector_names), index=labels.index)
    selector_metrics = _binary_metrics(labels, selector)
    candidate_lookup = {
        record["name"]: record for record in getattr(policy, "candidate_policies", [])
    }
    comparisons: list[dict[str, Any]] = []
    for idx, name in enumerate(baseline_names):
        record = candidate_lookup.get(name)
        if not record or "flagged_obs_names" not in record:
            comparisons.append(
                {
                    "baseline": name,
                    "status": "NOT_EVALUABLE",
                    "reason": "Exact baseline calls are absent from the policy.",
                }
            )
            continue
        baseline_names_set = set(record["flagged_obs_names"])
        baseline = pd.Series(labels.index.isin(baseline_names_set), index=labels.index)
        baseline_metrics = _binary_metrics(labels, baseline)
        gain = selector_metrics["low_quality_recall"] - baseline_metrics["low_quality_recall"]
        bootstrap = _grouped_bootstrap_difference(
            labels,
            selector,
            baseline,
            groups,
            n_bootstrap=n_bootstrap,
            seed=seed + idx,
        )
        comparisons.append(
            {
                "baseline": name,
                "status": bootstrap["status"],
                "baseline_metrics": baseline_metrics,
                "absolute_recall_gain": float(gain),
                "grouped_bootstrap": bootstrap,
                "passes_gain_gate": bool(
                    gain >= min_absolute_recall_gain
                    and bootstrap["status"] == "EVALUATED"
                    and bootstrap["ci95"][0] > 0
                ),
            }
        )

    gate = bool(
        selector_metrics["false_removal_rate"] <= max_keep_false_removal
        and comparisons
        and all(row.get("passes_gain_gate", False) for row in comparisons)
    )
    return {
        "status": "PASS" if gate else "FAIL",
        "uncertain_excluded": int((expert_truth.astype(str).str.upper() == "UNCERTAIN").sum()),
        "selector_metrics": selector_metrics,
        "comparisons": comparisons,
        "thresholds": {
            "max_keep_false_removal": max_keep_false_removal,
            "min_absolute_recall_gain": min_absolute_recall_gain,
            "grouped_bootstrap_ci_lower_must_exceed": 0.0,
        },
        "claim_boundary": (
            "PASS supports superiority only for the labeled datasets, tasks, and baselines."
            if gate
            else "Scientific superiority is not established."
        ),
    }


def evaluate_preprocess_policy_acceptance(
    results: pd.DataFrame,
    *,
    dataset_key: str = "dataset",
    candidate_key: str = "candidate",
    selected_key: str = "selected",
    utility_key: str = "preregistered_task_utility",
    biology_loss_key: str = "biology_loss",
    baseline_name: str = "standard_unintegrated",
    max_regret: float = 0.05,
    max_biology_loss: float = 0.02,
) -> dict[str, Any]:
    """Evaluate held-out regret and the simple-baseline fallback gate."""
    required = {
        dataset_key,
        candidate_key,
        selected_key,
        utility_key,
        biology_loss_key,
    }
    missing = sorted(required - set(results.columns))
    if missing:
        raise ValueError(f"Missing preprocessing acceptance columns: {missing}")

    rows: list[dict[str, Any]] = []
    for dataset, frame in results.groupby(dataset_key, observed=True):
        selected = frame[frame[selected_key].astype(bool)]
        baseline = frame[frame[candidate_key].astype(str) == baseline_name]
        if len(selected) != 1 or len(baseline) != 1:
            rows.append(
                {
                    "dataset": str(dataset),
                    "status": "NOT_EVALUABLE",
                    "reason": "Exactly one selected candidate and one simple baseline are required.",
                }
            )
            continue
        eligible = frame[frame[biology_loss_key].astype(float) <= max_biology_loss]
        if eligible.empty:
            rows.append(
                {
                    "dataset": str(dataset),
                    "status": "FAIL",
                    "reason": "No candidate satisfies the biology-loss guardrail.",
                }
            )
            continue
        selected_row = selected.iloc[0]
        baseline_row = baseline.iloc[0]
        best_utility = float(eligible[utility_key].max())
        selected_utility = float(selected_row[utility_key])
        regret = max(0.0, best_utility - selected_utility) / max(abs(best_utility), 1e-12)
        selected_is_baseline = str(selected_row[candidate_key]) == baseline_name
        complex_has_pareto_gain = bool(
            selected_utility > float(baseline_row[utility_key])
            and float(selected_row[biology_loss_key]) <= float(baseline_row[biology_loss_key])
        )
        simple_fallback_ok = selected_is_baseline or complex_has_pareto_gain
        passes = bool(
            float(selected_row[biology_loss_key]) <= max_biology_loss
            and regret <= max_regret
            and simple_fallback_ok
        )
        rows.append(
            {
                "dataset": str(dataset),
                "status": "PASS" if passes else "FAIL",
                "selected": str(selected_row[candidate_key]),
                "regret": float(regret),
                "biology_loss": float(selected_row[biology_loss_key]),
                "simple_fallback_or_pareto_gain": simple_fallback_ok,
            }
        )
    overall = bool(rows) and all(row["status"] == "PASS" for row in rows)
    return {
        "status": "PASS" if overall else "FAIL",
        "datasets": rows,
        "thresholds": {
            "max_regret": max_regret,
            "max_biology_loss": max_biology_loss,
        },
        "claim_boundary": (
            "PASS supports method selection only on the registered held-out tasks."
            if overall
            else "Preprocessing selector superiority is not established."
        ),
    }


def evaluate_real_project_ux_acceptance(
    records: pd.DataFrame,
    *,
    expected_projects: Sequence[str],
    min_config_reduction: float = 0.70,
) -> dict[str, Any]:
    """Evaluate the locked real-project usability contract without imputing records."""
    required = {
        "project",
        "legacy_config_fields",
        "current_config_fields",
        "manual_predicted_doublet_deletion",
        "manual_review_summary_edit",
        "schema_bypass",
        "project_specific_patch_count",
        "run_evidence_status",
    }
    missing_columns = sorted(required - set(records.columns))
    if missing_columns:
        return {
            "status": "BLOCKED",
            "reason": f"Missing UX evidence columns: {missing_columns}",
            "projects": [],
        }

    observed = set(records["project"].dropna().astype(str))
    missing_projects = sorted(set(map(str, expected_projects)) - observed)
    rows: list[dict[str, Any]] = []
    for project in expected_projects:
        frame = records[records["project"].astype(str) == str(project)]
        if len(frame) != 1:
            rows.append(
                {
                    "project": str(project),
                    "status": "BLOCKED",
                    "reason": "Exactly one executed UX record is required.",
                }
            )
            continue
        row = frame.iloc[0]
        try:
            legacy = int(row["legacy_config_fields"])
            current = int(row["current_config_fields"])
            patch_count = int(row["project_specific_patch_count"])
        except (TypeError, ValueError):
            rows.append(
                {
                    "project": str(project),
                    "status": "BLOCKED",
                    "reason": "Configuration and patch counts must be recorded as integers.",
                }
            )
            continue
        if legacy <= 0 or current < 0:
            rows.append(
                {
                    "project": str(project),
                    "status": "BLOCKED",
                    "reason": "legacy_config_fields must be positive and current count non-negative.",
                }
            )
            continue
        reduction = 1.0 - current / legacy

        def _is_false(value: Any) -> bool:
            return str(value).strip().lower() in {"false", "0", "no"}

        no_workarounds = all(
            _is_false(row[column])
            for column in (
                "manual_predicted_doublet_deletion",
                "manual_review_summary_edit",
                "schema_bypass",
            )
        )
        stable_evidence = str(row["run_evidence_status"]).strip().upper() in {
            "READY",
            "REVIEW",
        }
        passed = bool(
            reduction >= min_config_reduction
            and no_workarounds
            and patch_count == 0
            and stable_evidence
        )
        rows.append(
            {
                "project": str(project),
                "status": "PASS" if passed else "FAIL",
                "config_field_reduction": float(reduction),
                "no_manual_workarounds": no_workarounds,
                "project_specific_patch_count": patch_count,
                "run_evidence_status": str(row["run_evidence_status"]).strip().upper(),
            }
        )

    overall = not missing_projects and bool(rows) and all(row["status"] == "PASS" for row in rows)
    blocked = bool(missing_projects) or any(row["status"] == "BLOCKED" for row in rows)
    return {
        "status": "PASS" if overall else "BLOCKED" if blocked else "FAIL",
        "missing_projects": missing_projects,
        "projects": rows,
        "thresholds": {"min_config_field_reduction": min_config_reduction},
        "claim_boundary": (
            "PASS supports the maintained-path usability claim only for the recorded projects."
            if overall
            else "The real-project usability gate is not established."
        ),
    }


__all__ = [
    "evaluate_preprocess_policy_acceptance",
    "evaluate_qc_policy_acceptance",
    "evaluate_real_project_ux_acceptance",
]
