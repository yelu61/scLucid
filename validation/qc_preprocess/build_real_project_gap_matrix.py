"""Build a QC+preprocess gap matrix from local real-project Step1 notebooks.

The local Step1 notebooks are intentionally not tracked by git. This script
extracts their workflow blocks and maps recurring notebook logic to current
scLucid APIs, module gaps, and validation follow-up items.
"""

from __future__ import annotations

import ast
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_NOTEBOOKS = [
    ROOT / "Step1-QC_and_Preprocessing.ipynb",
    ROOT / "Step1-QC_and_Preprocessing_CT26.ipynb",
]
OUT_DIR = ROOT / "validation_outputs" / "qc_preprocess_real_project"


@dataclass(frozen=True)
class Rule:
    pattern: str
    block_name: str
    coverage: str
    gap: str
    owner: str
    priority: str


RULES = [
    Rule(
        pattern="run_iterative_qc",
        block_name="Reviewer-first QC",
        coverage="covered_by_stable_api",
        gap="Use validation to strengthen evidence; no new QC entrypoint needed.",
        owner="qc.workflow",
        priority="P1",
    ),
    Rule(
        pattern="summarize_qc_review_summary|qc_review|get\\(\"review_action_items\"",
        block_name="QC review checkpoint",
        coverage="covered_by_review_summary",
        gap="Ensure project-critical retention, doublet, ambient, and benchmark evidence is in the QC report.",
        owner="qc.trace/qc.reporting",
        priority="P1",
    ),
    Rule(
        pattern="run_preprocessing",
        block_name="Preprocess bootstrap",
        coverage="partially_workflow_owned",
        gap="Bootstrap remains valid for stepwise notebooks; compare with run_iterative_preprocessing evidence after execution.",
        owner="preprocess.workflow",
        priority="P1",
    ),
    Rule(
        pattern="run_iterative_preprocessing|RUN_CANONICAL_ITERATIVE_PREPROCESSING_COMPARISON|canonical_iterative_preprocess_summary",
        block_name="Canonical iterative preprocessing comparison",
        coverage="covered_by_stable_api",
        gap="Execute on real projects and compare review evidence against project-specific manual preprocessing.",
        owner="preprocess.workflow/preprocess.trace",
        priority="P0",
    ),
    Rule(
        pattern="recommend_analysis_parameters|recommendation_summary|workflow_recommendations",
        block_name="Recommendation advisor checkpoint",
        coverage="covered_by_recommendation_layer",
        gap="Compare advisor recommendations with project-selected QC/preprocess parameters.",
        owner="recommendation",
        priority="P0",
    ),
    Rule(
        pattern="find_hvgs|suggest_hvg_choice|select_and_audit_hvgs",
        block_name="Dual HVG strategy",
        coverage="low_level_api_only",
        gap="Compare notebook dual-HVG strategy with canonical iterative preprocessing HVG audit and reviewer evidence.",
        owner="preprocess.hvg/preprocess.workflow",
        priority="P0",
    ),
    Rule(
        pattern="evaluate_hvg_stability",
        block_name="HVG stability",
        coverage="low_level_api_only",
        gap="Verify HVG stability metrics are captured in preprocess review summaries and validation evidence tables.",
        owner="preprocess.hvg/preprocess.trace",
        priority="P0",
    ),
    Rule(
        pattern="diagnose_cell_cycle_regression|regress_out|scale_data|sc\\.tl\\.pca",
        block_name="Regression, scaling, PCA",
        coverage="low_level_api_only",
        gap="Verify regression/scaling/PCA decisions are explainable in preprocess review summaries.",
        owner="preprocess.workflow/preprocess.trace",
        priority="P0",
    ),
    Rule(
        pattern="run_embedding_pipeline",
        block_name="Diagnostic/final embedding",
        coverage="low_level_api_only",
        gap="Verify diagnostic and final embedding decisions are captured by iterative preprocessing summary.",
        owner="preprocess.neighbors/preprocess.workflow",
        priority="P0",
    ),
    Rule(
        pattern="decide_integration|batch_correction",
        block_name="Integration decision",
        coverage="low_level_api_only",
        gap="Verify integration decision, overcorrection risk, and final representation are captured as key evidence rows.",
        owner="preprocess.integrate/preprocess.workflow",
        priority="P0",
    ),
    Rule(
        pattern="finalize_manual_review_summary|validate_preprocess|summarize_preprocess",
        block_name="Final contract audit",
        coverage="manual_review_finalizer",
        gap="Use executed notebook outputs to build comparable review-summary completeness and key-decision evidence tables.",
        owner="preprocess.trace/preprocess.reporting",
        priority="P1",
    ),
]


def _call_name(node: ast.Call) -> str | None:
    func = node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        parts: list[str] = []
        current = func
        while isinstance(current, ast.Attribute):
            parts.append(current.attr)
            current = current.value
        if isinstance(current, ast.Name):
            parts.append(current.id)
            return ".".join(reversed(parts))
    return None


def _extract_calls(source: str) -> list[str]:
    calls: list[str] = []
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return calls
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = _call_name(node)
            if name:
                calls.append(name)
    return list(dict.fromkeys(calls))


def _nearest_markdown_heading(cells: list[dict], index: int) -> str:
    for j in range(index - 1, -1, -1):
        cell = cells[j]
        if cell.get("cell_type") != "markdown":
            continue
        text = "".join(cell.get("source", [])).strip()
        headings = [line.strip("# ").strip() for line in text.splitlines() if line.startswith("#")]
        if headings:
            return headings[-1]
        compact = re.sub(r"\s+", " ", text)
        if compact:
            return compact[:120]
    return "unlabeled"


def _classify(source: str, calls: Iterable[str]) -> list[Rule]:
    haystack = source + "\n" + "\n".join(calls)
    matched = [rule for rule in RULES if re.search(rule.pattern, haystack)]
    if matched:
        return matched
    return [
        Rule(
            pattern="",
            block_name="Project-specific or unsupported block",
            coverage="not_classified",
            gap="Review manually before productizing.",
            owner="manual_review",
            priority="P3",
        )
    ]


def build_rows(notebook_paths: Iterable[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in notebook_paths:
        if not path.exists():
            rows.append(
                {
                    "notebook": path.name,
                    "cell_index": "",
                    "heading": "",
                    "block_name": "Notebook missing",
                    "api_calls": "",
                    "coverage": "missing_input",
                    "gap": "Notebook not found in local workspace.",
                    "owner": "local_workspace",
                    "priority": "P2",
                }
            )
            continue
        nb = json.loads(path.read_text())
        cells = nb.get("cells", [])
        for idx, cell in enumerate(cells):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            calls = _extract_calls(source)
            relevant_calls = [
                call
                for call in calls
                if any(
                    token in call.lower()
                    for token in (
                        "qc",
                        "preprocess",
                        "hvg",
                        "pca",
                        "neighbor",
                        "umap",
                        "leiden",
                        "batch",
                        "integrat",
                        "scale",
                        "regress",
                        "embedding",
                    )
                )
            ]
            if not relevant_calls and not any(rule.pattern and re.search(rule.pattern, source) for rule in RULES):
                continue
            heading = _nearest_markdown_heading(cells, idx)
            for rule in _classify(source, relevant_calls):
                rows.append(
                    {
                        "notebook": path.name,
                        "cell_index": str(idx),
                        "heading": heading,
                        "block_name": rule.block_name,
                        "api_calls": ";".join(relevant_calls),
                        "coverage": rule.coverage,
                        "gap": rule.gap,
                        "owner": rule.owner,
                        "priority": rule.priority,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, str]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    matrix_path = OUT_DIR / "real_project_gap_matrix.tsv"
    fields = [
        "notebook",
        "cell_index",
        "heading",
        "block_name",
        "api_calls",
        "coverage",
        "gap",
        "owner",
        "priority",
    ]
    with matrix_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    by_block: dict[str, int] = {}
    by_priority: dict[str, int] = {}
    for row in rows:
        by_block[row["block_name"]] = by_block.get(row["block_name"], 0) + 1
        by_priority[row["priority"]] = by_priority.get(row["priority"], 0) + 1

    report_path = OUT_DIR / "qc_preprocess_evidence_report.md"
    with report_path.open("w") as handle:
        handle.write("# Real-project QC + preprocess gap report\n\n")
        handle.write("Generated from local Step1 notebooks without executing project data.\n\n")
        handle.write("## Priority counts\n\n")
        for key in sorted(by_priority):
            handle.write(f"- {key}: {by_priority[key]}\n")
        handle.write("\n## Block counts\n\n")
        for key, value in sorted(by_block.items()):
            handle.write(f"- {key}: {value}\n")
        handle.write("\n## Recommended next module work\n\n")
        handle.write("1. Execute the updated real-project notebooks and keep their sidecar JSON/TSV review outputs.\n")
        handle.write("2. Run `build_review_summary_evidence_tables.py` on final `.h5ad` files or result directories to build comparable evidence tables.\n")
        handle.write("3. Compare `scl.recommendation` advice, canonical iterative preprocessing evidence, and project-specific notebook decisions.\n")


def main() -> None:
    rows = build_rows(DEFAULT_NOTEBOOKS)
    write_outputs(rows)
    print(f"Wrote {len(rows)} rows to {OUT_DIR}")


if __name__ == "__main__":
    main()
