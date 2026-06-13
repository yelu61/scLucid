#!/usr/bin/env python3
"""Audit scLucid public API and regenerate FUNCTION_INVENTORY.md.

This script uses the stdlib ``ast`` module to scan ``src/scLucid`` without
importing the package, so it works even when optional dependencies are not
installed. It discovers public symbols from ``__all__`` and the custom
``_export()`` helper used by several subpackages, then updates the developer
API inventory.
"""

from __future__ import annotations

import argparse
import ast
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src" / "scLucid"
DOCS_DEV = REPO_ROOT / "docs" / "dev"
INVENTORY_PATH = DOCS_DEV / "FUNCTION_INVENTORY.md"

GEN_START = "<!-- AUTO-GENERATED INVENTORY START -->"
GEN_END = "<!-- AUTO-GENERATED INVENTORY END -->"

STABLE_KIND_ORDER = [
    "workflow_orchestrator",
    "config",
    "class",
    "function",
    "alias",
    "constant",
    "trace",
]
TRIAGE_KIND_ORDER = [
    "deprecated",
    "uncertain",
    "private_but_exposed",
]

KIND_LABELS = {
    "workflow_orchestrator": "Workflow Orchestrator",
    "config": "Config Class",
    "class": "Class",
    "function": "Function",
    "alias": "Alias",
    "constant": "Constant",
    "trace": "Trace / Contract",
    "deprecated": "Deprecated",
    "uncertain": "Uncertain",
    "private_but_exposed": "Private-but-Exposed",
}


@dataclass
class Symbol:
    name: str
    module_path: str
    source_file: str
    kind: str
    origin: Optional[str] = None
    is_alias: bool = False
    is_deprecated: bool = False
    is_private_but_exposed: bool = False
    is_optional: bool = False
    notes: Optional[str] = None


def _is_all_upper(name: str) -> bool:
    return name.isupper()


def _is_config_name(name: str) -> bool:
    return name.endswith("Config") or name.endswith("Config")


def _looks_like_trace(name: str) -> bool:
    return (
        "_TRACE_" in name
        or "_SCHEMA_VERSION" in name
        or "_REQUIRED_" in name
        or "_STABLE_" in name
        or "_CONTRACT" in name
        or "_REVIEW_" in name
        or name.startswith("validate_")
        or name.startswith("build_")
        and ("review" in name or "trace" in name or "maturity" in name or "contract" in name)
        or name.startswith("summarize_")
        and "review" in name
        or name.startswith("enrich_")
        and "review" in name
        or name.startswith("get_")
        and "contract" in name
    )


def _is_workflow_orchestrator(name: str, origin: Optional[str], module_path: str) -> bool:
    if not name.startswith("run_"):
        return False
    if origin and ("workflow" in origin or "workflow" in module_path):
        return True
    return name in {
        "run_pipeline",
        "run_standard_qc",
        "run_advanced_qc",
        "run_preprocessing",
        "run_standard_analysis",
        "run_custom_analysis",
        "run_annotation",
        "run_tumor_analysis",
    }


def _extract_constant_str(node: ast.AST) -> Optional[str]:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _extract_str_list(node: ast.AST) -> List[str]:
    values: List[str] = []
    if isinstance(node, (ast.List, ast.Tuple)):
        for elt in node.elts:
            s = _extract_constant_str(elt)
            if s is not None:
                values.append(s)
    return values


def _find_deprecated_functions(node: ast.AST) -> Set[str]:
    """Find function/class definitions that emit deprecation warnings."""
    deprecated: Set[str] = set()
    for child in ast.walk(node):
        if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for stmt in ast.walk(child):
            if not isinstance(stmt, ast.Call):
                continue
            # warnings.warn(..., FutureWarning/DeprecationWarning)
            if isinstance(stmt.func, ast.Attribute) and stmt.func.attr == "warn":
                if isinstance(stmt.func.value, ast.Name) and stmt.func.value.id == "warnings":
                    for arg in stmt.args:
                        if isinstance(arg, ast.Name) and arg.id in (
                            "FutureWarning",
                            "DeprecationWarning",
                        ):
                            deprecated.add(child.name)
                            break
                        if isinstance(arg, ast.Attribute) and arg.attr in (
                            "FutureWarning",
                            "DeprecationWarning",
                        ):
                            deprecated.add(child.name)
                            break
            # String literal containing "deprecated" anywhere in warn args
            for arg in stmt.args:
                s = _extract_constant_str(arg)
                if s and "deprecated" in s.lower():
                    deprecated.add(child.name)
                    break
    return deprecated


def _extract_all_names(node: ast.AST) -> List[str]:
    names: List[str] = []
    for stmt in getattr(node, "body", []):
        if isinstance(stmt, ast.Assign):
            for target in stmt.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    names.extend(_extract_str_list(stmt.value))
        elif isinstance(stmt, ast.AugAssign):
            if isinstance(stmt.target, ast.Name) and stmt.target.id == "__all__":
                names.extend(_extract_str_list(stmt.value))
        elif isinstance(stmt, ast.Expr):
            expr = stmt.value
            if isinstance(expr, ast.Call):
                func = expr.func
                if (
                    isinstance(func, ast.Attribute)
                    and func.attr == "append"
                    and isinstance(func.value, ast.Name)
                    and func.value.id == "__all__"
                ):
                    if expr.args:
                        s = _extract_constant_str(expr.args[0])
                        if s is not None:
                            names.append(s)
                if (
                    isinstance(func, ast.Attribute)
                    and func.attr == "extend"
                    and isinstance(func.value, ast.Name)
                    and func.value.id == "__all__"
                ):
                    names.extend(_extract_str_list(expr.args[0]))
    # Also walk for any __all__.append outside the top-level body (e.g., inside if blocks)
    for stmt in ast.walk(node):
        if isinstance(stmt, ast.Call):
            func = stmt.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "append"
                and isinstance(func.value, ast.Name)
                and func.value.id == "__all__"
            ):
                if stmt.args:
                    s = _extract_constant_str(stmt.args[0])
                    if s is not None:
                        names.append(s)
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "extend"
                and isinstance(func.value, ast.Name)
                and func.value.id == "__all__"
            ):
                names.extend(_extract_str_list(stmt.args[0]))
    return names


def _extract_export_calls(
    node: ast.AST,
) -> Tuple[List[Tuple[str, List[str], bool]], Dict[str, str]]:
    """Return list of (module, names, optional) from _export() calls.

    Also returns a mapping of name -> module for every exported name.
    """
    exports: List[Tuple[str, List[str], bool]] = []
    name_to_module: Dict[str, str] = {}
    for stmt in ast.walk(node):
        if not isinstance(stmt, ast.Call):
            continue
        func = stmt.func
        if not isinstance(func, ast.Name) or func.id != "_export":
            continue
        if len(stmt.args) < 2:
            continue
        module = _extract_constant_str(stmt.args[0])
        if module is None:
            continue
        names = _extract_str_list(stmt.args[1])
        optional = False
        for kw in stmt.keywords:
            if kw.arg == "optional" and isinstance(kw.value, ast.Constant):
                optional = bool(kw.value.value)
        exports.append((module, names, optional))
        for name in names:
            name_to_module[name] = module
    return exports, name_to_module


def _extract_simple_aliases(node: ast.AST) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for stmt in ast.walk(node):
        if not isinstance(stmt, ast.Assign):
            continue
        if len(stmt.targets) != 1:
            continue
        target = stmt.targets[0]
        if not isinstance(target, ast.Name):
            continue
        if isinstance(stmt.value, ast.Name):
            aliases[target.id] = stmt.value.id
    return aliases


def _extract_import_aliases(node: ast.AST) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for stmt in ast.walk(node):
        if not isinstance(stmt, ast.ImportFrom):
            continue
        module = stmt.module or ""
        for alias in stmt.names:
            if alias.asname:
                aliases[alias.asname] = f"{module}.{alias.name}" if module else alias.name
    return aliases


def _extract_import_optional_aliases(node: ast.AST) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for stmt in ast.walk(node):
        if not isinstance(stmt, ast.Assign):
            continue
        if len(stmt.targets) != 1 or not isinstance(stmt.targets[0], ast.Name):
            continue
        target = stmt.targets[0].id
        value = stmt.value
        if isinstance(value, ast.Call):
            if isinstance(value.func, ast.Name) and value.func.id == "_import_optional":
                if value.args:
                    s = _extract_constant_str(value.args[0])
                    if s is not None:
                        aliases[target] = s
    return aliases


def _extract_getattr_aliases(node: ast.AST) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    for stmt in ast.walk(node):
        if not isinstance(stmt, ast.Assign):
            continue
        if len(stmt.targets) != 1 or not isinstance(stmt.targets[0], ast.Name):
            continue
        target = stmt.targets[0].id
        value = stmt.value
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Attribute):
            if value.func.attr == "getattr":
                args = value.args
                if (
                    len(args) >= 2
                    and isinstance(args[1], ast.Constant)
                    and isinstance(args[1].value, str)
                ):
                    aliases[target] = args[1].value
    return aliases


def _classify_symbol(
    name: str,
    module_path: str,
    source_file: str,
    origin: Optional[str],
    is_alias: bool,
    is_deprecated: bool,
    is_private_but_exposed: bool,
    is_optional: bool,
    notes: Optional[str] = None,
) -> Symbol:
    if is_deprecated:
        kind = "deprecated"
    elif is_private_but_exposed:
        kind = "private_but_exposed"
    elif is_alias:
        kind = "alias"
    elif _is_workflow_orchestrator(name, origin, module_path):
        kind = "workflow_orchestrator"
    elif _is_config_name(name):
        kind = "config"
    elif _looks_like_trace(name):
        kind = "trace"
    elif _is_all_upper(name):
        kind = "constant"
    else:
        kind = "function"

    return Symbol(
        name=name,
        module_path=module_path,
        source_file=source_file,
        kind=kind,
        origin=origin,
        is_alias=is_alias,
        is_deprecated=is_deprecated,
        is_private_but_exposed=is_private_but_exposed,
        is_optional=is_optional,
        notes=notes,
    )


def _discover_subpackages() -> List[Tuple[str, Path]]:
    subpackages: List[Tuple[str, Path]] = []
    if not SRC_ROOT.exists():
        return subpackages
    for init_path in sorted(SRC_ROOT.rglob("__init__.py")):
        rel = init_path.relative_to(SRC_ROOT).parent
        if rel == Path("."):
            module_path = "scLucid"
        else:
            module_path = "scLucid." + ".".join(rel.parts)
        subpackages.append((module_path, init_path))
    return subpackages


def _parse_init(module_path: str, init_path: Path) -> List[Symbol]:
    text = init_path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        print(f"WARNING: syntax error in {init_path}: {exc}", file=sys.stderr)
        return []

    all_names = _extract_all_names(tree)
    deprecated_local = _find_deprecated_functions(tree)
    export_calls, name_to_export_module = _extract_export_calls(tree)
    simple_aliases = _extract_simple_aliases(tree)
    import_aliases = _extract_import_aliases(tree)
    optional_aliases = _extract_import_optional_aliases(tree)
    getattr_aliases = _extract_getattr_aliases(tree)

    # Combined alias resolution map: name -> descriptive origin
    alias_origins: Dict[str, str] = {}
    for name, target in simple_aliases.items():
        alias_origins[name] = target
    for name, target in import_aliases.items():
        alias_origins[name] = target
    for name, target in optional_aliases.items():
        alias_origins[name] = f"module {target}"
    for name, target in getattr_aliases.items():
        alias_origins[name] = f"module {target}"

    symbols: List[Symbol] = []
    seen: Set[str] = set()

    def add_symbol(
        name: str, origin: Optional[str], is_alias: bool, notes: Optional[str] = None
    ) -> None:
        if name in seen:
            return
        seen.add(name)
        is_deprecated = name in deprecated_local
        is_private_but_exposed = name.startswith("_")
        is_optional = name in name_to_export_module and name_to_export_module[name] in {
            mod for mod, _, opt in export_calls if opt
        }
        symbols.append(
            _classify_symbol(
                name=name,
                module_path=module_path,
                source_file=str(init_path.relative_to(REPO_ROOT)),
                origin=origin,
                is_alias=is_alias,
                is_deprecated=is_deprecated,
                is_private_but_exposed=is_private_but_exposed,
                is_optional=is_optional,
                notes=notes,
            )
        )

    # Process __all__ names
    for name in all_names:
        origin = None
        is_alias = False
        notes = None
        if name in alias_origins:
            origin = alias_origins[name]
            is_alias = True
        elif name in name_to_export_module:
            origin = name_to_export_module[name]
        elif module_path == "scLucid" and name in {
            "run_standard_qc",
            "run_advanced_qc",
            "run_preprocessing",
            "run_standard_analysis",
            "run_custom_analysis",
            "run_annotation",
            "characterize_clusters",
            "recommend_analysis_parameters",
            "run_tumor_analysis",
        }:
            origin = "dynamically resolved from submodule workflow"
        add_symbol(name, origin, is_alias, notes)

    # Process _export names that may not be in __all__ yet (safety net)
    for module, names, optional in export_calls:
        for name in names:
            if name in seen:
                continue
            origin = module
            notes = None
            if optional:
                notes = "optional"
            add_symbol(name, origin, False, notes)

    return symbols


def _group_symbols(symbols: List[Symbol]) -> Dict[str, List[Symbol]]:
    groups: Dict[str, List[Symbol]] = {}
    for sym in symbols:
        groups.setdefault(sym.module_path, []).append(sym)
    return groups


def _tags_for_symbol(sym: Symbol) -> List[str]:
    tags: List[str] = []
    if sym.is_alias:
        tags.append("[A]")
    if sym.kind == "config":
        tags.append("[C]")
    if sym.kind == "workflow_orchestrator":
        tags.append("[W]")
    if sym.kind == "trace":
        tags.append("[T]")
    if sym.is_deprecated:
        tags.append("[D]")
    if sym.is_private_but_exposed:
        tags.append("[P]")
    if sym.is_optional:
        tags.append("[O]")
    if sym.kind == "uncertain":
        tags.append("[?]")
    return tags


def _render_symbol_table(symbols: List[Symbol], include_kind_header: bool = True) -> List[str]:
    if not symbols:
        return ["*No symbols.*"]
    lines = [
        "| Symbol | Kind | Source | Notes |",
        "|--------|------|--------|-------|",
    ]
    for sym in symbols:
        tags = " ".join(_tags_for_symbol(sym))
        notes = sym.notes or ""
        if sym.origin:
            origin_note = f"from `{sym.origin}`" if " " not in sym.origin else sym.origin
            notes = f"{origin_note}; {notes}".strip("; ")
        if tags:
            notes = f"{tags} {notes}".strip()
        notes = notes.replace("|", "\\|")
        kind = KIND_LABELS.get(sym.kind, sym.kind)
        lines.append(f"| `{sym.name}` | {kind} | `{sym.source_file}` | {notes} |")
    return lines


def _render_subpackage_section(module_path: str, symbols: List[Symbol]) -> List[str]:
    stable = [s for s in symbols if s.kind not in TRIAGE_KIND_ORDER]
    triage = [s for s in symbols if s.kind in TRIAGE_KIND_ORDER]

    stable_by_kind: Dict[str, List[Symbol]] = {k: [] for k in STABLE_KIND_ORDER}
    for sym in stable:
        stable_by_kind.setdefault(sym.kind, []).append(sym)

    triage_by_kind: Dict[str, List[Symbol]] = {k: [] for k in TRIAGE_KIND_ORDER}
    for sym in triage:
        triage_by_kind.setdefault(sym.kind, []).append(sym)

    lines: List[str] = [f"## {module_path}", ""]

    # Stable table
    lines.append("### Stable APIs")
    lines.append("")
    has_stable = False
    for kind in STABLE_KIND_ORDER:
        group = stable_by_kind.get(kind, [])
        if not group:
            continue
        has_stable = True
        group.sort(key=lambda s: s.name.lower())
        lines.append(f"#### {KIND_LABELS.get(kind, kind)}")
        lines.append("")
        lines.extend(_render_symbol_table(group))
        lines.append("")
    if not has_stable:
        lines.append("*No stable public symbols.*")
        lines.append("")

    # Triage table
    lines.append("### Deprecated / Uncertain / Private-but-Exposed")
    lines.append("")
    has_triage = False
    for kind in TRIAGE_KIND_ORDER:
        group = triage_by_kind.get(kind, [])
        if not group:
            continue
        has_triage = True
        group.sort(key=lambda s: s.name.lower())
        lines.append(f"#### {KIND_LABELS.get(kind, kind)}")
        lines.append("")
        lines.extend(_render_symbol_table(group))
        lines.append("")
    if not has_triage:
        lines.append("*No flagged symbols.*")
        lines.append("")

    # Per-subpackage summary
    counts = dict.fromkeys(STABLE_KIND_ORDER + TRIAGE_KIND_ORDER, 0)
    for sym in symbols:
        counts[sym.kind] = counts.get(sym.kind, 0) + 1
    total = len(symbols)
    stable_total = sum(counts[k] for k in STABLE_KIND_ORDER)
    triage_total = sum(counts[k] for k in TRIAGE_KIND_ORDER)
    lines.append(
        f"**Summary:** {total} symbols ({stable_total} stable, {triage_total} flagged). "
        f"workflow={counts['workflow_orchestrator']}, config={counts['config']}, "
        f"class={counts['class']}, function={counts['function']}, alias={counts['alias']}, "
        f"constant={counts['constant']}, trace={counts['trace']}, "
        f"deprecated={counts['deprecated']}, uncertain={counts['uncertain']}, "
        f"private_but_exposed={counts['private_but_exposed']}."
    )
    lines.append("")
    return lines


def _generate_inventory(all_symbols: List[Symbol]) -> str:
    groups = _group_symbols(all_symbols)
    module_order = sorted(groups.keys(), key=lambda m: (m != "scLucid", m.lower()))

    lines: List[str] = [
        f"<!-- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by scripts/audit_public_api.py -->",
        f"<!-- Total public symbols: {len(all_symbols)} -->",
        "",
    ]

    totals = dict.fromkeys(STABLE_KIND_ORDER + TRIAGE_KIND_ORDER, 0)
    for sym in all_symbols:
        totals[sym.kind] = totals.get(sym.kind, 0) + 1

    for module_path in module_order:
        lines.extend(_render_subpackage_section(module_path, groups[module_path]))

    lines.append("## Global Summary")
    lines.append("")
    lines.append("| Kind | Count |")
    lines.append("|------|-------|")
    for kind in STABLE_KIND_ORDER + TRIAGE_KIND_ORDER:
        lines.append(f"| {KIND_LABELS.get(kind, kind)} | {totals[kind]} |")
    lines.append(f"| **Total** | **{len(all_symbols)}** |")
    lines.append("")

    return "\n".join(lines)


def _read_manual_sections(path: Path) -> Tuple[str, str]:
    if not path.exists():
        return "", ""
    content = path.read_text(encoding="utf-8")
    start_idx = content.find(GEN_START)
    end_idx = content.find(GEN_END)
    if start_idx == -1 or end_idx == -1:
        return content, ""
    header = content[:start_idx].strip()
    footer = content[end_idx + len(GEN_END) :].strip()
    return header, footer


def _write_inventory(path: Path, generated: str) -> None:
    DOCS_DEV.mkdir(parents=True, exist_ok=True)
    header, footer = _read_manual_sections(path)

    default_header = """# scLucid Public API Inventory

**Regenerate:** `python scripts/audit_public_api.py --write`

This document lists every public symbol in the scLucid API. The auto-generated
section below is maintained by `scripts/audit_public_api.py`. You may add notes
above the `AUTO-GENERATED` markers; they will be preserved across regenerations.

## Legend

| Tag | Meaning |
|-----|---------|
| `[A]` | Alias — points to another symbol or module |
| `[C]` | Config class |
| `[W]` | Workflow orchestrator |
| `[T]` | Trace / contract / schema constant |
| `[D]` | Deprecated |
| `[P]` | Private-but-exposed (starts with `_` but in `__all__`) |
| `[O]` | Optional / depends on extra dependencies |
| `[?]` | Uncertain / could not be traced cleanly |

## Maintainer Notes

- To add a symbol to the public API, add it to `__all__` in the subpackage `__init__.py`
  or register it via `_export()`.
- To deprecate a symbol, emit a `FutureWarning`/`DeprecationWarning`; the audit script
  will detect it automatically.
- To remove a symbol from the public API, remove it from `__all__` and run the audit
  script with `--write`.
"""

    if not header:
        header = default_header

    parts = [header, "", GEN_START, generated, GEN_END]
    if footer:
        parts.append("")
        parts.append(footer)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _check_inventory(path: Path, generated: str) -> int:
    if not path.exists():
        print(
            f"ERROR: {path.relative_to(REPO_ROOT)} does not exist. Run with --write.",
            file=sys.stderr,
        )
        return 1
    content = path.read_text(encoding="utf-8")
    start_idx = content.find(GEN_START)
    end_idx = content.find(GEN_END)
    if start_idx == -1 or end_idx == -1:
        print("ERROR: auto-generated delimiters not found in inventory file.", file=sys.stderr)
        return 1
    existing = content[start_idx + len(GEN_START) : end_idx].strip()
    if existing == generated.strip():
        print("OK: API inventory is up to date.")
        return 0
    print(
        "ERROR: API inventory is stale. Run 'python scripts/audit_public_api.py --write'",
        file=sys.stderr,
    )
    return 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit scLucid public API and manage FUNCTION_INVENTORY.md."
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Regenerate docs/dev/FUNCTION_INVENTORY.md.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit non-zero if the inventory is stale (for CI).",
    )
    args = parser.parse_args(argv)

    if not args.write and not args.check:
        parser.print_help()
        return 1

    subpackages = _discover_subpackages()
    all_symbols: List[Symbol] = []
    for module_path, init_path in subpackages:
        all_symbols.extend(_parse_init(module_path, init_path))

    generated = _generate_inventory(all_symbols)

    if args.write:
        _write_inventory(INVENTORY_PATH, generated)
        print(f"Wrote {INVENTORY_PATH.relative_to(REPO_ROOT)} ({len(all_symbols)} symbols).")
        return 0

    if args.check:
        return _check_inventory(INVENTORY_PATH, generated)

    return 0


if __name__ == "__main__":
    sys.exit(main())
