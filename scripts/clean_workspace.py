"""Remove local cache and generated workspace artifacts.

This script is intentionally conservative. By default it removes only files
that are deterministic local products and are already ignored by git.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_FILE_NAMES = {
    ".DS_Store",
    ".coverage",
    "coverage.xml",
}

DEFAULT_DIR_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    "htmlcov",
    "override",
}

OPTIONAL_OUTPUT_DIRS = {
    "results",
    "test_results",
    "old_results",
    "legacy_results",
    "out",
    "test_output",
    "test",
    "test_save",
    "from_dict",
    "from_json",
}


def _iter_targets(include_outputs: bool) -> list[Path]:
    targets: list[Path] = []
    ignored_dirs = {".git", ".venv", "venv", "env", "ENV"}
    dir_names = set(DEFAULT_DIR_NAMES)
    if include_outputs:
        dir_names.update(OPTIONAL_OUTPUT_DIRS)

    for path in ROOT.rglob("*"):
        if any(part in ignored_dirs for part in path.relative_to(ROOT).parts):
            continue
        if path.is_file() and path.name in DEFAULT_FILE_NAMES:
            targets.append(path)
        elif path.is_dir() and path.name in dir_names:
            targets.append(path)

    for name in dir_names:
        path = ROOT / name
        if path.exists() and path not in targets:
            targets.append(path)

    return sorted(targets, key=lambda item: (len(item.parts), str(item)))


def clean_workspace(*, include_outputs: bool = False, dry_run: bool = False) -> list[Path]:
    """Clean deterministic local artifacts and return removed targets."""
    targets = _iter_targets(include_outputs=include_outputs)
    for path in targets:
        if dry_run:
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink(missing_ok=True)
    return targets


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-outputs",
        action="store_true",
        help="Also remove ignored local output directories such as results/.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print targets without deleting them.",
    )
    args = parser.parse_args()

    targets = clean_workspace(
        include_outputs=args.include_outputs,
        dry_run=args.dry_run,
    )
    action = "Would remove" if args.dry_run else "Removed"
    for path in targets:
        print(path.relative_to(ROOT))
    print(f"{action} {len(targets)} local artifact(s).")


if __name__ == "__main__":
    main()
