"""Optional external ambient-RNA correction backends.

This module provides thin wrappers around Python-based ambient RNA correction
tools such as CellBender.  It is intentionally separate from :mod:`ambient` so
that scLucid's core QC path remains dependency-free.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from anndata import AnnData

from ..utils import sanitize_for_hdf5
from .ambient import (
    AMBIENT_CORRECTED_COUNTS_LAYER,
    infer_ambient_input_context,
    record_ambient_correction_status,
    record_ambient_layer_contract,
)

log = logging.getLogger(__name__)


def _check_cellbender_cli() -> Optional[str]:
    """Return the cellbender executable path if available, else None."""
    return shutil.which("cellbender")


def _check_cellbender_module() -> bool:
    """Return True if the cellbender Python package can be imported."""
    try:
        import cellbender  # noqa: F401

        return True
    except Exception:
        return False


def _estimate_expected_cells(total_barcodes: int) -> int:
    """Heuristic for CellBender's --expected-cells argument.

    Uses the knee of a crude barcode-rank approximation: cells above 10% of
    the maximum total counts.  This is intentionally conservative and should
    be reviewed by the user for real datasets.
    """
    return max(100, int(total_barcodes * 0.1))


def _estimate_total_droplets(total_barcodes: int, expected_cells: int) -> int:
    """Heuristic for CellBender's --total-droplets-included argument."""
    return max(expected_cells + 100, min(total_barcodes, expected_cells * 5))


def run_cellbender(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    expected_cells: Optional[int] = None,
    total_droplets_included: Optional[int] = None,
    fpr: float = 0.01,
    epochs: int = 150,
    remove_checkpoint: bool = True,
    record: bool = True,
    key_added: str = "ambient_correction_summary",
) -> Dict[str, Any]:
    """Run CellBender ``remove-background`` on an AnnData object.

    This is a thin CLI wrapper.  It writes a temporary ``.h5ad`` file, runs
    CellBender, reads the corrected matrix from the output, and stores it in
    ``adata.layers[output_layer]``.

    Parameters
    ----------
    adata
        AnnData object containing raw counts.  The matrix should be in
        ``adata.layers[layer]`` or ``adata.X``.
    output_layer
        Layer name for the corrected matrix.
    expected_cells
        Passed to CellBender ``--expected-cells``.  Estimated heuristically if
        not provided.
    total_droplets_included
        Passed to CellBender ``--total-droplets-included``.  Estimated
        heuristically if not provided.
    fpr
        CellBender ``--fpr`` value.
    epochs
        CellBender ``--epochs`` value.
    remove_checkpoint
        If True, delete the temporary input/output files after reading results.
    record
        Store the correction summary in ``adata.uns['sclucid']['qc']``.
    key_added
        Key used when ``record=True``.

    Returns:
    -------
    dict
        Correction summary with ``corrected``, ``output_layer``, ``backend``,
        ``command``, ``removed_counts``, ``risk_note``.

    Raises:
    ------
    RuntimeError
        If CellBender is not installed or the command fails.
    """
    if layer is not None and layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

    # Materialise as a CSR-ish float matrix for storage; h5ad handles both.
    adata.layers["_cellbender_input"] = matrix

    total_barcodes = int(adata.n_obs)
    expected = expected_cells or _estimate_expected_cells(total_barcodes)
    total_droplets = total_droplets_included or _estimate_total_droplets(
        total_barcodes, expected
    )

    with tempfile.TemporaryDirectory(prefix="sclucid_cellbender_") as tmpdir:
        tmp_path = Path(tmpdir)
        input_h5ad = tmp_path / "input.h5ad"
        output_h5ad = tmp_path / "output.h5ad"

        adata.write_h5ad(str(input_h5ad))

        cellbender_exe = _check_cellbender_cli()
        if cellbender_exe is None:
            raise RuntimeError(
                "CellBender CLI (cellbender) not found. "
                "Install with 'pip install cellbender' and ensure it is on PATH."
            )

        command = [
            cellbender_exe,
            "remove-background",
            "--input",
            str(input_h5ad),
            "--output",
            str(output_h5ad),
            "--expected-cells",
            str(expected),
            "--total-droplets-included",
            str(total_droplets),
            "--fpr",
            str(fpr),
            "--epochs",
            str(epochs),
        ]

        log.info(f"Running CellBender: {' '.join(command)}")
        try:
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"CellBender failed with exit code {exc.returncode}.\n{exc.output}"
            ) from exc

        corrected_adata = AnnData.read_h5ad(str(output_h5ad))
        if not corrected_adata.obs_names.equals(adata.obs_names):
            raise ValueError("CellBender output obs_names do not match input")
        if not corrected_adata.var_names.equals(adata.var_names):
            raise ValueError("CellBender output var_names do not match input")

        # CellBender stores the corrected counts in .X by default.
        corrected_X = corrected_adata.X
        adata.layers[output_layer] = corrected_X.copy()

        removed_counts = float(np.asarray(matrix.sum() - corrected_X.sum()))

        summary = {
            "corrected": True,
            "output_layer": output_layer,
            "backend": "cellbender",
            "method": "cellbender_remove_background",
            "command": " ".join(command),
            "n_cells": int(adata.n_obs),
            "expected_cells": expected,
            "total_droplets_included": total_droplets,
            "fpr": fpr,
            "epochs": epochs,
            "removed_counts": removed_counts,
            "removed_fraction": float(
                removed_counts / max(float(np.asarray(matrix.sum())), 1.0)
            ),
            "review_required": removed_counts > float(np.asarray(matrix.sum())) * 0.25,
            "risk_note": (
                "CellBender correction removed a large fraction of counts; inspect output."
                if removed_counts > float(np.asarray(matrix.sum())) * 0.25
                else "CellBender correction applied."
            ),
            "stdout_tail": "\n".join(result.stdout.splitlines()[-20:]),
        }

        if remove_checkpoint:
            adata.layers.pop("_cellbender_input", None)

        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                summary
            )
            record_ambient_correction_status(
                adata,
                corrected=True,
                backend="cellbender",
                output_layer=output_layer,
                details=summary,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(
                    adata,
                    layer=layer,
                    matrix_source="raw_like",
                ),
                correction_summary=summary,
                output_layer=output_layer,
            )
        return summary


def cellbender_available() -> bool:
    """Return True if CellBender can be invoked (CLI or Python module)."""
    return _check_cellbender_cli() is not None or _check_cellbender_module()


def correct_ambient_rna(
    adata: AnnData,
    *,
    method: str = "auto",
    backend: str = "auto",
    layer: Optional[str] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Unified ambient RNA correction entry point.

    Parameters
    ----------
    adata
        AnnData object with raw counts.
    method
        ``"linear"`` for the zero-dependency linear subtraction,
        ``"external"`` to attempt CellBender, or ``"auto"`` to choose based on
        ambient risk and backend availability.
    backend
        External backend to use.  ``"auto"`` tries CellBender.
    layer
        Layer to correct; defaults to ``adata.X``.
    output_layer
        Layer name for corrected counts.
    **kwargs
        Additional arguments forwarded to the selected backend or to
        :func:`correct_ambient_rna_linear`.

    Returns:
    -------
    dict
        Correction summary.  On backend failure, falls back to linear
        correction and adds ``backend_fallback`` to the summary.
    """
    # Avoid a circular import: ambient_backends is imported by ambient, but the
    # convenience entry point lives here and can import from ambient directly.
    from .ambient import correct_ambient_rna_linear, diagnose_ambient_rna

    if method == "auto":
        risk = diagnose_ambient_rna(adata, layer=layer)
        use_external = backend == "auto" and cellbender_available() and risk.get(
            "risk_level", "low"
        ) in {"moderate", "high"}
        method = "external" if use_external else "linear"

    if method == "external":
        if backend in ("auto", "cellbender"):
            if cellbender_available():
                try:
                    return run_cellbender(
                        adata,
                        layer=layer,
                        output_layer=output_layer,
                        **kwargs,
                    )
                except Exception as exc:
                    log.warning(f"CellBender failed: {exc}; falling back to linear correction.")
                    linear = correct_ambient_rna_linear(
                        adata,
                        layer=layer,
                        output_layer=output_layer,
                        **kwargs,
                    )
                    linear["backend_fallback"] = {
                        "backend": "cellbender",
                        "error": str(exc),
                    }
                    return linear
            else:
                log.warning(
                    "CellBender not available; falling back to linear ambient correction."
                )
        else:
            raise ValueError(f"Unknown ambient backend: {backend!r}")

    if method == "linear":
        return correct_ambient_rna_linear(
            adata,
            layer=layer,
            output_layer=output_layer,
            **kwargs,
        )

    raise ValueError(f"Unknown ambient correction method: {method!r}")
