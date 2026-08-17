"""Optional external ambient-RNA correction backends.

This module provides thin wrappers around Python-based and R-based ambient RNA
correction tools such as CellBender, SoupX, and DecontX. It is intentionally
separate from :mod:`ambient` so that scLucid's core QC path remains
dependency-free.

The backend registry pattern makes it easy to add new ambient RNA correction
methods without changing the unified :func:`correct_ambient_rna` entry point.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
import tempfile
from functools import cache
from pathlib import Path
from typing import Any, Callable, Dict, Optional

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

# ---------------------------------------------------------------------------
# Backend availability checks
# ---------------------------------------------------------------------------

_R_CAPABILITY_PROBE_TIMEOUT_SECONDS = 15.0
_R_CAPABILITY_PROBE_CODE = """
import sys

try:
    import rpy2.robjects  # noqa: F401
    package = sys.argv[1] if len(sys.argv) > 1 else ""
    if package:
        from rpy2.robjects.packages import importr
        importr(package)
except BaseException as exc:
    print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
    raise SystemExit(1)
"""


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


def cellbender_available() -> bool:
    """Return True if CellBender can be invoked (CLI or Python module)."""
    return _check_cellbender_cli() is not None or _check_cellbender_module()


@cache
def _probe_r_capability(package: Optional[str] = None) -> bool:
    """Probe rpy2/R in a child process so embedded-R failures stay isolated.

    Importing :mod:`rpy2.robjects` initializes embedded R and can terminate the
    interpreter when the local R installation is missing or ABI-incompatible.
    A Python ``try`` block cannot catch that signal, so capability discovery
    must never perform the import in scLucid's main process.
    """
    command = [sys.executable, "-c", _R_CAPABILITY_PROBE_CODE]
    if package:
        command.append(package)
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=_R_CAPABILITY_PROBE_TIMEOUT_SECONDS,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        log.debug("R capability probe failed for %r: %s", package or "rpy2", exc)
        return False
    if result.returncode != 0:
        log.debug(
            "R capability probe reported unavailable for %r (return code %s).",
            package or "rpy2",
            result.returncode,
        )
        return False
    return True


def _rpy2_available() -> bool:
    """Return True when rpy2 and embedded R initialize in an isolated process."""
    return _probe_r_capability()


def _r_package_available(name: str) -> bool:
    """Return True if an R package can be imported in an isolated process."""
    return _probe_r_capability(name)


def soupx_available() -> bool:
    """Return True if SoupX (R package) is available."""
    return _r_package_available("SoupX")


def decontx_available() -> bool:
    """Return True if celda/DecontX (R package) is available."""
    return _r_package_available("celda")


# ---------------------------------------------------------------------------
# CellBender backend
# ---------------------------------------------------------------------------


def _estimate_expected_cells(total_barcodes: int) -> int:
    """Heuristic for CellBender's --expected-cells argument."""
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
    """Run CellBender ``remove-background`` on an AnnData object."""
    if layer is not None and layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X

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


# ---------------------------------------------------------------------------
# SoupX backend (R bridge)
# ---------------------------------------------------------------------------


def run_soupx(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    clusters: Optional[str] = None,
    record: bool = True,
    key_added: str = "ambient_correction_summary",
) -> Dict[str, Any]:
    """Run SoupX ambient RNA correction via rpy2.

    SoupX works best when provided with a raw-feature matrix and preliminary
    clustering. If ``clusters`` is provided, it is used as the SoupX cluster
    argument; otherwise SoupX estimates clusters internally.
    """
    if not soupx_available():
        raise RuntimeError(
            "SoupX is not available. Install the R package with:\n"
            "  install.packages('SoupX')\n"
            "and ensure rpy2 is installed in Python."
        )

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    matrix = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    adata.layers["_soupx_input"] = matrix

    importr("SoupX")
    base = importr("base")

    with tempfile.TemporaryDirectory(prefix="sclucid_soupx_") as tmpdir:
        tmp_path = Path(tmpdir)
        input_h5ad = tmp_path / "input.h5ad"
        adata.write_h5ad(str(input_h5ad))

        ro.globalenv["input_h5ad"] = str(input_h5ad)
        ro.globalenv["output_h5ad"] = str(tmp_path / "output.h5ad")
        if clusters is not None and clusters in adata.obs.columns:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.globalenv["clusters"] = adata.obs[clusters].astype(str)
        else:
            ro.globalenv["clusters"] = ro.NULL

        r_code = """
        library(SoupX)
        library(anndata)
        # anndata R package is experimental; fall back to reticulate if needed.
        if (requireNamespace("anndata", quietly = TRUE)) {
          ad <- anndata::read_h5ad(input_h5ad)
        } else if (requireNamespace("reticulate", quietly = TRUE)) {
          ad <- reticulate::py_anndata(input_h5ad)
        } else {
          stop("SoupX backend requires the R anndata or reticulate package.")
        }
        sc <- SoupX::setDR(ad, clusters)
        sc <- SoupX::calculateContaminationFraction(sc)
        out <- SoupX::adjustCounts(sc)
        # Write corrected counts back via reticulate/anndata
        ad$layers[["ambient_corrected_counts"]] <- out
        anndata::write_h5ad(ad, output_h5ad)
        list(removed_counts = sum(ad$X - out), n_cells = nrow(out))
        """
        try:
            r_result = base.eval(ro.r(r_code))
        except Exception as exc:
            raise RuntimeError(f"SoupX correction failed: {exc}") from exc

        corrected_adata = AnnData.read_h5ad(str(tmp_path / "output.h5ad"))
        if output_layer not in corrected_adata.layers:
            raise RuntimeError(f"SoupX did not write layer '{output_layer}'")
        adata.layers[output_layer] = corrected_adata.layers[output_layer].copy()

        removed_counts = float(np.asarray(r_result.rx2("removed_counts")[0]))
        summary = {
            "corrected": True,
            "output_layer": output_layer,
            "backend": "soupx",
            "method": "soupx_adjust_counts",
            "n_cells": int(adata.n_obs),
            "clusters_used": clusters,
            "removed_counts": removed_counts,
            "removed_fraction": float(
                removed_counts / max(float(np.asarray(matrix.sum())), 1.0)
            ),
            "review_required": removed_counts > float(np.asarray(matrix.sum())) * 0.25,
            "risk_note": (
                "SoupX correction removed a large fraction of counts; inspect output."
                if removed_counts > float(np.asarray(matrix.sum())) * 0.25
                else "SoupX correction applied."
            ),
        }

        adata.layers.pop("_soupx_input", None)

        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                summary
            )
            record_ambient_correction_status(
                adata,
                corrected=True,
                backend="soupx",
                output_layer=output_layer,
                details=summary,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(
                    adata,
                    layer=layer,
                    matrix_source="filtered_like" if clusters else "raw_like",
                ),
                correction_summary=summary,
                output_layer=output_layer,
            )
        return summary


# ---------------------------------------------------------------------------
# DecontX backend (R bridge)
# ---------------------------------------------------------------------------


def run_decontx(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    batch_key: Optional[str] = None,
    max_iter: int = 500,
    record: bool = True,
    key_added: str = "ambient_correction_summary",
) -> Dict[str, Any]:
    """Run DecontX ambient RNA correction via rpy2.

    DecontX is designed for filtered-feature matrices where empty droplets have
    already been removed. It estimates ambient contamination using a
    variational Bayesian model.
    """
    if not decontx_available():
        raise RuntimeError(
            "DecontX (celda package) is not available. Install with:\n"
            "  BiocManager::install('celda')\n"
            "and ensure rpy2 is installed in Python."
        )

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    matrix = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    adata.layers["_decontx_input"] = matrix

    importr("celda")
    base = importr("base")

    with tempfile.TemporaryDirectory(prefix="sclucid_decontx_") as tmpdir:
        tmp_path = Path(tmpdir)
        input_h5ad = tmp_path / "input.h5ad"
        output_h5ad = tmp_path / "output.h5ad"
        adata.write_h5ad(str(input_h5ad))

        ro.globalenv["input_h5ad"] = str(input_h5ad)
        ro.globalenv["output_h5ad"] = str(output_h5ad)
        ro.globalenv["max_iter"] = max_iter
        if batch_key is not None and batch_key in adata.obs.columns:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.globalenv["batch"] = adata.obs[batch_key].astype(str)
        else:
            ro.globalenv["batch"] = ro.NULL

        r_code = """
        library(celda)
        library(anndata)
        if (requireNamespace("anndata", quietly = TRUE)) {
          ad <- anndata::read_h5ad(input_h5ad)
        } else if (requireNamespace("reticulate", quietly = TRUE)) {
          ad <- reticulate::py_anndata(input_h5ad)
        } else {
          stop("DecontX backend requires the R anndata or reticulate package.")
        }
        counts <- t(ad$X)
        decont <- celda::decontX(counts, batch = batch, maxIter = max_iter)
        corrected <- t(decont$decontXcounts)
        ad$layers[["ambient_corrected_counts"]] <- corrected
        anndata::write_h5ad(ad, output_h5ad)
        list(removed_counts = sum(ad$X - corrected), n_cells = nrow(corrected))
        """
        try:
            r_result = base.eval(ro.r(r_code))
        except Exception as exc:
            raise RuntimeError(f"DecontX correction failed: {exc}") from exc

        corrected_adata = AnnData.read_h5ad(str(output_h5ad))
        if output_layer not in corrected_adata.layers:
            raise RuntimeError(f"DecontX did not write layer '{output_layer}'")
        adata.layers[output_layer] = corrected_adata.layers[output_layer].copy()

        removed_counts = float(np.asarray(r_result.rx2("removed_counts")[0]))
        summary = {
            "corrected": True,
            "output_layer": output_layer,
            "backend": "decontx",
            "method": "decontx_variational",
            "n_cells": int(adata.n_obs),
            "batch_key": batch_key,
            "max_iter": max_iter,
            "removed_counts": removed_counts,
            "removed_fraction": float(
                removed_counts / max(float(np.asarray(matrix.sum())), 1.0)
            ),
            "review_required": removed_counts > float(np.asarray(matrix.sum())) * 0.25,
            "risk_note": (
                "DecontX correction removed a large fraction of counts; inspect output."
                if removed_counts > float(np.asarray(matrix.sum())) * 0.25
                else "DecontX correction applied."
            ),
        }

        adata.layers.pop("_decontx_input", None)

        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                summary
            )
            record_ambient_correction_status(
                adata,
                corrected=True,
                backend="decontx",
                output_layer=output_layer,
                details=summary,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(
                    adata,
                    layer=layer,
                    matrix_source="filtered_like",
                ),
                correction_summary=summary,
                output_layer=output_layer,
            )
        return summary


# ---------------------------------------------------------------------------
# Backend registry and unified entry point
# ---------------------------------------------------------------------------

_BackendFn = Callable[..., Dict[str, Any]]
_BACKEND_REGISTRY: Dict[str, Dict[str, Any]] = {
    "cellbender": {
        "run": run_cellbender,
        "available": cellbender_available,
        "matrix_types": {"raw_like"},
        "python": True,
        "r": False,
        "experimental": False,
        "auto_select": True,
        "recommended_for_filtered_matrix": False,
    },
    "soupx": {
        "run": run_soupx,
        "available": soupx_available,
        "matrix_types": {"raw_like", "filtered_like"},
        "python": False,
        "r": True,
        "experimental": True,
        "auto_select": False,
        "recommended_for_filtered_matrix": False,
    },
    "decontx": {
        "run": run_decontx,
        "available": decontx_available,
        "matrix_types": {"filtered_like"},
        "python": False,
        "r": True,
        "experimental": True,
        "auto_select": False,
        "recommended_for_filtered_matrix": False,
    },
}


def list_ambient_backends() -> Dict[str, Dict[str, Any]]:
    """Return ambient backend metadata with live availability flags."""
    return {
        name: {
            **meta,
            "available_now": bool(meta["available"]()),
        }
        for name, meta in _BACKEND_REGISTRY.items()
    }


def _choose_backend(
    matrix_type: str,
    backend: str,
    risk_level: str,
) -> Optional[str]:
    """Choose a concrete backend given matrix type and availability.

    Auto-selection only considers backends marked ``auto_select=True``.
    R-based backends (SoupX, DecontX) are kept available for explicit calls
    but are not chosen automatically because they require fragile rpy2/R setup.
    """
    if backend != "auto" and backend in _BACKEND_REGISTRY:
        chosen = backend
    elif backend == "auto":
        candidates = [
            name
            for name, meta in _BACKEND_REGISTRY.items()
            if meta.get("auto_select", False)
            and matrix_type in meta.get("matrix_types", set())
        ]
        # Preserve the historical preference order within auto-selectable backends.
        order = {"cellbender": 0}
        candidates = sorted(candidates, key=lambda name: order.get(name, 1))
        chosen = None
        for cand in candidates:
            if _BACKEND_REGISTRY[cand]["available"]():
                chosen = cand
                break
    else:
        raise ValueError(f"Unknown ambient backend: {backend!r}")

    if chosen is None:
        return None

    meta = _BACKEND_REGISTRY[chosen]
    if matrix_type not in meta["matrix_types"]:
        log.warning(
            "Backend '%s' is not designed for '%s' matrices; review before using.",
            chosen,
            matrix_type,
        )
    if meta.get("experimental", False):
        log.warning(
            "Backend '%s' is experimental and requires additional dependencies; "
            "results should be reviewed carefully.",
            chosen,
        )
    if not meta["available"]():
        return None
    return chosen


def _record_diagnostic_only_ambient_result(
    adata: AnnData,
    *,
    layer: Optional[str],
    output_layer: str,
    reason: str,
    backend_requested: str,
    key_added: str = "ambient_correction_summary",
) -> Dict[str, Any]:
    """Record that ambient correction was intentionally deferred for review."""
    from .ambient import diagnose_ambient_rna

    context = infer_ambient_input_context(adata, layer=layer)
    diagnostic = diagnose_ambient_rna(adata, layer=layer)
    summary = {
        "corrected": False,
        "output_layer": output_layer,
        "backend": None,
        "backend_requested": backend_requested,
        "method": "diagnostic_only",
        "reason": reason,
        "matrix_type": context.get("matrix_type"),
        "diagnostic": diagnostic,
        "review_required": True,
        "risk_note": (
            "Filtered matrices do not contain the empty-droplet background needed "
            "for CellBender-style modeling. Prefer external DecontX/SoupX-like "
            "correction and register the result, or explicitly request method='linear' "
            "for conservative background subtraction."
        ),
        "recommendation": (
            "Run a validated filtered-matrix ambient workflow externally when correction "
            "is needed, then call register_external_ambient_result."
        ),
    }
    if key_added:
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
            summary
        )
    record_ambient_correction_status(
        adata,
        corrected=False,
        backend="diagnostic_only",
        output_layer=None,
        details=summary,
    )
    record_ambient_layer_contract(
        adata,
        input_context=context,
        correction_summary=summary,
        output_layer=output_layer,
    )
    return summary


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
        AnnData object with raw or filtered counts.
    method
        ``"linear"`` for the zero-dependency linear subtraction,
        ``"external"`` to attempt a registered backend, or ``"auto"`` to choose
        based on ambient risk and backend availability.
    backend
        External backend to use. ``"auto"`` picks a suitable backend based on the
        matrix type (raw vs filtered). Registered backends: ``cellbender``,
        ``soupx``, ``decontx``.
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
        Correction summary. On backend failure, falls back to linear correction
        and adds ``backend_fallback`` to the summary.
    """
    from .ambient import correct_ambient_rna_linear, diagnose_ambient_rna

    context = infer_ambient_input_context(adata, layer=layer)
    matrix_type = context.get("matrix_type", "filtered_like")

    if method == "auto":
        risk = diagnose_ambient_rna(adata, layer=layer)
        if matrix_type == "filtered_like":
            chosen = (
                _choose_backend(matrix_type, backend, risk.get("risk_level", "low"))
                if risk.get("risk_level", "low") in {"moderate", "high"}
                else None
            )
            if chosen is None:
                return _record_diagnostic_only_ambient_result(
                    adata,
                    layer=layer,
                    output_layer=output_layer,
                    reason="filtered_matrix_requires_external_or_explicit_linear_correction",
                    backend_requested=backend,
                )
            method = "external"
            backend = chosen
        else:
            use_external = risk.get("risk_level", "low") in {"moderate", "high"}
            method = "external" if use_external else "linear"

    if method == "external":
        chosen = _choose_backend(matrix_type, backend, "moderate")
        if chosen is None:
            if matrix_type == "filtered_like":
                return _record_diagnostic_only_ambient_result(
                    adata,
                    layer=layer,
                    output_layer=output_layer,
                    reason="no_filtered_matrix_backend_available",
                    backend_requested=backend,
                )
            log.warning(
                "No suitable ambient backend available for %s matrix; "
                "falling back to linear correction.",
                matrix_type,
            )
            return correct_ambient_rna_linear(
                adata,
                layer=layer,
                output_layer=output_layer,
                **kwargs,
            )

        try:
            return _BACKEND_REGISTRY[chosen]["run"](
                adata,
                layer=layer,
                output_layer=output_layer,
                **kwargs,
            )
        except Exception as exc:
            log.warning(f"{chosen} failed: {exc}; falling back to linear correction.")
            linear = correct_ambient_rna_linear(
                adata,
                layer=layer,
                output_layer=output_layer,
                **kwargs,
            )
            linear["backend_fallback"] = {
                "backend": chosen,
                "error": str(exc),
            }
            return linear

    if method == "linear":
        return correct_ambient_rna_linear(
            adata,
            layer=layer,
            output_layer=output_layer,
            **kwargs,
        )

    raise ValueError(f"Unknown ambient correction method: {method!r}")
