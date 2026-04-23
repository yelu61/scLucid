"""
Backend abstraction layer for preprocessing operations.

This module provides a plugin architecture that allows swapping the underlying
implementation of preprocessing operations (e.g., scanpy vs rapids-singlecell).
"""

import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Union

from anndata import AnnData

log = logging.getLogger(__name__)


class PreprocessingBackend(ABC):
    """
    Abstract base class for preprocessing backends.

    Implementations of this class provide the actual computation for
    preprocessing operations. The default implementation uses scanpy,
    but alternative backends (e.g., GPU-accelerated) can be implemented.
    """

    name: str = "abstract"

    @abstractmethod
    def normalize_total(
        self,
        adata: AnnData,
        target_sum: float = 1e4,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
        **kwargs,
    ) -> None:
        """Normalize counts per cell."""
        pass

    @abstractmethod
    def log1p(self, adata: AnnData, **kwargs) -> None:
        """Log-transform data."""
        pass

    @abstractmethod
    def scale(self, adata: AnnData, max_value: Optional[float] = None, **kwargs) -> None:
        """Scale data to unit variance and zero mean."""
        pass

    @abstractmethod
    def pca(self, adata: AnnData, n_comps: int = 50, **kwargs) -> None:
        """Run PCA."""
        pass

    @abstractmethod
    def neighbors(
        self,
        adata: AnnData,
        n_neighbors: int = 15,
        n_pcs: Optional[int] = None,
        use_rep: Optional[str] = None,
        **kwargs,
    ) -> None:
        """Compute nearest neighbors."""
        pass

    @abstractmethod
    def umap(self, adata: AnnData, **kwargs) -> None:
        """Run UMAP."""
        pass

    @abstractmethod
    def highly_variable_genes(
        self,
        adata: AnnData,
        n_top_genes: Optional[int] = None,
        flavor: str = "seurat_v3",
        batch_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        """Identify highly variable genes."""
        pass

    @abstractmethod
    def regress_out(self, adata: AnnData, keys: List[str], **kwargs) -> None:
        """Regress out covariates."""
        pass


class ScanpyBackend(PreprocessingBackend):
    """
    Default scanpy-based preprocessing backend.
    """

    name = "scanpy"

    def __init__(self):
        import scanpy as sc

        self.sc = sc

    def normalize_total(
        self,
        adata: AnnData,
        target_sum: float = 1e4,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
        **kwargs,
    ) -> None:
        self.sc.pp.normalize_total(
            adata,
            target_sum=target_sum,
            exclude_highly_expressed=exclude_highly_expressed,
            max_fraction=max_fraction,
            **kwargs,
        )

    def log1p(self, adata: AnnData, **kwargs) -> None:
        self.sc.pp.log1p(adata, **kwargs)

    def scale(self, adata: AnnData, max_value: Optional[float] = None, **kwargs) -> None:
        self.sc.pp.scale(adata, max_value=max_value, **kwargs)

    def pca(self, adata: AnnData, n_comps: int = 50, **kwargs) -> None:
        self.sc.tl.pca(adata, n_comps=n_comps, **kwargs)

    def neighbors(
        self,
        adata: AnnData,
        n_neighbors: int = 15,
        n_pcs: Optional[int] = None,
        use_rep: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep, **kwargs)

    def umap(self, adata: AnnData, **kwargs) -> None:
        self.sc.tl.umap(adata, **kwargs)

    def highly_variable_genes(
        self,
        adata: AnnData,
        n_top_genes: Optional[int] = None,
        flavor: str = "seurat_v3",
        batch_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes, flavor=flavor, batch_key=batch_key, **kwargs
        )

    def regress_out(self, adata: AnnData, keys: List[str], **kwargs) -> None:
        self.sc.pp.regress_out(adata, keys, **kwargs)


class RapidsBackend(PreprocessingBackend):
    """
    GPU-accelerated backend using rapids-singlecell.

    .. warning::
        This backend is **experimental**. The rapids-singlecell ecosystem is
        small and the wrapper may not be fully tested on your platform.  For
        production work, prefer ``ScanpyBackend``.

    Requires rapids-singlecell to be installed.
    """

    name = "rapids"
    _experimental = True

    def __init__(self):
        try:
            import rapids_singlecell as rsc

            self.rsc = rsc
        except ImportError:
            raise ImportError(
                "rapids-singlecell not installed. " "Install with: pip install rapids-singlecell"
            )

    def normalize_total(
        self,
        adata: AnnData,
        target_sum: float = 1e4,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
        **kwargs,
    ) -> None:
        # rapids-singlecell API may differ, this is illustrative
        self.rsc.pp.normalize_total(adata, target_sum=target_sum, **kwargs)

    def log1p(self, adata: AnnData, **kwargs) -> None:
        self.rsc.pp.log1p(adata, **kwargs)

    def scale(self, adata: AnnData, max_value: Optional[float] = None, **kwargs) -> None:
        self.rsc.pp.scale(adata, max_value=max_value, **kwargs)

    def pca(self, adata: AnnData, n_comps: int = 50, **kwargs) -> None:
        self.rsc.tl.pca(adata, n_comps=n_comps, **kwargs)

    def neighbors(
        self,
        adata: AnnData,
        n_neighbors: int = 15,
        n_pcs: Optional[int] = None,
        use_rep: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.rsc.pp.neighbors(
            adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep, **kwargs
        )

    def umap(self, adata: AnnData, **kwargs) -> None:
        self.rsc.tl.umap(adata, **kwargs)

    def highly_variable_genes(
        self,
        adata: AnnData,
        n_top_genes: Optional[int] = None,
        flavor: str = "seurat_v3",
        batch_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.rsc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes, flavor=flavor, batch_key=batch_key, **kwargs
        )

    def regress_out(self, adata: AnnData, keys: List[str], **kwargs) -> None:
        # rapids-singlecell may not support regress_out
        log.warning("rapids-singlecell may not support regress_out, falling back to scanpy")
        # Fall back to scanpy for unsupported operations
        fallback = ScanpyBackend()
        fallback.regress_out(adata, keys, **kwargs)


# Global backend instance
_backend: Optional[PreprocessingBackend] = None


def get_backend() -> PreprocessingBackend:
    """Get the current preprocessing backend."""
    global _backend
    if _backend is None:
        _backend = ScanpyBackend()  # Default to scanpy
    return _backend


def set_backend(backend: Union[str, PreprocessingBackend]) -> None:
    """
    Set the preprocessing backend.

    Args:
        backend: Either a string ("scanpy", "rapids") or a PreprocessingBackend instance.

    Example:
        >>> # Use default scanpy backend
        >>> set_backend("scanpy")

        >>> # Use custom backend
        >>> set_backend(MyCustomBackend())
    """
    global _backend

    if isinstance(backend, str):
        if backend == "scanpy":
            _backend = ScanpyBackend()
        elif backend == "rapids":
            log.warning(
                "RapidsBackend is experimental. "
                "Consider using ScanpyBackend for production analyses."
            )
            _backend = RapidsBackend()
        else:
            raise ValueError(f"Unknown backend: {backend}. Use 'scanpy' or 'rapids'.")
    elif isinstance(backend, PreprocessingBackend):
        _backend = backend
    else:
        raise TypeError(f"Backend must be string or PreprocessingBackend, got {type(backend)}")

    log.info(f"Set preprocessing backend to: {_backend.name}")


def list_available_backends() -> Dict[str, Union[bool, str]]:
    """
    List available backends and their installation status.

    Returns:
        Dict mapping backend name to availability. Values are:
        - ``True``: available and recommended
        - ``False``: not installed
        - ``"experimental"``: installed but not recommended for production
    """
    available: Dict[str, Union[bool, str]] = {"scanpy": True}  # Always available

    # Check rapids-singlecell
    try:
        import rapids_singlecell  # noqa: F401

        available["rapids"] = "experimental"
    except ImportError:
        available["rapids"] = False

    return available


__all__ = [
    "PreprocessingBackend",
    "ScanpyBackend",
    "RapidsBackend",
    "get_backend",
    "set_backend",
    "list_available_backends",
]
