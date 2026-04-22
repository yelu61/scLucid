"""Unit tests for scLucid.preprocess.backend."""

import pytest

from scLucid.preprocess.backend import (
    ScanpyBackend,
    RapidsBackend,
    get_backend,
    set_backend,
    list_available_backends,
    PreprocessingBackend,
)


@pytest.mark.unit
class TestBackendBasics:
    """Tests for backend abstraction basics."""

    def test_get_backend_returns_scanpy_by_default(self):
        backend = get_backend()
        assert isinstance(backend, ScanpyBackend)
        assert backend.name == "scanpy"

    def test_set_backend_with_string_scanpy(self):
        set_backend("scanpy")
        backend = get_backend()
        assert isinstance(backend, ScanpyBackend)
        assert backend.name == "scanpy"

    def test_set_backend_with_instance(self):
        instance = ScanpyBackend()
        set_backend(instance)
        backend = get_backend()
        assert backend is instance

    def test_set_backend_invalid_string(self):
        with pytest.raises(ValueError, match="Unknown backend"):
            set_backend("nonexistent_backend")

    def test_list_available_backends(self):
        backends = list_available_backends()
        assert "scanpy" in backends
        # rapids may or may not be available depending on environment


@pytest.mark.unit
class TestScanpyBackend:
    """Tests for ScanpyBackend implementation."""

    def test_methods_exist(self):
        backend = ScanpyBackend()
        required_methods = [
            "normalize_total",
            "log1p",
            "scale",
            "pca",
            "neighbors",
            "umap",
            "highly_variable_genes",
            "regress_out",
        ]
        for method in required_methods:
            assert hasattr(backend, method)
            assert callable(getattr(backend, method))

    def test_normalize_total_runs(self, minimal_adata):
        backend = ScanpyBackend()
        adata = minimal_adata.copy()
        backend.normalize_total(adata, target_sum=1e4)
        # normalize_total modifies adata in place
        assert adata.X is not None

    def test_scale_runs(self, minimal_adata):
        backend = ScanpyBackend()
        adata = minimal_adata.copy()
        backend.scale(adata, max_value=10)
        assert adata.X is not None

    def test_log1p_runs(self, minimal_adata):
        backend = ScanpyBackend()
        adata = minimal_adata.copy()
        backend.log1p(adata)
        assert adata.X is not None


@pytest.mark.unit
class TestRapidsBackend:
    """Tests for RapidsBackend implementation (conditional on availability)."""

    def test_init_requires_rapids(self):
        try:
            backend = RapidsBackend()
            # If we get here, rapids is installed
            assert backend.name == "rapids"
        except ImportError:
            pytest.skip("rapids-singlecell not installed")

    def test_methods_exist_if_available(self):
        try:
            backend = RapidsBackend()
        except ImportError:
            pytest.skip("rapids-singlecell not installed")

        required_methods = [
            "normalize_total",
            "log1p",
            "scale",
            "pca",
            "neighbors",
            "umap",
            "highly_variable_genes",
            "regress_out",
        ]
        for method in required_methods:
            assert hasattr(backend, method)
            assert callable(getattr(backend, method))


@pytest.mark.unit
class TestBackendAbstract:
    """Tests for PreprocessingBackend abstract base class."""

    def test_cannot_instantiate_abstract(self):
        with pytest.raises(TypeError):
            PreprocessingBackend()

    def test_custom_backend_subclass(self):
        class DummyBackend(PreprocessingBackend):
            name = "dummy"

            def normalize_total(self, adata, **kwargs):
                pass

            def log1p(self, adata, **kwargs):
                pass

            def scale(self, adata, **kwargs):
                pass

            def pca(self, adata, **kwargs):
                pass

            def neighbors(self, adata, **kwargs):
                pass

            def umap(self, adata, **kwargs):
                pass

            def highly_variable_genes(self, adata, **kwargs):
                pass

            def regress_out(self, adata, **kwargs):
                pass

        backend = DummyBackend()
        assert backend.name == "dummy"
        set_backend(backend)
        assert get_backend() is backend

        # Reset to scanpy for other tests
        set_backend("scanpy")
