Installation Guide
==================

Requirements
------------

**Python Version**: 3.9 or higher

**Core Dependencies**:
- scanpy >= 1.10
- anndata >= 0.10
- pandas >= 1.5
- numpy >= 1.23
- scipy >= 1.10
- pydantic >= 2.0

Quick Install
-------------

Install from PyPI (recommended for most users)::

    pip install sclucid

Install with all optional dependencies::

    pip install sclucid[all]

Install with specific extras::

    # Analysis tools (scrublet, celltypist, gseapy, etc.)
    pip install sclucid[analysis]

    # External tools (scVelo, CellPhoneDB, etc.)
    pip install sclucid[tools]

    # Development version
    pip install sclucid[dev]

Install from Source
-------------------

For the latest development version::

    git clone https://github.com/yourusername/scLucid.git
    cd scLucid
    pip install -e .

Install with extras::

    pip install -e ".[all]"

Verify Installation
-------------------

After installation, verify that scLucid is working::

    python -c "import scLucid; print(scLucid.__version__)"

Run the test suite (optional)::

    pytest tests/

Optional Dependencies
---------------------

Analysis Extras
~~~~~~~~~~~~~~~~
- **scrublet**: Doublet detection
- **celltypist**: Automated cell type annotation
- **gseapy**: Gene set enrichment analysis
- **cosg**: Fast marker gene detection
- **hdbscan**: Density-based clustering

Tools Extras
~~~~~~~~~~~~
- **scvelo**: RNA velocity analysis
- **rpy2**: R bridge for CellPhoneDB, SCCODA
- **harmonypy**, **scanorama**, **scvi-tools**: Batch correction methods
- **triku**: Alternative HVG selection
- **squidpy**: Spatial analysis

Development Setup
------------------

For contributors::

    # Clone repository
    git clone https://github.com/yourusername/scLucid.git
    cd scLucid

    # Create virtual environment
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate

    # Install development dependencies
    pip install -e ".[dev]"

    # Install pre-commit hooks (optional)
    pre-commit install

Troubleshooting
---------------

Pydantic Import Error
~~~~~~~~~~~~~~~~~~~~~

If you see ``ModuleNotFoundError: No module named 'pydantic'``::

    pip install pydantic>=2.0

Scanpy Import Error
~~~~~~~~~~~~~~~~~~~

If scanpy fails to import with library loading errors (macOS)::

    # Reinstall scipy with conda
    conda install -c conda-forge scipy

    # Or use micromamba
    micromamba install -c conda-forge scipy

Memory Issues
~~~~~~~~~~~~~

For large datasets, enable low-memory mode::

    from scLucid.config import set_config
    set_config(low_memory_mode=True, chunk_size=1000)

Next Steps
----------

- :doc:`quickstart` - Get started with a basic workflow
- :doc:`../notebooks/01_quality_control` - Interactive QC tutorial
- :doc:`../notebooks/02_preprocessing` - Interactive preprocessing tutorial
- :doc:`../notebooks/03_clustering_annotation` - Clustering and annotation tutorial
- :doc:`best_practices` - Recommended practices for scLucid usage
