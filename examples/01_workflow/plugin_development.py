"""Plugin extension example for scLucid.

This example shows how to define, register, and directly execute custom
``AnalysisStep`` implementations. The registry is an extension mechanism for
creating step instances; it is not automatically wired into
``scLucid.analysis.run_custom_analysis``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from pydantic import Field

from scLucid.base_config import SclucidBaseConfig
from scLucid.base_interfaces import AnalysisStep, AnalysisStepFactory


class CustomQCConfig(SclucidBaseConfig):
    """Configuration for a simple custom QC step."""

    min_genes: int = Field(default=200, ge=0)
    max_mt_percent: float = Field(default=20.0, ge=0, le=100)


class HighStringencyQC(AnalysisStep):
    """Filter cells using stricter gene-count and mitochondrial thresholds."""

    def __init__(self, config: Optional[CustomQCConfig] = None):
        super().__init__(config or CustomQCConfig())

    def validate_input(self, adata: AnnData) -> bool:
        if adata.n_obs == 0:
            raise ValueError("Input AnnData has no cells.")
        if "pct_counts_mt" not in adata.obs:
            raise ValueError("Run QC metrics first so adata.obs['pct_counts_mt'] exists.")
        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        self.validate_input(adata)
        before = int(adata.n_obs)
        filtered = adata[
            (adata.obs["n_genes_by_counts"] >= self.config.min_genes)
            & (adata.obs["pct_counts_mt"].astype(float) < self.config.max_mt_percent)
        ].copy()
        self._results = {
            "n_cells_before": before,
            "n_cells_after": int(filtered.n_obs),
            "min_genes": self.config.min_genes,
            "max_mt_percent": self.config.max_mt_percent,
        }
        return filtered


class CustomAnnotatorConfig(SclucidBaseConfig):
    """Configuration for a toy custom annotation step."""

    reference_path: Optional[str] = Field(default=None)
    similarity_threshold: float = Field(default=0.8, ge=0, le=1)
    key_added: str = Field(default="custom_cell_type")


class MyCustomAnnotator(AnalysisStep):
    """Toy annotator that demonstrates where a real model would write labels."""

    def __init__(self, config: Optional[CustomAnnotatorConfig] = None):
        super().__init__(config or CustomAnnotatorConfig())

    def validate_input(self, adata: AnnData) -> bool:
        if "X_pca" not in adata.obsm:
            raise ValueError("Run PCA first so adata.obsm['X_pca'] exists.")
        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        self.validate_input(adata)
        cell_types = np.array(["T cells", "B cells", "NK cells", "Monocytes"])
        rng = np.random.default_rng(kwargs.get("random_state", 42))
        adata.obs[self.config.key_added] = rng.choice(cell_types, size=adata.n_obs)
        adata.obs[f"{self.config.key_added}_confidence"] = rng.random(adata.n_obs)
        self._results = {
            "method": "toy_custom_annotation",
            "key_added": self.config.key_added,
            "n_cell_types": int(len(cell_types)),
        }
        return adata


def register_custom_plugins() -> None:
    """Register custom steps with the factory."""

    AnalysisStepFactory.register("high_stringency_qc", HighStringencyQC)
    AnalysisStepFactory.register("my_annotator", MyCustomAnnotator)


def prepare_demo_data() -> AnnData:
    """Create a small demo AnnData with QC metrics and PCA."""

    adata = sc.datasets.pbmc3k()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    return adata


def example_usage() -> AnnData:
    """Run registered plugin steps directly."""

    register_custom_plugins()
    adata = prepare_demo_data()

    qc_step = AnalysisStepFactory.create(
        "high_stringency_qc",
        config=CustomQCConfig(min_genes=300, max_mt_percent=15.0),
    )
    annotator = AnalysisStepFactory.create(
        "my_annotator",
        config=CustomAnnotatorConfig(key_added="custom_cell_type"),
    )

    adata = qc_step.run(adata)
    adata = annotator.run(adata, random_state=42)

    print("Registered steps:", AnalysisStepFactory.list_steps())
    print("QC summary:", qc_step.get_summary())
    print("Annotation summary:", annotator.get_summary())
    return adata


if __name__ == "__main__":
    print("Plugin extension example")
    print("This example registers custom steps and executes them directly.")
    print("For package workflows, use supported workflow entrypoints unless a")
    print("registry bridge is explicitly implemented.")
    example_usage()
