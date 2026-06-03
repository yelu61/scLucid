"""Tests for enrichment analysis helpers."""

import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.analysis.config import EnrichmentConfig
from scLucid.analysis.differential_expression import enrichment
from scLucid.analysis.differential_expression.enrichment import run_enrichment


class _FakeEnrichResult:
    def __init__(self):
        self.results = pd.DataFrame(
            {
                "Term": ["Pathway"],
                "P-value": [0.01],
                "Adjusted P-value": [0.02],
                "Genes": ["GeneA;GeneB"],
            }
        )


def test_run_enrichment_partitions_pseudobulk_group_and_contrast(monkeypatch, tmp_path):
    adata = AnnData(X=np.ones((4, 3)))
    adata.var_names = ["GeneA", "GeneB", "GeneC"]
    marker_df = pd.DataFrame(
        {
            "names": ["GeneA", "GeneB", "GeneA", "GeneB", "GeneA", "GeneB", "GeneA", "GeneB"],
            "logfoldchanges": [2, 1, -1, -2, 1.5, 0.5, -0.5, -1.5],
            "scores": [4, 3, -3, -4, 2, 1, -1, -2],
            "group": ["T", "T", "T", "T", "B", "B", "B", "B"],
            "contrast": ["B_vs_A", "B_vs_A", "A_vs_B", "A_vs_B"] * 2,
        }
    )
    adata.uns["sclucid"] = {"analysis": {"de": {"pb_de": marker_df}}}
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("Pathway\tna\tGeneA\tGeneB\n")

    calls = []

    def fake_enrich(gene_list, gene_sets, background, outdir, cutoff):
        calls.append(tuple(gene_list))
        return _FakeEnrichResult()

    monkeypatch.setattr(enrichment.gp, "enrich", fake_enrich)
    config = EnrichmentConfig(
        de_key="pb_de",
        mode="offline",
        method="ora",
        custom_gene_sets=str(gmt),
        min_genes_for_ora=1,
        n_top_genes_ora=2,
        cutoff_pval=1.0,
    )

    result = run_enrichment(adata, groupby="cell_type", config=config)

    assert set(result) == {"T|B_vs_A", "T|A_vs_B", "B|B_vs_A", "B|A_vs_B"}
    assert len(calls) == 4
