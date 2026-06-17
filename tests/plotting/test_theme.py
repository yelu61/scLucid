"""Tests for plotting theme helpers."""

import pandas as pd
from anndata import AnnData

from scLucid.plotting.theme import build_obs_palette


def test_build_obs_palette_uses_obs_values_not_column_names():
    adata = AnnData(
        obs=pd.DataFrame(
            {
                "sampleID": ["S1", "S2", "S1"],
                "group": ["ctrl", "treat", "ctrl"],
            },
            index=["c1", "c2", "c3"],
        )
    )

    palette = build_obs_palette(
        adata,
        ["sampleID", "group"],
        color_maps={
            "samples": {"S1": "#111111", "S2": "#222222"},
            "groups": {"ctrl": "#333333", "treat": "#444444"},
        },
    )

    assert palette == {
        "S1": "#111111",
        "S2": "#222222",
        "ctrl": "#333333",
        "treat": "#444444",
    }
    assert "sampleID" not in palette
    assert "group" not in palette
