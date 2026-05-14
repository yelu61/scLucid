# Publication-Quality Figure Examples

Each script in this directory generates **one** standalone publication-quality
figure from synthetic single-cell data. They are designed to be:

- **Self-contained** — no external datasets, no internet, no GPUs
- **Fast** — each script runs in under 30 seconds
- **Editable** — outputs are PDFs with embedded TrueType fonts
  (`pdf.fonttype=42`) so every text label can be edited in Illustrator
- **Themed** — uses scLucid's Nature theme by default; swap to
  `"science"` to retarget another journal

Use these scripts as starting points when adapting figures for your own
manuscript: copy the script, swap in your AnnData, tweak colors / titles.

## Scripts

| Script | Figure type | Demonstrates |
|---|---|---|
| `01_umap_annotation.py` | UMAP scatter colored by cell type | `scl.pl.plot_embedding` + Nature theme |
| `02_marker_heatmap.py` | Marker gene heatmap (cluster × gene) | `scl.pl.plot_marker_heatmap` + cluster ordering |
| `03_volcano_de.py` | Differential expression volcano plot | `scl.pl.plot_volcano` + significance highlighting |
| `04_cnv_heatmap.py` | CNV profile heatmap by cell group | Tumor module + chromosome-ordered heatmap |

## Running

From the repository root:

```bash
/path/to/python examples/04_publication_figures/01_umap_annotation.py
```

Each script writes its output PDF to `results/publication_figures/`.

## Customising For Your Journal

Each script ends with a single line that calls `apply_theme("nature")`. To
target another journal, change it to one of:

- `"nature"` — Arial sans-serif, ticks-out, no top/right spines
- `"science"` — Arial sans-serif with light grid, ticks-out

Custom themes can be created by following `src/scLucid/plotting/theme.py`.

## Notes

- All scripts use synthetic data so they can run anywhere and in CI. For your
  own analysis, replace the `_make_synthetic_data()` call at the top of each
  script with your real data loading step.
- Fonts: PDFs use TrueType (`pdf.fonttype=42`). To verify embedded fonts:
  `pdffonts your_figure.pdf` — every font row should show `TrueType`.
- Figure size defaults to `(6, 5)` inches — adjust the `figsize` kwarg to
  match your journal's column width specifications (e.g. Nature single
  column ≈ 89 mm = 3.5 in, double column ≈ 183 mm = 7.2 in).
