# QC Summary

- **Cells before**: 1500
- **Cells after**: 1330
- **Genes**: 19882
- **Threshold mode**: hierarchical
- **Strategy**: tumor_aware
- **Overall confidence**: 0.5
- **Tissue type**: pancreatic_tumor

## Filtering

- **Criteria used**: outlier_min_genes
- **Removed cells**: 170
- **Removed fraction**: 0.11333333333333333

## Concerns

- Missing metrics detected (n_genes, n_counts, pct_counts_mt); attempting automatic QC metric calculation

## Warnings

- None

## Tumor-aware Flags

```json
{
  "tissue_type": "pancreatic_tumor",
  "tumor_aware_enabled": true,
  "high_mt_population_flagged": false,
  "mean_pct_counts_mt": 4.994709185954394,
  "fraction_mt_above_10pct": 0.104,
  "mean_pct_counts_ribo": 24.275450177066222,
  "note": "Tumor-aware QC active: elevated mitochondrial content is flagged rather than hard-filtered. Review thresholds manually."
}
```