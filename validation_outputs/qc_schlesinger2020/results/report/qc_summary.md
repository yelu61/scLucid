# QC Summary

- **Cells before**: 6499
- **Cells after**: 5889
- **Genes**: 19882
- **Threshold mode**: hierarchical
- **Strategy**: tumor_aware
- **Overall confidence**: 0.5
- **Tissue type**: tumor

## Filtering

- **Criteria used**: outlier_min_genes
- **Removed cells**: 610
- **Removed fraction**: 0.09386059393752885

## Concerns

- Missing metrics detected (n_genes, n_counts, pct_counts_mt); attempting automatic QC metric calculation

## Warnings

- Tumor-aware QC: outlier_mt excluded from filtering criteria.

## Tumor-aware Flags

```json
{
  "tissue_type": "tumor",
  "tumor_aware_enabled": true,
  "high_mt_population_flagged": false,
  "mean_pct_counts_mt": 5.065386918172671,
  "fraction_mt_above_10pct": 0.09924603785197723,
  "mean_pct_counts_ribo": 24.203867839266703,
  "note": "Tumor-aware QC active: elevated mitochondrial content is flagged rather than hard-filtered. Review thresholds manually."
}
```