# QC Summary

- **Cells before**: 9621
- **Cells after**: 9593
- **Genes**: 49056
- **Threshold mode**: hierarchical
- **Strategy**: tumor_aware
- **Overall confidence**: 0.5
- **Tissue type**: tumor

## Filtering

- **Criteria used**: outlier_min_genes
- **Removed cells**: 28
- **Removed fraction**: 0.0029103003845754078

## Concerns

- Missing metrics detected (n_genes, n_counts, pct_counts_mt); attempting automatic QC metric calculation

## Warnings

- Tumor-aware QC: outlier_mt excluded from filtering criteria.

## Tumor-aware Flags

```json
{
  "tissue_type": "tumor",
  "tumor_aware_enabled": true,
  "high_mt_population_flagged": true,
  "mean_pct_counts_mt": 15.291096996204471,
  "fraction_mt_above_10pct": 0.3027751792952916,
  "mean_pct_counts_ribo": 18.427964376590342,
  "note": "Tumor-aware QC active: elevated mitochondrial content is flagged rather than hard-filtered. Review thresholds manually."
}
```