# QC Validation: schlesinger2020.pdac.h5ad

- **Dataset**: schlesinger2020.pdac.h5ad
- **Cells before/after**: 6499 / 5892
- **Threshold mode**: hierarchical
- **Strategy**: tumor_aware
- **Tumor-aware enabled**: True

## Tumor-aware flags

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

## Filtering summary

```json
{
  "initial_cells": 6499,
  "final_cells": "5892",
  "removed_cells": "607",
  "removed_fraction": 0.09339898445914756,
  "criteria_used": [
    "outlier_min_genes"
  ],
  "combination_logic": "any",
  "criteria_counts": {
    "outlier_min_genes": "607"
  },
  "config": {
    "save_dir": null,
    "verbose": true,
    "plot": true,
    "report": true,
    "criteria_to_filter": [
      "outlier_min_genes"
    ],
    "combination_logic": "any",
    "custom_logic_expr": null,
    "min_criteria_for_removal": 2,
    "metadata_filters": null
  }
}
```

## Warnings

- Tumor-aware QC: outlier_mt excluded from filtering criteria.
