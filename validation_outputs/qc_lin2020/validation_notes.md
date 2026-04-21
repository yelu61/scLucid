# QC Validation: lin2020.pdac.h5ad

- **Dataset**: lin2020.pdac.h5ad
- **Cells before**: 9621
- **Cells after (pooled)**: 8613
- **Cells after (hierarchical)**: 9593

## Suspect sample comparison

### GSM4679533
```json
{
  "pooled": {},
  "hierarchical": {
    "n_genes_by_counts": {
      "lower": 0.0,
      "upper": 2269.2789653981736,
      "method": "hierarchical",
      "shrinkage_factor": 0.11286681715575622,
      "n_cells": 786
    },
    "total_counts": {
      "lower": 0.0,
      "upper": 13274.637840458874,
      "method": "hierarchical",
      "shrinkage_factor": 0.11286681715575622,
      "n_cells": 786
    },
    "pct_counts_mt": {
      "lower": 19.932206208667402,
      "upper": 100.0,
      "method": "hierarchical",
      "shrinkage_factor": 0.11286681715575622,
      "n_cells": 786
    }
  }
}
```

### GSM4679535
```json
{
  "pooled": {},
  "hierarchical": {
    "n_genes_by_counts": {
      "lower": 0.0,
      "upper": 1810.0912021760757,
      "method": "hierarchical",
      "shrinkage_factor": 0.08880994671403197,
      "n_cells": 1026
    },
    "total_counts": {
      "lower": 0.0,
      "upper": 7041.327268970044,
      "method": "hierarchical",
      "shrinkage_factor": 0.08880994671403197,
      "n_cells": 1026
    },
    "pct_counts_mt": {
      "lower": 0.0,
      "upper": 53.20158644628591,
      "method": "hierarchical",
      "shrinkage_factor": 0.08880994671403197,
      "n_cells": 1026
    }
  }
}
```

## Pooled warnings

- Tumor-aware QC: outlier_mt excluded from filtering criteria.

## Hierarchical warnings

- Tumor-aware QC: outlier_mt excluded from filtering criteria.
