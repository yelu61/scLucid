# QC Review Summary

## Recommendation Summary

- **Strategy**: aggressive
- **Confidence**: 0.5
- **Data Quality Score**: 100.0
- **Concerns**:
  - Missing metrics detected (pct_counts_mt); attempting automatic QC metric calculation

| Parameter | Recommended | User Provided |
|-----------|-------------|---------------|
| min_genes | 643 | 300 |
| pc_mt | None | 20.0 |
| min_counts | None | None |
| doublet_threshold | 0.5 | None |

## Decision Table

| Parameter | Applied | Source | Filter Enabled | Method | Confidence |
|-----------|---------|--------|----------------|--------|------------|
| min_genes | 300 | user_override | True | GMM + Bootstrap | 0.5 |
| max_genes | None | disabled_or_not_available | False | None | None |
| n_counts | 1503 | recommendation | False | log-normal distribution + bootstrap | 1.0 |
| max_counts | None | disabled_or_not_available | False | None | None |
| max_mt_percent | 20.0 | user_override | True | distribution fitting (beta) | 0.5 |
| max_hb_percent | 20.0 | default_or_config | False | None | None |
| doublet_threshold | None | disabled_or_not_available | True | no_doublet_scores | 0.0 |
| nmads | 5.0 | default_or_config | False | None | None |

## Output Health

- **Status**: ok
- **Cells**: 2523
- **Genes**: 32738

## Benchmark Summary

- **Profile**: PBMC / immune suspension (pbmc)
- **Status**: pass
- **Retention rate**: 0.9344444444444444
- **Marker fidelity**: 0.9786796153336377

## Applied Thresholds

| Parameter | Value |
|-----------|-------|
| min_genes | 300 |
| max_genes | None |
| min_counts | 1503 |
| max_counts | None |
| pc_mt | 20.0 |
| pc_hb | 20.0 |
| nmads | 5.0 |

## User Overrides

- **Overrides detected**: True
- **Details**:
  - min_genes: recommended=643, user=300
  - max_mt_percent: recommended=3.1, user=20.0

## Sample-Level Thresholds

- **Mode**: hierarchical
- **Samples with thresholds**: 4

```json
{
  "sample3": {
    "n_genes_by_counts": {
      "lower": 414.8,
      "upper": 1299.49
    },
    "total_counts": {
      "lower": 735.76,
      "upper": 4067.66
    },
    "pct_counts_mt": {
      "lower": 0.11,
      "upper": 4.3
    }
  },
  "sample4": {
    "n_genes_by_counts": {
      "lower": 365.74,
      "upper": 1336.52
    },
    "total_counts": {
      "lower": 430.9,
      "upper": 4335.37
    },
    "pct_counts_mt": {
      "lower": 0.46,
      "upper": 3.96
    }
  },
  "sample1": {
    "n_genes_by_counts": {
      "lower": 392.29,
      "upper": 1288.72
    },
    "total_counts": {
      "lower": 644.72,
      "upper": 4018.76
    },
    "pct_counts_mt": {
      "lower": 0.26,
      "upper": 4.15
    }
  },
  "sample2": {
    "n_genes_by_counts": {
      "lower": 360.39,
      "upper": 1318.73
    },
    "total_counts": {
      "lower": 472.28,
      "upper": 4233.33
    },
    "pct_counts_mt": {
      "lower": 0.38,
      "upper": 4.1
    }
  }
}
```

## Tumor-Aware QC

- **Enabled**: False

## Filtering Results

- **Initial cells**: 2700
- **Final cells**: 2523
- **Removed**: 177 (6.6%)
- **Criteria used**: ['outlier_min_genes', 'outlier_mt']
