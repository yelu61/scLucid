# QC Review Summary

## Recommendation Summary

- **Strategy**: tumor_aware
- **Confidence**: 0.5
- **Data Quality Score**: 100.0
- **Concerns**:
  - Missing metrics detected (n_genes, n_counts, pct_counts_mt); attempting automatic QC metric calculation

| Parameter | Recommended | User Provided |
|-----------|-------------|---------------|
| min_genes | 369 | 200 |
| pc_mt | None | 20.0 |
| min_counts | None | None |
| doublet_threshold | 0.5 | None |

## Decision Table

| Parameter | Applied | Source | Filter Enabled | Method | Confidence |
|-----------|---------|--------|----------------|--------|------------|
| min_genes | 369 | recommendation | True | GMM + Bootstrap | 0.5 |
| max_genes | None | disabled_or_not_available | False | None | None |
| n_counts | 885 | recommendation | False | log-normal distribution + bootstrap | 1.0 |
| max_counts | None | disabled_or_not_available | False | None | None |
| max_mt_percent | 10.8 | recommendation | False | distribution fitting (beta) | 0.5 |
| max_hb_percent | 20.0 | default_or_config | False | None | None |
| doublet_threshold | None | disabled_or_not_available | True | no_doublet_scores | 0.0 |
| nmads | 5.0 | default_or_config | False | None | None |

## Output Health

- **Status**: ok
- **Cells**: 1330
- **Genes**: 19882

## Benchmark Summary

- **Profile**: Tumor tissue (tumor)
- **Status**: pass
- **Retention rate**: 0.8866666666666667
- **Marker fidelity**: 0.9192560862697707

## Applied Thresholds

| Parameter | Value |
|-----------|-------|
| min_genes | 369 |
| max_genes | None |
| min_counts | 885 |
| max_counts | None |
| pc_mt | 10.8 |
| pc_hb | 20.0 |
| nmads | 5.0 |

## User Overrides

- **Overrides detected**: False

## Sample-Level Thresholds

- **Mode**: hierarchical
- **Samples with thresholds**: 0


## Tumor-Aware QC

- **Enabled**: True
- Tumor-aware QC is active: elevated mitochondrial content is flagged rather than hard-filtered.
- Mitochondrial outlier filtering was disabled for this tumor dataset.

## Filtering Results

- **Initial cells**: 1500
- **Final cells**: 1330
- **Removed**: 170 (11.3%)
- **Criteria used**: ['outlier_min_genes']
