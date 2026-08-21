# Controlled Preprocess validation

## Question and estimand

The controlled mixology benchmark asks whether the preprocessing policy chosen
by scLucid remains close to the best registered candidate when the measurement
protocol is unseen. The primary estimand is mean balanced cell-line identity
accuracy across three leave-one-protocol-out folds. The experimental unit is
protocol, not cell.

The registered candidates are:

- standard log normalization, batch-aware HVG, and unintegrated PCA;
- analytic Pearson residuals with batch-aware HVG;
- standard log normalization with multinomial deviance features;
- Pearson residuals with multinomial deviance features.

Feature selection and PCA fitting use only the two training protocols in each
fold. The held-out protocol is used only for evaluation. The primary task
utility is known cell-line identity accuracy; batch mixing, cross-protocol
identity purity, identity silhouette, clustering agreement, and graph seed
stability remain separate diagnostics.

## Run the controlled benchmark

```bash
python validation/preprocess/run_mixology_preprocess_benchmark.py \
  --input data/public_mixology.h5ad \
  --output-dir validation_outputs/current/preprocess_mixology \
  --n-top-genes 2000 \
  --n-pcs 30 \
  --contract validation/qc_preprocess/acceptance_contract.json
```

The runner calls the maintained
`recommend_preprocess_policy -> apply_preprocess_policy` path and verifies the
counts, full-gene normalized, discovery-feature, and PCA representation
contracts. It then evaluates the actual selected policy rather than a
hard-coded surrogate.

The locked configuration currently gives:

- selected standard baseline held-out accuracy: 0.9835;
- best eligible candidate accuracy: 0.9883;
- selected-policy regret: 0.0049, below the 0.05 limit;
- selected-policy biology loss: 0;
- sensitivity regret: 0.0016 with 1,000 genes/20 PCs and 0.0026 with 3,000
  genes/50 PCs.

Harmony improves batch-mixing similarity by 0.1803 and identity silhouette by
0.0806. Its graph seed stability decreases by 0.0215, slightly beyond the
locked 0.02 loss tolerance. It therefore does not Pareto-dominate the
unintegrated baseline in this benchmark.

These observations support retaining the simple baseline provisionally for
this controlled mixture. They do not establish universal superiority or tumor
project benefit.

## Real-project usability gate

Initialize the three-project evidence record with:

```bash
python validation/qc_preprocess/run_real_project_ux_acceptance.py \
  --output-dir validation_outputs/current/real_project_ux \
  --contract validation/qc_preprocess/acceptance_contract.json
```

The generated TSV records legacy and current user-edited configuration fields,
manual doublet deletion, manual review-summary edits, schema bypasses,
project-specific patches, and stable `RunEvidence`. Blank or incomplete rows
remain `BLOCKED`; the runner never imputes a usability success.

After the maintained four-action path has been executed in 202604JJH,
202507LPJ, and 202603AK112, rerun with `--input` pointing to the completed TSV.
Passing requires at least 70% fewer edited fields, no listed manual workaround,
no project-specific patch, and a stable evidence record for every project.

## Claim boundary

- Supported: controlled cross-protocol candidate regret and representation
  contract results for public mixology.
- Exploratory: whether the same method ranking and integration trade-off holds
  in heterogeneous tumor samples.
- Unsupported: universal preprocessing superiority, automatic integration, or
  release of QC/Preprocess as scientific `CORE`.
