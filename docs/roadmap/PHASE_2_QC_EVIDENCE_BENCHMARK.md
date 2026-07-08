# Phase 2: QC Evidence Benchmark

**Recommended duration**: 1-2 个月  
**Primary output**: QC 证据链 benchmark、Figure 2 source data、claim scorecard
**Main claim**: scLucid QC 比固定阈值流程更可审计，并在肿瘤场景中更少误伤关键生物信号。当前证据支持系统性证据框架和若干明确优势方向，但不应表述为所有场景下的全面准确性优越。

## 目标

证明 scLucid 的 QC 不只是“自动过滤细胞”，而是能形成 reviewable decision workflow：

- 阈值不是黑箱。
- tumor-aware warning 能避免误删真实肿瘤状态。
- doublet evidence 能整合算法、谱系共表达和外部证据。
- QC 输出对后续 preprocess/analysis 有明确影响说明。
- QC-to-preprocess handoff 能明确推荐 counts layer、保留
  `qc_decision` / `qc_remove` / `qc_review_required`，并说明哪些细胞应进入
  review 或 sensitivity-only 记录。

本 phase 的 vision task 是把 QC 从“过滤坏细胞”推进为“生物学风险归因”：
高 MT%、低检测基因数、ambient、doublet-like、stress-high 或 hypoxia-like
细胞不应只被二值删除，而应进入可复查的风险、状态和敏感性记录。

## 准备

### 数据集

至少准备：

- PBMC normal baseline：1-2 个。
- PDAC tumor：1 个。
- 第二癌种 tumor：1-2 个。
- 低质量或复杂样本：高 mitochondrial、高 doublet、强 batch、样本量不均衡。

每个数据集需要 metadata：

- sample key。
- known cell type, 如果有。
- tumor/normal label, 如果有。
- doublet ground truth 或 proxy, 如果有 hashing/genotype/synthetic doublets。

### 对照流程

必须建立：

- Scanpy fixed-threshold QC baseline。
- Seurat fixed-threshold QC baseline。
- scLucid adaptive QC。
- scLucid tumor-aware QC。
- doublet baseline：Scrublet、scDblFinder/pyscdblfinder、scLucid heuristic。

## 具体步骤

### Step 1. 固定 QC 指标面板

每个数据集统一输出：

- cells before/after QC。
- retention rate。
- sample-level retention bias。
- cell-type retention bias。
- `n_genes_by_counts` distribution。
- `total_counts` distribution。
- `pct_counts_mt` distribution。
- top gene fraction。
- doublet score and predicted doublet rate。
- marker fidelity before/after QC。
- **count mixture model for `n_genes_by_counts`**：fitted family (NB/ZINB/Poisson), AIC, fallback status。
- **ambient RNA correction residual score**：linear and/or CellBender backend。
- **MT% review band**：hard threshold and review-band lower bound for tumor datasets。

### Step 2. 建立 threshold decision table

每个 QC 结果必须输出：

- recommended threshold。
- applied threshold。
- source：user / default / recommendation / count_mixture / bimodal_gmm / sample_aware。
- rationale。
- confidence。
- affected cells。
- risk note。
- **model evidence**：count-model AIC, MT component separation, stratum baselines。
- **handoff readiness**：recommended preprocess counts layer,
  review/sensitivity cell fractions, safe-to-continue flags, blockers, and
  downstream handling requirements.

标准字段：

| Field | Meaning |
|-------|---------|
| `parameter` | QC 参数 |
| `recommended` | 推荐值 |
| `applied` | 实际值 |
| `source` | 推荐来源 |
| `confidence` | 可信度 |
| `evidence` | 证据摘要 |
| `review_required` | 是否需要人工复查 |

### Step 3. Doublet evidence benchmark

重点比较：

- Scrublet。
- scDblFinder / pyscdblfinder。
- scLucid lineage co-expression heuristic。
- algorithm-only vs algorithm-plus-heuristic fusion, including candidate
  `algorithm_weight` values.
- external evidence, 如果有。

输出：

- algorithm score distribution。
- algorithm overlap。
- doublet calls by sample。
- doublet calls by cluster/cell type。
- heterotypic vs homotypic risk decomposition。
- top doublet evidence profiles。

完成标准：

- 缺少 optional dependency 时 workflow 有清晰 fallback。
- scDblFinder 专属 `dbr` 和通用 expected rate 行为明确。
- raw-count guard 能阻止 normalized data 误用。
- doublet evidence report 不只给二值标签。
- Python/R scDblFinder parity 和 disagreement group 进入报告。
- Kang demuxlet 作为 external genotype-based reference，而不是全局金标准。

### Step 4. Tumor-aware QC validation

重点回答：

- 高 mitochondrial 是否可能代表肿瘤状态、缺氧、应激或技术坏细胞？
- 高 cell cycle 是否是 proliferating tumor，而不是过滤目标？
- epithelial/stromal/immune marker 是否在 QC 后被保留？
- malignant-like clusters 是否被过度删除？
- **代谢重编程状态（OXPHOS / glycolysis / mt biogenesis）是否在高-MT 细胞中保留？**

输出：

- tumor cell preservation metric。
- marker fidelity。
- retained tumor program score。
- **MT% threshold scorecard**：fixed 20%、fixed 5%、bimodal GMM、sample-aware、multicomponent 的 program retention 对比。
- QC effect on tumor purity / immune composition / stromal composition。
- "what if retained" sensitivity summary for high-risk cell bands, when data support it。
- warning examples。

### Step 5. 对照固定阈值流程

每个数据集比较：

- Scanpy/Seurat fixed threshold。
- scLucid recommendation。
- scLucid user override。
- **scLucid count-adaptive threshold for `n_genes_by_counts`。**
- **scLucid sample-aware / multicomponent MT% threshold for tumor datasets。**

核心判断不是“保留越多越好”，而是：

- 是否减少样本偏置。
- 是否保护关键 marker/program。
- 是否识别低质量风险。
- 是否生成可审计理由。
- **统计模型是否更适合数据类型（count mixture vs GMM/percentile）。**

## 完成标准

Phase 2 通过条件：

- 至少 5 个数据集完成 QC benchmark。当前本地 inventory 已有 8 个 h5ad
  validation datasets，其中 7 个进入 threshold/tumor/doublet/ambient QC
  evidence 路径。
- 每个 benchmark 输出都有 machine-readable reviewer/source table。
- 标准 QC review summary 包含 `qc_handoff_readiness`，可直接被 preprocess
  和 Figure evidence package 消费。
- 至少 2 个 tumor case 显示 scLucid QC 相比固定阈值更少误伤关键生物信号。当前 PDAC/NSCLC/CRC 通过 marker/program retention proxy 支撑该方向。
- doublet evidence 有算法对照、fallback、Python/R parity、threshold
  calibration、algorithm_weight recommendation 和解释表。
- Figure 2 source-data package 完成，并区分 supported / partial /
  contract-only claims。

## 交付物

- `validation/qc/` benchmark scripts。
- `validation_outputs/qc_*` 汇总表。
- `validation_outputs/qc_figure2_package/figure2_qc_source_data.tsv`。
- `validation_outputs/qc_figure2_package/qc_claim_scorecard.tsv`。
- `validation_outputs/qc_figure2_package/qc_evidence_report.md`。
- Figure 2 数据表和绘图脚本。
- QC audit report 示例。
- QC-to-preprocess handoff source rows for the combined QC/preprocess Figure X
  package and real-project acceptance review.

## 推荐主图

Figure 2A: QC workflow and evidence table。  
Figure 2B: retention and marker fidelity across datasets。  
Figure 2C: tumor-aware QC preserves malignant/TME programs。  
Figure 2D: doublet evidence overlap, calibration, parity, and risk decomposition。
Figure 2E: reviewer-facing tumor QC decision narrative。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 无 ground truth | 很难说更准 | 使用 marker fidelity、retention bias、synthetic doublets |
| scDblFinder Python port 与 R 不一致 | doublet calls 差异大 | 已有 Python/R parity table；继续报告 disagreement donor/sample/cell type |
| tumor-aware 被认为主观 | warning 太像建议 | 所有 warning 配 evidence key 和 review status |
| 只赢固定阈值太弱 | 审稿人说基线简单 | 加入 Scanpy recommended practices 和 Seurat baseline |
| Ambient claim 过强 | CellBender tiny 只能证明接口 | 标记为 contract-only，补 full raw 10x + SoupX/CellBender reference |
| count mixture 稳定性 | NB/ZINB 退化或 fallback 到 GMM | 增加 equidispersion 预筛、log-space 优化、Poisson 平局 |
| ambient backend 可选依赖 | CellBender 未安装 | 自动 fallback 到 linear correction，并在 summary 标注 backend |
| MT% 代谢异质性 | 单阈值误删高 OXPHOS / 增殖肿瘤细胞 | 提供 sample-aware baseline、review band、multicomponent GMM |
