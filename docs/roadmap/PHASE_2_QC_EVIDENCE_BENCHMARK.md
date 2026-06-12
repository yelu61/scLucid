# Phase 2: QC Evidence Benchmark

**Recommended duration**: 1-2 个月  
**Primary output**: QC 证据链 benchmark 和 Figure 2  
**Main claim**: scLucid QC 比固定阈值流程更可审计，并在肿瘤场景中更少误伤关键生物信号。

## 目标

证明 scLucid 的 QC 不只是“自动过滤细胞”，而是能形成 reviewable decision workflow：

- 阈值不是黑箱。
- tumor-aware warning 能避免误删真实肿瘤状态。
- doublet evidence 能整合算法、谱系共表达和外部证据。
- QC 输出对后续 preprocess/analysis 有明确影响说明。

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

### Step 2. 建立 threshold decision table

每个 QC 结果必须输出：

- recommended threshold。
- applied threshold。
- source：user / default / recommendation。
- rationale。
- confidence。
- affected cells。
- risk note。

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

### Step 4. Tumor-aware QC validation

重点回答：

- 高 mitochondrial 是否可能代表肿瘤状态、缺氧、应激或技术坏细胞？
- 高 cell cycle 是否是 proliferating tumor，而不是过滤目标？
- epithelial/stromal/immune marker 是否在 QC 后被保留？
- malignant-like clusters 是否被过度删除？

输出：

- tumor cell preservation metric。
- marker fidelity。
- retained tumor program score。
- warning examples。

### Step 5. 对照固定阈值流程

每个数据集比较：

- Scanpy/Seurat fixed threshold。
- scLucid recommendation。
- scLucid user override。

核心判断不是“保留越多越好”，而是：

- 是否减少样本偏置。
- 是否保护关键 marker/program。
- 是否识别低质量风险。
- 是否生成可审计理由。

## 完成标准

Phase 2 通过条件：

- 至少 5 个数据集完成 QC benchmark。
- 每个数据集都有 machine-readable QC review summary。
- 至少 2 个 tumor case 显示 scLucid QC 相比固定阈值更少误伤关键生物信号。
- doublet evidence 有算法对照、fallback 和解释表。
- Figure 2 草图完成。

## 交付物

- `validation/qc/` benchmark scripts。
- `validation_outputs/qc_*` 汇总表。
- `docs/VALIDATION_QC_EVIDENCE.md` 或等价报告。
- Figure 2 数据表和绘图脚本。
- QC audit report 示例。

## 推荐主图

Figure 2A: QC workflow and evidence table。  
Figure 2B: retention and marker fidelity across datasets。  
Figure 2C: tumor-aware QC preserves malignant/TME programs。  
Figure 2D: doublet evidence overlap and risk decomposition。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 无 ground truth | 很难说更准 | 使用 marker fidelity、retention bias、synthetic doublets |
| scDblFinder Python port 与 R 不一致 | doublet calls 差异大 | Phase 5 做 parity benchmark |
| tumor-aware 被认为主观 | warning 太像建议 | 所有 warning 配 evidence key 和 review status |
| 只赢固定阈值太弱 | 审稿人说基线简单 | 加入 Scanpy recommended practices 和 Seurat baseline |

