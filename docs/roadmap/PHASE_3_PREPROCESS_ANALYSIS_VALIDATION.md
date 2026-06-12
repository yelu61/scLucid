# Phase 3: Preprocess + Analysis 稳定性验证

**Recommended duration**: 1-2 个月  
**Primary output**: Figure 3 and workflow stability benchmark  
**Main claim**: scLucid 的预处理和分析不是只生成漂亮 UMAP，而是保留生物信号、记录推断边界并提供可复查的 annotation evidence。

## 目标

建立 preprocess 和 analysis 的分层验证，证明 scLucid 在真实项目中能提供：

- 稳定的数据层语义。
- 保守的 batch correction 建议。
- 可解释的 clustering resolution。
- marker/reference/negative-marker 组成的 annotation evidence。
- 明确的 DE/proportion inference semantics。

## 准备

### 数据集

使用 Phase 2 中完成 QC 的数据，至少包括：

- PBMC baseline。
- PDAC tumor。
- 第二癌种 tumor。
- 一个 batch-heavy 或 multi-patient 数据集。

### 对照流程

- Scanpy standard workflow。
- Seurat standard workflow。
- scVI/scANVI, 如果数据规模和依赖允许。
- CellTypist/reference annotation, 如果适合。

### 基准资源

- marker panels。
- known cell type labels。
- expected compartments：
  - epithelial/malignant。
  - T/NK。
  - B/plasma。
  - myeloid。
  - CAF/stromal。
  - endothelial。

## 具体步骤

### Step 1. Preprocess layer contract validation

每个 workflow 后检查：

- `layers["counts"]` 是否保留 raw counts。
- `adata.X` 当前语义是否记录。
- `adata.raw` 是否设置且语义清楚。
- normalized/log data 保存位置。
- PCA/neighbors/UMAP key 是否符合命名。
- highly variable genes 是否记录参数和版本。

输出：

- data contract table。
- h5ad round-trip test。
- invalid input failure examples。

### Step 2. HVG 与 marker preservation benchmark

比较不同流程对 marker/program genes 的保留：

- standard HVG。
- tumor-aware retained genes。
- marker/resource genes。
- pathway/module genes。

指标：

- HVG overlap。
- marker inclusion rate。
- tumor program gene retention。
- downstream cluster marker recovery。

### Step 3. Batch correction diagnostic benchmark

scLucid 不应默认“强行校正”，而应先诊断：

- batch mixing。
- biological separation。
- patient/tumor structure。
- overcorrection risk。

对照：

- no correction。
- Harmony/BBKNN/scVI, 视项目支持。
- Seurat integration。

输出：

- batch mixing score。
- biological conservation score。
- marker fidelity after correction。
- warning and recommendation summary。

### Step 4. Clustering stability

对 resolution、neighbors、random seed 做稳定性评估：

- ARI/NMI between resolutions。
- marker specificity per cluster。
- rare cell population preservation。
- cluster-level QC flags。
- post-hoc doublet/low-quality cluster review。

完成标准：

- scLucid 能解释为什么推荐某个 resolution。
- 不把 artifact cluster 直接当作真实 cell type。
- cluster review 有 action items。

### Step 5. Annotation evidence benchmark

每个 annotation 需要证据表：

- positive markers。
- negative markers。
- reference/celltypist evidence。
- artifact/stress evidence。
- confidence。
- review status。
- unresolved/suspect label。

对照：

- manual marker annotation。
- CellTypist/reference。
- Seurat label transfer, 如果可行。

指标：

- major lineage agreement。
- marker consistency。
- ambiguous cluster rate。
- annotation confidence calibration。

### Step 6. DE/proportion inference safety

输出必须区分：

- exploratory cell-level。
- descriptive sample-level。
- valid sample-level inference。

检查：

- cell-level DE 是否有 warning。
- pseudobulk 是否使用 sample/patient replicate。
- proportion analysis 是否区分 descriptive 和 inferential。
- plots 是否避免过度宣称。

## 完成标准

Phase 3 通过条件：

- 至少 3 个数据集完成 preprocess + analysis benchmark。
- 每个数据集有 data contract report。
- annotation evidence table 可以人工复查。
- batch correction 建议有 evidence，而不是默认开关。
- DE/proportion 输出含 inference level。
- Figure 3 草图完成。

## 交付物

- preprocess contract tests。
- analysis contract tests。
- `validation/preprocess_analysis/` scripts。
- marker preservation summary。
- clustering stability report。
- annotation evidence report。
- Figure 3 绘图脚本。

## 推荐主图

Figure 3A: preprocessing contract and evidence flow。  
Figure 3B: marker/program preservation across workflows。  
Figure 3C: clustering stability and annotation confidence。  
Figure 3D: inference semantics guardrails for DE/proportion。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| UMAP 视觉差异难量化 | 结果像主观比较 | 用 marker fidelity、cluster stability、batch metrics |
| annotation ground truth 不可靠 | label agreement 误导 | 分 major lineage 和 subtype/state 两层评估 |
| batch correction 争议大 | 不同工具目标不同 | 强调 diagnostic recommendation 而非单一最优 |
| DE 结果被过度解释 | 审稿人质疑统计 | inference semantics 作为核心防线 |

