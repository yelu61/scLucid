# Phase 4: Tumor Interpretation 与 Case Studies

**Recommended duration**: 2-3 个月  
**Primary output**: Figure 4 and tumor case-study evidence  
**Main claim**: scLucid 的核心特色是肿瘤单细胞解释系统，而不是通用单细胞工具集合。

## 目标

把 scLucid 的肿瘤特色变成可投稿的生物学叙事：

- malignant cell evidence。
- CNV/malignancy signature。
- tumor programs。
- TME states。
- CAF/myeloid/immune exhaustion。
- spatial/bulk concordance。
- ecotype/cellular ecosystem summaries。

Phase 4 的重点不是发明所有肿瘤算法，而是把高水平肿瘤单细胞论文中常见分析内容整理为可复用、可审计、可验证的 interpretation workflow。

## 准备

### 癌种选择

建议至少 3 个 case study：

1. PDAC：CAF/myeloid/TME 复杂，适合做首个 tumor benchmark。
2. NSCLC 或 CRC：免疫/TME 较丰富，公开数据多。
3. BRCA/HCC/melanoma/GBM 任选一个：展示泛化能力。

每个 case study 准备：

- raw or processed h5ad。
- sample/patient metadata。
- published cell type annotations, 如果有。
- paper-reported marker/programs。
- clinical/response metadata, 如果有。
- bulk/spatial paired data, 如果有。

### marker/program 资源

至少准备模块：

- epithelial/malignant markers。
- proliferation。
- EMT。
- hypoxia。
- interferon/inflammation。
- cytotoxicity。
- exhaustion。
- suppressive myeloid。
- inflammatory myeloid。
- CAF subtypes。
- endothelial states。

## 具体步骤

### Step 1. Malignancy evidence workflow

对 epithelial/tumor-like cells 输出：

- marker-based malignancy evidence。
- CNV evidence, 如果有。
- malignancy signature score。
- tumor/normal epithelial distinction。
- confidence。
- unresolved/suspect calls。

标准输出：

| Field | Meaning |
|-------|---------|
| `malignancy_call` | malignant / non_malignant / suspect / unresolved |
| `confidence` | 证据强度 |
| `evidence_sources` | marker, CNV, signature, manual |
| `review_required` | 是否需要人工确认 |
| `limitations` | 缺失 CNV 或 reference 等限制 |

### Step 2. TME state interpretation

对主要 compartment 输出 state/program：

- T/NK：cytotoxicity、exhaustion、IFN response。
- Myeloid：inflammatory、suppressive、TAM-like、DC-like。
- CAF：myCAF、iCAF、apCAF 或文献定义变体。
- Endothelial：angiogenic、lymphatic、vascular。
- Plasma/B：plasma、germinal-center-like、memory-like, 视数据支持。

每个 state 必须有：

- gene set 来源。
- score distribution。
- sample-level abundance。
- marker heatmap。
- review note。

### Step 3. Ecotype / cellular ecosystem summary

目标不是复制某个 ecotype 工具，而是引入可解释生态型概念：

- sample-level TME composition。
- malignant program + stromal/immune state co-occurrence。
- cellular community / ecosystem score。
- patient-level clustering。
- ecotype heatmap。

输出：

- sample x program matrix。
- sample x cell-state abundance matrix。
- ecosystem summary table。
- ecotype label, 如果证据足够。
- exploratory status 标注。

### Step 4. Bulk/spatial evidence linkage

使用 `tools.bulk` 和 `tools.spatial` 作为 evidence modules：

- bulk expression 与 single-cell program concordance。
- pseudobulk DE 与 bulk DE concordance。
- deconvolution 与 single-cell composition concordance。
- spatial SVG / niche 与 single-cell cell-state program concordance。

所有 bulk/spatial 输出需要标注：

- exploratory。
- descriptive。
- sample-level。
- spatial evidence。

### Step 5. Case-study narrative

每个癌种输出一页 case-study summary：

- 数据集介绍。
- QC 与保留情况。
- major compartments。
- malignancy evidence。
- TME state map。
- ecosystem summary。
- 与原论文/已知生物学一致点。
- scLucid 额外提供的 evidence/audit 价值。

## 完成标准

Phase 4 通过条件：

- 至少 3 个 tumor case studies。
- 每个 case study 有 malignancy/TME/ecotype summary。
- 至少 1 个 case 使用 bulk 或 spatial evidence module。
- 肿瘤解释输出不依赖单一 UMAP。
- 所有 tumor calls 都有 confidence 和 limitations。
- Figure 4 草图完成。

## 交付物

- `validation/tumor_case_studies/` scripts/notebooks。
- tumor program gene sets。
- tumor interpretation review summary。
- case-study markdown report。
- Figure 4 绘图脚本。

## 推荐主图

Figure 4A: tumor interpretation workflow。  
Figure 4B: malignancy evidence and CNV/signature support。  
Figure 4C: TME state landscape across patients。  
Figure 4D: ecotype/cellular ecosystem summary。  
Figure 4E: bulk/spatial concordance example。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 肿瘤解释太散 | 功能很多但无主张 | 聚焦 malignancy、TME states、ecotype summary |
| gene sets 主观 | 审稿人质疑来源 | 每个 gene set 记录来源和版本 |
| ecotype 太像过度概念化 | 数据支持不足 | 默认标注 exploratory，强调 summary 而非 definitive subtype |
| bulk/spatial 证据薄 | paired data 少 | 作为 optional evidence，不做主结论基础 |

