# Phase 5: Tumor Ecosystem Modeling

**Recommended duration**: 2-3 个月  
**Primary output**: sample-level ecosystem feature matrix, ecotype prototype, Figure 5 candidate  
**Main claim**: scLucid 的长期学术壁垒来自从 cell-level annotation 走向 sample-level / patient-level tumor ecosystem interpretation。

## 目标

把 Phase 4 的 tumor state interpretation 进一步聚合为：

- cell communities。
- sample-level tumor states。
- microenvironment archetypes。
- ecotype-style summaries。
- patient stratification。
- optional spatial/bulk/clinical support evidence。

这一 phase 的重点不是做一个通用 clustering 工具，而是建立可审计的 tumor ecosystem interpretation object。

## 准备

### 输入数据

至少需要 2-3 个 tumor case studies，其中至少一个具备多病人或多样本结构：

- cell type / state annotation。
- malignancy calls。
- tumor program scores。
- TME state scores。
- sample/patient metadata。
- optional spatial zones。
- optional bulk deconvolution or pseudobulk data。
- optional clinical response/outcome metadata。

### Ecosystem acceptance datasets

优先建立三个 acceptance routes：

1. **PDAC ecosystem**
   - 目标：CAF/myeloid/TME-rich ecosystem interpretation。
   - 重点：malignant-stromal-immune co-abundance、CAF states、suppressive myeloid patterns。
2. **NSCLC ecosystem**
   - 目标：immune infiltration/exclusion and exhaustion archetypes。
   - 重点：T/NK exhaustion、cytotoxicity、myeloid suppression、therapy-relevant immune context。
3. **CRC ecosystem**
   - 目标：epithelial state, inflammation, stromal/immune co-abundance。
   - 重点：tumor epithelial programs、inflammatory TME、stromal remodeling。

### 前置条件

- Phase 1-3 核心 workflow 和 review summary 稳定。
- Phase 4 已能输出 malignancy、program、TME state evidence。
- 所有输入 features 都有 evidence source、confidence 和 limitations。

## 具体步骤

### Step 1. 定义 ecosystem feature schema

样本或区域级 feature matrix 至少包括：

- cell type proportions。
- malignant program averages。
- TME state abundance。
- immune exhaustion/cytotoxicity balance。
- myeloid suppressive/inflammatory balance。
- CAF state abundance。
- malignancy confidence summary。
- optional spatial niche/zone scores。
- optional bulk concordance scores。
- optional clinical metadata。

输出建议：

```text
adata.uns["sclucid"]["tumor"]["ecosystem_features"]
adata.uns["sclucid"]["tumor"]["ecosystem_review_summary"]
```

### Step 2. Cell community construction

从 cell-level 结果构建 community evidence：

- cluster/state composition。
- neighborhood or co-abundance pattern。
- sample-level co-occurrence matrix。
- patient-level composition similarity。

注意：

- community 不等于 spatial neighborhood，除非有 spatial evidence。
- 没有空间数据时应称为 co-abundance community 或 ecosystem component。

### Step 3. Ecotype / archetype discovery

可选方法：

- NMF。
- consensus clustering。
- hierarchical clustering。
- graph community detection。
- manually defined archetype score。

每个 ecotype/archetype 输出：

- assignment。
- confidence。
- stability。
- dominant programs。
- dominant cell states。
- sample membership。
- limitations。
- exploratory/sample-level inference status。

### Step 4. Stability and validation

至少评估：

- resampling stability。
- feature sensitivity。
- cancer-type specificity。
- sample-size limitation。
- agreement with known paper subtype, if available。
- association with clinical metadata, if available。
- program/state co-occurrence consistency。
- patient-level separability when metadata support it。
- optional spatial/bulk/clinical support evidence concordance。

### Step 5. Visualization templates

优先做：

- ecosystem heatmap。
- sample-level archetype map。
- ecotype composition stacked bar。
- program/state contribution plot。
- patient stratification summary。
- optional spatial region overlay。

## 完成标准

Phase 5 通过条件：

- 至少 1 个 tumor dataset 生成完整 ecosystem feature matrix。
- 至少 1 个 case study 生成 ecotype/archetype prototype。
- 所有 ecotype 输出标注 exploratory 或 sample-level inference status。
- 有 stability summary。
- 有 review summary 和 limitations。
- 不把 ecotype 当作未经验证的临床 subtype。

## 交付物

- `validation/tumor_ecosystem/` scripts/notebooks。
- ecosystem feature schema。
- ecotype/archetype prototype。
- stability report。
- ecosystem visualization templates。
- manuscript figure source data。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 样本量不足 | ecotype 不稳定 | 默认 exploratory，强调 prototype |
| 概念过度包装 | 被审稿人认为过度解释 | 输出 limitations 和 stability |
| spatial 证据缺失 | community 被误解为空间结构 | 明确 co-abundance vs spatial niche |
| 临床关联过度 | 小样本也谈预测 | 只谈 association，不谈 prediction |
