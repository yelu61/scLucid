# Phase 7: Support Evidence Modules And R/Python Parity

**Recommended duration**: 1-2 个月起，持续维护  
**Primary output**: validated support evidence adapters and parity matrices  
**Main claim**: bulk、spatial、clinical 和 R/Python parity 是肿瘤解释的支持证据，不是 scLucid 的主产品方向。

## 目标

把外部证据模块统一纳入 scLucid audit contract：

- bulk evidence。
- spatial evidence。
- clinical evidence。
- selected R/Python parity。

同时明确边界，避免 scLucid 变成 bulk platform、spatial platform 或 R 包重写项目。

## Support Layer 1: Bulk Evidence

### 目标

支持：

- pseudobulk validation。
- bulk deconvolution。
- bulk-pseudobulk concordance。
- sample-level DE support。

### 完成标准

- bulk 输出有 inference level。
- 与 single-cell interpretation 的连接明确。
- 不把 scLucid 扩展为完整 bulk RNA-seq 平台。

## Support Layer 2: Spatial Evidence

### 目标

Spatial 只作为 tumor interpretation evidence：

- spatial organization of tumor programs。
- TME infiltration/exclusion。
- tumor-stroma boundaries。
- spatial niches。
- region-level ecosystem support。

### 明确不做

- general spatial analysis platform。
- image segmentation platform。
- broad spatial alignment suite。
- full deep spatial wrapper ecosystem。

### 完成标准

- spatial 输出标注为 spatial evidence。
- spatial niche/community 与 non-spatial co-abundance community 区分清楚。
- 至少一个 tumor case study 使用 spatial evidence 支持解释。

## Support Layer 3: Clinical Evidence

### 目标

连接 sample-level ecosystem features 和 clinical metadata：

- response association。
- outcome association。
- subtype context。
- translational interpretation。

### 完成标准

- clinical 输出只在 sample-level 分析。
- 区分 association 和 prediction。
- 小样本不做预测性宣称。

## Support Layer 4: Selected R/Python Parity

### 方法优先级

Priority 1:

- scDblFinder / pyscdblfinder doublet evidence。
- Scrublet compatibility。
- CNV/malignancy evidence consumption。
- pseudobulk/sample-level DE evidence。

Priority 2:

- reference annotation evidence。
- program/pathway scoring parity。
- selected deconvolution evidence。

Defer:

- broad trajectory parity。
- general cell-cell communication wrappers。
- deep spatial wrappers。
- methods not tied to tumor interpretation evidence。

## 具体步骤

### Step 1. 定义接入等级

每个方法标记为：

- `wrapped`。
- `ported`。
- `approximated`。
- `consumed`。

### Step 2. 建立 parity matrix

每个方法必须记录：

- reference implementation。
- Python implementation。
- input contract。
- output contract。
- validation data。
- acceptance metrics。
- limitations。
- license status。

### Step 3. Synthetic + real validation

至少：

- one synthetic test。
- one real-data smoke test。
- missing dependency fallback test。

### Step 4. Evidence contract integration

输出进入：

- `adata.obs` or `adata.uns`。
- method evidence table。
- review summary。
- limitation and dependency metadata。

## 完成标准

Phase 7 通过条件：

- bulk/spatial/clinical/R parity 都被定义为 support evidence。
- scDblFinder parity 或 documented pyscdblfinder behavior 完成。
- 每个接入方法有 validation 或明确 experimental 标记。
- 没有 unvalidated method 被推荐为默认结论。
- Spatial 边界在 docs 和 README 中一致。

## 交付物

- support evidence matrix。
- R/Python parity matrix。
- method evidence schemas。
- validation notebooks。
- optional dependency documentation。
- spatial evidence boundary documentation。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 支持模块喧宾夺主 | 项目变成 multi-omics platform | support evidence guardrails |
| spatial 功能越做越宽 | 审稿人按 spatial platform 要求评估 | 只服务 tumor interpretation |
| R parity 无限扩展 | 维护成本爆炸 | selected high-value only |
| clinical 过度宣称 | 小样本谈预测 | association only unless properly validated |

