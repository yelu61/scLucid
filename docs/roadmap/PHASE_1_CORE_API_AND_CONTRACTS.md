# Phase 1: 核心 API 与 Evidence Contract 冻结

**Recommended duration**: 2-4 周  
**Primary output**: v0.1 core workflow contract  
**Main figure support**: Figure 1 framework overview and workflow spine

## 目标

冻结 scLucid 的核心工作流边界，保证 `qc -> preprocess -> analysis -> tumor interpretation` 是稳定主线，`tools.bulk`、`tools.spatial` 和 R/Python parity 作为外部证据模块接入，而不是把项目扩散成大而全工具箱。

Phase 1 完成后，新用户应该能用稳定 API 跑通标准流程；高级用户应该能在 `adata.uns["sclucid"]` 中复查每一步关键决策。

Phase 1 还需要把 scLucid 的差异化 vision 压成稳定 contract：
每个核心步骤都应记录 context、decision、rationale、risk、evidence、
limitation 和 review action。这里先冻结字段和语义，不要求一次性完成
所有未来交互界面或高级模型。

## 准备

### 数据准备

- PBMC baseline：`data/pbmc3k.h5ad`。
- 肿瘤 baseline：`data/lin2020.pdac.h5ad` 或等价 PDAC 数据。
- 第二肿瘤候选：NSCLC / CRC / BRCA / HCC 中至少准备一个公开数据集。
- 每个数据集准备 metadata mapping：
  - sample key。
  - patient key, 如果存在。
  - known cell type annotation, 如果存在。
  - cancer type / tissue / species。

### 工程准备

- 确认核心 import 轻量：
  - `import scLucid as scl` 不应强依赖 heavy optional packages。
  - 缺少 scrublet、pyscdblfinder、squidpy、R bridge 时不能导致核心包不可用。
- 确认 workflow API 命名：
  - `scl.qc.run_standard_qc`
  - `scl.pp.run_preprocessing` 或项目当前等价入口。
  - `scl.analysis.run_standard_analysis`
  - `scl.tumor.malignancy.run_malignancy_interpretation`
- 确认 `adata.uns["sclucid"]` contract 不再频繁变化。

## 具体步骤

### Step 1. 定义核心 workflow contract

为每个核心模块列出：

- 输入要求：
  - `adata.X` 语义。
  - `layers["counts"]` 是否必须。
  - `obs` 必需列。
  - `var_names` 基因命名要求。
- 输出要求：
  - 新增 `obs`、`var`、`obsm`、`uns` 字段。
  - 字段命名。
  - 是否可写入 h5ad。
- 审计要求：
  - 参数记录。
  - 决策理由。
  - warning。
  - inference level。
  - review action items。

同时统一跨模块语义：

- `claim_level`：结果可支持的科学声明强度。
- `inference_level`：cell-level、sample-level、descriptive、exploratory 等推断边界。
- `evidence_level`：证据来源和成熟度，例如 validated_core、curated、heuristic、exploratory、unavailable。
- `confidence`：模块内部对当前输出可靠性的紧凑判断。
- `limitations`：缺失参考、样本量不足、optional dependency 不可用、confounding 等限制。

建议产出表格：

| Module | Required input | Required output | Evidence key | Failure mode |
|--------|----------------|-----------------|--------------|--------------|
| QC | counts, sample key | QC metrics, flags | `sclucid.qc.review_summary` | missing counts |
| Preprocess | counts layer, QC flags | PCA, neighbors, UMAP | `sclucid.preprocess.review_summary` | invalid layer semantics |
| Analysis | graph, embeddings | clusters, markers, annotations | `sclucid.analysis.review_summary` | weak marker evidence |
| Tumor | annotations, marker/CNV evidence | malignancy/TME calls | `sclucid.tumor.review_summary` | insufficient tumor context |

### Step 2. 固定 public API

整理用户最常用入口：

- beginner workflow。
- composable API。
- advanced audit notebook。

每个入口需要：

- 最小示例。
- 参数说明。
- 失败时错误信息。
- 输出位置说明。

不稳定或实验性 API 需要标记：

- `experimental`
- `requires_optional_dependency`
- `not_yet_benchmark_validated`

### Step 3. 完成 core smoke tests

增加或更新测试：

- import smoke test。
- PBMC mini workflow smoke test。
- tumor mini workflow smoke test。
- missing optional dependency test。
- h5ad round-trip test。

### Step 4. 建立文档入口

README 必须回答：

- scLucid 是什么。
- scLucid 不是什么。
- 核心 workflow 是什么。
- 证据契约在哪里。
- 如何跑最小例子。

`docs/README.md` 必须指向：

- 战略总纲。
- phase playbook。
- data usage。
- bulk/spatial design。
- marker resource。

## 完成标准

Phase 1 通过条件：

- `pytest tests/qc tests/preprocess tests/analysis -q --no-cov` 核心通过，或记录明确的非阻塞失败。
- `ruff check` 通过关键源码和新增测试。
- PBMC 和 PDAC 最小 workflow 可复现。
- 每个核心模块有稳定 review summary。
- README 中有框架图和定位说明。
- `docs/roadmap/` 文档可作为执行索引。

## 交付物

- `docs/roadmap/PHASE_1_CORE_API_AND_CONTRACTS.md`
- API contract checklist。
- PBMC minimal workflow notebook/script。
- PDAC minimal workflow notebook/script。
- README framework figure。
- CI/smoke test 更新。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| API 过早冻结 | 后续发现字段不够 | 使用 schema_version 和 migration note |
| optional dependency 污染核心 | import 失败 | lazy import + clear fallback |
| 功能边界扩散 | bulk/spatial/tumor 混入核心 | 保持 tools/tumor 分层 |
| 文档与代码漂移 | README 写法不能跑 | 文档示例纳入 smoke test |
