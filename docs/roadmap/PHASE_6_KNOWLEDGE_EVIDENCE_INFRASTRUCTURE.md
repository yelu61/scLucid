# Phase 6: Knowledge & Evidence Infrastructure

**Recommended duration**: 持续推进，P0 后启动系统化建设  
**Primary output**: source-aware knowledge schemas and evidence routing  
**Main claim**: scLucid 的解释能力来自可审计的知识基础设施，而不是散落的 marker list 或黑箱 annotation。

## 目标

统一管理 scLucid 的知识资产：

- Marker Manager。
- Gene Set Manager。
- Atlas Manager。
- Literature Evidence Manager。
- Ontology Manager。
- Therapy Knowledge Manager。
- LLM Evidence Manager。

这些知识不应只是文件资源，而应成为 evidence contract 的一部分。

Phase 6 承接全项目的 evidence ontology：同一套 source、confidence、
evidence_level、claim_level、inference_level、limitations 和 review_status
字段应该能服务 annotation、program scoring、tumor interpretation、bulk/spatial
support evidence，以及未来的 report/agent interface。

## 准备

### 资源盘点

整理当前资源：

- TOML marker resources。
- gene set JSON/GMT。
- marker manager views。
- artifact/stress/cell-cycle markers。
- tumor markers。
- therapy signatures, 如果已有。
- LLM/CellTypist/reference evidence outputs。

### 统一元数据字段

每个 knowledge item 至少需要：

| Field | Meaning |
|-------|---------|
| `name` | resource or evidence item name |
| `type` | marker, gene_set, atlas, literature, ontology, therapy, llm |
| `species` | human/mouse/etc. |
| `tissue` | tissue context |
| `cancer_type` | optional cancer context |
| `source` | citation, package, atlas, curator |
| `version` | resource version |
| `evidence_level` | canonical, curated, atlas-derived, exploratory, suggested |
| `limitations` | scope limits |
| `review_status` | reviewed, needs_review, experimental |

## 具体步骤

### Step 1. 定义 knowledge infrastructure routing

资源按任务路由：

- lineage annotation。
- subtype annotation。
- state annotation。
- artifact review。
- program scoring。
- tumor interpretation。
- therapy interpretation。
- ecosystem modeling。

### Step 2. Marker and gene-set governance

建立规则：

- concise marker TOML 用于 annotation。
- broad program/pathway gene sets 放入 JSON/GMT。
- tumor-state signatures 必须记录来源。
- negative markers 与 artifact markers 必须可查询。
- mouse parity 后于 human route。

### Step 3. Atlas/reference evidence

将 CellTypist、SingleR-like、atlas reference、Seurat transfer 等输出统一成 evidence table：

- label。
- score。
- reference source。
- reference version。
- confidence。
- disagreement with marker evidence。
- review action。

### Step 4. Literature and ontology evidence

为未来高水平论文和 agent layer 准备：

- literature claim schema。
- ontology term mapping。
- tumor-state naming conventions。
- therapy signature provenance。
- analysis-intent terms, such as rare-cell discovery, tumor microenvironment
  interpretation, trajectory analysis, therapy response, and ecosystem
  modeling。
- claim-strength terms that downstream reports can use to avoid overstating
  exploratory results。

### Step 5. LLM evidence interface guardrails

LLM 输出只能作为：

- suggestion。
- explanation draft。
- literature summary。
- review prompt。

LLM 输出不能作为：

- final annotation。
- final malignancy call。
- validated ecotype。
- clinical prediction。

## 完成标准

Phase 6 通过条件：

- marker/gene-set/atlas/literature/ontology/therapy/LLM evidence 有统一 schema。
- knowledge infrastructure 可按任务路由。
- review summary 能引用 knowledge source。
- LLM evidence 有 provenance 和 limitations。
- 至少 annotation、program scoring、tumor interpretation 三类任务使用 knowledge routing。

## 交付物

- `docs/KNOWLEDGE_EVIDENCE_INFRASTRUCTURE.md`。
- resource schemas。
- knowledge routing tests。
- updated marker/gene-set governance。
- ontology and naming mapping draft。
- LLM evidence guardrail documentation。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 资源无限膨胀 | 维护困难 | 先做 schema 和 routing，不急于大量收录 |
| marker/gene set 混用 | annotation 和 program scoring 混乱 | concise marker vs broad signature 分层 |
| LLM 产生幻觉 | 未验证解释进入结论 | LLM 只能作为 evidence suggestion |
| ontology 过度工程化 | 花很多时间做命名 | 先支持高频 cell/tumor terms |
