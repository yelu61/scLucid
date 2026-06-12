# Phase 6: Release、Manuscript 与投稿包

**Recommended duration**: 2-3 个月  
**Primary output**: manuscript-ready release package  
**Main claim**: scLucid 已达到可复现、可审计、可验证的投稿级工具/框架标准。

## 目标

把前 5 个 phase 的工程、验证和生物学案例收敛成投稿资产：

- release tag。
- Zenodo DOI。
- reproducible benchmark repository。
- manuscript figures。
- supplementary tables。
- tutorial notebooks。
- documentation website。
- submission-ready manuscript。

## 准备

### 论文目标选择

按结果强度选择：

- 如果 validation 系统性强、framework 贡献突出：Nature Computational Science。
- 如果 genomics software resource + 多癌种 case study 强：Genome Biology。
- 如果生物学 case study 更强：Nature Communications。
- 如果方法学创新和性能比较极强：Nature Methods 预投稿咨询。

### 资料准备

- 所有 datasets accession。
- 所有 software versions。
- 所有 random seeds。
- 所有 figure source data。
- 所有 benchmark scripts。
- package release notes。

## 具体步骤

### Step 1. 冻结投稿版本

创建 release branch 或 tag 前完成：

- 版本号。
- changelog。
- dependency matrix。
- optional extras。
- API docs。
- deprecated APIs 清理。
- known limitations。

建议版本：

- `v0.2.0` for preprint/submission candidate。
- `v1.0.0` 只在 API 和 validation 稳定后使用。

### Step 2. Reproducible benchmark package

目录建议：

```text
validation/
  datasets/
  qc/
  preprocess_analysis/
  tumor_case_studies/
  r_python_parity/
  figures/
  tables/
```

每个 benchmark 需要：

- config file。
- run script。
- output summary。
- figure source data。
- environment note。

### Step 3. 主图与补充图

主图建议：

- Figure 1: scLucid framework overview。
- Figure 2: QC evidence benchmark。
- Figure 3: preprocess/analysis stability。
- Figure 4: tumor interpretation case studies。
- Figure 5: engineering, reproducibility, and scalability。

补充图：

- dataset summary。
- QC threshold distributions。
- doublet method parity。
- marker fidelity by dataset。
- annotation confusion/agreement。
- runtime/memory。
- failure mode examples。

### Step 4. Manuscript skeleton

建议结构：

1. Introduction
   - single-cell tumor workflows are powerful but hard to audit。
   - current tools are excellent but fragmented。
   - scLucid introduces evidence-audited interpretation。
2. Results
   - framework overview。
   - QC benchmark。
   - preprocess/analysis validation。
   - tumor case studies。
   - engineering/reproducibility。
3. Discussion
   - scope and limitations。
   - not replacing specialist tools。
   - future bulk/spatial/multi-omics evidence。
4. Methods
   - workflow contracts。
   - benchmark design。
   - datasets。
   - statistics。
   - implementation。

### Step 5. 投稿前审计

检查清单：

- 主图是否可一键重建。
- figure source data 是否完整。
- 所有 claims 是否有数据支持。
- 所有 benchmark 是否有 baseline。
- 所有 optional methods 是否标注依赖。
- 所有 inference claims 是否不过度。
- GitHub README 是否能让审稿人快速理解。
- notebooks 是否能在干净环境运行。

### Step 6. Preprint and submission

推荐流程：

1. 内部 frozen release。
2. bioRxiv/arXiv preprint, 视团队策略。
3. Nature Computational Science / Genome Biology / Nature Communications 预投稿咨询。
4. 根据编辑反馈选择正式投稿。
5. 准备 reviewer response scaffold。

## 完成标准

Phase 6 通过条件：

- release tag 创建。
- DOI 可引用。
- documentation online。
- benchmark scripts 可运行。
- manuscript figures 完成。
- supplementary tables 完成。
- paper draft 可送合作者审阅。

## 交付物

- release tag。
- Zenodo DOI。
- manuscript draft。
- figure source data。
- supplementary tables。
- reproducibility package。
- submission cover letter draft。

## Cover letter 核心信息

建议强调：

- scLucid solves an auditability and interpretation gap in tumor single-cell workflows。
- It complements, rather than replaces, Scanpy/Seurat/specialist tools。
- It provides evidence contracts, inference semantics, and tumor-focused interpretation。
- It is validated across layered benchmarks and tumor case studies。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 工具包感太强 | 编辑认为 novelty 不够 | 强调 evidence-audited framework 和 tumor-specific validation |
| case study 太弱 | 像软件说明书 | 至少 3 个 tumor datasets，提炼 biological insight |
| benchmark 不够公平 | 审稿人质疑 | baseline scripts 开源，参数透明 |
| 复现困难 | 审稿人跑不动 | 最小数据、Docker/conda、cached outputs |

