# scLucid 投稿倒逼式 Phase Execution Playbook

> **ARCHIVED / SUPERSEDED (2026-08-14)**
>
> This point-in-time execution playbook is retained for planning provenance.
> The canonical roadmap entry is now `docs/roadmap/index.md`; current
> implementation status is defined by code, tests, and
> `docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md`.

**Status**: 当前执行手册
**Updated**: 2026-06-30
**Goal**: 按 Nature Methods 的方法学标准倒逼建设 scLucid，并以 Nature Computational Science / Genome Biology / Nature Communications 为主攻投稿目标。

> Documentation note: this directory is the phase-level execution playbook.
> It records where scLucid is going and how claims should be validated. For
> the current user-facing module contract, start with
> `docs/user/module_features_and_plan.md` and
> `docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md`.

## 总策略

scLucid 的投稿主张不能是“又一个单细胞工具包”，而应是：

> scLucid is a diagnostic-first, evidence-audited framework for tumor single-cell interpretation, turning routine workflows into reproducible, reviewable, and biologically validated decision systems.

因此每个 phase 都必须服务三件事：

1. **方法学可信**：有清楚的输入输出契约、审计证据和失败模式。
2. **生物学有用**：能保护并解释肿瘤单细胞研究中的关键生物信号。
3. **投稿可证**：能形成 figure、benchmark table、case study 或 supplementary material。

## Phase 文档

| Phase | 文件 | 核心问题 | 推荐时长 |
|-------|------|----------|----------|
| Phase 1 | `PHASE_1_CORE_API_AND_CONTRACTS.md` | scLucid 的核心 API 和证据契约是否稳定？ | 2-4 周 |
| Phase 2 | `PHASE_2_QC_EVIDENCE_BENCHMARK.md` | QC 是否真的比固定阈值流程更安全、更可审计？ | 1-2 个月 |
| Phase 3 | `PHASE_3_PREPROCESS_ANALYSIS_VALIDATION.md` | preprocess/analysis 是否稳定保留生物信号？ | 1-2 个月 |
| Phase 4 | `PHASE_4_TUMOR_INTERPRETATION_CASE_STUDIES.md` | malignancy、tumor programs、TME states 是否形成可复查的肿瘤解释？ | 2-3 个月 |
| Phase 5 | `PHASE_5_TUMOR_ECOSYSTEM_MODELING.md` | 是否能从 cell-level 结果推进到 sample-level ecotype / ecosystem interpretation？ | 2-3 个月 |
| Phase 6 | `PHASE_6_KNOWLEDGE_EVIDENCE_INFRASTRUCTURE.md` | marker、gene set、atlas、literature、ontology、therapy、LLM evidence 是否形成统一知识基础设施？ | 持续推进 |
| Phase 7 | `PHASE_7_SUPPORT_EVIDENCE_AND_R_PYTHON_PARITY.md` | bulk/spatial/clinical/R parity 是否作为支持证据，而不是主线平台？ | 1-2 个月起 |
| Phase 8 | `PHASE_8_RELEASE_MANUSCRIPT_AND_SUBMISSION.md` | 是否达到投稿级 reproducible research package？ | 2-3 个月 |

## 与战略优先级的对应关系

| Strategic Priority | Roadmap Phases |
|--------------------|----------------|
| **P0: Core + Audit + Validation** | Phase 1, Phase 2, Phase 3 |
| **P1: Tumor Interpretation + Ecosystem + Scale** | Phase 4, Phase 5, engineering items across all phases |
| **P2: Knowledge Infrastructure + Support Evidence + Agent Interfaces** | Phase 6, Phase 7, long-term agent notes |

Agent-assisted interpretation is intentionally not a near-term implementation
phase. It is a long-term interface layer over mature evidence bundles after the
P0/P1 contracts are stable.

## Vision 到任务的转译

scLucid 的差异化 vision 可以写成一句话：

> 从“正确执行单细胞分析”升级为“在特定生物医学情境下，做出可审计、可证伪、可验证的科学解释决策”。

这个 vision 不新增平行路线图，而是压入现有 phase：

| Vision 方向 | Roadmap 落点 | 近期任务 |
|-------------|--------------|----------|
| 情境驱动参数和解释 | Phase 1, Phase 3, Phase 4 | 在 review summary 中稳定记录 `AnalysisContext`、分析目标、推荐理由、风险和人工复核点。 |
| QC 从过滤变成状态/风险归因 | Phase 2 | 比较固定阈值、scLucid adaptive、tumor-aware 等策略对 marker fidelity、tumor purity、immune/TME composition 的影响。 |
| 统一 claim/evidence/inference ontology | Phase 1, Phase 3, Phase 6, Phase 7 | 收敛 `claim_level`、`inference_level`、`evidence_level`、confidence、source、limitation 字段和文档定义。 |
| 肿瘤解释系统而非通用工具箱 | Phase 4, Phase 5 | 让 malignancy、CNV、tumor program、TME state、therapy signal、ecotype 输出都带 evidence source、confidence、limitations。 |
| 支持证据桥接 | Phase 6, Phase 7 | bulk、spatial、clinical、marker、atlas、literature 作为支持证据进入统一 evidence contract，而不是扩成新主产品线。 |
| 可产品化的 review surface | Phase 8 | 审计报告、Methods-ready 摘要、figure source data、release package 都必须体现 decision/rationale/risk/evidence/limitation。 |

Deferred ideas such as interactive review apps, natural-language interfaces,
foundation-model adapters, and fully declarative pipeline DSLs should remain
post-P0/P1 interface or adapter work until the evidence contracts and validation
figures are stable.

## 投稿目标与闸门

### 主攻目标

- **Genome Biology**：适合 genomics/post-genomics 方法和软件工具，要求有真实生物学应用与清楚资源价值。
- **Nature Communications**：适合方法 + 肿瘤 case study 都较强的工作。
- **Nature Computational Science**：适合强调 evidence-audited computational framework、benchmark 和工程可复现性的版本。

### 倒逼标准

按 Nature Methods 标准要求自己：

- 与现有方法系统比较，而不是只展示 scLucid 跑通。
- 有强验证，包括 synthetic、public real data、tumor case study、failure mode。
- 有重要生物问题应用，而不是工具截图。
- 有足够技术细节，方便别人立即应用和复现。

## 全局完成标准

投稿前至少满足：

- 3 个核心 workflow 可稳定复现：PBMC baseline、PDAC tumor、第二癌种 tumor。
- 6-10 个公开数据集完成 benchmark。
- 所有主图由 scripts/notebooks 一键重建。
- 所有关键方法输出进入 `adata.uns["sclucid"]` evidence contract。
- README、API docs、tutorial notebooks、release tag、Zenodo DOI 准备完成。
- 论文主张可用数据支持，而不是概念性宣称。
- Spatial、bulk、clinical 和 R/Python parity 保持为 support evidence，不变成新的主产品方向。
