# scLucid 投稿倒逼式 Phase Execution Playbook

**Status**: 当前执行手册  
**Updated**: 2026-06-12  
**Goal**: 按 Nature Methods 的方法学标准倒逼建设 scLucid，并以 Nature Computational Science / Genome Biology / Nature Communications 为主攻投稿目标。

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
| Phase 4 | `PHASE_4_TUMOR_INTERPRETATION_CASE_STUDIES.md` | scLucid 的肿瘤特色是否形成高水平研究叙事？ | 2-3 个月 |
| Phase 5 | `PHASE_5_R_PYTHON_PARITY_AND_METHOD_PORTS.md` | 成熟 R 方法 Python 化是否准确、必要、可验证？ | 1-2 个月 |
| Phase 6 | `PHASE_6_RELEASE_MANUSCRIPT_AND_SUBMISSION.md` | 是否达到投稿级 reproducible research package？ | 2-3 个月 |

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

