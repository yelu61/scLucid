# scLucid 五方向战略梳理与实施计划

**Status**: 当前战略总纲  
**Updated**: 2026-06-12  
**Scope**: 围绕 scLucid 作为 tumor single-cell interpretation system 的核心定位，整理后续 5 个推荐方向、实施路径、交付物和验证标准。

**Execution Playbook**: 投稿倒逼式的逐阶段执行细则见 `docs/roadmap/README.md`。本文档保留战略方向，roadmap 文档负责每个 phase 的准备、步骤、验收标准、交付物和风险控制。

## 0. 总定位

scLucid 的差异化不应是“比 Scanpy/Seurat 多几个函数”，也不应是“复制 OmicVerse 式大而全生态”。更准确的定位是：

> scLucid 是一个 diagnostic-first、audit-ready、tumor-focused 的单细胞研究解释系统。它让肿瘤单细胞分析更 lucid：假设清楚、证据可追溯、推断边界明确、结果适合项目复查和论文级解释。

核心边界：

- `qc + preprocess + analysis` 是核心单细胞工作流主干。
- `tumor` 是肿瘤解释层，消费 analysis 产物并整合 CNV、malignancy、TME、therapy、heterogeneity、ecosystem 等证据。
- `tools.bulk` 和 `tools.spatial` 是外部证据增强模块，服务于肿瘤单细胞解释，不追求泛 multi-omics 全覆盖。
- R/Python parity 只做精选成熟方法的 Python 接入或复现，并必须进入 scLucid 的 evidence/audit contract。
- 工程化目标不是“功能最多”，而是轻依赖、可扩展、可验证、可复现、能处理真实大数据集。

## 1. 方向一：核心 QC -> Preprocess -> Analysis 工作流打磨

### 目标

把 `qc + preprocess + analysis` 打磨成 scLucid 的稳定核心路径，使其可以在真实项目中提供可靠 handoff、清晰审计记录和保守可解释的结果。

这不是继续堆方法，而是让每个核心阶段都做到：

- 输入输出契约清楚。
- 默认路径轻依赖且稳定。
- 每个关键决策有诊断依据。
- 所有重要参数、警告和证据写入 `adata.uns["sclucid"]`。
- 结果能被 audit report、notebook 和人工 review 复查。

### 实施重点

#### 1.1 QC 模块

当前优势：

- adaptive threshold 思路已经成型。
- doublet heuristic、Scrublet optional evidence、tumor-aware warnings 已有基础。
- review summary 和 audit report 方向明确。

下一步：

- 梳理 QC 输出 contract：
  - `adata.obs` 必需列。
  - `adata.uns["sclucid"]["qc"]` 必需字段。
  - filtering 前后统计。
  - threshold rationale。
- 强化 tumor-aware QC：
  - 高 mitochondrial/stress 是否是坏细胞还是肿瘤状态。
  - cell cycle 高表达是否提示增殖肿瘤群。
  - doublet 的 tumor-normal 混合风险。
- 建立 QC acceptance tests：
  - PBMC baseline。
  - PDAC tumor sample。
  - 高线粒体、高 doublet、低质量边界样本。

交付物：

- `docs/source/qc_preprocess_maturity.rst` 更新。
- `tests/integration/test_qc_contract.py`。
- `scripts/run_qc_acceptance.py`。
- QC audit report 示例。

#### 1.2 Preprocess 模块

当前优势：

- counts/layers contract 已有基础。
- normalization、HVG、PCA、neighbors、UMAP 路径完整。
- batch correction 已经倾向显式 opt-in。

下一步：

- 明确 preprocessing handoff contract：
  - `layers["counts"]` 是否保留。
  - `X` 当前代表什么。
  - `raw` 是否设置。
  - PCA/neighbors/UMAP key 命名。
- 强化 HVG 解释：
  - 标准 HVG。
  - tumor-aware retained genes。
  - marker/program genes 是否被保留。
- 对 batch correction 做诊断优先：
  - 不默认校正。
  - 先评估 batch mixing 与 biological separation。
  - 明确过度校正风险。

交付物：

- `docs/source/data_contracts.rst` 增加 preprocessing contract。
- HVG stability / tumor marker retention 测试。
- `scripts/run_preprocess_acceptance.py`。

#### 1.3 Analysis 模块

当前优势：

- clustering review、marker evidence、annotation evidence、consensus、posthoc QC review 已经形成闭环。
- marker manager 和资源视图是 scLucid 的关键差异化。

下一步：

- 继续把 analysis 做成第二个 benchmark module。
- 强化 annotation evidence table：
  - marker-based evidence。
  - reference/CellTypist evidence。
  - negative markers。
  - artifact/stress/cell-cycle flags。
  - confidence 与 review status。
- DE/proportion 输出统一：
  - `inference_level`。
  - `valid_for_publication_inference`。
  - sample-level vs exploratory cell-level。
- 明确 first-pass annotation 策略：
  - 先 major lineage / compartment。
  - subtype/state 依赖 subset analysis 或用户显式请求。

交付物：

- `tests/integration/test_analysis_contract.py`。
- `scripts/run_analysis_acceptance.py` 扩展。
- analysis audit report 示例。

### 成功标准

- PBMC 和 PDAC golden path 可重复运行。
- 核心输出在 `adata.uns["sclucid"]` 中可追溯。
- 新用户能用 workflow API 跑通，专家能用 advanced notebook 复查每个决策。
- 不因 optional heavy dependency 缺失导致 core import 失败。

## 2. 方向二：分层验证，而不是只证明跑通

### 目标

建立 scLucid 与 Scanpy/Seurat/成熟专用工具的分层验证框架。scLucid 不应直接宣称“更好更准”，而应证明：

- 在哪些任务上结果等价。
- 在哪些决策上更可解释。
- 在哪些肿瘤场景中更安全。
- 在哪些工作流中更省人工、更适合审计。

### 分层验证框架

#### Level 1: Execution Parity

问题：同样的数据，scLucid 是否能稳定产生 Scanpy/Seurat 可比的基础结果？

指标：

- cell retention rate。
- detected genes / counts distributions。
- HVG overlap。
- PCA/UMAP neighborhood preservation。
- clustering resolution 可比性。
- major cell type annotation 一致性。

交付物：

- `validation_outputs/parity_scanpy/`。
- Scanpy baseline script。
- Seurat baseline notebook 或 R script。

#### Level 2: Decision Quality

问题：scLucid 的自动建议是否更容易解释和复查？

指标：

- QC threshold rationale 是否完整。
- batch correction 是否有诊断依据。
- clustering resolution 是否有 stability/marker evidence。
- annotation 是否有多证据表。

交付物：

- human-review checklist。
- decision audit score。
- `docs/VALIDATION_DECISION_QUALITY.md`。

#### Level 3: Biological Concordance

问题：scLucid 的结果是否符合已知生物学和文献？

指标：

- PBMC major lineage 期望结构。
- PDAC malignant epithelial / CAF / myeloid / T/NK / endothelial compartment 是否合理。
- canonical marker expression 是否匹配。
- tumor program / TME signature 是否符合文献。

交付物：

- `docs/VALIDATION_BIOLOGICAL_CONCORDANCE.md`。
- curated expected marker panels。
- dataset-specific acceptance reports。

#### Level 4: Inference Safety

问题：scLucid 是否减少常见统计误用？

重点：

- cell-level DE 不冒充 biological replicate inference。
- pseudobulk / sample-level DE 明确标注。
- proportion analysis 区分 descriptive、exploratory、sample_level。
- bulk/spatial 结果默认保守标注。

交付物：

- inference semantics test suite。
- example notebook: bad practice vs scLucid guarded workflow。

#### Level 5: Usability And Engineering

问题：scLucid 是否真的更方便、更可复现？

指标：

- 完成同一 workflow 的代码行数。
- 需要人工记录的参数数量。
- audit report 完整性。
- runtime / memory。
- failure message clarity。

交付物：

- `docs/VALIDATION_USABILITY_ENGINEERING.md`。
- benchmark table。

### 成功标准

- README 中的差异化声明都有验证材料支持。
- 每个 benchmark dataset 都有 manifest、final h5ad、figures、audit report、acceptance summary。
- 对 Scanpy/Seurat 的比较措辞严谨：better / equivalent / weaker 都记录。

## 3. 方向三：精选 R/Python Parity 与成熟方法接入

### 目标

打通 Python/R 语言壁垒，但不做“为了 port 而 port”。只选择肿瘤单细胞研究中高频、成熟、Python 生态缺口明显、可验证的方法。

### 方法选择标准

纳入前必须回答：

1. 这个方法是否在高水平肿瘤单细胞研究中常见？
2. Python 生态是否已有可靠替代？
3. 结果是否能进入 scLucid 的 evidence/audit contract？
4. 是否有 benchmark 或可复现预期结果？
5. 依赖是否能隔离为 optional extra？
6. license 是否允许包装或复现？

### 优先级

#### P0：直接影响推断安全和肿瘤解释

- pseudobulk/sample-level DE：
  - DESeq2。
  - edgeR。
  - limma-voom。
  - muscat。
- CNV / malignancy evidence：
  - inferCNV-like evidence。
  - CopyKAT-like calls。
- doublet evidence：
  - ScDblFinder parity or wrapper。
- bulk deconvolution：
  - BayesPrism/DWLS 当前已内置 Python 版本。
  - 后续只做精选补充。

#### P1：高频研究模块

- cell-cell communication：
  - CellChat。
  - NicheNet。
  - CellPhoneDB。
- compositional analysis：
  - scCODA。
  - propeller / Milo-like ideas。
- pathway/program scoring：
  - AUCell/GSVA/ssGSEA 等可验证路线。

#### P2：复杂但有潜力

- trajectory / lineage：
  - Monocle3 parity。
  - Slingshot-like workflows。
- spatial methods：
  - Tangram。
  - cell2location。
  - RCTD。
  - STAGATE / STAligner 等 deep spatial methods 延后。

### 实施方式

三类实现模式：

1. **Optional wrapper**
   - 调用原工具或成熟 Python/R bridge。
   - 记录版本、参数、命令、输入输出。
2. **Clean-room reimplementation**
   - 只针对公开公式、通用统计或基础算法。
   - 不复制 GPL 或不兼容 license 代码。
3. **Parity adapter**
   - 将外部结果标准化为 scLucid evidence table。
   - 不重写方法本身。

### 验证标准

- 每个 port/wrapper 必须有 parity matrix。
- 至少一个 synthetic test。
- 至少一个 real-data smoke test。
- 记录与原工具的差异：
  - 输入限制。
  - 输出字段差异。
  - 不支持功能。
  - 推荐使用场景。

交付物：

- `docs/source/r_parity.rst` 持续更新。
- `docs/R_PARITY_SELECTION_MATRIX.md`。
- `tests/tools/test_r_parity_contracts.py`。
- method-specific validation notebooks。

## 4. 方向四：加强肿瘤单细胞高水平研究概念、图表和解释层

### 目标

scLucid 的核心特色应该是 tumor research system，而不仅是标准 scRNA workflow。需要系统引入肿瘤单细胞论文中常见的研究对象、概念、可视化和解释模板。

### 核心概念层

#### 4.1 Malignant vs Non-malignant

功能：

- malignant epithelial identification。
- CNV evidence integration。
- tumor marker evidence。
- normal epithelial vs malignant epithelial separation。
- suspect/unresolved status。

输出：

- malignancy evidence table。
- malignant call confidence。
- reason strings。
- review-required flags。

#### 4.2 Tumor Programs

常见 programs：

- proliferation。
- EMT。
- hypoxia。
- interferon response。
- stress。
- apoptosis。
- stemness。
- antigen presentation。
- metabolic rewiring。

输出：

- program score matrix。
- tumor program heatmap。
- program-by-cluster summary。
- program vs clinical response association。

#### 4.3 TME Composition And States

重点细胞群：

- T/NK exhaustion and cytotoxicity。
- macrophage / monocyte / DC states。
- neutrophil states。
- CAF subtypes。
- endothelial states。
- B/plasma states。

输出：

- compartment summary。
- cell type proportion shift。
- state score plots。
- patient-level TME profile。

#### 4.4 Tumor Ecosystem / Ecotype

生态型是 scLucid 可以形成特色的高级方向。

定义建议：

> Ecosystem/ecotype 是由 malignant programs、immune states、stromal states、spatial zones、clinical phenotype 共同构成的样本或区域级肿瘤微环境模式。

实施路线：

- 从 cell-level annotation 聚合到 sample-level / region-level。
- 输入：
  - malignant program scores。
  - TME proportions。
  - immune exhaustion / myeloid / CAF states。
  - CNV/malignancy confidence。
  - optional spatial zones。
  - optional clinical metadata。
- 方法：
  - NMF / clustering / consensus clustering。
  - stability analysis。
  - marker/program interpretation。
- 输出：
  - ecotype assignment。
  - ecotype signature。
  - ecotype composition heatmap。
  - ecotype clinical association。
  - review summary。

建议模块：

- `src/scLucid/tumor/microenvironment/ecosystem.py`
- `src/scLucid/tumor/heterogeneity/programs.py`
- `src/scLucid/plotting/tumor_ecosystem.py`

#### 4.5 Spatial Tumor Concepts

结合 `tools.spatial`：

- tumor-stroma boundary。
- immune infiltration score。
- spatial niches。
- tissue zones。
- malignant program spatial distribution。
- ICI response signature spatial enrichment。

输出：

- spatial niche map。
- boundary score。
- infiltration gradient。
- zone-by-program heatmap。

### 可视化模板

优先添加：

- tumor compartment UMAP。
- malignancy evidence heatmap。
- CNV score violin / heatmap。
- TME composition stacked bar。
- patient-level TME heatmap。
- tumor program dotplot。
- ecotype composition heatmap。
- ecotype Sankey/alluvial。
- spatial niche plot。
- tumor-stroma boundary plot。
- bulk-pseudobulk concordance scatter。

### 成功标准

- 用户可以从一个 tumor AnnData 直接生成：
  - major annotation。
  - malignancy evidence。
  - TME composition。
  - tumor programs。
  - optional ecotype summary。
  - audit report。
- 每个高级概念都有“输入、方法、输出、限制、推断语义”说明。
- 不把 exploratory programs/ecotypes 过度包装成 causal conclusion。

## 5. 方向五：工程化、性能和大数据能力

### 目标

scLucid 要在真实大数据项目中可靠、可复现、方便交付。工程化不是附属项，而是 scLucid 能否成为研究系统的关键。

### 工程原则

- Core import must be lightweight。
- Optional dependencies must be isolated。
- Sparse matrices should stay sparse whenever possible。
- Large data workflows should avoid unnecessary `.toarray()`。
- Every workflow should be reproducible from config。
- Errors should teach the user what to fix。
- Reports and manifests should make project handoff easy。

### 实施重点

#### 5.1 Dependency And Import Hygiene

- core package 不依赖 heavy spatial/deep/R packages。
- `tools.bulk` / `tools.spatial` 可以 import，但 heavy backend 函数内 lazy import。
- optional extras：
  - `analysis`
  - `bulk`
  - `spatial`
  - `tools`
  - `spatial-deep`
  - future `r-parity`

交付物：

- import smoke tests。
- optional dependency missing tests。

#### 5.2 Sparse And Memory Awareness

重点排查：

- `.toarray()` 使用点。
- 大矩阵复制。
- DataFrame 转换。
- all-gene operations。

策略：

- 对大数据默认使用 sparse-aware 计算。
- 对需要 dense 的函数增加 warning 或 chunking。
- 对 plotting 函数限制默认 features。

交付物：

- `docs/PERFORMANCE_MEMORY_GUIDE.md`。
- memory benchmark scripts。

#### 5.3 Runtime Benchmarks

基线：

- Scanpy standard pipeline。
- scLucid workflow pipeline。

数据：

- PBMC3k。
- 10k / 50k / 100k simulated or real subset。
- PDAC tumor dataset。

指标：

- runtime。
- peak memory。
- output artifact size。
- failure rate。

交付物：

- `benchmarks/` scripts。
- `docs/PERFORMANCE_BENCHMARKS.md`。

#### 5.4 User Experience

目标：

- 新手可用 workflow API。
- 项目分析者可用 simple API。
- 高级用户可用 advanced notebooks。

重点：

- clear error messages。
- progress display。
- manifest files。
- HTML audit report。
- compact project outputs。

交付物：

- updated examples。
- workflow templates。
- report templates。

## 6. 总体优先级建议

### P0：当前最应该做

1. 固化 `qc -> preprocess -> analysis` handoff contract。
2. 建立 PBMC + PDAC real-data acceptance gates。
3. 完成 analysis evidence/annotation review loop。
4. 明确 README 和 docs 中的 scLucid 定位。
5. 修正 bulk/spatial/tools 的 namespace 与 optional dependency hygiene。

### P1：下一阶段主攻

1. 分层验证框架。
2. tumor interpretation API 和 plots。
3. selected R/Python parity matrix。
4. bulk-pseudobulk-spatial concordance。
5. memory/runtime benchmark。

### P2：长期增强

1. ecotype/ecosystem 高级解释层。
2. spatial deep wrappers。
3. second/third tumor type acceptance workflows。
4. clinical response association。
5. polished publication figure templates。

## 7. 推荐里程碑

### Milestone A：Core Lucid Workflow

时间建议：2-4 周。

交付：

- PBMC golden path。
- PDAC golden path。
- QC/preprocess/analysis contracts。
- audit report examples。
- README 定位更新。

验收：

- 核心 workflow 可重复运行。
- audit summary 完整。
- optional heavy deps 缺失不影响 core import。

### Milestone B：Layered Validation

时间建议：4-6 周。

交付：

- Scanpy baseline comparison。
- Seurat baseline comparison。
- decision quality checklist。
- inference safety test suite。

验收：

- 至少 2 个数据集完成分层验证。
- 有明确 better/equivalent/weaker 记录。

### Milestone C：Tumor Interpretation System

时间建议：6-10 周。

交付：

- malignancy evidence table。
- tumor program scoring。
- TME composition summary。
- therapy/ICI signature skeleton。
- tumor-specific plots。

验收：

- PDAC workflow 可以输出完整 tumor interpretation report。

### Milestone D：Evidence Modules And R Parity

时间建议：持续推进。

交付：

- bulk evidence module validation。
- spatial evidence module validation。
- selected R parity matrix。
- method-specific validation notebooks。

验收：

- 每个接入方法都有 parity/validation 说明。
- 没有 unvalidated method 被推荐为默认结论。

### Milestone E：Scale And Productization

时间建议：长期。

交付：

- memory/runtime benchmarks。
- large dataset guide。
- project templates。
- polished reports。

验收：

- 真实大数据集可控运行。
- 项目交付物结构稳定。

## 8. 一句话战略

scLucid 的下一阶段不应追求“功能最多”，而应追求：

> 在肿瘤单细胞研究中，把最常见、最容易误判、最需要复查的分析环节，做成诊断优先、证据可追溯、推断边界清楚、工程上可靠的 lucid workflow。
