# Phase 5: R/Python Parity 与成熟方法精选移植

**Recommended duration**: 1-2 个月  
**Primary output**: parity matrices, method validation notebooks, supplementary methods  
**Main claim**: scLucid 不是盲目重写 R 生态，而是精选肿瘤单细胞高价值方法，并把它们纳入 Python workflow 和 evidence contract。

## 目标

打通 Python/R 语言壁垒，但保持克制：

- 只迁移或包装高价值方法。
- 每个方法必须有 parity validation。
- 每个输出必须进入 scLucid evidence/audit contract。
- 不能因为某个 R/Python optional dependency 缺失破坏核心 workflow。

## 方法优先级

### Priority 1: Doublet evidence

- scDblFinder / pyscdblfinder。
- Scrublet compatibility。
- lineage co-expression heuristic。
- external hashing/genotype/manual evidence。

目标：

- 与 R scDblFinder 在 synthetic 和 real data 上比较。
- 记录差异来源。
- 完成 fallback 和 raw-count guard。

### Priority 2: Annotation evidence

- SingleR-like reference annotation。
- Azimuth/Seurat label transfer-like evidence。
- CellTypist integration。

目标：

- 不把 reference label 当 ground truth。
- 作为 annotation evidence 与 marker evidence 合并。

### Priority 3: Malignancy/CNV evidence

- inferCNV-like output consumption。
- CopyKAT-like output consumption。
- lightweight CNV scoring。

目标：

- 先支持结果接入和 evidence storage。
- 再考虑 Python-native approximation。

### Priority 4: Program/pathway scoring

- AUCell/UCell/GSVA-like scoring。
- module score parity。

目标：

- 肿瘤 program 和 TME state scoring 稳定。
- 与成熟实现比较 rank consistency。

## 准备

每个方法建立 parity sheet：

| Field | Requirement |
|-------|-------------|
| method name | R/Python method |
| task | doublet, annotation, CNV, scoring |
| input contract | counts/log expression/metadata |
| output contract | obs/var/uns fields |
| reference implementation | R package/version |
| Python implementation | wrapper/port/native |
| validation data | synthetic and real |
| acceptance metric | AUC, agreement, correlation, rank consistency |
| fallback | missing dependency behavior |

## 具体步骤

### Step 1. 定义方法接入等级

每个方法标记为：

- `wrapped`: 调用外部成熟实现。
- `ported`: Python 复现核心算法。
- `approximated`: Python 近似，不声称完全等价。
- `consumed`: 读取外部结果作为 evidence。

README/API 文档必须明确等级，避免误导。

### Step 2. 输入输出契约

每个方法必须说明：

- 是否需要 raw counts。
- 是否支持 sparse matrix。
- 是否支持 multi-sample。
- 是否需要 sample/patient key。
- 是否需要 species/tissue/cancer context。
- 输出列名。
- 输出 inference/evidence level。

### Step 3. Synthetic validation

为每个方法准备可控 synthetic data：

- doublet：known synthetic doublets。
- annotation：known marker-positive populations。
- CNV：simulated chromosome arm shift。
- scoring：known spiked gene sets。

输出：

- expected vs observed。
- failure cases。
- runtime/memory。

### Step 4. Real-data parity validation

每个方法至少 2 个真实数据集：

- PBMC or normal baseline。
- tumor dataset。

比较：

- label agreement。
- score correlation。
- top-k overlap。
- AUC/ARI/NMI。
- discordant examples with interpretation。

### Step 5. Evidence contract integration

输出必须进入：

- `adata.obs` for per-cell labels/scores。
- `adata.uns["sclucid"][module]["method_evidence"]`。
- review summary。
- limitations。
- dependency/version metadata。

## 完成标准

Phase 5 通过条件：

- scDblFinder parity 完成并文档化。
- 至少 2 个其他高价值方法完成 validation matrix。
- 所有 optional dependency 缺失时有 graceful fallback。
- 没有 direct copy GPL 代码风险。
- 每个方法的 claims 与接入等级一致。

## 交付物

- `docs/R_PYTHON_PARITY_MATRIX.md` 或更新 `docs/source/r_parity.rst`。
- parity notebooks。
- synthetic tests。
- real-data validation reports。
- method evidence schemas。
- supplementary method paragraphs。

## 风险与处理

| 风险 | 表现 | 处理 |
|------|------|------|
| 盲目 Python 化 | 维护成本爆炸 | 只做 selected high-value methods |
| 与 R 结果不一致 | 审稿人质疑准确性 | 明确 wrapped/ported/approximated/consumed 等级 |
| GPL 污染 | 复制代码 | clean-room 实现或只做 wrapper/consumer |
| optional dependency 太重 | 用户安装困难 | extras + lazy import + fallback |

