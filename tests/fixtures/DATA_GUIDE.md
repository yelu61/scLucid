# scLucid 测试数据规范

## 数据来源总览

| 类型 | 名称 | 规模 | 用途 | 标记 |
|------|------|------|------|------|
| 合成 | `tiny_adata` | 8 cells x 6 genes | smoke/导入测试 | pytest fixture |
| 合成 | `minimal_adata` | 100 cells x 200 genes | 通用单元测试 | pytest fixture |
| 合成 | `qc_test_adata` | 500 cells, 含异常值 | QC过滤测试 | pytest fixture |
| 合成 | `integration_test_adata` | 1000 cells, 3 batches | 批次校正测试 | pytest fixture |
| 合成 | `doublet_test_adata` | 500 cells, 5%双细胞 | 双细胞检测测试 | pytest fixture |
| 真实 | `pbmc3k.h5ad` | 2700 cells, PBMC, 4 samples | 非肿瘤pipeline验证 | 本地文件 |
| 真实 | `schlesinger2020.pdac` | 6499 cells, 单样本胰腺癌 | 单样本肿瘤验证 | 本地文件 |
| 真实 | `lin2020.pdac` | 9621 cells, 10样本胰腺癌 | 多样本批次校正验证 | 本地文件 |

## 数据分配矩阵

```
                    tiny  minimal  qc_test  integ_test  doublet  pbmc3k  pdac_s  pdac_m
smoke tests          ●
config tests                 ●
QC metrics                   ●        ●                                    ●       ●       ●
QC filtering                          ●                                     ●       ●       ●
QC workflow                  ●        ●                                    ●       ●       ●
Preprocess normalize         ●        ●                                    ●       ●       ●
Preprocess HVG               ●        ●                                    ●       ●       ●
Preprocess integrate                                   ●                     ●               ●
Preprocess workflow          ●        ●                                    ●       ●       ●
Analysis clustering          ●        ●                                    ●       ●       ●
Analysis DE                  ●        ●                                    ●       ●       ●
Analysis workflow            ●        ●                                    ●       ●       ●
Tools (pyCellChat etc)       ●
Tumor module                                                    (todo)     (todo)  (todo)  (todo)
Integration (synthetic)      ●        ●
Integration (real)                                                              ●       ●       ●
```

## Fixture 定义位置

| Fixture | 定义位置 | 作用域 |
|---------|----------|--------|
| `synthetic_generator` | `tests/fixtures/synthetic_data.py` | function |
| `minimal_adata` | `tests/fixtures/synthetic_data.py` | function |
| `qc_test_adata` | `tests/fixtures/synthetic_data.py` | function |
| `integration_test_adata` | `tests/fixtures/synthetic_data.py` | function |
| `doublet_test_adata` | `tests/fixtures/synthetic_data.py` | function |
| `tiny_adata` | `tests/conftest.py` | function |
| `temp_output_dir` | `tests/conftest.py` | function |
| `test_data_dir` | `tests/conftest.py` | session |

## 真实数据加载规范

Integration 测试通过 `_load_subset()` 加载真实数据并子采样：
- 目的：控制运行时间（完整pipeline在数千细胞上需要数分钟）
- 默认子采样：400-500 cells
- 必须确保 `counts` layer 存在（PDAC数据的raw counts在 `.X`，需要复制到 `layers['counts']`）

## 运行命令速查

```bash
# 仅合成数据快速测试 (~30秒)
pytest -m "not slow and not integration"

# 含真实数据的集成测试 (~1分钟)
pytest tests/integration -m "slow"

# 全部测试 (~2分钟)
pytest

# 特定模块
pytest tests/qc
pytest tests/preprocess
pytest tests/analysis
```
