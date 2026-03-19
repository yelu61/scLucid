## 1. Core Module Migration

### 1.1 QC Module
- [x] 1.1.1 Backup old `src/scLucid/qc/config.py`
- [x] 1.1.2 Rename `config_v2.py` to `config.py`
- [x] 1.1.3 Remove `.validate()` calls in `workflow.py` and `workflow_v2.py`
- [x] 1.1.4 Replace `asdict(config)` with `config.to_dict()` in:
  - [x] `filtering.py` (2 locations)
  - [x] `doublet.py` (1 location)
  - [x] `metrics.py` (1 location)
- [x] 1.1.5 Remove `from dataclasses import asdict` imports
- [x] 1.1.6 Update `__init__.py` exports (already correct)

### 1.2 Preprocess Module
- [x] 1.2.1 Backup old `src/scLucid/preprocess/config.py`
- [x] 1.2.2 Rename `config_v2.py` to `config.py`
- [x] 1.2.3 Remove `.validate()` calls in workflow files
- [x] 1.2.4 Replace `asdict(config)` with `config.to_dict()` in:
  - [x] `workflow.py` (1 location)
  - [x] `workflow_v2.py` (1 location)
  - [x] `hvg.py` (1 location)
  - [x] `scale.py` (1 location)
  - [x] `normalize.py` (1 location)
  - [x] `integrate.py` (1 location)
- [x] 1.2.5 Update `dataclasses.replace()` to use `from dataclasses import replace`
- [x] 1.2.6 Update `__init__.py` exports

### 1.3 Analysis Module
- [x] 1.3.1 Backup old `src/scLucid/analysis/config.py`
- [x] 1.3.2 Rename `config_v2.py` to `config.py`
- [x] 1.3.3 Fix `de_enrichment.py` import:
  - [x] Change `from .config import BaseConfig` to `from ..base_config import SclucidBaseConfig`
  - [x] Update type hint from `BaseConfig` to `SclucidBaseConfig`
- [x] 1.3.4 Update `__init__.py` exports (already correct)

## 2. Tools Module Migration

### 2.1 BayesPrism Rewrite
- [x] 2.1.1 Rewrite `src/scLucid/tools/pyBayesPrism/config.py`:
  - [x] Convert `@dataclass` to Pydantic `SclucidBaseConfig` inheritance
  - [x] Add `Field()` with validation and descriptions
  - [x] Implement `@model_validator(mode="after")` for parameter validation
  - [x] Preserve all existing fields and behavior
- [x] 2.1.2 Remove `self.config.validate()` call in `core.py`
- [x] 2.1.3 Update `__init__.py` exports (already correct)

## 3. Global Configuration Migration

### 3.1 GlobalConfig Rewrite
- [x] 3.1.1 Rewrite `src/scLucid/config.py`:
  - [x] Convert `@dataclass` to Pydantic `SclucidBaseConfig` inheritance
  - [x] Move `__post_init__` logic to `@model_validator(mode="after")`
  - [x] Preserve singleton pattern with threading lock
  - [x] Keep all custom methods (`set()`, `reset()`, `_setup_logging()`)
  - [x] Preserve context manager (`config_context()`)

## 4. Testing

### 4.1 Test Updates
- [x] 4.1.1 Rename `tests/test_config_v2.py` to `tests/test_configs.py`
- [x] 4.1.2 Update all imports from `.config_v2` to `.config`
- [x] 4.1.3 Fix `test_cache_config_validation` to expect Pydantic `ValidationError`
- [x] 4.1.4 Run tests: `pytest tests/test_configs.py -v` (39 passed)
- [x] 4.1.5 Run QC tests: `pytest tests/qc/ -v` (51 passed)

## 5. Documentation

### 5.1 Update User Docs
- [x] 5.1.1 Update `CLAUDE.md`:
  - [x] Replace dataclass references with Pydantic
  - [x] Add "Configuration System (Pydantic)" section
  - [x] Document Pydantic features and examples
- [ ] 5.1.2 Update `README.md` with breaking changes notice
- [ ] 5.1.3 Create migration guide for users (optional)

### 5.2 Cleanup
- [ ] 5.2.1 Delete backup files after verification:
  - [ ] `src/scLucid/qc/config.py.old`
  - [ ] `src/scLucid/preprocess/config.py.old`
  - [ ] `src/scLucid/analysis/config.py.old`
- [ ] 5.2.2 Update CHANGELOG.md with breaking changes
