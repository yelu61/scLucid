# Change: Migrate to Pydantic Configuration System

## Why

The project previously maintained a dual configuration system:
- **Old system**: `@dataclass` configs requiring manual validation
- **New system**: Pydantic configs with automatic validation

This dual state caused:
1. Increased maintenance burden
2. Inconsistent user API
3. Underutilized Pydantic features (validation, serialization, schemas)

## What Changes

- **BREAKING**: Remove all dataclass-based configuration classes
- Replace with Pydantic-based configs inheriting from `SclucidBaseConfig`
- Remove manual `.validate()` calls (Pydantic auto-validates on instantiation)
- Replace `asdict()` with `config.to_dict()`
- Update all imports across QC, Preprocess, Analysis, Tools, and Global modules

## Impact

### Affected specs
- `qc` - QC workflow, thresholds, doublet detection configs
- `preprocess` - Normalization, HVG, scaling, integration configs
- `analysis` - Clustering, annotation, DE, enrichment configs
- `tools` - BayesPrism deconvolution configs
- `global` - Runtime configuration (logging, plotting, cache)

### Affected code
- `src/scLucid/qc/config.py` - Replaced with Pydantic version
- `src/scLucid/preprocess/config.py` - Replaced with Pydantic version
- `src/scLucid/analysis/config.py` - Replaced with Pydantic version
- `src/scLucid/tools/pyBayesPrism/config.py` - Rewritten with Pydantic
- `src/scLucid/config.py` - Migrated to Pydantic
- 15+ files updated to remove `.validate()` and `asdict()` calls
- Test files updated to use Pydantic validation error types

### Migration notes
Users upgrading from previous version will need to:
- Ensure `pydantic` is installed
- Update any code calling `config.validate()` (no longer needed)
- Update any code using `asdict(config)` to use `config.to_dict()`

### Benefits
- ✅ Automatic type validation and conversion
- ✅ Built-in JSON serialization/deserialization
- ✅ Clearer error messages
- ✅ Auto-generated JSON Schema for documentation
- ✅ Unified configuration API across all modules
