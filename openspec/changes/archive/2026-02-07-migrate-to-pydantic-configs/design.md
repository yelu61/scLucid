## Context

scLucid has historically used Python's `@dataclass` for configuration classes. While functional, this approach required:
- Manual validation logic (custom `.validate()` methods)
- Manual serialization (`asdict()` from dataclasses module)
- No built-in type coercion
- Limited error messages

The project had already created Pydantic v2 config files (`config_v2.py`) but they were not being used, creating a dual system.

### Constraints
- **Breaking change acceptable**: User approved direct replacement strategy
- Must maintain existing API where possible
- Preserve singleton pattern for GlobalConfig
- Keep all custom methods and behaviors

## Goals / Non-Goals

### Goals
- Unified Pydantic-based configuration system across all modules
- Automatic validation on config instantiation
- Built-in serialization methods (`to_dict()`, `to_json()`)
- Self-documenting configs with field descriptions
- Auto-generated JSON Schema for documentation

### Non-Goals
- Backward compatibility with dataclass configs (explicit breaking change)
- Adding new configuration fields (migration only)
- Performance optimization (Pydantic v2 already performant)

## Decisions

### Decision 1: Base Configuration Class
**Choice**: All configs inherit from `SclucidBaseConfig` or `WorkflowConfigBase`

**Rationale**:
- Provides common fields (`save_dir`, `verbose`, `plot`, `n_jobs`, `random_state`)
- Built-in serialization methods (`to_dict()`, `to_json()`, `from_json()`)
- Consistent API across all modules

**Alternatives considered**:
- Direct Pydantic `BaseModel` inheritance → Rejected: Loses common functionality
- Multiple base classes → Rejected: Unnecessary complexity

### Decision 2: Validation Strategy
**Choice**: Remove all `.validate()` calls, rely on Pydantic's auto-validation

**Rationale**:
- Pydantic validates on instantiation automatically
- `@field_validator()` and `@model_validator()` for custom logic
- Simpler user API (no manual validation needed)

**Code pattern**:
```python
# OLD (dataclass):
config = QCThresholds(min_genes=200)
config.validate()  # Manual

# NEW (Pydantic):
config = QCThresholds(min_genes=200)  # Auto-validated
# Raises ValidationError if invalid
```

### Decision 3: Serialization API
**Choice**: Use Pydantic's `.model_dump()` wrapped as `.to_dict()`

**Rationale**:
- `.to_dict()` is already defined in `SclucidBaseConfig`
- Consistent with existing v2 configs
- User-facing API remains familiar

**Implementation**:
```python
# OLD:
from dataclasses import asdict
config_dict = asdict(config)

# NEW:
config_dict = config.to_dict()  # Calls model_dump() internally
```

### Decision 4: Error Handling in Tests
**Choice**: Update tests to expect `pydantic.ValidationError`

**Rationale**:
- Pydantic raises detailed validation errors
- Different from `ValueError` used by old dataclass validators
- Tests need to match new error type and message format

**Example**:
```python
# OLD:
with pytest.raises(ValueError, match="must be non-negative"):
    config = CacheConfig(max_cache_age_days=-1)

# NEW:
from pydantic import ValidationError as PydanticValidationError
with pytest.raises(PydanticValidationError, match="max_cache_age_days"):
    config = CacheConfig(max_cache_age_days=-1)
```

## Risks / Trade-offs

### Risk 1: Breaking User Code
**Impact**: Users with custom config classes or calling `.validate()` will break

**Mitigation**:
- Clear breaking change documentation
- Migration guide in CLAUDE.md
- Comprehensive test coverage

### Risk 2: Pydantic Dependency
**Impact**: New required dependency

**Mitigation**:
- Pydantic is lightweight and widely adopted
- Already used in v2 configs (just not activated)
- Add to `requirements.txt` or `pyproject.toml`

### Risk 3: Learning Curve
**Impact**: Users unfamiliar with Pydantic may struggle

**Mitigation**:
- Extensive examples in CLAUDE.md
- Pydantic error messages are self-explanatory
- User API remains simple (`ConfigClass(field=value)`)

## Migration Plan

### Phase 1: Core Modules (Completed)
1. QC: Replace config_v2 → config, update imports
2. Preprocess: Replace config_v2 → config, update imports
3. Analysis: Replace config_v2 → config, fix BaseConfig import

### Phase 2: Tools Module (Completed)
1. BayesPrism: Rewrite config.py to use Pydantic
2. Remove manual validation in core.py

### Phase 3: Global Config (Completed)
1. Rewrite GlobalConfig with Pydantic
2. Preserve singleton pattern and threading
3. Keep all custom methods

### Phase 4: Testing (Completed)
1. Rename test_config_v2 → test_configs
2. Update imports
3. Fix validation error expectations
4. Run full test suite

### Phase 5: Documentation (In Progress)
1. Update CLAUDE.md ✅
2. Update README.md
3. Add migration guide (optional)
4. Delete backup files after verification

## Open Questions

None - migration is complete and tested.

### Verification Results
- ✅ All 39 config tests pass
- ✅ All 51 QC tests pass
- ✅ Integration tests pass
- ✅ Manual verification successful
