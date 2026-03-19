## MODIFIED Requirements

### Requirement: Configuration Validation
QC configuration classes MUST automatically validate on instantiation using Pydantic.

**Before**: Dataclass configs required manual `.validate()` calls
**After**: Pydantic configs validate automatically

#### Scenario: Invalid threshold raises error
- **WHEN** creating `QCThresholds(min_genes=-1)`
- **THEN** system raises `pydantic.ValidationError` immediately

#### Scenario: Valid config succeeds
- **WHEN** creating `QCThresholds(min_genes=200, pc_mt=20.0)`
- **THEN** config object is created successfully

### Requirement: Configuration Serialization
QC configuration MUST support serialization to dictionary and JSON formats.

**Before**: Used `asdict(config)` from dataclasses module
**After**: Use `config.to_dict()` method (Pydantic's `model_dump()`)

#### Scenario: Serialize to dict
- **WHEN** calling `config.to_dict()`
- **THEN** returns Python dictionary with all config fields

#### Scenario: Serialize to JSON
- **WHEN** calling `config.to_json()`
- **THEN** returns JSON string representation

#### Scenario: Save to file
- **WHEN** calling `config.to_json_file("path/to/config.json")`
- **THEN** saves config to specified file path

### Requirement: Configuration Classes
All QC configs MUST inherit from `SclucidBaseConfig` and use Pydantic `Field()` for validation.

#### Scenario: Create workflow config
- **WHEN** instantiating `QCWorkflowConfig()`
- **THEN** config has all default values set and validated

#### Scenario: Create thresholds config
- **WHEN** instantiating `QCThresholds(min_genes=200)`
- **THEN** config validates that `min_genes >= 0`

### Requirement: Field Descriptions
All configuration fields MUST include human-readable descriptions via Pydantic `Field()`.

#### Scenario: Inspect field description
- **WHEN** accessing `help(QCThresholds.pc_mt)`
- **THEN** displays "Mitochondrial percentage threshold"
