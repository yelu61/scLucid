## MODIFIED Requirements

### Requirement: Preprocessing Configuration Validation
Preprocessing configuration classes MUST automatically validate on instantiation using Pydantic.

**Before**: Dataclass configs required manual `.validate()` calls
**After**: Pydantic configs validate automatically

#### Scenario: Invalid target_sum raises error
- **WHEN** creating `NormalizationConfig(target_sum=-1)`
- **THEN** system raises `pydantic.ValidationError` with `target_sum` field error

#### Scenario: Valid normalization config
- **WHEN** creating `NormalizationConfig(target_sum=1e4)`
- **THEN** config is created successfully

### Requirement: Configuration Serialization
Preprocessing configs MUST support `to_dict()` method for serialization.

**Before**: Used `asdict(config)` from dataclasses
**After**: Use `config.to_dict()` (Pydantic's `model_dump()`)

#### Scenario: Serialize workflow config
- **WHEN** workflow stores config to `adata.uns`
- **THEN** calls `config.to_dict()` to get dictionary representation

#### Scenario: Store config in AnnData
- **WHEN** preprocessing workflow completes
- **THEN** config is serialized via `config.to_dict()` and stored in `adata.uns["sclucid"]["preprocess"]["workflow_config"]`

### Requirement: Configuration Replacement
Existing `dataclasses.replace()` calls MUST be updated to use `from dataclasses import replace`.

#### Scenario: Clone config with changes
- **WHEN** code calls `replace(config, field=new_value)`
- **THEN** returns new config instance with specified field changed

### Requirement: Config Classes Inherit Base
All preprocessing configs MUST inherit from `SclucidBaseConfig` or `WorkflowConfigBase`.

#### Scenario: Workflow config base
- **WHEN** inspecting `WorkflowConfig` class
- **THEN** it inherits from `WorkflowConfigBase` with `n_jobs`, `random_state` fields

#### Scenario: Sub-config base
- **WHEN** inspecting `NormalizationConfig` class
- **THEN** it inherits from `SclucidBaseConfig` with `save_dir`, `verbose` fields
