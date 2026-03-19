## MODIFIED Requirements

### Requirement: Global Configuration Validation
GlobalConfig MUST automatically validate settings on instantiation using Pydantic.

**Before**: Validation happened in `__post_init__` with manual checks
**After**: Pydantic validates via `@model_validator(mode="after")`

#### Scenario: Invalid n_jobs value
- **WHEN** creating `GlobalConfig(n_jobs=0)` (invalid, must be -1 or >= 1)
- **THEN** validation warning is issued (warning, not error, for backwards compatibility)

#### Scenario: Invalid verbosity
- **WHEN** creating `GlobalConfig(verbosity=5)` (invalid, must be 0, 1, or 2)
- **THEN** warning is issued indicating valid values

### Requirement: Logging Setup
Logging configuration MUST be set up automatically after config validation.

**Before**: Called in `__post_init__`
**After**: Called in `@model_validator(mode="after")`'s `setup_and_validate()`

#### Scenario: Auto-logging setup
- **WHEN** creating `GlobalConfig(verbosity=1, log_file="app.log")`
- **THEN** logging is configured at INFO level with file handler

### Requirement: Path Type Conversion
String paths MUST be automatically converted to `Path` objects.

#### Scenario: Convert string log_file
- **WHEN** creating `GlobalConfig(log_file="/path/to/log.txt")`
- **THEN** string is converted to `Path("/path/to/log.txt")`

#### Scenario: Convert string cache_dir
- **WHEN** creating `GlobalConfig(cache_dir="~/.cache")`
- **THEN** string is converted to `Path` object

### Requirement: Singleton Pattern Preservation
GlobalConfig MUST preserve its singleton pattern with threading lock.

**Before**: Used with `@dataclass` + module-level `_config` instance
**After**: Same pattern, but with Pydantic `GlobalConfig` class

#### Scenario: Get global instance
- **WHEN** calling `get_config()`
- **THEN** returns the singleton `_config` instance

#### Scenario: Thread-safe updates
- **WHEN** multiple threads call `set_config(**kwargs)`
- **THEN** `_config_lock` ensures thread-safe updates

#### Scenario: Config context manager
- **WHEN** using `with config_context(n_jobs=4):`
- **THEN** config is temporarily changed within context, restored after

### Requirement: Custom Methods
GlobalConfig MUST preserve all custom methods from dataclass version.

#### Scenario: Set configuration
- **WHEN** calling `config.set(n_jobs=4, verbosity=2)`
- **THEN** updates fields and re-validates

#### Scenario: Reset configuration
- **WHEN** calling `config.reset()`
- **THEN** resets all fields to default values

#### Scenario: Config updates figure params
- **WHEN** calling `config.set(figure_dpi=300)`
- **THEN** also calls `set_figure_params()` to update matplotlib settings
