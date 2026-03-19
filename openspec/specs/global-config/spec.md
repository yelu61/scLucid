# Global Configuration Specification

## Purpose

Provide runtime configuration system for scLucid, including logging, plotting settings, computational parameters, cache management, and resource paths.

## Requirements

### Requirement: Global Configuration Instance
The system SHALL provide a singleton global configuration instance.

#### Scenario: Get global config
- **WHEN** calling `get_config()`
- **THEN** returns singleton `_config` instance of `GlobalConfig`

#### Scenario: Thread-safe access
- **WHEN** multiple threads access `get_config()`
- **THEN** all threads receive the same instance safely

### Requirement: Configuration Updates
The system SHALL allow updating global configuration parameters.

#### Scenario: Set configuration
- **GIVEN** current global config
- **WHEN** calling `set_config(n_jobs=4, verbosity=2)`
- **THEN** system:
  - Updates n_jobs to 4
  - Updates verbosity to 2
  - Re-validates configuration
  - Reconfigures logging if verbosity changed

#### Scenario: Invalid n_jobs warning
- **GIVEN** `set_config(n_jobs=0)` (invalid value)
- **WHEN** setting configuration
- **THEN** system:
  - Issues warning (n_jobs=0 is invalid)
  - Does not modify value (user must set correctly)

#### Scenario: Reset to defaults
- **WHEN** calling `reset_config()`
- **THEN** system resets all fields to default values

### Requirement: Temporary Configuration Context
The system SHALL support temporary configuration changes via context manager.

#### Scenario: Use config context
- **GIVEN** global config with `n_jobs=-1`
- **WHEN** using:
  ```python
  with config_context(n_jobs=4, verbosity=2):
      # Some operation
  ```
- **THEN** system:
  - Temporarily sets n_jobs=4, verbosity=2 within context
  - Restores original values after context exits
  - Ensures thread-safety with lock

#### Scenario: Exception in context
- **GIVEN** config context with exception raised inside
- **WHEN** exception occurs
- **THEN** system still restores original configuration

### Requirement: Logging Configuration
The system SHALL automatically configure logging based on verbosity.

#### Scenario: Configure logging on init
- **GIVEN** `GlobalConfig(verbosity=1, log_file="analysis.log")`
- **WHEN** creating config
- **THEN** system:
  - Sets up root logger at INFO level
  - Creates console handler with formatter
  - Creates file handler if log_file specified
  - Log format: `[timestamp] LEVEL - name - message`

#### Scenario: Verbosity levels
- **GIVEN** different verbosity values
- **WHEN** creating config:
  - `verbosity=0` → WARNING level
  - `verbosity=1` → INFO level
  - `verbosity=2` → DEBUG level
- **THEN** logging level matches verbosity

#### Scenario: Update logging
- **GIVEN** existing config with `verbosity=1`
- **WHEN** calling `set_config(verbosity=2)`
- **THEN** system reconfigures logging at DEBUG level

#### Scenario: Log file directory creation
- **GIVEN** `log_file="/new/path/logs/analysis.log"`
- **WHEN** creating config
- **THEN** system creates parent directories if needed

### Requirement: Path Type Conversion
The system SHALL automatically convert string paths to Path objects.

#### Scenario: Convert log_file path
- **GIVEN** `GlobalConfig(log_file="/path/to/log.txt")`
- **WHEN** creating config
- **THEN** string is converted to `Path("/path/to/log.txt")`

#### Scenario: Convert cache_dir path
- **GIVEN** `GlobalConfig(cache_dir="~/.cache/sclucid")`
- **WHEN** creating config
- **THEN** string is converted to `Path` object

#### Scenario: Convert marker_db_path
- **GIVEN** `GlobalConfig(marker_db_path="/data/markers.toml")`
- **WHEN** creating config
- **THEN** path is converted to `Path` object

### Requirement: Computational Settings
The system SHALL provide configurable computational parameters.

#### Scenario: Configure parallel jobs
- **GIVEN** `GlobalConfig(n_jobs=4, backend="loky")`
- **WHEN** using parallel operations
- **THEN** system:
  - Uses 4 CPU cores
  - Uses loky backend for joblib

#### Scenario: Use all cores
- **GIVEN** `GlobalConfig(n_jobs=-1)`
- **WHEN** running parallel operations
- **THEN** system uses all available CPU cores

#### Scenario: Configure random state
- **GIVEN** `GlobalConfig(random_state=42)`
- **WHEN** running stochastic operations
- **THEN** system uses seed 42 for reproducibility

### Requirement: Plotting Settings
The system SHALL provide global plotting configuration.

#### Scenario: Configure figure parameters
- **GIVEN** `GlobalConfig(figure_dpi=300, figure_format="svg", color_palette="tab20")`
- **WHEN** calling `set_figure_params()` (imported from settings)
- **THEN** system:
  - Sets matplotlib DPI to 300
  - Sets default figure format to SVG
  - Uses tab20 color palette

#### Scenario: Auto-update figure params on config change
- **GIVEN** existing config
- **WHEN** calling `set_config(figure_dpi=150, plot_theme="dark")`
- **THEN** system:
  - Updates config values
  - Calls `set_figure_params()` to apply changes

#### Scenario: Academic font style
- **GIVEN** `GlobalConfig(font_style="nature")`
- **WHEN** creating figures
- **THEN** system uses Nature journal font style

### Requirement: Cache Configuration
The system SHALL provide cache management settings.

#### Scenario: Enable caching
- **GIVEN** `GlobalConfig(use_cache=True, cache_dir="~/.sclucid/cache")`
- **WHEN** running cacheable operations
- **THEN** system caches results in specified directory

#### Scenario: Disable caching
- **GIVEN** `GlobalConfig(use_cache=False)`
- **WHEN** running operations
- **THEN** system does not use caching

### Requirement: Memory Management
The system SHALL provide memory-efficient processing options.

#### Scenario: Low memory mode
- **GIVEN** `GlobalConfig(low_memory_mode=True)`
- **WHEN** processing large datasets
- **THEN** system:
  - Processes data in chunks
  - Uses `chunk_size` parameter

#### Scenario: Configure chunk size
- **GIVEN** `GlobalConfig(chunk_size=1000)`
- **WHEN** chunking operations
- **THEN** system processes 1000 items per chunk

### Requirement: Species-Specific Settings
The system SHALL provide default species configuration.

#### Scenario: Set default species
- **GIVEN** `GlobalConfig(default_species="human")`
- **WHEN** loading species-specific resources (e.g., cell cycle genes)
- **THEN** system uses human as default species

### Requirement: Resource Paths
The system SHALL provide configurable paths to resource files.

#### Scenario: Configure marker database path
- **GIVEN** `GlobalConfig(marker_db_path="/custom/path/markers.toml")`
- **WHEN** loading marker database
- **THEN** system loads from custom path instead of default

#### Scenario: Configure gene set path
- **GIVEN** `GlobalConfig(gene_set_path="/data/genesets/")`
- **WHEN** loading gene sets
- **THEN** system uses custom gene set directory

### Requirement: Configuration Validation
The system SHALL validate configuration parameters.

#### Scenario: Validate n_jobs
- **GIVEN** `n_jobs=-2` (invalid)
- **WHEN** creating config
- **THEN** system issues warning (must be -1 or >= 1)

#### Scenario: Validate verbosity
- **GIVEN** `verbosity=5` (invalid)
- **WHEN** creating config
- **THEN** system issues warning (must be 0, 1, or 2)

### Requirement: Configuration Serialization
GlobalConfig MUST support serialization.

#### Scenario: Serialize to dict
- **GIVEN** a `GlobalConfig` instance
- **WHEN** calling `config.to_dict()`
- **THEN** returns dictionary with all configuration fields

#### Scenario: Serialize to JSON
- **GIVEN** a `GlobalConfig` instance
- **WHEN** calling `config.to_json()`
- **THEN** returns JSON string representation

#### Scenario: Save to file
- **GIVEN** a `GlobalConfig` instance
- **WHEN** calling `config.to_json_file("global_config.json")`
- **THEN** saves configuration to JSON file

#### Scenario: Load from file
- **GIVEN** a saved config file
- **WHEN** calling `GlobalConfig.from_json_file("global_config.json")`
- **THEN** loads and validates configuration

### Requirement: Extra Fields Handling
The system SHALL ignore extra fields (Pydantic `extra="ignore"`).

#### Scenario: Extra fields are ignored
- **GIVEN** `GlobalConfig(invalid_field="value")`
- **WHEN** creating config
- **THEN** system:
  - Ignores unknown field
  - Does not raise error
  - Continues with valid fields only
