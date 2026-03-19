# Extensibility and Plugin Architecture Specification

## Purpose

Define a plugin-based architecture that allows developers to extend scLucid's functionality without modifying core code, using abstract base classes and factory patterns.

## Background

scLucid has grown significantly and needs to support:
1. Custom analysis algorithms developed by users
2. Third-party integrations (new annotation methods, scoring methods, etc.)
3. Experimental features without destabilizing core code
4. Community contributions without requiring core maintainer review for each change

## Requirements

### Requirement: Abstract Base Classes for Analysis Steps

The system SHALL provide abstract base classes that define interfaces for all major analysis components.

#### Scenario: Create custom QC filter
- **GIVEN** `AnalysisStep` abstract base class
- **WHEN** developer creates class `CustomQC(AnalysisStep)`
- **THEN** the class MUST implement:
  - `validate_input(adata) -> bool`
  - `run(adata, **kwargs) -> AnnData`
  - `get_summary() -> dict`

#### Scenario: Create custom cell annotator
- **GIVEN** `CellAnnotator` abstract base class
- **WHEN** developer creates class `MyAnnotator(CellAnnotator)`
- **THEN** the class MUST implement:
  - `predict(adata, reference, **kwargs) -> AnnData`
  - `get_confidence(adata) -> pd.Series`

#### Scenario: Attempt to use incomplete implementation
- **GIVEN** class `IncompleteStep(AnalysisStep)` missing required methods
- **WHEN** user tries to instantiate or run it
- **THEN** system raises `NotImplementedError` or `TypeError`

### Requirement: Plugin Registration Factory

The system SHALL provide a factory pattern for dynamic plugin registration and instantiation.

#### Scenario: Register custom plugin
- **GIVEN** custom analysis class `MyCustomAnalysis(AnalysisStep)`
- **WHEN** calling `AnalysisStepFactory.register('my_plugin', MyCustomAnalysis)`
- **THEN** system:
  - Validates class is subclass of AnalysisStep
  - Stores plugin in registry
  - Makes it available for instantiation

#### Scenario: Create plugin instance
- **GIVEN** registered plugin 'my_plugin'
- **WHEN** calling `AnalysisStepFactory.create('my_plugin', config=config)`
- **THEN** system:
  - Retrieves plugin class from registry
  - Instantiates with provided parameters
  - Returns plugin instance

#### Scenario: List available plugins
- **GIVEN** multiple registered plugins
- **WHEN** calling `AnalysisStepFactory.list_steps()`
- **THEN** system returns list of all registered plugin names

#### Scenario: Register non-compliant plugin
- **GIVEN** class `NotAnAnalysisStep` that doesn't inherit from AnalysisStep
- **WHEN** calling `AnalysisStepFactory.register('invalid', NotAnAnalysisStep)`
- **THEN** system raises `TypeError` (must be subclass of AnalysisStep)

#### Scenario: Create unregistered plugin
- **GIVEN** factory with plugins 'plugin_a' and 'plugin_b'
- **WHEN** calling `AnalysisStepFactory.create('plugin_c')`
- **THEN** system raises `KeyError` (unknown plugin name)

### Requirement: Plugin Configuration Integration

The system SHALL allow plugins to use Pydantic-based configuration like core components.

#### Scenario: Plugin with custom config
- **GIVEN** plugin class `MyPlugin(AnalysisStep)` and `MyPluginConfig(SclucidBaseConfig)`
- **WHEN** instantiating plugin with `MyPlugin(config=MyPluginConfig(param=value))`
- **THEN** system:
  - Validates config parameters
  - Provides type hints and field descriptions
  - Allows JSON serialization/deserialization

#### Scenario: Config validation in plugins
- **GIVEN** `MyPluginConfig(min_value: int = Field(ge=0, le=100))`
- **WHEN** user creates `MyPluginConfig(min_value=-1)`
- **THEN** system raises `ValidationError` (min_value must be between 0 and 100)

### Requirement: Plugin Discovery and Loading

The system SHALL support automatic plugin discovery from specified directories.

#### Scenario: Load plugins from directory
- **GIVEN** directory `/path/to/plugins/` containing Python files with plugin classes
- **WHEN** calling `load_plugins_from_directory('/path/to/plugins/')`
- **THEN** system:
  - Scans directory for `.py` files
  - Imports each module
  - Auto-registers classes inheriting from AnalysisStep
  - Returns number of plugins loaded

#### Scenario: Prevent plugin name conflicts
- **GIVEN** plugin 'my_filter' already registered
- **WHEN** attempting to register another plugin with same name
- **THEN** system raises `ValueError` (plugin name already exists) or overwrites with warning

### Requirement: Backward Compatibility

The system SHALL maintain backward compatibility with existing code that doesn't use plugins.

#### Scenario: Use core functions without plugins
- **GIVEN** user doesn't register any custom plugins
- **WHEN** using standard analysis functions (run_standard_qc, etc.)
- **THEN** all functions work exactly as before

#### Scenario: Gradual migration to plugins
- **GIVEN** existing code using direct function calls
- **WHEN** gradually refactoring to use plugin interface
- **THEN** both approaches can coexist without breaking changes

## Architecture

### Abstract Base Classes Hierarchy

```
scLucid.base_interfaces
├── AnalysisStep (ABC)
│   └── Base for all analysis steps
│
├── QCFilter (ABC)
│   └── Base for QC filtering operations
│
├── CellAnnotator (ABC)
│   └── Base for annotation methods
│
├── ScoringMethod (ABC)
│   └── Base for scoring algorithms
│
├── PlottingBackend (ABC)
│   └── Base for visualization backends
│
└── ProportionMethod (ABC)
    └── Base for proportion analysis methods
```

### Factory Pattern

```
AnalysisStepFactory
├── register(name, class)
├── create(name, **kwargs)
└── list_steps()
```

### Plugin Examples

1. **Custom QC Filter**: `HighStringencyQC` - Stricter filtering thresholds
2. **Custom Annotator**: `DatabaseAnnotator` - Annotation from external database
3. **Custom Scorer**: `PathwayScorer` - Score by pathway activity
4. **Custom Plotter**: `InteractivePlotter` - Interactive visualizations

## Benefits

### For Developers
- **No core modification needed**: Extend functionality without touching scLucid core
- **Sandboxed development**: Test plugins independently
- **Clear interface**: Abstract base classes define contract
- **Type safety**: Pydantic configs provide validation

### For Users
- **Custom workflows**: Mix core and custom components
- **Experimental features**: Try new algorithms safely
- **Community sharing**: Contribute and use plugins from others

### For Maintainers
- **Stable core**: Core code changes less frequently
- **Easier review**: Plugin PRs are isolated
- **Flexibility**: Accept/reject plugins without affecting core

## Non-Functional Requirements

### Performance
- Plugin instantiation overhead: < 10ms
- Factory lookup: < 1ms
- Plugin execution: Comparable to core functions

### Security
- Plugins run in same Python process (no sandboxing)
- No resource limits (users must trust plugins they install)

### Documentation
- All abstract base classes MUST have docstrings
- Plugin examples provided in `examples/plugin_development_example.py`
- Naming conventions documented in `docs/NAMING_CONVENTIONS.md`

## Migration Path

### Phase 1: Abstract Base Classes (Current)
- ✅ Create abstract base classes in `scLucid.base_interfaces`
- ✅ Export in top-level `__init__.py`
- ✅ Provide examples and documentation

### Phase 2: Factory Pattern (Current)
- ✅ Implement `AnalysisStepFactory`
- ✅ Support dynamic registration
- ✅ Provide examples

### Phase 3: Plugin Discovery (Future)
- ⏳ Implement automatic plugin loading
- ⏳ Create plugin registry (optional)
- ⏳ Support plugin versioning

### Phase 4: Core Refactoring (Future)
- ⏳ Gradually refactor core code to use abstract base classes
- ⏳ Maintain backward compatibility
- ⏳ Deprecate old interfaces gracefully

## Success Metrics

1. **Usability**: Developer can create and register plugin in < 50 lines of code
2. **Safety**: Invalid plugins fail fast with clear error messages
3. **Performance**: No performance degradation vs core functions
4. **Adoption**: At least 3 community plugins created within 6 months
5. **Quality**: All plugins pass type checking and linters

## Examples

See `examples/plugin_development_example.py` for complete working examples of:
- Creating custom QC filter
- Creating custom annotator
- Registering and using plugins
- Integrating plugins into workflows
