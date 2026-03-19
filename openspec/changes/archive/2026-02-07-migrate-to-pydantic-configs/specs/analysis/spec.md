## MODIFIED Requirements

### Requirement: Analysis Configuration Base Class
Analysis configs MUST import `SclucidBaseConfig` from base module, not local `BaseConfig`.

**Before**: `from .config import BaseConfig`
**After**: `from ..base_config import SclucidBaseConfig`

#### Scenario: Import base config
- **WHEN** analysis modules need base config class
- **THEN** they import from `..base_config` (parent directory)

#### Scenario: Type hints use SclucidBaseConfig
- **WHEN** function signatures accept config parameters
- **THEN** they use `SclucidBaseConfig` as type hint, not `BaseConfig`

### Requirement: Configuration Validation
Analysis configs MUST automatically validate using Pydantic.

#### Scenario: Invalid clustering resolution
- **WHEN** creating `ClusteringConfig(resolution=-1.0)`
- **THEN** raises `ValidationError` (resolution must be > 0)

#### Scenario: Invalid annotation confidence
- **WHEN** creating `AnnotationConfig(min_confidence=1.5)`
- **THEN** raises `ValidationError` (min_confidence must be <= 1.0)

### Requirement: Configuration Inheritance
All analysis config classes MUST inherit from `SclucidBaseConfig`.

#### Scenario: Clustering config inheritance
- **WHEN** instantiating `ClusteringConfig(method="leiden")`
- **THEN** config has base fields (`save_dir`, `verbose`) plus clustering-specific fields

#### Scenario: Enrichment config inheritance
- **WHEN** instantiating `EnrichmentConfig(method="gsea")`
- **THEN** config inherits from `SclucidBaseConfig` with validation enabled
