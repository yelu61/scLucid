## ADDED Requirements

### Requirement: Pydantic-Based BayesPrism Configuration
BayesPrism tools MUST use Pydantic configuration classes with automatic validation.

#### Scenario: Create Prism config
- **WHEN** instantiating `PrismConfig(n_iter=200, n_chains=4)`
- **THEN** config validates that `n_iter > 0` and `n_chains > 0`

#### Scenario: Invalid burnin parameter
- **WHEN** creating `PrismConfig(n_iter=100, burnin=150)`
- **THEN** raises `ValidationError` (burnin must be < n_iter)

### Requirement: Gibbs Control Initialization
Gibbs sampler control parameters MUST be initialized with defaults if not provided.

#### Scenario: Default gibbs_control
- **WHEN** creating `PrismConfig()` without `gibbs_control`
- **THEN** system auto-initializes with `{'chain_length': n_iter, 'burn_in': burnin, 'thinning': 1, 'verbose': False}`

#### Scenario: Custom gibbs_control
- **WHEN** creating `PrismConfig(gibbs_control={'thinning': 5})`
- **THEN** custom thinning value is used, other defaults still apply

### Requirement: Field Validation with Pydantic
All BayesPrism config fields MUST use Pydantic `Field()` for validation and descriptions.

#### Scenario: Field with constraints
- **WHEN** defining `n_iter: int = Field(default=100, gt=0, description="...")`
- **THEN** Pydantic validates `n_iter > 0` on instantiation

#### Scenario: Optional fields
- **WHEN** defining `key: Optional[str] = Field(default=None, description="...")`
- **THEN** config accepts `None` or string values

### Requirement: Input Type Validation
Reference config MUST validate `input_type` against allowed values.

#### Scenario: Valid input type
- **WHEN** creating `ReferenceConfig(input_type="count.matrix")`
- **THEN** config is created successfully

#### Scenario: Invalid input type
- **WHEN** creating `ReferenceConfig(input_type="invalid_type")`
- **THEN** raises `ValidationError` with allowed values: ["count.matrix", "GEP"]
