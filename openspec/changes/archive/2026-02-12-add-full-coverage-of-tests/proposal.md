# Change: Add comprehensive test coverage across scLucid modules

## Why
The current test surface does not fully cover core public functions across `qc`, `preprocess`, `analysis`, `tools`, `utils`, and global configuration behavior. This increases regression risk when module structures, imports, and optional dependency boundaries change.

A recent set of import regressions showed that package-level breakage can reach users before feature-level tests fail. We need a structured, maintainable test module strategy that validates both API behavior and package import integrity.

## What Changes
- Add a consolidated testing capability and test architecture for full-package coverage of public functions.
- Expand unit tests for config validation, helpers, and deterministic utility behavior.
- Add integration tests for core workflows using lightweight synthetic `AnnData` fixtures.
- Add and enforce import-smoke tests for top-level and submodule public APIs.
- Add optional-dependency boundary tests to ensure graceful degradation when extras are absent.
- Standardize test layout, markers, fixtures, and CI execution gates.

## Impact
- Affected specs: testing, development-model
- Affected code: `tests/`, `pytest.ini`, selected modules under `src/scLucid/` where testability seams are needed
- Risk level: medium
