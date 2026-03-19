## 1. Specification
- [x] 1.1 Define testing capability deltas under `openspec/changes/add-full-coverage-of-tests/specs/testing/spec.md`
- [x] 1.2 Align with existing requirements in `openspec/specs/development-model/spec.md`
- [x] 1.3 Validate proposal: `openspec validate add-full-coverage-of-tests --strict`

## 2. Test Architecture
- [x] 2.1 Standardize test directory structure by module (`tests/qc`, `tests/preprocess`, `tests/analysis`, `tests/tools`, `tests/utils`, `tests/config`)
- [x] 2.2 Add shared fixtures for minimal synthetic `AnnData` objects and reusable mock data
- [x] 2.3 Add/standardize pytest markers (`unit`, `integration`, `slow`, `optional`)

## 3. Coverage Expansion
- [x] 3.1 Add unit tests for configuration models and validation paths
- [x] 3.2 Add unit tests for pure helpers and storage/validation utilities
- [x] 3.3 Add integration tests for representative QC, preprocessing, and analysis workflows
- [x] 3.4 Add tests for optional dependency fallbacks and lazy import behavior
- [x] 3.5 Add import-smoke tests for `scLucid`, `scLucid.qc`, `scLucid.preprocess`, `scLucid.analysis`, `scLucid.tools`, `scLucid.config`

## 4. CI and Quality Gates
- [x] 4.1 Define baseline test commands for PR validation
- [x] 4.2 Ensure changed modules require targeted tests before merge
- [x] 4.3 Add syntax/import checks for modified Python files in CI

## 5. Verification
- [x] 5.1 Run targeted pytest suites for all touched modules
- [x] 5.2 Run import smoke suite in a minimal-dependency environment profile
- [x] 5.3 Re-run `openspec validate add-full-coverage-of-tests --strict`
