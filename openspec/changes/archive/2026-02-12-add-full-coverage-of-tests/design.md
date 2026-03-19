## Context
scLucid has broad module coverage and both core and optional dependencies. Existing tests are uneven across modules, and package-level import regressions can bypass function-level tests.

## Goals / Non-Goals
- Goals:
  - Ensure core public imports are continuously validated.
  - Increase function-level test coverage across all primary package modules.
  - Keep tests fast by default while preserving deeper integration coverage.
- Non-Goals:
  - Achieve exhaustive integration coverage for every optional backend in one change.
  - Redesign core workflows solely for test convenience.

## Decisions
- Decision: Split tests into `unit`, `integration`, and `optional` layers with shared fixtures.
- Rationale: Balances speed and reliability while making optional-backend behavior explicit.
- Alternatives considered:
  - Single monolithic test suite: rejected due to poor signal and long runtimes.
  - Only unit tests: rejected because workflow and import regressions would remain under-tested.

- Decision: Treat import-smoke checks as first-class regression gates.
- Rationale: scLucid package entry points are heavily used and prone to breakage when module layout changes.

## Risks / Trade-offs
- Broader tests increase CI runtime -> Mitigate with marker-based split and targeted execution in PRs.
- Optional-backend tests may be flaky across environments -> Mitigate with conditional markers and deterministic mocks where possible.

## Migration Plan
1. Define fixtures and marker conventions.
2. Add high-value tests for public APIs and core workflows.
3. Add optional dependency and import-smoke regressions.
4. Enforce baseline gates in CI and document local commands.
5. Iterate to close remaining module gaps.

## Open Questions
- Which exact coverage threshold should block merge (for total and per-module)?
- Which optional backends should run in default CI vs scheduled/nightly jobs?
