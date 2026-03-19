## 1. Specification
- [ ] 1.1 Confirm affected capability specs (`openspec/specs/*/spec.md`)
- [ ] 1.2 Add/modify delta specs under `openspec/changes/<change-id>/specs/`
- [ ] 1.3 Validate proposal: `openspec validate <change-id> --strict`

## 2. Implementation
- [ ] 2.1 Implement code changes in `src/scLucid/...`
- [ ] 2.2 Add/update tests in `tests/...`
- [ ] 2.3 Add import-smoke coverage for affected public modules when module layout changes

## 3. Documentation
- [ ] 3.1 Update user docs in `docs/source/...` for behavior/API changes
- [ ] 3.2 Update examples in `examples/...` if workflows changed

## 4. Verification
- [ ] 4.1 Run targeted pytest for changed modules
- [ ] 4.2 Run import smoke checks for `scLucid` and affected submodules
- [ ] 4.3 Re-run `openspec validate <change-id> --strict`
