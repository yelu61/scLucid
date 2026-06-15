# scLucid Improvement Log

This document captures concrete problems, gaps, and enhancement ideas that emerge
from real-world use of scLucid in analysis projects. It is the bridge between
**project-driven discovery** and **scLucid source-code polishing**.

## How to Use This Log

1. **During analysis** (in your project directory):
   - Do not modify scLucid source code.
   - Record each issue using the template below.
   - Attach only minimal, structured context (not full AnnData objects or long
     notebooks).

2. **Sync to this log** (in scLucid directory):
   - Periodically move high-priority items into this file.
   - Group related issues into themes.

3. **Source-code polishing** (in scLucid directory):
   - Pick 1–3 related items per session.
   - Use Plan mode for behavior changes; direct edits for obvious bugs.
   - Add or update tests.
   - Commit with a reference to the log entry.

4. **Verify in project** (back in project directory):
   - Install the updated scLucid (`pip install -e /path/to/scLucid`).
   - Run the analysis cell or script that exposed the issue.
   - Update the log entry status.

## Relationship to Other Docs

- `API_TRIAGE.md` — tracks uncertain *public API symbols* (promote/demote/remove).
- `CORE_WORKFLOW_COVERAGE_MATRIX.md` — tracks *workflow-stage coverage*.
- `IMPROVEMENT_LOG.md` (this file) — tracks *behavioral improvements, bugs, and
  usability issues* discovered during real analyses.

## Entry Template

```markdown
### Ixxx: Short title

| Field | Value |
|-------|-------|
| **Module** | qc / preprocess / analysis / tumor / tools / utils / plotting |
| **Function(s)** | `run_standard_qc`, `suggest_qc_thresholds`, ... |
| **Discovered in** | project name / dataset / date |
| **Severity** | blocker / high / medium / low |
| **Status** | reported / triaged / in_progress / resolved / wontfix |
| **Phenomenon** | What went wrong or felt awkward? |
| **Expected behavior** | What should happen instead? |
| **Minimal repro** | 5–15 lines of code that trigger the issue |
| **Key output / error** | One error message, one number, or one small table |
| **Proposed fix** | Optional rough direction |
| **Resolution commit** | Filled when fixed |
| **Verification** | Filled when verified in a real project |
```

## Active Items

*(Populate as you discover issues during analysis.)*

## Resolved Items

*(Move items here once fixed and verified.)*

## Theme Backlog

Use this section to group related issues that should be tackled together in a
single polishing session.

### QC robustness

- *no entries yet*

### Preprocess defaults

- *no entries yet*

### Analysis / annotation usability

- *no entries yet*

### Tumor interpretation edge cases

- *no entries yet*

### Documentation / discoverability

- *no entries yet*

## Polishing Session Checklist

When you sit down to polish scLucid:

- [ ] Review active items and pick 1–3 related ones.
- [ ] Confirm each item has a minimal repro or clear expected behavior.
- [ ] Implement the fix.
- [ ] Add/update unit or integration tests.
- [ ] Update relevant docs (docstrings, this log, coverage matrix if needed).
- [ ] Run `python scripts/audit_public_api.py --write` if public API changed.
- [ ] Run `pytest` for affected modules.
- [ ] Commit with `Refs Ixxx` in the message.
- [ ] Verify in the project that reported the issue.
- [ ] Move the item to **Resolved Items** with resolution commit and verification.
