## Summary

<!-- 1-3 sentences describing what changed and why. -->

## Type of change

- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change (alters a contract, public API, or saved-file format)
- [ ] Documentation
- [ ] Test / infrastructure only

## Which scLucid goal does this advance?

- [ ] Efficiency / accuracy / flexibility
- [ ] Real-world EDA integration
- [ ] Automated parameter selection
- [ ] Traceability ("有据可循")
- [ ] R-package Python port (which: ___)
- [ ] Tumor-specific capability
- [ ] Publication-quality visualization

## Test plan

- [ ] `pytest` clean (full suite, ~4 min)
- [ ] Added tests for new code paths
- [ ] Manually verified on real data (specify dataset):
- [ ] `import scLucid` produces zero `ImportWarning`
- [ ] If touching workflows: PBMC golden path still runs

```
<paste command + key output excerpts>
```

## Contract changes

- [ ] None
- [ ] Changed AnnData layout — updated `docs/source/data_contracts.rst`
- [ ] Bumped `SCHEMA_VERSION` in `src/scLucid/utils/contracts.py`
- [ ] Updated review_summary keys

## Documentation

- [ ] Docstrings cover new public API
- [ ] Sphinx docs build (`make -C docs html`)
- [ ] Updated relevant `.rst` files / examples / README

## Reviewer notes

<!-- Anything specific you want feedback on, known limitations, or follow-ups
that intentionally fall outside this PR. -->
