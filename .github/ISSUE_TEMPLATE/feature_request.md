---
name: Feature request
about: Suggest a new capability, API, or improvement
title: '[FEAT] '
labels: enhancement
---

## Problem

What workflow friction or missing capability is this aimed at? Concrete example
is helpful — e.g. *"I'm analyzing PDAC data and currently have to drop to R to
run BayesPrism deconvolution; I'd like a Python-native path."*

## Proposed solution

What the new API/behavior would look like, ideally with a code sketch:

```python
adata = scl.tl.deconvolve_bulk(
    adata_ref,
    bulk_data,
    method="DWLS",
    dampen_factor=1.0,
)
```

## Which user layer does this serve?

- [ ] Workflow (`scl.run_pipeline()`-style one-liner)
- [ ] Simple API (`scl.pp.*`, `scl.qc.*`, etc.)
- [ ] Advanced notebook user
- [ ] Plugin / extension author

## Which of the project's six core goals does this serve?

(Pick one or more; see the project vision in the README / `docs/`.)

- [ ] Efficiency / accuracy / flexibility
- [ ] Real-world EDA integration
- [ ] Automated parameter selection
- [ ] Traceability ("有据可循")
- [ ] R-package Python port
- [ ] Tumor-specific capability
- [ ] Publication-quality visualization

## Alternatives considered

Other tools, existing scLucid paths, or design alternatives you weighed.

## Willing to contribute?

- [ ] I can submit a PR
- [ ] I can help review
- [ ] I can test on real data once implemented
- [ ] I'm only proposing — would like a maintainer to drive it
