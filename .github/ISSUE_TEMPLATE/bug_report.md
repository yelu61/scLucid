---
name: Bug report
about: Report a bug, regression, or unexpected behavior in scLucid
title: '[BUG] '
labels: bug
---

## Summary

A clear, one-sentence description of what went wrong.

## Environment

- scLucid version: `import scLucid; print(scLucid.__version__)`
- Python version:
- OS:
- Install command: `pip install sclucid[...]` / `pip install -e .[...]`
- Key dependency versions: `pip freeze | grep -E 'scanpy|anndata|scrublet|celltypist|pydeseq2'`

## Dataset profile

- Shape: `adata.shape =`
- Species: human / mouse / other
- Tissue / cancer type:
- Single-sample or multi-sample:
- Counts layer present: yes / no
- (If tumor) primary / metastatic / treated:

## Steps to reproduce

```python
import scanpy as sc
import scLucid as scl

adata = sc.read_h5ad(...)
# minimal example that triggers the bug
adata = scl.run_pipeline(adata, ...)
```

## Expected behavior

What you expected to happen.

## Actual behavior

What actually happened, including the full traceback if applicable:

```
<paste traceback here>
```

## Additional context

Workarounds you've tried, related issues, or anything else that helps.
