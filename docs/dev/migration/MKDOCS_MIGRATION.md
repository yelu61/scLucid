# MkDocs Migration

## Decision

scLucid is moving the maintained documentation site from Sphinx/RST to
MkDocs Material/Markdown.

The reason is practical: most project documentation is now design guidance,
workflow contracts, reviewer tables, roadmaps, and audit notes. Markdown fits
that style better than RST, while `mkdocstrings-python` can still generate API
reference pages from the package.

## Migration Policy

- The legacy Sphinx tree has been moved out of the active documentation path to
  `docs/archive/sphinx_source_legacy/`.
- Treat `mkdocs.yml`, `docs/index.md`, `docs/user/*.md`, and `docs/api/*.md`
  as the new primary documentation path.
- New user-facing documentation should be written in Markdown under `docs/`.
- Sphinx-specific files should not be updated for current docs; they are kept
  only as historical fallback during the transition window.

## Current Validation

The initial MkDocs migration has been validated locally with:

```bash
/Users/luye/micromamba/envs/scrna-env/bin/python -m mkdocs build --strict --site-dir /tmp/sclucid-mkdocs-site
```

The build completed successfully. Material for MkDocs prints an upstream notice
about the future MkDocs 2.0 direction, but it does not fail the build.

## Follow-Up Cleanup

- Decide whether to include roadmap phase pages in the main navigation or keep
  them reachable from `docs/roadmap/index.md` only.
- Delete `docs/archive/sphinx_source_legacy/` in a dedicated cleanup commit
  only if the project no longer needs a fallback copy of the old Sphinx source.
- If Quarto reports are added later, keep them under a separate report or
  validation-report path instead of mixing them into the package API docs.
