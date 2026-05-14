# Minimal Makefile for Sphinx documentation.
#
# Usage:
#     make html      build HTML into build/html
#     make clean     remove built docs
#     make linkcheck verify external URLs
#
# Set SCLUCID_PYTHON if you want to use a specific interpreter
# (defaults to whatever sphinx-build resolves to on PATH).

SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help html clean linkcheck

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(O)

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(O)
	@echo
	@echo "Build finished. HTML pages are in $(BUILDDIR)/html."

clean:
	rm -rf "$(BUILDDIR)"

linkcheck:
	$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(O)
