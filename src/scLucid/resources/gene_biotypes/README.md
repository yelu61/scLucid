Gene biotype resources for scLucid.

Resolution order used by `scLucid.preprocess.load_gene_biotypes()`:

1. Bundled package resource under `scLucid/resources/gene_biotypes/`
2. User-local cache under `~/.sclucid/gene_annotations/`
3. One-time download from Ensembl BioMart, then cache locally

Important:
- scLucid does not write downloaded files into installed package resources.
- To ship fully offline biotype support in a release, include a reference file here,
  for example `human_ensembl_latest.csv.gz`.
- Required columns: `gene_name`, `biotype`
- Optional columns: `gene_id`, `chromosome`, `start`, `end`

Current bundled references:
- `human_reference_latest.csv.gz`: extracted from official GENCODE human gene annotation GTF
- `mouse_reference_latest.csv.gz`: extracted from official GENCODE mouse gene annotation GTF
- `human_ensembl_latest.csv.gz`: extracted from official GENCODE human gene annotation GTF
- `mouse_ensembl_latest.csv.gz`: extracted from official GENCODE mouse gene annotation GTF

Note:
- The preferred filenames are now `*_reference_latest.csv.gz`.
- The historical `*_ensembl_latest.csv.gz` filenames are retained as compatibility aliases.
- The bundled tables themselves can be sourced from either official Ensembl BioMart exports
  or official GENCODE GTF releases, as long as the required columns are preserved.
