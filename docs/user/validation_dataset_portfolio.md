# Validation Dataset Portfolio

## Current conclusion

The recommended portfolio is now explicit for QC, Preprocess, and the current
Analysis scope (annotation, clustering/identity preservation, pseudobulk DE,
and composition). Dataset registration is not scientific validation:

- `LOCAL_READY` or `DOWNLOADABLE` means that acquisition is possible.
- `PASS` requires an executed, versioned endpoint with the preregistered metric.
- missing metadata, missing artifacts, or inaccessible truth fail closed.
- LARRY and sci-Plex are future advanced-analysis rows and do not count toward
  current Analysis readiness.

The machine-readable source of truth is
`validation/dataset_evidence_registry.json`. It contains the full download
record, license/redistribution boundary, required metadata, validation
objective, endpoint mapping, and limitations for every row. The executed
evidence bindings are separate in `validation/evidence_run_index.json`.

## Preregistered endpoint families

| Module | Primary endpoint | Locked acceptance boundary |
|---|---|---|
| QC | Cell filtering | expert KEEP false removal <=2%; absolute low-quality recall gain >=5 percentage points against both global and per-sample MAD baselines; grouped-bootstrap 95% CI lower bound >0 |
| QC | Sample failure | catastrophic-sample recall 100%; locked high-quality false-block rate 0% |
| QC | Doublets | AUPRC regret <=5% versus the best registered method; ECE <=0.10 when probabilities exist; no unconfirmed automatic deletion |
| QC | Cell calling | true-cell recall >=95%; low-RNA-cell recall >=90%; empty-droplet FDR <=1% |
| QC | Ambient | contamination reduction >=50%; native-marker and identity loss each <=2% |
| Preprocess | Input/representation contract | exact count preservation and representation provenance 100%; zero consumer-contract violations or review mutations |
| Preprocess | Normalization selection | held-out utility regret <=5%; biology loss <=2%; shifted-log baseline required; zero silent fallback |
| Preprocess | Feature selection | utility regret <=5%; rare-class and program loss each <=2%; three feature-count variants pass; zero protected-marker injection |
| Preprocess | Selector generalization | held-out regret <=5%; biology loss <=2%; selected candidate consistent across three preregistered parameter variants |
| Preprocess | Graph stability | graph-seed and neighbor-identity loss each <=2% relative to baseline; partition-stability regret <=5% |
| Preprocess | Integration need/confounding | unintegrated baseline required; false `READY` under confounding 0%; unnecessary integration <=5%; Cramér's V >=0.7 blocks automatic integration |
| Preprocess | Integration Pareto | complex method must Pareto-improve; biology, rare-population, program, and graph-stability loss each <=2% |
| Preprocess | Identity/rare preservation | utility regret <=5%; rare-class recall loss and absolute abundance bias each <=2% |
| Preprocess | Tumor structure/programs | lineage, program-correlation, and patient-structure loss each <=2% relative to the simple unintegrated baseline |
| Preprocess | Policy execution | review mutation 0; policy/apply agreement, repeat agreement, count preservation, and provenance each 100% |
| Preprocess | DecisionCard UX | edited fields reduced >=70%; zero critical errors, workarounds, or project patches; RunEvidence 100% |
| Preprocess | Scalability | zero dense expansion or failures; repeat policy agreement 100%; runtime/peak memory and reference hardware reported |
| Analysis | Annotation | macro-F1 and balanced-accuracy regret <=5%; unknown detection must be reported |
| Analysis | Pseudobulk DE | empirical FDR <=0.06; sign concordance >=95% for detectable truth; donor/sample is the experimental unit |
| Analysis | Composition | MAE <=0.05; absolute bias <=0.02; nominal 95% interval coverage >=90% |

These thresholds are version `1.2.0`. Changing them after examining a new
benchmark requires a registry version bump and a written reason.

## Dataset rows

`Raw` below distinguishes raw reads, raw/unfiltered counts, and processed
objects. For GEO rows, public availability follows NCBI repository terms; it is
not represented as a dataset-specific open-source license. NCBI explicitly
notes that it does not transfer all submitter rights, so redistribution remains
under review even when analysis access is public.

| Dataset / accession | Modules; truth | Raw / processed availability | License boundary | Required metadata | Main endpoint IDs | Current evidence |
|---|---|---|---|---|---|---|
| Lin PDAC / [GSE154778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778) | QC, PP; RV | SRA reads; GEO counts/processed; local H5AD | GEO repository terms; redistribution review | sample, condition, chemistry, primary/metastatic, blinded sample/cell labels | profile, selector, catastrophic sample, damage, iterative review, rare/tumor retention, tumor structure | **BLOCKED**: expert labels incomplete; GSM4679533 remains a failure control |
| Schlesinger human PDAC / [GSE141017, GSM4293555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4293555) | QC, PP; RV | SRA reads; GEO count/processed; local H5AD | GEO repository terms | sample, species, tissue, expert cell-region label | `qc_rare_population_preservation`, `pp_tumor_structure_preservation` | **NOT RUN**; single sample only |
| Moncada multi-patient PDAC / [GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672) | QC, PP, Analysis; SS/RV | SRA reads; GEO counts/processed; local H5AD | GEO repository terms | patient, sample, condition, cell type/subtype, assay | rare-population retention, tumor structure, annotation | **NOT RUN**; author labels are silver standard |
| Zilionis NSCLC / [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465) | QC, PP, Analysis; SS/RV | human reads restricted; processed human matrix available; mouse reads available | GEO repository terms | patient, sample, tumor/blood, emulsion, cell type, species | `qc_rare_population_preservation`, `pp_tumor_structure_preservation`, `an_annotation_accuracy` | **NOT RUN** for locked endpoints |
| Lee CRC / [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) | QC, PP, Analysis; SS/RV | reads unavailable for privacy; processed UMI/annotation available | GEO repository terms | patient, sample, tumor/normal, cell type, pairing | retention, integration Pareto, annotation, composition | **NOT RUN** for locked endpoints |
| Kang/demuxlet / [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583) | QC, PP, Analysis; GT-A/B | SRA reads and GEO counts/processed; local H5AD | GEO repository terms | donor, condition, capture, demuxlet class, cell type, pairing | `qc_doublet_calibration`, `pp_identity_preservation`, `an_pseudobulk_de` | doublet **REVIEW** (AUPRC 0.605); donor-level execution contract passed, but FDR/power endpoint not run |
| Cell Hashing / [GSE108313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313) | QC, Analysis; GT-A | GEO candidate RNA matrix plus raw HTO counts acquired; complete unfiltered RNA droplets unavailable | GEO repository terms | capture, HTO class, sample hash, RNA barcode, identity | orthogonal doublet calibration, annotation | **ACQUIRED / PERFORMANCE NOT RUN**; HTO-derived labels do not support complete RNA cell-calling claims |
| Hashed mouse kidney / [GSE140262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140262) | QC; GT-A | SRA reads; GEO RNA/CMO matrices and labels | GEO repository terms | capture, CMO hash, doublet truth, cell type, kidney source | doublet calibration, rare-population preservation | **NOT DOWNLOADED / NOT RUN** |
| 10x human-mouse 6k / [10x dataset](https://www.10xgenomics.com/datasets/6-k-1-1-mixture-of-fresh-frozen-human-hek-293-t-and-mouse-nih-3-t-3-cells-2-standard-1-2-0) | QC, PP; GT-A/B | raw and filtered matrices acquired and checksum-verified; FASTQ intentionally not duplicated | **CC BY 4.0** | barcode, human/mouse UMI, filtered call, raw droplets, library | cell calling, doublet, ambient, identity preservation, PP scalability | **ACQUIRED / PERFORMANCE NOT RUN**; 737,280 droplets and 6,806 vendor-called barcodes prepared |
| SampleQC simulations + 172-sample QC metrics / [PMC9912498](https://pmc.ncbi.nlm.nih.gov/articles/PMC9912498/) | QC; SIM/SS | full expression unavailable; anonymized QC metrics and simulation code available | article/data CC BY 4.0 unless credited otherwise; code MIT | sample, QC group, simulated truth, QC metrics, rare QC population | profile, selector, sample failure, damage, iterative review, rare-population preservation | **NOT DOWNLOADED / NOT RUN** |
| Microscopy-labelled mESC capture sites / [E-MTAB-2600](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600) ([ENA PRJEB6455](https://www.ebi.ac.uk/ena/browser/view/PRJEB6455)) | QC; GT-A/SS | ENA reads; 960 source capture labels; local 168-run processed reconstruction with 147 explicit labels | mixed-source file rights under review; do not bundle | capture site, microscopy label, library, run mapping | legacy cell calling, damaged-cell classification | **ACQUIRED**; damage endpoint **FAIL** (0% damaged-cell REMOVE recall at 0% intact false removal); cell calling **BLOCKED** because no scLucid call output/low-RNA truth |
| Deterministic controlled QC truth / local generator | QC; SIM/ENG | locally generated counts and truth | scLucid code/output rights kept separate | library, assay, droplet class, policy label, lineage, ambient and damage fractions | input, mechanism controls, policy execution, scalability | engineering evidence only; never external scientific truth |
| scMixology / [GSE118767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767) | PP, Analysis; GT-A/B | SRA/subseries reads; GEO/CellBench matrices; local H5AD | GEO repository terms | protocol, mixture, identity, cell/RNA mixture, expected proportion | representation, normalization, features, selector, graph, integration need/Pareto, identity, policy execution, annotation, composition | selector binding **PASS**; integration binding **PASS_BASELINE**; other bindings **NOT RUN** |
| Mereu multicenter protocols / [GSE133544](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133544) | PP, Analysis; GT-B | SRA reads; GEO counts/processed | GEO repository terms | center, protocol, species, cell line, spike-in, replicate | normalization, features, selector, integration Pareto, identity, annotation, composition | **NOT DOWNLOADED / NOT RUN** |
| Ding cross-technology / [GSE132044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044) | PP, Analysis; GT-B/SS | SRA reads; GEO/SCP objects | GEO repository terms | sample type, technology, experiment day, replicate, species, cell type | normalization, features, selector, integration Pareto, identity, annotation | **NOT DOWNLOADED / NOT RUN** |
| scIB pancreas / [Figshare 25953868](https://figshare.com/articles/dataset/scIB_pancreas_dataset/25953868) | PP, Analysis; SS | full raw-count and processed H5AD objects | **CC BY 4.0** | batch, study, cell type, reference/query, donor if available | graph stability, integration need/Pareto, identity, annotation | **NOT DOWNLOADED / NOT RUN** |
| Zheng FACS mixtures / `DuoClustering2018` | PP, Analysis; GT-A | ExperimentHub count/processed objects; source FASTQ requires review | package/source rights must be reviewed separately | FACS cell type, mixture, expected proportion, barcode | identity, annotation, composition | **NOT DOWNLOADED / NOT RUN** |
| Baron pancreas / [GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133) | PP, Analysis; SS | SRA reads; GEO matrix; local H5AD | GEO repository terms | donor, sample, cell type, species | integration, identity, annotation | local object exists; locked endpoints **NOT RUN** |
| Kinker pan-cancer lines / [GSE157220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157220) | PP, Analysis; GT-B/SS | SRA reads; GEO/SCP542 processed | GEO/SCP source review | cell line, cancer type, library, program, low-quality/doublet label | features, selector, identity, tumor programs, PP scalability, annotation | **NOT DOWNLOADED / NOT RUN** |
| scDesign3 / [GitHub](https://github.com/SONGDONGYUAN1994/scDesign3) | QC, PP, Analysis; SIM | generated counts and truth | MIT software; generated data inherit source constraints | seed, version, source fingerprint, all simulated truths | QC stress tests, PP selector, integration need/confounding and Pareto, DE and composition | **NOT GENERATED / NOT RUN** |
| muscat simulations / [Bioconductor](https://bioconductor.org/packages/release/bioc/html/muscat.html) | Analysis; SIM | generated counts/truth | GPL-3 software; generated data inherit reference constraints | sample, group, cluster, gene truth class, true logFC, seed | `an_pseudobulk_de` | **NOT GENERATED / NOT RUN** |
| Squair DE compendium / [Zenodo 5048449](https://doi.org/10.5281/zenodo.5048449) | Analysis; GT-A/SS | benchmark objects available; original reads vary by source | record and source-dataset rights require verification before bundling | dataset, replicate, condition, cell type, matched bulk truth, null/signal | `an_pseudobulk_de` | **NOT DOWNLOADED / NOT RUN** |
| Three real projects / local controlled access | QC, PP, Analysis; RV | project-controlled | no redistribution without approval | project, sample, experimental unit, condition, batch/pairing, expert review, RunEvidence | all QC and all 12 PP heads plus current Analysis endpoints | **BLOCKED**: acceptance records incomplete |
| Deterministic Preprocess contract fixture / local generator | PP; SIM/ENG | transient synthetic UMI input; versioned benchmark report retained | scLucid code/output rights kept separate | seed, shape, count fingerprint, design keys, consumer, public API, check results | PP representation and policy execution contracts; scalability smoke input only | **CONTRACT PASS, NOT PERFORMANCE**; 120 cells × 160 genes cannot establish scale or scientific benefit |
| 10x PBMC 3k / [10x dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0) | QC, PP; ENG | FASTQ, raw matrix and Cell Ranger outputs available; local curated H5AD | **CC BY 4.0** | counts, species, single-donor status | PP representation, policy execution, and scalability contracts only | **ENGINEERING ONLY**; cannot support a scientific superiority claim |
| CellBender tiny fixture / tutorial-derived | QC; ENG | local derived counts | source rights not yet verified | barcode, rank, likely cell/empty, counts | none; contract tests only | **ENGINEERING ONLY** |
| LARRY / [GSE140802](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140802) | Analysis; GT-A/B | SRA reads; GEO clone/count/metadata | GEO repository terms | clone, time, condition, state, replicate | trajectory future gate | **P2; excluded from current CORE gate** |
| sci-Plex / [GSE139944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139944) | Analysis; GT-B | SRA reads; GEO matrices/annotations | GEO repository terms | cell line, compound, dose, well, replicate, hash | perturbation-program future gate | **P2; excluded from current CORE gate** |

## Required P0 release portfolios

- QC: Lin failure control, Kang, Cell Hashing, 10x HGMM, SampleQC
  simulation/QC metrics, E-MTAB-2600/PRJEB6455 microscopy labels, the deterministic
  engineering truth suite, and the three real projects.
- Preprocess engineering contracts: the deterministic controlled contract
  fixture, PBMC3k, scMixology, and the three real projects. Controlled
  technical generalization: scMixology, Mereu, Ding,
  scIB pancreas, 10x HGMM, and scDesign3. Tumor/scale generalization: Lin,
  Moncada, Zilionis, Lee, Kinker, and the three real projects.
- Current Analysis: Kang, scMixology, muscat, Squair, and the three real
  projects.

QC and Preprocess release are evaluated from their own exact
`endpoint × required dataset` combinations in `required_endpoint_portfolio`,
not by requiring every declared endpoint on every dataset. The complete
Preprocess mapping is documented in
[Full Preprocess Validation Protocol](preprocess_full_validation_protocol.md).
Endpoints belonging to another module do not cross-block it.
`NOT_EVALUABLE`, simulation-only, and contract-only evidence do not count as a
scientific pass. The existing scMixology selector result therefore passes only
that single binding and cannot promote the other 11 Preprocess heads.
