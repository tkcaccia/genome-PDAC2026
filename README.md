# genome-PDAC2026

Code for the figures, tables and supporting analyses in the current pancreatic
cancer multi-omics study, reorganized from the collaborators' supplied source.
The current analysis distinguishes the primary PDAC cohort from the exploratory
periampullary cohort.

**Code only.** No patient data, expression matrices, variant calls, result tables,
figures, manuscript, supplementary document or reviewer correspondence is included.
No analysis, plotting, package installation or reference download was performed
during this reorganization.

## Find the Code

- [Figure and table index](docs/figures-and-tables.md): start with a figure number, panel or table.
- [Analysis workflow](docs/workflow.md): where each input comes from and which scripts produce it.
- [Private inputs and configuration](docs/inputs.md): file formats, directory expectations and sample ordering.
- [Source coverage](docs/source-coverage.md): supplied code versus missing original producers.
- [Source manifest](docs/source_manifest.json): original filenames, SHA-256 checksums and export changes.
- [Website](https://tkcaccia.github.io/genome-PDAC2026/): the same documentation on GitHub Pages.

## Organization

| Folder | Contents |
| --- | --- |
| `scripts/wes/` | DeepSomatic/VEP parsing, gene-level filtering, mutation summaries, KRAS evidence, MSI launcher |
| `scripts/rnaseq/` | Salmon count/TPM assembly, tumour-normal responses, native R DESeq2, PCA and enrichment |
| `scripts/tme/` | quanTIseq, ESTIMATE, MCP-counter, EPIC, xCell, Hallmark ssGSEA and TME group assignment |
| `scripts/proteomics/` | Supplied abundance-matrix and normalization helpers |
| `scripts/mofa/` | MOFA fitting, reconstruction summaries, candidate scoring and sensitivity analyses |
| `scripts/tcga/` | Public TCGA inputs, joint RNA PCA, mutation frequencies, pathway themes and modality availability |
| `scripts/periampullary/` | Separate PDAC/periampullary comparisons and KRAS plotting |
| `scripts/figures/` | Final supplied vector rendering and PDF conversion helpers |
| `scripts/tables/` | KRAS/MSI and RNA pathway table-formatting helpers |
| `scripts/shared/`, `pdac_runtime.py` | Private configuration readers |
| `config/` | Empty example configuration only |
| `env/` | Declared dependencies, not a newly installed or verified environment |
| `tools/` | Static repository checks that never import the analysis scripts |

Original filenames are retained where practical so the code can be traced back
to the supplied bundle. Older figure numbers inside filenames are mapped to the
current numbering in the index. These are not additional biological analyses.

## Scope and Reproducibility

This is an organized source release, not a newly validated end-to-end pipeline.
The scripts retain the collaborators' statistical choices and expected private
relative paths. Machine locations, patient lists, exclusion scenarios and
embedded patient results have been removed from the public export or moved to
private configuration. The original supplied files remain untouched.

Native DESeq2 is the current differential-expression implementation. Previous
custom differential-expression fallbacks, earlier manuscript workflows, unrelated
Sarek/COSMIC/rnafusion versions, report builders and duplicate plot revisions are
not part of this branch's working tree. They remain recoverable through Git
history, not in a confusing second "old versions" code directory.

Several original proteomics figure/table producers were not supplied. The index
identifies them explicitly. PDF-to-SVG conversion preserves an existing figure;
it does not reproduce its statistical analysis or recover editable plot data.

## Before Any Future Execution

Read [the input guide](docs/inputs.md), populate a configuration **outside this
repository**, and prepare a new private output directory. Do not run all files
alphabetically: some fit models, some only format existing results, and several
create directories or read data at import time. There is deliberately no automatic
`run_all.sh` in this source release.

For a source-only check, `python tools/check_repository.py` reads code and
documentation without importing or executing analysis scripts. GitHub Actions
uses only static checks; it does not run any study analysis.
