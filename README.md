# genome-PDAC2026

Analysis code and methods documentation for the PDAC2026 pancreatic cancer
multi-omics study. This repository makes the analytical data flow, statistical
methods, parameter choices and figure/table generation transparent.

The study integrates whole-exome sequencing, RNA sequencing and proteomics in
pancreatic ductal adenocarcinoma (PDAC). PDAC and periampullary cases are analysed
as separate cohorts.

**Code only.** No patient data, expression matrices, variant calls, result tables,
figures, manuscript, supplementary document or reviewer correspondence is included.

## Find the Code

- [Figure and table index](docs/figures-and-tables.md): start with a figure number, panel or table.
- [Analysis workflow](docs/workflow.md): where each input comes from and which scripts produce it.
- [Private inputs and configuration](docs/inputs.md): file formats, directory expectations and sample ordering.
- [Reproducibility](docs/source-coverage.md): implementation coverage and requirements for reproducing each analysis.
- [Code manifest](docs/source_manifest.json): script identifiers and SHA-256 checksums.
- [Website](https://tkcaccia.github.io/genome-PDAC2026/): the same documentation on GitHub Pages.

## Organization

| Folder | Contents |
| --- | --- |
| `scripts/wes/` | DeepSomatic/VEP parsing, gene-level filtering, mutation summaries, KRAS evidence, MSI launcher |
| `scripts/rnaseq/` | Salmon count/TPM assembly, tumour-normal responses, native R DESeq2, PCA and enrichment |
| `scripts/tme/` | quanTIseq, ESTIMATE, MCP-counter, EPIC, xCell, Hallmark ssGSEA and TME group assignment |
| `scripts/proteomics/` | Protein abundance-matrix assembly and normalization |
| `scripts/mofa/` | MOFA fitting, reconstruction summaries, candidate scoring and sensitivity analyses |
| `scripts/tcga/` | Public TCGA inputs, joint RNA PCA, mutation frequencies, pathway themes and modality availability |
| `scripts/periampullary/` | Separate PDAC/periampullary comparisons and KRAS plotting |
| `scripts/figures/` | Vector figure rendering and PDF conversion |
| `scripts/tables/` | KRAS/MSI and RNA pathway table-formatting helpers |
| `scripts/shared/`, `pdac_runtime.py` | Private configuration readers |
| `config/` | Empty example configuration only |
| `env/` | Python and R dependencies |
| `tools/` | Static repository checks that never import the analysis scripts |

The figure/table index links each panel to its input files, computational
producer and rendering function. It is the reference for matching publication
panels to script and output filenames.

## Scope and Reproducibility

The documentation distinguishes statistical analysis from presentation. Paired
DESeq2 tests tumour-normal gene expression; deconvolution estimates immune and
stromal features; ssGSEA scores tumour programmes; and MOFA integrates the
molecular layers. Their inputs, transformations and interpretation are documented
in the workflow.

Patient membership, private paths and result-dependent annotations are read from
configuration outside the repository. [Reproducibility](docs/source-coverage.md)
identifies analyses with computational implementations here and panels requiring
additional source material. PDF conversion is a presentation operation, not
statistical reproduction.

## Configuration and Checks

Read [the input guide](docs/inputs.md), populate a configuration **outside this
repository**, and prepare a new private output directory. Do not run all files
alphabetically: some fit models, some only format existing results, and several
create directories or read data at import time. Follow the documented dependencies
rather than treating the directory listing as an execution order.

For a source-only check, `python tools/check_repository.py` reads code and
documentation without importing or executing analysis scripts. GitHub Actions
uses only static checks; it does not run any study analysis.
