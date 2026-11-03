---
layout: default
title: Reproducibility
---

# Reproducibility

[Overview](index.md) | [Figure and table index](figures-and-tables.md) | [Workflow](workflow.md)

Reproducing an analysis requires its code, software environment, reference files,
input data, sample matching and parameters. This repository provides public code
and method descriptions; patient-level inputs and results remain private.

## Implementation Coverage

| Component | Repository implementation | Additional requirements |
| --- | --- | --- |
| Downstream WES and KRAS | VEP-annotated VCF parsing, variant filtering, cohort summaries and evidence extraction | Annotated VCFs and sample configuration; raw-read alignment/calling/annotation job definitions are not included |
| MSI | Paired MSIsensor-pro launcher, raw-summary aggregation and table formatting | Matched BAMs/indexes, executable version and microsatellite-list provenance; full-cohort table assembly extends beyond the subset formatter |
| RNA differential expression | Salmon count import, paired R DESeq2, response-matrix assembly, PCA and enrichment | Gene quantifications, identifier mapping and the software/reference versions used for the analysis |
| Immune/stromal inference | Five inference methods, paired tests and EPIC convergence sensitivity | Gene-symbol TPM, matched metadata, method references and fit diagnostics |
| Programme and phenotype scoring | Named Hallmark sets, ssGSEA, feature standardisation, composites and median-based grouping | Selected GMT, mapped-gene coverage and complete tumour score matrices |
| Protein abundance and normalization | Protein-report parsing and integration/QC normalization | Quantified protein report and aligned sample metadata |
| Protein differential abundance, extraction-method figures and selected protein enrichment tables | Exact statistical/plotting producers are not included | Analysis scripts, numeric source tables, search/modification settings and extraction-method replicate definitions |
| MOFA | Model fitting, factor/loadings export, candidate scoring, reconstruction and sensitivity calculations | Aligned molecular matrices, identifier maps, model files and configured exclusions |
| TCGA context | Public matrix preparation, joint PCA, mutation frequencies, pathway themes and ancestry/modality summaries | Public expression/mutation/protein files and clinical metadata with ancestry labels |
| Clinical summary | No table-generation implementation included | Private clinical metadata and table-assembly code |
| Figure presentation | Vector rendering and PDF conversion | Numeric plotting inputs, or the source PDF for conversion-only panels |

The [figure/table index](figures-and-tables.md) gives panel-specific inputs and
identifies conversion-only figures. Converting a PDF into SVG preserves its
contents; it does not recover numeric data or reproduce the underlying statistics.

## Analytical Records

Record the following with each execution in the private analysis environment:

- Input and reference checksums, genome build, annotation release and sample ordering.
- Software and container versions, commands, thresholds and execution logs.
- Gene-mapping coverage, tested-gene sets and enrichment background definitions.
- Pairing, convergence flags, exclusions, multiple-testing families and score definitions.
- Model inputs, preprocessing, diagnostic outputs and the tables used for plotting.

The [code manifest](source_manifest.json) identifies script filenames and SHA-256
checksums. File hashes establish code identity, not the validity of a biological
conclusion.

## Private Configuration

Private input/output locations, cohort membership, KRAS comparison sets, EPIC
eligibility, exclusion scenarios and result-dependent annotations are read from
configuration outside the repository. The public template contains no patient
identifiers or results. See [the input guide](inputs.md).

## Automated Checks

`tools/check_repository.py` checks source parsing, configuration keys, file types,
privacy patterns, code checksums and documentation links. GitHub Actions runs
these static checks and shell-syntax validation without importing analysis
modules. These checks are distinct from executing or statistically validating
an analysis.
