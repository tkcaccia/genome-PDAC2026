---
layout: default
title: Source coverage
---

# Source Coverage

[Overview](index.md) | [Figure and table index](figures-and-tables.md)

The source release selects current collaborator code from `Original Source and
Code`, the supplied figure bundle, its immune/TME update, and the additional WES
helpers in `All_figures/scripts`. Identical or superseded plot versions are not
published as competing implementations.

## Export Changes

- Scripts are grouped by scientific purpose; original filenames and SHA-256 checksums are retained in [the source manifest](source_manifest.json).
- Machine-specific absolute paths and actual cohort membership are replaced by private configuration readers. Population order, EPIC eligibility, KRAS groups and exclusion scenarios remain distinct inputs.
- Patient-specific MSI fallback values are removed. Without an original file the public code does not invent or recover its result from a literal embedded in the script.
- Result-dependent MOFA programme annotations and patient label offsets are private configuration, not public patient results.
- The mixed TCGA plotting helper contributes only its current mutation-frequency function. Its older mutation/RNA-candidate overlap analysis is excluded.
- The MSI table caption is aligned with the supplied categorization function: <5%, 5 to <20%, and >=20%, with exploratory interpretation retained.
- Shared import paths are updated to the new folder structure. Headers explain whether a file analyses data or only renders existing results.

Statistical algorithms are otherwise preserved. A source-only organization does
not establish numerical equivalence of a future rerun. No analysis functions,
plotting functions, R workflows, Slurm jobs, package installers or data downloaders
were executed during this update.

## Missing Original Producers

The following were not recoverable as exact generating code from the supplied
files: the primary/larger-cohort proteomics differential-analysis tables and their
volcano producers; extraction-method CV/peptide/protein/hydrophobicity/pI plots;
the additional peptide/protein supplementary panel; the clinical characteristics
table generator; proteomics search-parameter tables; and the original selected
proteomics enrichment table producers.

The protein-matrix, normalization and exploratory MOFA helpers are useful code,
but do not substitute for those missing original analyses. Likewise, a public
TCGA downloader does not provide the separate ancestry metadata automatically.

The figure bundle's `convert_figures.py` only wraps supplied PDF contents into
SVG/PDF outputs. Raster-backed source artwork remains raster-backed. Conversion
does not prove that underlying statistics were reproduced.

The source index records presentation discrepancies visible in the supplied code,
including volcano colouring versus formal DEG filtering, original versus current
supplement numbering, the subset MSI formatter versus the full manuscript table,
and enrichment top-five selection. These are documented, not silently reanalysed.

## What Was Retired

The current branch no longer contains earlier manuscript/report builders,
reviewer correspondence, custom differential-expression fallbacks, unrelated
historical Sarek/COSMIC/rnafusion workflows, example GMT substitutions or earlier
phenotype classifiers. Git history preserves the preceding repository state;
there is no legacy directory in the current source layout.

## Checks

`tools/check_repository.py` checks source parsing, configuration keys, provenance
hashes, file types, private-literal patterns and local documentation links. It
does not import the analysis modules or evaluate their results. The GitHub
workflow has no execution step for the study code, even with `--help`.
