# RNA fusion analysis

This folder documents the patient-data-safe RNA fusion workflow used for the
PDAC project. It runs the official `nf-core/rnafusion` release `4.1.0`, pinned
to Git commit `eabd3f2d60d3e70d5e80d6c62b81d7d3473c36d1`.

No FASTQ files, samplesheets, sample identifiers, result tables, figures,
execution logs, or private filesystem paths belong in this repository.

## What this analysis uses

RNA fusion calling starts from the original paired-end RNA-seq FASTQ files. It
does not use WES variants, normalized expression values, or the tumour
phenotype assignments. Existing STAR `ReadsPerGene.out.tab` files are used
only to verify library strandedness when the samplesheet is prepared.

The configured tools are:

- Arriba for STAR-based fusion calling.
- FusionCatcher for an independent fusion-calling method.
- STAR-Fusion for an additional STAR-based caller.
- Salmon for transcript quantification used by the workflow.

`--no_cosmic` prevents `fusionreport` from requiring a licensed COSMIC data
download. COSMIC annotation, if licensed data are available locally, is a
separate downstream annotation step and must not be confused with fusion
calling.

## Files in this folder

- `prepare_samplesheet_from_fastqs.sh` identifies one original R1/R2 pair per
  sample and infers strandedness from STAR gene-count columns.
- `sync_references.sh` downloads and verifies the official public reference
  bundle from the nf-core AWS bucket.
- `rnafusion.config` limits Nextflow to one task, four CPUs, and at most 29 GB
  RAM for a workstation with 31 GiB physical memory.
- `run_rnafusion.sh` runs one samplesheet directly.
- `run_sample_queue.sh` runs one sample at a time, checks free RAM and disk,
  resumes interrupted work, and cleans only work from successful runs.
- `launch_example.sh` is a fully commented command template.
- `validate_samplesheet.py` checks the samplesheet structure.

## 1. Clone and pin the pipeline

```bash
git clone --branch 4.1.0 --depth 1 \
  https://github.com/nf-core/rnafusion.git /path/to/nf-core-rnafusion-4.1.0
git -C /path/to/nf-core-rnafusion-4.1.0 rev-parse HEAD
```

The second command should report
`eabd3f2d60d3e70d5e80d6c62b81d7d3473c36d1`.

## 2. Download and verify references

The complete public reference bundle is approximately 131 GiB. Choose a
destination with substantially more free space because containers, work files,
and outputs also require space.

```bash
./sync_references.sh \
  /large-disk/references/rnafusion-4.1.0 \
  /large-disk/logs/rnafusion-references \
  /large-disk/manifests/rnafusion_reference_manifest.tsv
```

The script compares the local file count and total bytes with a fresh AWS S3
inventory. The queue refuses to start unless the manifest status is
`verified`.

## 3. Build the nf-core samplesheet

Expected FASTQ layout:

```text
/fastq-root/SAMPLE_A/SAMPLE_A_R1_001.fastq.gz
/fastq-root/SAMPLE_A/SAMPLE_A_R2_001.fastq.gz
```

Expected STAR count layout:

```text
/star-count-root/SAMPLE_A/*ReadsPerGene.out.tab
```

Run:

```bash
./prepare_samplesheet_from_fastqs.sh \
  /fastq-root \
  /star-count-root \
  /private-config/rnafusion_samplesheet.csv
```

The script deliberately matches only `*_R1_001.fastq.gz` and
`*_R2_001.fastq.gz`; files whose names contain an additional trimming suffix
are excluded so the same library is not analysed twice. For STAR gene counts,
column 2 is unstranded, while columns 3 and 4 represent the two stranded count
orientations. The larger of columns 3 and 4 determines `forward` or `reverse`.
The PDAC samples all had a larger reverse-strand total and were therefore
entered as `reverse`. The accompanying `input_check.tsv` records the evidence
and fails if a sample has a missing or duplicated pair.

## 4. Run safely

For a single samplesheet:

```bash
./run_rnafusion.sh \
  --pipeline-root /path/to/nf-core-rnafusion-4.1.0 \
  --references /large-disk/references/rnafusion-4.1.0 \
  --samplesheet /private-config/rnafusion_samplesheet.csv \
  --outdir /large-disk/results/rnafusion \
  --work-dir /large-disk/work/rnafusion \
  --tools arriba,fusioncatcher,starfusion,salmon \
  -resume
```

For the low-RAM sequential mode used in this project, edit the paths in
`launch_example.sh` and run it. The queue requires at least 24 GiB available
RAM and 150 GiB free disk before starting each sample. If another desktop
application consumes too much RAM, the queue waits; it does not kill the
application or force Nextflow to start. Failed work is retained for `-resume`.
Only a successful sample whose execution report exists receives a completion
marker and has its Nextflow work cleaned.

The 29 GB cap is specific to a 31 GiB workstation. Lower it on a smaller host,
or use a scheduler/cloud profile on a larger system. Do not restore the former
48 GB request on a 31 GiB workstation because it can exhaust RAM and swap.

## 5. Verify completeness

A sample is complete only when all of the following are present:

1. Nextflow exits with status zero.
2. `pipeline_info/execution_report.html` exists and is non-empty.
3. `pipeline_info/execution_trace.tsv` shows completed requested modules.
4. The queue writes the sample status as `success` and creates `.completed`.

Caller-specific result files must then be inspected and merged using explicit
gene-pair, breakpoint, caller-support, read-support, blacklist, and biological
relevance criteria. A candidate is exploratory until independently validated;
the pipeline output alone is not proof of a true fusion.

## Relationship to the other RNA analyses

Differential expression, GSVA/ssGSEA programme scoring, and immune/stromal
deconvolution are separate downstream analyses. They use normalized or log
expression matrices. RNA fusion calling instead uses raw paired FASTQs. The
two result layers may be integrated later, for example to ask whether a fusion
is expressed, but one is not an input for the other.

Official documentation: <https://nf-co.re/rnafusion/4.1.0/>
