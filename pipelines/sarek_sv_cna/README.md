## Sarek SV/CNA

This folder preserves the code configured to add somatic structural-variant and
copy-number analysis to the PDAC2026 matched tumour-normal WES cohort.

The archived original tumour-normal configuration enabled:

- `mutect2`
- `strelka`
- `vep`

That configuration requested somatic single-nucleotide variant (SNV),
insertion/deletion (indel) calling and annotation, but did not request:

- `manta` for somatic structural variants
- `ascat` for allele-specific copy number

The scripts in this folder run a dedicated Sarek retry with:

- `--tools manta,ascat`
- the same tumor-normal samplesheet
- the same capture BED
- the same shared work directory for cache reuse
- a separate output directory so existing SNV/indel outputs remain untouched

ASCAT on WES requires WES-specific loci, allele, GC, and replication-timing reference bundles. The remote launcher downloads the official hg38 WES bundles published by the ASCAT authors and then rewrites the loci files to `chr`-based format because the cohort CRAMs are aligned against a `chr`-based GRCh38 reference.

Interpretation limits:

- This is a WES-compatible SV/CNA workflow, not a WGS-grade structural-variation workflow.
- ASCAT WES outputs are useful for purity/ploidy, allele-specific copy-number and gene-level CNA summaries at captured loci.
- Manta WES outputs are useful for candidate SV and breakend prioritization in covered regions, but exome capture gives lower and less uniform breakpoint sensitivity than whole-genome sequencing.
- Interchromosomal translocations, CNV burden, SV burden and circos links should therefore be interpreted as exploratory summaries requiring orthogonal review before driver-level conclusions.

Primary references:

- [ASCAT GitHub](https://github.com/VanLoo-lab/ascat)
- [ASCAT WES reference README](https://github.com/VanLoo-lab/ascat/blob/master/ReferenceFiles/WES/README.md)
- [nf-core/sarek usage](https://nf-co.re/sarek/usage)

Required private configuration:

- `PDAC2026_RESOURCES_ENV`: reference environment file used by the Sarek wrapper.
- `PDAC_SV_CNA_SAMPLESHEET`: matched tumour-normal Sarek samplesheet.
- `PDAC_SV_CNA_INTERVALS`: exome capture BED.
- `PDAC_SV_CNA_OUTDIR`: restricted result directory.
- `PDAC_SV_CNA_WORKDIR`: restricted Nextflow work directory.
- `PDAC_ASCAT_REFERENCE_ROOT`: ASCAT WES reference directory.

No machine-specific drive address is embedded in the public launcher.

Launch flow:

1. Run [`launch_remote.sh`](launch_remote.sh) locally after setting `PDAC2026_SSH_HELPER` to a private SSH wrapper, or invoke [`run_sarek_sv_cna_remote.sh`](run_sarek_sv_cna_remote.sh) directly on the analysis workstation.
2. It uses the shared SSH helper to execute [`run_sarek_sv_cna_remote.sh`](run_sarek_sv_cna_remote.sh) on the remote Ubuntu workstation.
3. The remote script backgrounds the job under `nohup`, downloads/prepares ASCAT references if needed, and starts Sarek with `-resume`.

## Downstream Summaries

[`summarize_interchromosomal_translocations.py`](summarize_interchromosomal_translocations.py) summarizes Manta breakend records and reports interchromosomal translocations per tumor-normal pair.

[`make_patient_circos.py`](make_patient_circos.py) builds one circos-style PNG per patient when the required Sarek source outputs are available. Each plot includes all GRCh38 autosomes plus chrX/chrY, cytobands and centromeres, a chromosome-level PASS SNV load ring, 5 Mb PASS SNV density bars, ASCAT total-copy-number segments, Manta intrachromosomal SV links, highlighted interchromosomal translocations, PDAC driver-gene markers, matched tumor-normal RNA expression change, and likely biallelic-inactivation flags. Generated PNGs are result files and are intentionally not committed to GitHub.

Because the circos plots use WES-derived Manta and ASCAT calls, absence of a link or segment should not be interpreted as absence of a genome-wide SV/CNA event.

[`run_enhanced_circos_remote.sh`](run_enhanced_circos_remote.sh) provides the reusable command used to download a missing UCSC cytoband reference and regenerate enhanced circos plots. All private input/output paths are supplied with environment variables.

The enhanced circos run uses these additional inputs:

- UCSC hg38 cytobands supplied through `CYTOBAND`.
- GENCODE v46 gene coordinates supplied through `GTF`.
- STAR `ReadsPerGene.out.tab` root supplied through `STAR_ROOT`.
- Mutect2 and VEP roots supplied through `MUTECT2_ROOT` and `VEP_ROOT`.

Audit boundary: the August 2026 connected storage did not contain the detailed ASCAT and Manta source files needed to regenerate these plots. Existing plot code is preserved, but figure/event claims require restoration of those primary inputs.
