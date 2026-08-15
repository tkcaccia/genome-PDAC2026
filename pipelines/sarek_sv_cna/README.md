## Sarek SV/CNA

This folder contains the code used to add the missing somatic structural-variant and copy-number analysis to the PDAC 2026 tumor-normal WES cohort after the initial production run.

The original tumor-normal launch only enabled:

- `mutect2`
- `strelka`
- `vep`

That produced somatic SNV/indel calls and annotation, but not:

- `manta` for somatic structural variants
- `ascat` for allele-specific copy number

The scripts in this folder run a dedicated Sarek retry with:

- `--tools manta,ascat`
- the same tumor-normal samplesheet
- the same capture BED
- the same shared work directory for cache reuse
- a separate output directory so the finished SNV/indel results stay untouched

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

Expected remote outputs:

- `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna`
- `/media/user/New_Volume3/Lion/PDAC/SEQ_refs/ascat_wes_hg38`
- `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/circos_plots_fixed`
- `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/circos_plots_enhanced`

Launch flow:

1. Run [launch_remote.sh](/Users/stefano/Documents/SEQ/PDAC2026/sarek_sv_cna/launch_remote.sh) locally.
2. It uses the shared SSH helper to execute [run_sarek_sv_cna_remote.sh](/Users/stefano/Documents/SEQ/PDAC2026/sarek_sv_cna/run_sarek_sv_cna_remote.sh) on the remote Ubuntu workstation.
3. The remote script backgrounds the job under `nohup`, downloads/prepares ASCAT references if needed, and starts Sarek with `-resume`.

## Downstream Summaries

[summarize_interchromosomal_translocations.py](/Users/stefano/Documents/SEQ/PDAC2026/sarek_sv_cna/summarize_interchromosomal_translocations.py) summarizes Manta breakend records and reports interchromosomal translocations per tumor-normal pair.

[make_patient_circos.py](/Users/stefano/Documents/SEQ/PDAC2026/sarek_sv_cna/make_patient_circos.py) builds one circos-style PNG per patient from the completed Sarek outputs. Each plot includes all GRCh38 autosomes plus chrX/chrY, cytobands and centromeres, a chromosome-level PASS SNV load ring, 5 Mb PASS SNV density bars, ASCAT total-copy-number segments, Manta intrachromosomal SV links, highlighted interchromosomal translocations, PDAC driver-gene markers, matched tumor-normal RNA expression change, and likely biallelic-inactivation flags. Generated PNGs are result files and are intentionally not committed to GitHub.

Because the circos plots use WES-derived Manta and ASCAT calls, absence of a link or segment should not be interpreted as absence of a genome-wide SV/CNA event.

[run_enhanced_circos_remote.sh](/Users/stefano/Documents/SEQ/PDAC2026/sarek_sv_cna/run_enhanced_circos_remote.sh) records the exact remote command used to download the missing UCSC cytoband reference and regenerate the enhanced circos plots.

The enhanced circos run uses these additional inputs:

- UCSC hg38 cytobands downloaded to `/media/user/SEQ/refs/ucsc/hg38.cytoBand.txt`
- GENCODE v46 gene coordinates from `/media/user/SEQ/refs/annotation/gencode.v46.primary_assembly.annotation.gtf`
- STAR `ReadsPerGene.out.tab` files from `/home/user/PDAC_SEQ_native_results/rnafusion_pdac/star`
- Mutect2 VEP annotations from `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal/annotation/mutect2`
