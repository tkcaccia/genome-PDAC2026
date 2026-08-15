#!/usr/bin/env bash
/media/user/SEQ/scripts/run_sarek.sh \
  --mode germline \
  --samplesheet /media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \
  --intervals /media/user/SEQ/refs/optional/PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed \
  --outdir /media/user/PDAC_SEQ_analysis/results/sarek_germline \
  -resume
