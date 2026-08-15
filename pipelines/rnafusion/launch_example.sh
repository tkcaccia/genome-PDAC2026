#!/usr/bin/env bash
/media/user/SEQ/scripts/run_rnafusion.sh \
  --samplesheet /media/user/SEQ/samplesheets/rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv \
  --tools arriba,fusioncatcher,salmon \
  --outdir /media/user/SEQ/results/rnafusion_pdac \
  -resume
