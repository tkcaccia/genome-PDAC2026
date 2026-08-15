#!/usr/bin/env bash
set -Eeuo pipefail

run_root="/media/user/SEQ/logs/pdac_production_$(date +%Y%m%d_%H%M%S)_resume"
mkdir -p "$run_root"

cat > "$run_root/queue.sh" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail
export PATH="/home/user/.local/bin:\$PATH"
run_root="$run_root"
log() {
  printf '[%s] %s\n' "\$(date '+%F %T')" "\$*" | tee -a "\$run_root/master.log"
}
run_stage() {
  local name="\$1"
  shift
  log "Starting \$name"
  if "\$@" >"\$run_root/\${name}.log" 2>&1; then
    log "Completed \$name"
  else
    local rc=\$?
    log "Failed \$name (exit \$rc)"
    exit \$rc
  fi
}
run_stage rnaseq \
  /media/user/SEQ/scripts/run_rnaseq.sh \
  --samplesheet /media/user/SEQ/samplesheets/rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv \
  --outdir /media/user/PDAC_SEQ_analysis/results/rnaseq_expression \
  -resume
run_stage sarek_germline \
  /media/user/SEQ/scripts/run_sarek.sh \
  --mode germline \
  --samplesheet /media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \
  --intervals /media/user/SEQ/refs/optional/PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed \
  --outdir /media/user/PDAC_SEQ_analysis/results/sarek_germline \
  -resume
run_stage sarek_tumor_normal \
  /media/user/SEQ/scripts/run_sarek.sh \
  --mode tumor-normal \
  --samplesheet /media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \
  --intervals /media/user/SEQ/refs/optional/PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed \
  --outdir /media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal \
  -resume
run_stage rnafusion \
  /media/user/SEQ/scripts/run_rnafusion.sh \
  --samplesheet /media/user/SEQ/samplesheets/rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv \
  --outdir /media/user/SEQ/results/rnafusion_pdac \
  -resume
EOF

chmod +x "$run_root/queue.sh"
nohup "$run_root/queue.sh" > "$run_root/nohup.out" 2>&1 &
pid=$!

printf 'RUN_ROOT=%s\n' "$run_root"
printf 'PID=%s\n' "$pid"
ps -fp "$pid" || true

echo '== FSTAB CHECK =='
sudo -k -S grep -n "New Volume3\|C0D20D80D20D7C42\|sda2" /etc/fstab || true
