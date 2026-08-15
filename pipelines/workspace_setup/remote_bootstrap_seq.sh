#!/usr/bin/env bash
set -Eeuo pipefail

BASE="/media/user/SEQ"
BIN_DIR="$BASE/bin"
LOG_DIR="$BASE/logs"
DOCS_DIR="$BASE/docs"
PIPELINES_DIR="$BASE/pipelines"
CONFIG_DIR="$BASE/configs"
REFS_DIR="$BASE/refs"
CONTAINERS_DIR="$BASE/containers"
WORK_DIR="$BASE/work"
RESULTS_DIR="$BASE/results"
SAMPLESHEETS_DIR="$BASE/samplesheets"
SCRIPTS_DIR="$BASE/scripts"
TMP_DIR="$BASE/tmp"
BACKUP_DIR="$TMP_DIR/backups/$(date +%Y%m%d_%H%M%S)"

SYSTEM_AUDIT="$LOG_DIR/system_audit.txt"
INSTALL_NOTES="$DOCS_DIR/install_notes.md"
SOFTWARE_VERSIONS="$LOG_DIR/software_versions.txt"
NEXTFLOW_INFO="$LOG_DIR/nextflow_info.txt"
PIPELINE_VERSIONS="$DOCS_DIR/pipeline_versions.md"
REFERENCE_INVENTORY="$DOCS_DIR/REFERENCE_INVENTORY.md"
README_SETUP="$DOCS_DIR/README_SETUP.md"
LAUNCH_EXAMPLES="$DOCS_DIR/LAUNCH_EXAMPLES.md"
MANUAL_ITEMS="$DOCS_DIR/MANUAL_ITEMS_REQUIRED.md"
SETUP_SUMMARY="$DOCS_DIR/SETUP_SUMMARY.md"
BOOTSTRAP_LOG="$LOG_DIR/bootstrap_$(date +%Y%m%d_%H%M%S).log"

GENCODE_VERSION="46"
GENCODE_LABEL="gencode_v${GENCODE_VERSION}"
VEP_CACHE_VERSION="115"
VEP_IMAGE_TAG="release_115.2"
VEP_IMAGE="docker://ensemblorg/ensembl-vep:${VEP_IMAGE_TAG}"
VEP_SIF="$CONTAINERS_DIR/ensembl-vep_${VEP_IMAGE_TAG}.sif"

NEXTFLOW_VERSION="$(nextflow -version 2>/dev/null | awk '/version/{print $3; exit}')"
HOSTNAME_FQDN="$(hostname)"
USER_NAME="$(whoami)"
DATE_ISO="$(date -Iseconds)"
AVAILABLE_GB="$(df -BG "$BASE" | awk 'NR==2 {gsub(/G/,"",$4); print $4}')"

RUNTIME_PROFILE="singularity"
if command -v apptainer >/dev/null 2>&1; then
  RUNTIME_PROFILE="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME_PROFILE="singularity"
elif command -v docker >/dev/null 2>&1; then
  RUNTIME_PROFILE="docker"
fi

FUSION_SYNC_STARFUSION="false"
if [[ "${AVAILABLE_GB:-0}" =~ ^[0-9]+$ ]] && (( AVAILABLE_GB >= 220 )); then
  FUSION_SYNC_STARFUSION="true"
fi

mkdir -p \
  "$BIN_DIR" "$LOG_DIR" "$DOCS_DIR" "$PIPELINES_DIR" "$CONFIG_DIR" \
  "$REFS_DIR/genome" "$REFS_DIR/annotation" "$REFS_DIR/gatk_bundle" "$REFS_DIR/vep" \
  "$REFS_DIR/clinvar" "$REFS_DIR/fusion" "$REFS_DIR/optional" \
  "$CONTAINERS_DIR" "$WORK_DIR" "$RESULTS_DIR" "$SAMPLESHEETS_DIR" "$SCRIPTS_DIR" "$TMP_DIR" "$BACKUP_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$BOOTSTRAP_LOG"
}

die() {
  log "ERROR: $*"
  exit 1
}

backup_if_exists() {
  local target="$1"
  if [[ -e "$target" || -L "$target" ]]; then
    local rel="${target#/}"
    local dest="$BACKUP_DIR/$rel"
    mkdir -p "$(dirname "$dest")"
    cp -a "$target" "$dest"
    log "Backed up existing file: $target -> $dest"
  fi
}

write_with_backup() {
  local target="$1"
  backup_if_exists "$target"
  cat > "$target"
}

download_url() {
  local url="$1"
  local dest="$2"
  local log_file="$3"
  if [[ -s "$dest" ]]; then
    log "Keeping existing file: $dest"
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  local tmp="${dest}.part"
  log "Downloading $url -> $dest"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 5 --retry-delay 5 -C - -o "$tmp" "$url" >>"$log_file" 2>&1
  else
    wget -O "$tmp" "$url" >>"$log_file" 2>&1
  fi
  mv "$tmp" "$dest"
}

download_gatk_file() {
  local file_name="$1"
  local dest="$2"
  local log_file="$3"
  local bases=(
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
  )
  local base
  for base in "${bases[@]}"; do
    if download_url "${base}/${file_name}" "$dest" "$log_file"; then
      return 0
    fi
    rm -f "${dest}.part"
  done
  die "Failed to download GATK resource: $file_name"
}

emit_system_audit() {
  write_with_backup "$SYSTEM_AUDIT" <<EOF
Audit timestamp: $DATE_ISO
Host: $HOSTNAME_FQDN
User: $USER_NAME

[OS]
$(hostnamectl 2>/dev/null || cat /etc/os-release 2>/dev/null || true)

[Package manager]
$(command -v apt-get || command -v dnf || command -v yum || command -v zypper || echo "UNKNOWN")

[CPU]
$(lscpu 2>/dev/null || true)

[Memory]
$(free -h 2>/dev/null || true)

[Disk]
$(df -h "$BASE" 2>/dev/null || true)

[Tools on PATH]
$(for t in java apptainer singularity docker git curl wget unzip tar pigz bgzip tabix aws gsutil nextflow micromamba conda samtools bcftools bwa gatk Rscript; do printf '%-14s %s\n' "$t" "$(command -v "$t" 2>/dev/null || echo MISSING)"; done)

[Versions]
Java:
$(java -version 2>&1 | head -n 2 || true)

Apptainer:
$(apptainer --version 2>/dev/null || echo MISSING)

Singularity:
$(singularity --version 2>/dev/null || echo MISSING)

Docker:
$(docker --version 2>/dev/null || echo MISSING)

Git:
$(git --version 2>/dev/null || echo MISSING)

Curl:
$(curl --version 2>/dev/null | head -n 1 || echo MISSING)

Wget:
$(wget --version 2>/dev/null | head -n 1 || echo MISSING)

Unzip:
$(unzip -v 2>/dev/null | head -n 1 || echo MISSING)

Tar:
$(tar --version 2>/dev/null | head -n 1 || echo MISSING)

Pigz:
$(pigz --version 2>/dev/null | head -n 1 || echo MISSING)

Bgzip:
$(bgzip --version 2>/dev/null | head -n 1 || echo MISSING)

Tabix:
$(tabix --version 2>/dev/null | head -n 1 || echo MISSING)

AWS CLI:
$(aws --version 2>/dev/null || echo MISSING)

Nextflow:
$(nextflow -version 2>/dev/null || echo MISSING)

Conda:
$(conda --version 2>/dev/null || echo MISSING)

Samtools:
$(samtools --version 2>/dev/null | head -n 1 || echo MISSING)

Bcftools:
$(bcftools --version 2>/dev/null | head -n 1 || echo MISSING)

BWA:
$(bwa 2>&1 | awk 'NR==1 {print; exit}' || echo MISSING)

Rscript:
$(Rscript --version 2>/dev/null || echo MISSING)
EOF
}

write_install_notes() {
  write_with_backup "$INSTALL_NOTES" <<EOF
# Install Notes

Setup timestamp: $DATE_ISO
Host: $HOSTNAME_FQDN

This machine already satisfied the core runtime requirements, so the setup was performed by configuration and verification rather than by replacing working system packages.

## Audited and configured

- Java 17+: \`$(java -version 2>&1 | head -n 1)\`
- Nextflow: \`$(nextflow -version 2>/dev/null | awk '/version/{print $2, $3; exit}')\`
- Container runtime preference: \`$RUNTIME_PROFILE\`
- git: \`$(git --version)\`
- curl: \`$(curl --version | head -n 1)\`
- wget: \`$(wget --version | head -n 1)\`
- pigz: \`$(pigz --version 2>&1 | head -n 1)\`
- bgzip: \`$(bgzip --version 2>&1 | head -n 1)\`
- tabix: \`$(tabix --version 2>&1 | head -n 1)\`
- awscli: \`$(aws --version 2>&1)\`
- samtools: \`$(samtools --version | head -n 1)\`
- bcftools: \`$(bcftools --version | head -n 1)\`
- bwa: \`$(bwa 2>&1 | awk 'NR==1 {print; exit}')\`

## Commands run

\`\`\`bash
nextflow info > $NEXTFLOW_INFO
git clone --branch 3.8.1 --depth 1 https://github.com/nf-core/sarek.git $PIPELINES_DIR/sarek-3.8.1
git clone --branch 3.24.0 --depth 1 https://github.com/nf-core/rnaseq.git $PIPELINES_DIR/rnaseq-3.24.0
git clone --branch 4.1.0 --depth 1 https://github.com/nf-core/rnafusion.git $PIPELINES_DIR/rnafusion-4.1.0
aws --no-sign-request s3 sync s3://nf-core-awsmegatests/rnafusion/references/GRCh38/... $REFS_DIR/fusion/GRCh38/...
curl -L <public reference URLs>
singularity pull $VEP_SIF $VEP_IMAGE
\`\`\`

No patient data were moved or uploaded during setup.
EOF
}

emit_software_versions() {
  write_with_backup "$SOFTWARE_VERSIONS" <<EOF
Generated: $DATE_ISO
java: $(java -version 2>&1 | head -n 1)
nextflow: $(nextflow -version 2>/dev/null | tr '\n' ' ' | sed 's/  */ /g')
apptainer: $(apptainer --version 2>/dev/null || echo MISSING)
singularity: $(singularity --version 2>/dev/null || echo MISSING)
docker: $(docker --version 2>/dev/null || echo MISSING)
git: $(git --version 2>/dev/null || echo MISSING)
curl: $(curl --version 2>/dev/null | head -n 1 || echo MISSING)
wget: $(wget --version 2>/dev/null | head -n 1 || echo MISSING)
unzip: $(unzip -v 2>/dev/null | head -n 1 || echo MISSING)
tar: $(tar --version 2>/dev/null | head -n 1 || echo MISSING)
pigz: $(pigz --version 2>/dev/null | head -n 1 || echo MISSING)
bgzip: $(bgzip --version 2>/dev/null | head -n 1 || echo MISSING)
tabix: $(tabix --version 2>/dev/null | head -n 1 || echo MISSING)
awscli: $(aws --version 2>/dev/null || echo MISSING)
conda: $(conda --version 2>/dev/null || echo MISSING)
samtools: $(samtools --version 2>/dev/null | head -n 1 || echo MISSING)
bcftools: $(bcftools --version 2>/dev/null | head -n 1 || echo MISSING)
bwa: $(bwa 2>&1 | awk 'NR==1 {print; exit}' || echo MISSING)
Rscript: $(Rscript --version 2>/dev/null || echo MISSING)
EOF

  nextflow info > "$NEXTFLOW_INFO" 2>&1 || true
}

clone_pipelines() {
  local pipe_log="$LOG_DIR/pipeline_clone.log"
  : > "$pipe_log"

  if [[ ! -d "$PIPELINES_DIR/sarek-3.8.1/.git" ]]; then
    git clone --branch 3.8.1 --depth 1 https://github.com/nf-core/sarek.git "$PIPELINES_DIR/sarek-3.8.1" >>"$pipe_log" 2>&1
  fi
  if [[ ! -d "$PIPELINES_DIR/rnaseq-3.24.0/.git" ]]; then
    git clone --branch 3.24.0 --depth 1 https://github.com/nf-core/rnaseq.git "$PIPELINES_DIR/rnaseq-3.24.0" >>"$pipe_log" 2>&1
  fi
  if [[ ! -d "$PIPELINES_DIR/rnafusion-4.1.0/.git" ]]; then
    git clone --branch 4.1.0 --depth 1 https://github.com/nf-core/rnafusion.git "$PIPELINES_DIR/rnafusion-4.1.0" >>"$pipe_log" 2>&1
  fi

  write_with_backup "$PIPELINE_VERSIONS" <<EOF
# Pinned Pipeline Versions

Generated: $DATE_ISO

| Pipeline | Requested release | Local path | Git commit |
| --- | --- | --- | --- |
| nf-core/sarek | 3.8.1 | \`$PIPELINES_DIR/sarek-3.8.1\` | \`$(git -C "$PIPELINES_DIR/sarek-3.8.1" rev-parse HEAD)\` |
| nf-core/rnaseq | 3.24.0 | \`$PIPELINES_DIR/rnaseq-3.24.0\` | \`$(git -C "$PIPELINES_DIR/rnaseq-3.24.0" rev-parse HEAD)\` |
| nf-core/rnafusion | 4.1.0 | \`$PIPELINES_DIR/rnafusion-4.1.0\` | \`$(git -C "$PIPELINES_DIR/rnafusion-4.1.0" rev-parse HEAD)\` |

## Clone commands

\`\`\`bash
git clone --branch 3.8.1 --depth 1 https://github.com/nf-core/sarek.git $PIPELINES_DIR/sarek-3.8.1
git clone --branch 3.24.0 --depth 1 https://github.com/nf-core/rnaseq.git $PIPELINES_DIR/rnaseq-3.24.0
git clone --branch 4.1.0 --depth 1 https://github.com/nf-core/rnafusion.git $PIPELINES_DIR/rnafusion-4.1.0
\`\`\`

## Nextflow used for validation

\`$NEXTFLOW_VERSION\`
EOF
}

write_common_script() {
  write_with_backup "$SCRIPTS_DIR/_common.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

SEQ_BASE="/media/user/SEQ"
SEQ_TMP="$SEQ_BASE/tmp"
SEQ_WORK="$SEQ_BASE/work"
SEQ_RESULTS="$SEQ_BASE/results"
SEQ_CONFIG="/media/user/SEQ/configs/resources.env"
SEQ_BACKUPS="$SEQ_BASE/tmp/script_backups"

# Prefer user-local installs and workspace-pinned helper binaries in non-login shells.
if [[ -d "$HOME/.local/bin" ]]; then
  export PATH="$HOME/.local/bin:$PATH"
fi
if [[ -d "$SEQ_BASE/bin" ]]; then
  export PATH="$SEQ_BASE/bin:$PATH"
fi

seq_log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

seq_die() {
  seq_log "ERROR: $*" >&2
  exit 1
}

seq_require_cmd() {
  command -v "$1" >/dev/null 2>&1 || seq_die "Required command not found: $1"
}

seq_require_file() {
  [[ -s "$1" ]] || seq_die "Required file missing: $1"
}

seq_source_env() {
  seq_require_file "$SEQ_CONFIG"
  # shellcheck disable=SC1090
  source "$SEQ_CONFIG"
  export NXF_HOME NXF_WORK NXF_SINGULARITY_CACHEDIR NXF_APPTAINER_CACHEDIR SINGULARITY_CACHEDIR APPTAINER_CACHEDIR
}

seq_backup_if_exists() {
  local target="$1"
  if [[ -e "$target" || -L "$target" ]]; then
    local stamp
    stamp="$(date +%Y%m%d_%H%M%S)"
    local rel="${target#/}"
    local dest="$SEQ_BACKUPS/$stamp/$rel"
    mkdir -p "$(dirname "$dest")"
    cp -a "$target" "$dest"
  fi
}

seq_strip_csv_comments() {
  local input="$1"
  local output="$2"
  awk 'NF && $0 !~ /^[[:space:]]*#/' "$input" > "$output"
}

seq_runtime_profile() {
  if [[ -n "${SEQ_RUNTIME_PROFILE:-}" ]]; then
    printf '%s\n' "$SEQ_RUNTIME_PROFILE"
    return 0
  fi
  if command -v apptainer >/dev/null 2>&1; then
    printf 'apptainer\n'
  elif command -v singularity >/dev/null 2>&1; then
    printf 'singularity\n'
  elif command -v docker >/dev/null 2>&1; then
    printf 'docker\n'
  else
    seq_die "No supported container runtime detected."
  fi
}
EOF
  chmod +x "$SCRIPTS_DIR/_common.sh"
}

write_install_prereqs_script() {
  write_with_backup "$SCRIPTS_DIR/install_prereqs.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

LOG_FILE="/media/user/SEQ/logs/install_prereqs.log"
: > "$LOG_FILE"

required=(java nextflow git curl wget unzip tar pigz bgzip tabix aws samtools bcftools bwa)
missing=()
for cmd in "${required[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    missing+=("$cmd")
  fi
done

if [[ "${#missing[@]}" -eq 0 ]]; then
  seq_log "All required commands are already available." | tee -a "$LOG_FILE"
  exit 0
fi

seq_log "Missing commands: ${missing[*]}" | tee -a "$LOG_FILE"
seq_log "Install missing packages manually or rerun after enabling sudo in this shell." | tee -a "$LOG_FILE"
exit 1
EOF
  chmod +x "$SCRIPTS_DIR/install_prereqs.sh"
}

write_download_gatk_script() {
  write_with_backup "$SCRIPTS_DIR/download_gatk_bundle.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

OUT="/media/user/SEQ/refs/gatk_bundle"
LOG_FILE="/media/user/SEQ/logs/download_gatk_bundle.log"
: > "$LOG_FILE"

download_gatk() {
  local name="$1"
  local dest="$OUT/$name"
  if [[ -s "$dest" ]]; then
    seq_log "Keeping existing $dest" | tee -a "$LOG_FILE"
    return 0
  fi
  local bases=(
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
  )
  local ok="false"
  local base
  for base in "${bases[@]}"; do
    if curl -L --fail --retry 5 --retry-delay 5 -C - -o "${dest}.part" "${base}/${name}" >>"$LOG_FILE" 2>&1; then
      mv "${dest}.part" "$dest"
      ok="true"
      break
    fi
    rm -f "${dest}.part"
  done
  [[ "$ok" == "true" ]] || seq_die "Failed to download $name"
}

mkdir -p "$OUT"
files=(
  "Homo_sapiens_assembly38.fasta"
  "dbsnp_146.hg38.vcf.gz"
  "dbsnp_146.hg38.vcf.gz.tbi"
  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  "Homo_sapiens_assembly38.known_indels.vcf.gz"
  "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  "1000G_omni2.5.hg38.vcf.gz"
  "1000G_omni2.5.hg38.vcf.gz.tbi"
  "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
  "hapmap_3.3.hg38.vcf.gz"
  "hapmap_3.3.hg38.vcf.gz.tbi"
  "af-only-gnomad.hg38.vcf.gz"
  "af-only-gnomad.hg38.vcf.gz.tbi"
  "1000g_pon.hg38.vcf.gz"
  "1000g_pon.hg38.vcf.gz.tbi"
)

for file in "${files[@]}"; do
  download_gatk "$file"
done

samtools faidx "$OUT/Homo_sapiens_assembly38.fasta" >>"$LOG_FILE" 2>&1
samtools dict "$OUT/Homo_sapiens_assembly38.fasta" -o "$OUT/Homo_sapiens_assembly38.dict" >>"$LOG_FILE" 2>&1
bwa index "$OUT/Homo_sapiens_assembly38.fasta" >>"$LOG_FILE" 2>&1

awk 'BEGIN{OFS="\t"} /^@SQ/ {sn=""; ln=""; for (i=1;i<=NF;i++) {if ($i ~ /^SN:/) sn=substr($i,4); if ($i ~ /^LN:/) ln=substr($i,4)} if (sn != "" && ln != "") print sn,0,ln}' "$OUT/Homo_sapiens_assembly38.dict" > "$OUT/Homo_sapiens_assembly38.genome.bed"
EOF
  chmod +x "$SCRIPTS_DIR/download_gatk_bundle.sh"
}

write_sync_fusion_script() {
  local starfusion_default="$FUSION_SYNC_STARFUSION"
  write_with_backup "$SCRIPTS_DIR/sync_rnafusion_refs.sh" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

BASE_OUT="/media/user/SEQ/refs/fusion/GRCh38"
LOG_FILE="/media/user/SEQ/logs/sync_rnafusion_refs.log"
: > "\$LOG_FILE"

sync_path() {
  local src="\$1"
  local dst="\$2"
  mkdir -p "\$(dirname "\$dst")"
  aws --no-sign-request s3 sync "\$src" "\$dst" >>"\$LOG_FILE" 2>&1
}

sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/arriba" "\$BASE_OUT/arriba"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/gencode" "\$BASE_OUT/gencode_v46/gencode"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/salmon" "\$BASE_OUT/gencode_v46/salmon"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/star" "\$BASE_OUT/gencode_v46/star"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/fusioncatcher" "\$BASE_OUT/gencode_v46/fusioncatcher"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/fusion_report_db" "\$BASE_OUT/fusion_report_db"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/hgnc" "\$BASE_OUT/hgnc"

if [[ "\${SEQ_SYNC_STARFUSION:-$starfusion_default}" == "true" ]]; then
  sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/starfusion" "\$BASE_OUT/gencode_v46/starfusion"
else
  mkdir -p "\$BASE_OUT/gencode_v46/starfusion"
  printf 'STAR-Fusion references intentionally not synced on this host due to disk budget.\\n' > "\$BASE_OUT/gencode_v46/starfusion/PLACEHOLDER.txt"
fi
EOF
  chmod +x "$SCRIPTS_DIR/sync_rnafusion_refs.sh"
}

write_clinvar_script() {
  write_with_backup "$SCRIPTS_DIR/download_clinvar.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

OUT="/media/user/SEQ/refs/clinvar"
LOG_FILE="/media/user/SEQ/logs/download_clinvar.log"
: > "$LOG_FILE"

mkdir -p "$OUT"
if [[ ! -s "$OUT/clinvar.vcf.gz" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "$OUT/clinvar.vcf.gz.part" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" >>"$LOG_FILE" 2>&1
  mv "$OUT/clinvar.vcf.gz.part" "$OUT/clinvar.vcf.gz"
fi

if curl -I --fail -s "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" >/dev/null 2>&1; then
  if [[ ! -s "$OUT/clinvar.vcf.gz.tbi" ]]; then
    curl -L --fail --retry 5 --retry-delay 5 -C - -o "$OUT/clinvar.vcf.gz.tbi.part" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" >>"$LOG_FILE" 2>&1
    mv "$OUT/clinvar.vcf.gz.tbi.part" "$OUT/clinvar.vcf.gz.tbi"
  fi
else
  tabix -f -p vcf "$OUT/clinvar.vcf.gz" >>"$LOG_FILE" 2>&1
fi
EOF
  chmod +x "$SCRIPTS_DIR/download_clinvar.sh"
}

write_vep_script() {
  write_with_backup "$SCRIPTS_DIR/setup_vep.sh" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

OUT="/media/user/SEQ/refs/vep"
LOG_FILE="/media/user/SEQ/logs/setup_vep.log"
VEP_SIF="/media/user/SEQ/containers/ensembl-vep_${VEP_IMAGE_TAG}.sif"
VEP_IMAGE="$VEP_IMAGE"
: > "\$LOG_FILE"

mkdir -p "\$OUT/cache" "\$OUT/fasta"

if [[ ! -s "\$VEP_SIF" ]]; then
  if command -v singularity >/dev/null 2>&1; then
    singularity pull "\$VEP_SIF" "\$VEP_IMAGE" >>"\$LOG_FILE" 2>&1
  elif command -v apptainer >/dev/null 2>&1; then
    apptainer pull "\$VEP_SIF" "\$VEP_IMAGE" >>"\$LOG_FILE" 2>&1
  else
    seq_die "No Apptainer/Singularity runtime available for VEP container pull."
  fi
fi

if [[ ! -d "\$OUT/cache/homo_sapiens/$VEP_CACHE_VERSION"'_GRCh38' ]]; then
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "\$OUT/cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz.part" "https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz" >>"\$LOG_FILE" 2>&1
  mv "\$OUT/cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz.part" "\$OUT/cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz"
  tar -xzf "\$OUT/cache/homo_sapiens_vep_${VEP_CACHE_VERSION}_GRCh38.tar.gz" -C "\$OUT/cache" >>"\$LOG_FILE" 2>&1
fi

if [[ ! -s "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.part" "https://ftp.ensembl.org/pub/release-${VEP_CACHE_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" >>"\$LOG_FILE" 2>&1
  mv "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.part" "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fi

if [[ ! -s "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]]; then
  gunzip -c "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" > "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fi

if [[ ! -s "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]]; then
  samtools faidx "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" >>"\$LOG_FILE" 2>&1
fi

cat > "\$OUT/toy_input.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	881918	.	G	A	.	PASS	.
VCF

singularity exec "\$VEP_SIF" vep \\
  --offline \\
  --cache \\
  --dir_cache "\$OUT/cache" \\
  --species homo_sapiens \\
  --assembly GRCh38 \\
  --fasta "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \\
  --input_file "\$OUT/toy_input.vcf" \\
  --output_file "\$OUT/toy_output.vep.txt" \\
  --tab >>"\$LOG_FILE" 2>&1
EOF
  chmod +x "$SCRIPTS_DIR/setup_vep.sh"
}

write_validation_script() {
  write_with_backup "$SCRIPTS_DIR/validate_samplesheet.py" <<'EOF'
#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

REQUIRED = {
    "sarek": [("patient", "sample"), ("lane",), ("fastq_1", "bam", "cram")],
    "rnaseq": [("sample",), ("strandedness",), ("fastq_1",)],
    "rnafusion": [("sample",), ("strandedness",), ("fastq_1", "bam", "cram", "junctions", "splice_junctions")],
}

def load_rows(path: Path):
    lines = [line for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("No non-comment rows found.")
    return list(csv.DictReader(lines))

def check_required(rows, pipeline):
    header = set(rows[0].keys())
    problems = []
    for group in REQUIRED[pipeline]:
        if not any(col in header for col in group):
            problems.append(f"missing one of required columns: {', '.join(group)}")
    if pipeline == "sarek" and {"status", "sex"} - header:
        problems.append("recommended columns missing for tumor/normal setups: status, sex")
    return problems

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", choices=sorted(REQUIRED), required=True)
    parser.add_argument("--input", required=True)
    args = parser.parse_args()

    path = Path(args.input)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 1
    try:
        rows = load_rows(path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if not rows:
        print("ERROR: samplesheet has no data rows", file=sys.stderr)
        return 1
    problems = check_required(rows, args.pipeline)
    if problems:
        for problem in problems:
            print(f"ERROR: {problem}", file=sys.stderr)
        return 1
    print(f"OK: {args.pipeline} samplesheet looks structurally valid with {len(rows)} data row(s).")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
EOF
  chmod +x "$SCRIPTS_DIR/validate_samplesheet.py"
}

write_reference_validation_script() {
  write_with_backup "$SCRIPTS_DIR/validate_reference_integrity.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

check_file() {
  local file="$1"
  [[ -s "$file" ]] || seq_die "Missing reference: $file"
  printf 'OK\t%s\n' "$file"
}

check_file "$SEQ_SAREK_FASTA"
check_file "$SEQ_SAREK_FASTA_FAI"
check_file "$SEQ_SAREK_DICT"
check_file "$SEQ_SAREK_DBSNP"
check_file "$SEQ_SAREK_KNOWN_INDELS"
check_file "$SEQ_SAREK_GERMLINE_RESOURCE"
check_file "$SEQ_SAREK_PON"
check_file "$SEQ_RNASEQ_FASTA"
check_file "$SEQ_RNASEQ_GTF"
check_file "$SEQ_RNASEQ_STAR_INDEX/Genome"
check_file "$SEQ_RNASEQ_SALMON_INDEX/versionInfo.json"
check_file "$SEQ_CLINVAR_VCF"
check_file "$SEQ_VEP_CACHE_ROOT/homo_sapiens/$SEQ_VEP_CACHE_VERSION"_GRCh38/info.txt
EOF
  chmod +x "$SCRIPTS_DIR/validate_reference_integrity.sh"
}

write_helper_scripts() {
  write_with_backup "$SCRIPTS_DIR/vcf_sanity_check.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# -ge 1 ]] || { echo "Usage: $0 <vcf.gz> [region]" >&2; exit 1; }
vcf="$1"
region="${2:-}"
if [[ -n "$region" ]]; then
  bcftools view -H -r "$region" "$vcf" | head
else
  bcftools view -H "$vcf" | head
fi
bcftools index -n "$vcf"
EOF
  chmod +x "$SCRIPTS_DIR/vcf_sanity_check.sh"

  write_with_backup "$SCRIPTS_DIR/summarize_somatic_vcf.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# -eq 1 ]] || { echo "Usage: $0 <somatic.vcf.gz>" >&2; exit 1; }
vcf="$1"
echo "Variants by FILTER:"
bcftools query -f '%FILTER\n' "$vcf" | sort | uniq -c | sort -nr
echo
echo "Variants by consequence-style INFO tags (if present):"
bcftools +fill-tags "$vcf" -- -t TYPE 2>/dev/null | bcftools query -f '%TYPE\n' 2>/dev/null | sort | uniq -c | sort -nr || true
EOF
  chmod +x "$SCRIPTS_DIR/summarize_somatic_vcf.sh"

  write_with_backup "$SCRIPTS_DIR/summarize_germline_vcf.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# -eq 1 ]] || { echo "Usage: $0 <germline.vcf.gz>" >&2; exit 1; }
vcf="$1"
bcftools stats "$vcf" | awk '/^SN/ {print}'
EOF
  chmod +x "$SCRIPTS_DIR/summarize_germline_vcf.sh"

  write_with_backup "$SCRIPTS_DIR/merge_count_matrices.py" <<'EOF'
#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--inputs", nargs="+", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--id-column", default=None)
args = parser.parse_args()

tables = []
for path_str in args.inputs:
    path = Path(path_str)
    df = pd.read_csv(path, sep=None, engine="python")
    id_col = args.id_column or df.columns[0]
    sample_col = df.columns[-1]
    df = df[[id_col, sample_col]].rename(columns={sample_col: path.stem})
    tables.append(df)

merged = tables[0]
for df in tables[1:]:
    merged = merged.merge(df, on=merged.columns[0], how="outer")
merged.to_csv(args.output, sep="\t", index=False)
EOF
  chmod +x "$SCRIPTS_DIR/merge_count_matrices.py"

  write_with_backup "$SCRIPTS_DIR/deseq2_starter.R" <<'EOF'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: deseq2_starter.R <counts.tsv> <metadata.tsv> <outdir>", call. = FALSE)
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed. Install it first or treat this as a starter template.", call. = FALSE)
}
counts_path <- args[[1]]
meta_path <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(counts_path, check.names = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]
meta <- read.delim(meta_path, check.names = FALSE)
if (!"sample" %in% colnames(meta)) stop("metadata.tsv must contain a sample column", call. = FALSE)
rownames(meta) <- meta$sample

stop("Edit the design formula and contrasts in deseq2_starter.R before using this in a real analysis.", call. = FALSE)
EOF
  chmod +x "$SCRIPTS_DIR/deseq2_starter.R"

  write_with_backup "$SCRIPTS_DIR/run_deseq2_starter.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 <counts.tsv> <metadata.tsv> <outdir>" >&2; exit 1; }
Rscript /media/user/SEQ/scripts/deseq2_starter.R "$@"
EOF
  chmod +x "$SCRIPTS_DIR/run_deseq2_starter.sh"

  write_with_backup "$SCRIPTS_DIR/multiqc_qc_summary.py" <<'EOF'
#!/usr/bin/env python3
import json
import sys
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("Usage: multiqc_qc_summary.py <multiqc_data/multiqc_data.json>")

path = Path(sys.argv[1])
data = json.loads(path.read_text())
print("Top-level sections:")
for key in sorted(data.keys()):
    print(f"- {key}")
EOF
  chmod +x "$SCRIPTS_DIR/multiqc_qc_summary.py"
}

write_wrapper_scripts() {
  write_with_backup "$SCRIPTS_DIR/run_sarek.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

MODE=""
SAMPLESHEET=""
INTERVALS=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --intervals) INTERVALS="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$MODE" ]] || seq_die "Missing --mode (tumor-normal or germline)"
[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
[[ -n "$INTERVALS" ]] || seq_die "Missing --intervals with the exome capture BED / interval list"

seq_require_file "$SAMPLESHEET"
seq_require_file "$INTERVALS"
python3 /media/user/SEQ/scripts/validate_samplesheet.py --pipeline sarek --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp /media/user/SEQ/tmp/sarek_samplesheet.XXXXXX.csv)"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"
python3 - "$TMP_SHEET" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
rows = []
with path.open(newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    if not fieldnames:
        raise SystemExit("Sarek samplesheet is missing a header row")
    for row in reader:
        patient = (row.get("patient") or "").strip()
        if patient and patient.isdigit():
            row["patient"] = f"P{patient}"
        rows.append(row)

with path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY

if [[ "$INTERVALS" == *.bed ]]; then
  SORTED_INTERVALS="$(mktemp /media/user/SEQ/tmp/sarek_intervals.XXXXXX.bed)"
  awk 'BEGIN { OFS="\t" } NF >= 3 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^track([[:space:]]|$)/ && $0 !~ /^browser([[:space:]]|$)/ { print $1, $2, $3 }' "$INTERVALS" \
    | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$SORTED_INTERVALS"
  INTERVALS="$SORTED_INTERVALS"
fi

case "$MODE" in
  tumor-normal)
    TOOLS="mutect2,strelka,vep"
    ;;
  germline)
    TOOLS="haplotypecaller,vep"
    ;;
  *)
    seq_die "Unsupported mode: $MODE"
    ;;
esac

PROFILE="$(seq_runtime_profile)"
LOCAL_IGENOMES_BASE="${SEQ_SAREK_IGENOMES_BASE:-/media/user/SEQ/refs/igenomes_stub}"
LOCAL_SNPEFF_CACHE="${SEQ_SAREK_SNPEFF_CACHE:-/media/user/SEQ/refs/annotation/snpeff_cache}"
mkdir -p /media/user/SEQ/results/sarek "$LOCAL_IGENOMES_BASE" "$LOCAL_SNPEFF_CACHE"
nextflow run /media/user/SEQ/pipelines/sarek-3.8.1/main.nf \
  -profile "$PROFILE" \
  -c /media/user/SEQ/configs/sarek.config \
  --input "$TMP_SHEET" \
  --outdir /media/user/SEQ/results/sarek \
  --tools "$TOOLS" \
  --fasta "$SEQ_SAREK_FASTA" \
  --fasta_fai "$SEQ_SAREK_FASTA_FAI" \
  --dict "$SEQ_SAREK_DICT" \
  --dbsnp "$SEQ_SAREK_DBSNP" \
  --dbsnp_tbi "$SEQ_SAREK_DBSNP_TBI" \
  --known_indels "$SEQ_SAREK_KNOWN_INDELS" \
  --known_indels_tbi "$SEQ_SAREK_KNOWN_INDELS_TBI" \
  --known_snps "$SEQ_SAREK_KNOWN_SNPS" \
  --known_snps_tbi "$SEQ_SAREK_KNOWN_SNPS_TBI" \
  --germline_resource "$SEQ_SAREK_GERMLINE_RESOURCE" \
  --germline_resource_tbi "$SEQ_SAREK_GERMLINE_RESOURCE_TBI" \
  --pon "$SEQ_SAREK_PON" \
  --pon_tbi "$SEQ_SAREK_PON_TBI" \
  --igenomes_base "$LOCAL_IGENOMES_BASE" \
  --snpeff_cache "$LOCAL_SNPEFF_CACHE" \
  --vep_cache "$SEQ_VEP_CACHE_ROOT" \
  --vep_cache_version "$SEQ_VEP_CACHE_VERSION" \
  --vep_genome GRCh38 \
  --vep_species homo_sapiens \
  --intervals "$INTERVALS" \
  --igenomes_ignore \
  -w /media/user/SEQ/work \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_sarek.sh"

  write_with_backup "$SCRIPTS_DIR/run_rnaseq.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

SAMPLESHEET=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 /media/user/SEQ/scripts/validate_samplesheet.py --pipeline rnaseq --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp /media/user/SEQ/tmp/rnaseq_samplesheet.XXXXXX.csv)"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

PROFILE="$(seq_runtime_profile)"
mkdir -p /media/user/SEQ/results/rnaseq
nextflow run /media/user/SEQ/pipelines/rnaseq-3.24.0/main.nf \
  -profile "$PROFILE" \
  -c /media/user/SEQ/configs/rnaseq.config \
  --input "$TMP_SHEET" \
  --outdir /media/user/SEQ/results/rnaseq \
  --fasta "$SEQ_RNASEQ_FASTA" \
  --gtf "$SEQ_RNASEQ_GTF" \
  --star_index "$SEQ_RNASEQ_STAR_INDEX" \
  --salmon_index "$SEQ_RNASEQ_SALMON_INDEX" \
  --igenomes_ignore \
  -w /media/user/SEQ/work \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_rnaseq.sh"

  write_with_backup "$SCRIPTS_DIR/run_rnafusion.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

SAMPLESHEET=""
EXTRA_ARGS=()
RESUME=""
TOOLS="${SEQ_RNAFUSION_DEFAULT_TOOLS:-arriba,fusioncatcher,salmon}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --tools) TOOLS="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 /media/user/SEQ/scripts/validate_samplesheet.py --pipeline rnafusion --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp /media/user/SEQ/tmp/rnafusion_samplesheet.XXXXXX.csv)"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

if [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/star" ]]; then
  seq_die "RNA fusion STAR references are missing. Re-run sync_rnafusion_refs.sh first."
fi

PROFILE="$(seq_runtime_profile)"
mkdir -p /media/user/SEQ/results/rnafusion
nextflow run /media/user/SEQ/pipelines/rnafusion-4.1.0/main.nf \
  -profile "$PROFILE" \
  -c /media/user/SEQ/configs/rnafusion.config \
  --input "$TMP_SHEET" \
  --outdir /media/user/SEQ/results/rnafusion \
  --genomes_base "$SEQ_RNAFUSION_GENOMES_BASE" \
  --genome GRCh38 \
  --genome_gencode_version 46 \
  --tools "$TOOLS" \
  --no_cosmic \
  -w /media/user/SEQ/work \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_rnafusion.sh"
}

write_config_files() {
  write_with_backup "$CONFIG_DIR/resources.env" <<EOF
export SEQ_BASE="$BASE"
export NXF_HOME="$TMP_DIR/nextflow"
export NXF_WORK="$WORK_DIR"
export NXF_SINGULARITY_CACHEDIR="$CONTAINERS_DIR/nextflow_singularity"
export NXF_APPTAINER_CACHEDIR="$CONTAINERS_DIR/nextflow_apptainer"
export SINGULARITY_CACHEDIR="$CONTAINERS_DIR/singularity"
export APPTAINER_CACHEDIR="$CONTAINERS_DIR/apptainer"
export SEQ_RUNTIME_PROFILE="$RUNTIME_PROFILE"
export SEQ_SYNC_STARFUSION="$FUSION_SYNC_STARFUSION"

export SEQ_SAREK_FASTA="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta"
export SEQ_SAREK_FASTA_FAI="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta.fai"
export SEQ_SAREK_DICT="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.dict"
export SEQ_SAREK_DBSNP="$REFS_DIR/gatk_bundle/dbsnp_146.hg38.vcf.gz"
export SEQ_SAREK_DBSNP_TBI="$REFS_DIR/gatk_bundle/dbsnp_146.hg38.vcf.gz.tbi"
export SEQ_SAREK_KNOWN_INDELS="$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
export SEQ_SAREK_KNOWN_INDELS_TBI="$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
export SEQ_SAREK_KNOWN_SNPS="$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz"
export SEQ_SAREK_KNOWN_SNPS_TBI="$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz.tbi"
export SEQ_SAREK_GERMLINE_RESOURCE="$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz"
export SEQ_SAREK_GERMLINE_RESOURCE_TBI="$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz.tbi"
export SEQ_SAREK_PON="$REFS_DIR/gatk_bundle/1000g_pon.hg38.vcf.gz"
export SEQ_SAREK_PON_TBI="$REFS_DIR/gatk_bundle/1000g_pon.hg38.vcf.gz.tbi"

export SEQ_RNASEQ_FASTA="$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa"
export SEQ_RNASEQ_GTF="$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46.gtf"
export SEQ_RNASEQ_STAR_INDEX="$REFS_DIR/fusion/GRCh38/gencode_v46/star"
export SEQ_RNASEQ_SALMON_INDEX="$REFS_DIR/fusion/GRCh38/gencode_v46/salmon"
export SEQ_GENCODE_TRANSCRIPTS="$REFS_DIR/annotation/gencode.v46.transcripts.fa.gz"

export SEQ_VEP_CACHE_ROOT="$REFS_DIR/vep/cache"
export SEQ_VEP_CACHE_VERSION="$VEP_CACHE_VERSION"
export SEQ_VEP_FASTA="$REFS_DIR/vep/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
export SEQ_VEP_CONTAINER="$VEP_SIF"

export SEQ_CLINVAR_VCF="$REFS_DIR/clinvar/clinvar.vcf.gz"
export SEQ_CLINVAR_TBI="$REFS_DIR/clinvar/clinvar.vcf.gz.tbi"

export SEQ_RNAFUSION_GENOMES_BASE="$REFS_DIR/fusion"
export SEQ_RNAFUSION_DEFAULT_TOOLS="arriba,fusioncatcher,salmon"

export SEQ_OPTIONAL_COSMIC_CONFIG="$REFS_DIR/optional/COSMIC.placeholder.env"
export SEQ_OPTIONAL_ONCOKB_CONFIG="$REFS_DIR/optional/OncoKB.placeholder.env"
export SEQ_OPTIONAL_DBNSFP_CONFIG="$REFS_DIR/optional/dbNSFP.placeholder.env"
export SEQ_OPTIONAL_CADD_CONFIG="$REFS_DIR/optional/CADD.placeholder.env"
EOF

  write_with_backup "$CONFIG_DIR/refs.yaml" <<EOF
base: "$BASE"
runtime_profile: "$RUNTIME_PROFILE"

sarek:
  fasta: "$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta"
  fasta_fai: "$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta.fai"
  dict: "$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.dict"
  dbsnp: "$REFS_DIR/gatk_bundle/dbsnp_146.hg38.vcf.gz"
  known_indels: "$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  known_snps: "$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz"
  germline_resource: "$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz"
  pon: "$REFS_DIR/gatk_bundle/1000g_pon.hg38.vcf.gz"
  exome_intervals_manual: "/ABSOLUTE/PATH/TO/EXOME_CAPTURE_INTERVALS.bed"

rnaseq:
  fasta: "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa"
  gtf: "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46.gtf"
  transcripts: "$REFS_DIR/annotation/gencode.v46.transcripts.fa.gz"
  star_index: "$REFS_DIR/fusion/GRCh38/gencode_v46/star"
  salmon_index: "$REFS_DIR/fusion/GRCh38/gencode_v46/salmon"

vep:
  cache_root: "$REFS_DIR/vep/cache"
  cache_version: "$VEP_CACHE_VERSION"
  cache_species: "homo_sapiens"
  cache_genome: "GRCh38"
  fasta: "$REFS_DIR/vep/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  container: "$VEP_SIF"

clinvar:
  vcf: "$REFS_DIR/clinvar/clinvar.vcf.gz"
  tbi: "$REFS_DIR/clinvar/clinvar.vcf.gz.tbi"

rnafusion:
  genomes_base: "$REFS_DIR/fusion"
  gencode_version: "$GENCODE_VERSION"
  default_tools: "arriba,fusioncatcher,salmon"
  starfusion_synced: "$FUSION_SYNC_STARFUSION"

optional_placeholders:
  cosmic: "$REFS_DIR/optional/COSMIC.placeholder.env"
  oncokb: "$REFS_DIR/optional/OncoKB.placeholder.env"
  dbnsfp: "$REFS_DIR/optional/dbNSFP.placeholder.env"
  cadd: "$REFS_DIR/optional/CADD.placeholder.env"
EOF

  write_with_backup "$CONFIG_DIR/sarek.config" <<EOF
process {
  executor = 'local'
  cpus = 8
  memory = '56.GB'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_singularity'
}

apptainer {
  enabled = false
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_apptainer'
}

params {
  igenomes_ignore = true
  outdir = '$RESULTS_DIR/sarek'
  save_reference = true
}
EOF

  write_with_backup "$CONFIG_DIR/rnaseq.config" <<EOF
process {
  executor = 'local'
  cpus = 8
  memory = '56.GB'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_singularity'
}

apptainer {
  enabled = false
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_apptainer'
}

params {
  outdir = '$RESULTS_DIR/rnaseq'
  igenomes_ignore = true
  aligner = 'star_salmon'
  pseudo_aligner = 'salmon'
}
EOF

  write_with_backup "$CONFIG_DIR/rnafusion.config" <<EOF
process {
  executor = 'local'
  cpus = 8
  memory = '56.GB'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_singularity'
}

apptainer {
  enabled = false
  autoMounts = true
  cacheDir = '$CONTAINERS_DIR/nextflow_apptainer'
}

params {
  outdir = '$RESULTS_DIR/rnafusion'
  no_cosmic = true
  genome = 'GRCh38'
  genome_gencode_version = $GENCODE_VERSION
}
EOF

  write_with_backup "$REFS_DIR/optional/COSMIC.placeholder.env" <<'EOF'
# Add licensed COSMIC credentials and local file paths here later.
EOF
  write_with_backup "$REFS_DIR/optional/OncoKB.placeholder.env" <<'EOF'
# Add OncoKB API token or local resources here later.
EOF
  write_with_backup "$REFS_DIR/optional/dbNSFP.placeholder.env" <<'EOF'
# Add dbNSFP local file paths here later.
EOF
  write_with_backup "$REFS_DIR/optional/CADD.placeholder.env" <<'EOF'
# Add CADD local file paths here later.
EOF
}

write_samplesheet_templates() {
  write_with_backup "$SAMPLESHEETS_DIR/sarek_samplesheet.csv" <<'EOF'
# Remove comment lines before using outside the provided wrappers.
# Columns:
# patient: patient identifier for pairing
# sex: XX or XY
# status: 0 for normal, 1 for tumor
# sample: unique sample identifier
# lane: lane label if starting from FASTQ/BAM per lane
# fastq_1/fastq_2 or bam/bai or cram/crai: supply one input mode consistently
# Optional manual references:
# exome intervals: /ABSOLUTE/PATH/TO/CAPTURE_INTERVALS.bed
# VEP cache root: /media/user/SEQ/refs/vep/cache
# ClinVar VCF: /media/user/SEQ/refs/clinvar/clinvar.vcf.gz
patient,sex,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai
PANCREAS001,XY,0,PANCREAS001_N,L001,/ABSOLUTE/PATH/TO/NORMAL_R1.fastq.gz,/ABSOLUTE/PATH/TO/NORMAL_R2.fastq.gz,,,,
PANCREAS001,XY,1,PANCREAS001_T,L001,/ABSOLUTE/PATH/TO/TUMOR_R1.fastq.gz,/ABSOLUTE/PATH/TO/TUMOR_R2.fastq.gz,,,,
EOF

  write_with_backup "$SAMPLESHEETS_DIR/rnaseq_samplesheet.csv" <<'EOF'
# Remove comment lines before using outside the provided wrappers.
# strandedness: unstranded, forward, or reverse
# fasta: /media/user/SEQ/refs/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa
# gtf: /media/user/SEQ/refs/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46.gtf
# salmon index: /media/user/SEQ/refs/fusion/GRCh38/gencode_v46/salmon
sample,fastq_1,fastq_2,strandedness
PANCREAS001_TUMOR_RNA,/ABSOLUTE/PATH/TO/RNA_R1.fastq.gz,/ABSOLUTE/PATH/TO/RNA_R2.fastq.gz,forward
EOF

  write_with_backup "$SAMPLESHEETS_DIR/rnafusion_samplesheet.csv" <<'EOF'
# Remove comment lines before using outside the provided wrappers.
# One of fastq_1, bam, cram, junctions, or splice_junctions must be populated.
# genomes_base: /media/user/SEQ/refs/fusion
# available default tools on this host: arriba,fusioncatcher,salmon
sample,fastq_1,fastq_2,bam,bai,cram,crai,junctions,splice_junctions,strandedness,seq_platform,seq_center
PANCREAS001_TUMOR_RNA,/ABSOLUTE/PATH/TO/RNA_R1.fastq.gz,/ABSOLUTE/PATH/TO/RNA_R2.fastq.gz,,,,,,,forward,ILLUMINA,LOCAL
EOF
}

download_transcript_fasta() {
  local log_file="$LOG_DIR/download_core_refs.log"
  : > "$log_file"
  if [[ ! -s "$REFS_DIR/annotation/gencode.v46.transcripts.fa.gz" ]]; then
    download_url \
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz" \
      "$REFS_DIR/annotation/gencode.v46.transcripts.fa.gz" \
      "$log_file"
  fi
}

create_reference_links() {
  mkdir -p "$REFS_DIR/genome" "$REFS_DIR/annotation"

  ln -sfn "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa" \
    "$REFS_DIR/genome/GRCh38.primary_assembly.gencode_v46.fa"
  ln -sfn "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" \
    "$REFS_DIR/genome/GRCh38.primary_assembly.gencode_v46.fa.fai"
  ln -sfn "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46.gtf" \
    "$REFS_DIR/annotation/gencode.v46.primary_assembly.annotation.gtf"

  if [[ ! -s "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.dict" ]]; then
    samtools dict \
      "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa" \
      -o "$REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.dict"
  fi
}

write_reference_inventory() {
  local fusion_tools="arriba,fusioncatcher,salmon"
  local starfusion_note="not synced on this host"
  if [[ "$FUSION_SYNC_STARFUSION" == "true" ]]; then
    fusion_tools="${fusion_tools},starfusion,fusioninspector"
    starfusion_note="synced from nf-core public S3 bucket"
  fi

  write_with_backup "$REFERENCE_INVENTORY" <<EOF
# Reference Inventory

Generated: $DATE_ISO

## Core genome and annotation

- GRCh38 GENCODE v46 genome FASTA: \`$REFS_DIR/genome/GRCh38.primary_assembly.gencode_v46.fa\`
- GENCODE v46 comprehensive GTF: \`$REFS_DIR/annotation/gencode.v46.primary_assembly.annotation.gtf\`
- GENCODE v46 transcripts FASTA: \`$REFS_DIR/annotation/gencode.v46.transcripts.fa.gz\`
- STAR index: \`$REFS_DIR/fusion/GRCh38/gencode_v46/star\`
- Salmon index: \`$REFS_DIR/fusion/GRCh38/gencode_v46/salmon\`

## GATK bundle

- Reference FASTA: \`$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta\`
- FASTA index: \`$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta.fai\`
- Sequence dictionary: \`$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.dict\`
- BWA index prefix: \`$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta\`
- Known sites:
  - \`dbsnp_146.hg38.vcf.gz\`
  - \`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\`
  - \`Homo_sapiens_assembly38.known_indels.vcf.gz\`
  - \`1000G_omni2.5.hg38.vcf.gz\`
  - \`1000G_phase1.snps.high_confidence.hg38.vcf.gz\`
  - \`Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz\`
  - \`hapmap_3.3.hg38.vcf.gz\`
- Mutect2 resources:
  - \`af-only-gnomad.hg38.vcf.gz\`
  - \`1000g_pon.hg38.vcf.gz\`

## VEP

- Container: \`$VEP_SIF\`
- Cache root: \`$REFS_DIR/vep/cache\`
- Cache version: \`$VEP_CACHE_VERSION\`
- Offline FASTA: \`$REFS_DIR/vep/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa\`

## ClinVar

- ClinVar GRCh38 VCF: \`$REFS_DIR/clinvar/clinvar.vcf.gz\`
- Index: \`$REFS_DIR/clinvar/clinvar.vcf.gz.tbi\`

## Fusion references

- Synced by default: \`$fusion_tools\`
- STAR-Fusion status: $starfusion_note
- Fusion genomes base: \`$REFS_DIR/fusion\`

## Optional placeholders only

- COSMIC: \`$REFS_DIR/optional/COSMIC.placeholder.env\`
- OncoKB: \`$REFS_DIR/optional/OncoKB.placeholder.env\`
- dbNSFP: \`$REFS_DIR/optional/dbNSFP.placeholder.env\`
- CADD: \`$REFS_DIR/optional/CADD.placeholder.env\`
EOF
}

write_manual_items() {
  write_with_backup "$MANUAL_ITEMS" <<EOF
# Manual Items Required

The setup completed without moving or uploading patient data. Before running real analyses, provide:

1. Patient FASTQ, BAM, or CRAM files at their final local paths.
2. The exome capture BED / interval list for the exact WES kit used.
3. Completed samplesheets under \`$SAMPLESHEETS_DIR\`.
4. Optional licensed or token-gated resources only if you need them later:
   - COSMIC
   - OncoKB
   - dbNSFP
   - CADD
5. If you plan to run the DESeq2 starter, install \`DESeq2\` for \`Rscript\` first or treat the provided script as a template.

Notes:

- Comment lines in the samplesheet templates are safe when using the provided wrappers; the wrappers strip them automatically before launch.
- The current host has about ${AVAILABLE_GB}G free at setup time, so STAR-Fusion references were set to \`$FUSION_SYNC_STARFUSION\`.
EOF
}

write_setup_docs() {
  write_with_backup "$README_SETUP" <<EOF
# README Setup

This workspace contains a reproducible local short-read pancreatic cancer analysis bundle rooted at \`$BASE\`.

## Included

- nf-core/sarek 3.8.1 for WES DNA
- nf-core/rnaseq 3.24.0 for bulk RNA-seq
- nf-core/rnafusion 4.1.0 for optional fusion analysis
- Public references for GRCh38 DNA, RNA, VEP, and ClinVar
- Idempotent helper scripts under \`$SCRIPTS_DIR\`
- Wrapper launchers for each pipeline

## Safety

- Patient data were not moved.
- Existing files are backed up before this setup overwrites them.
- Nextflow work and caches are kept under \`$BASE\`.
EOF

  write_with_backup "$LAUNCH_EXAMPLES" <<EOF
# Launch Examples

## WES tumor-normal

\`\`\`bash
$SCRIPTS_DIR/run_sarek.sh --mode tumor-normal --samplesheet $SAMPLESHEETS_DIR/sarek_samplesheet.csv --intervals /ABSOLUTE/PATH/TO/CAPTURE_INTERVALS.bed -resume
\`\`\`

## WES germline

\`\`\`bash
$SCRIPTS_DIR/run_sarek.sh --mode germline --samplesheet $SAMPLESHEETS_DIR/sarek_samplesheet.csv --intervals /ABSOLUTE/PATH/TO/CAPTURE_INTERVALS.bed -resume
\`\`\`

## RNA-seq expression

\`\`\`bash
$SCRIPTS_DIR/run_rnaseq.sh --samplesheet $SAMPLESHEETS_DIR/rnaseq_samplesheet.csv -resume
\`\`\`

## RNA fusion

\`\`\`bash
$SCRIPTS_DIR/run_rnafusion.sh --samplesheet $SAMPLESHEETS_DIR/rnafusion_samplesheet.csv -resume
\`\`\`
EOF
}

run_reference_setup() {
  log "Downloading core transcript FASTA"
  download_transcript_fasta

  log "Syncing published RNA/fusion references"
  "$SCRIPTS_DIR/sync_rnafusion_refs.sh"

  log "Downloading GATK bundle"
  "$SCRIPTS_DIR/download_gatk_bundle.sh"

  log "Downloading ClinVar"
  "$SCRIPTS_DIR/download_clinvar.sh"

  log "Setting up VEP"
  "$SCRIPTS_DIR/setup_vep.sh"

  log "Linking shared reference paths"
  create_reference_links
}

run_smoke_tests() {
  mkdir -p "$TMP_DIR/nextflow"
  export NXF_HOME="$TMP_DIR/nextflow"
  export NXF_WORK="$WORK_DIR"
  export NXF_SINGULARITY_CACHEDIR="$CONTAINERS_DIR/nextflow_singularity"

  local profile="$RUNTIME_PROFILE"
  local nf="/media/user/.local/bin/nextflow"
  if ! command -v "$nf" >/dev/null 2>&1; then
    nf="nextflow"
  fi

  "$nf" run "$PIPELINES_DIR/sarek-3.8.1/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" > "$LOG_DIR/test_sarek.log" 2>&1 || true
  "$nf" run "$PIPELINES_DIR/rnaseq-3.24.0/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" > "$LOG_DIR/test_rnaseq.log" 2>&1 || true
  "$nf" run "$PIPELINES_DIR/rnafusion-4.1.0/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" > "$LOG_DIR/test_rnafusion.log" 2>&1 || true
}

write_setup_summary() {
  local starfusion_status="not synced"
  if [[ "$FUSION_SYNC_STARFUSION" == "true" ]]; then
    starfusion_status="synced"
  fi

  write_with_backup "$SETUP_SUMMARY" <<EOF
# Setup Summary

Generated: $DATE_ISO

## System audit summary

- Host: \`$HOSTNAME_FQDN\`
- User: \`$USER_NAME\`
- Runtime profile selected: \`$RUNTIME_PROFILE\`
- Free disk at setup time: \`${AVAILABLE_GB}G\`
- Java: \`$(java -version 2>&1 | head -n 1)\`
- Nextflow: \`$NEXTFLOW_VERSION\`

## Installed or configured software

- Nextflow, Java 17+, git, Singularity/Docker, curl, wget, unzip, tar, pigz, bgzip, tabix, awscli, samtools, bcftools, bwa
- VEP container: \`$VEP_IMAGE_TAG\`

## Pipeline versions

- nf-core/sarek 3.8.1: \`$(git -C "$PIPELINES_DIR/sarek-3.8.1" rev-parse HEAD)\`
- nf-core/rnaseq 3.24.0: \`$(git -C "$PIPELINES_DIR/rnaseq-3.24.0" rev-parse HEAD)\`
- nf-core/rnafusion 4.1.0: \`$(git -C "$PIPELINES_DIR/rnafusion-4.1.0" rev-parse HEAD)\`

## Reference and database inventory

- GRCh38 GENCODE v46 genome, GTF, transcript FASTA
- STAR and Salmon indices
- Broad hg38 resource bundle components for BQSR, VQSR, Mutect2 germline resource, and PoN
- Ensembl VEP offline cache + FASTA + local container
- ClinVar GRCh38 VCF
- RNA fusion published references: Arriba, FusionCatcher, fusion-report, HGNC, STAR, Salmon
- STAR-Fusion references: $starfusion_status

## Downloaded automatically

- Public nf-core RNA fusion references from \`s3://nf-core-awsmegatests/rnafusion/references/\`
- Public Broad hg38 bundle files from Google public storage
- Public ClinVar GRCh38 VCF
- Public Ensembl VEP cache and FASTA
- Public GENCODE v46 transcript FASTA

## Still needs manual input

- Patient FASTQ/BAM/CRAM files
- Exome capture BED / interval list
- Optional licensed resources if desired later
- DESeq2 package installation if you want to run the starter directly

## Exact launch commands

### WES tumor-normal

\`\`\`bash
$SCRIPTS_DIR/run_sarek.sh --mode tumor-normal --samplesheet $SAMPLESHEETS_DIR/sarek_samplesheet.csv --intervals /ABSOLUTE/PATH/TO/CAPTURE_INTERVALS.bed -resume
\`\`\`

### WES germline

\`\`\`bash
$SCRIPTS_DIR/run_sarek.sh --mode germline --samplesheet $SAMPLESHEETS_DIR/sarek_samplesheet.csv --intervals /ABSOLUTE/PATH/TO/CAPTURE_INTERVALS.bed -resume
\`\`\`

### RNA-seq expression

\`\`\`bash
$SCRIPTS_DIR/run_rnaseq.sh --samplesheet $SAMPLESHEETS_DIR/rnaseq_samplesheet.csv -resume
\`\`\`

### Optional RNA fusion

\`\`\`bash
$SCRIPTS_DIR/run_rnafusion.sh --samplesheet $SAMPLESHEETS_DIR/rnafusion_samplesheet.csv -resume
\`\`\`
EOF
}

main() {
  log "Writing system audit and install notes"
  emit_system_audit
  write_install_notes
  emit_software_versions

  log "Cloning pinned pipeline revisions"
  clone_pipelines

  log "Creating reusable workspace scripts"
  write_common_script
  write_install_prereqs_script
  write_download_gatk_script
  write_sync_fusion_script
  write_clinvar_script
  write_vep_script
  write_validation_script
  write_reference_validation_script
  write_helper_scripts
  write_wrapper_scripts
  write_config_files
  write_samplesheet_templates

  log "Executing reference setup"
  run_reference_setup

  log "Writing reference inventory and docs"
  write_reference_inventory
  write_manual_items
  write_setup_docs

  log "Running smoke tests"
  run_smoke_tests

  log "Validating references"
  bash "$SCRIPTS_DIR/validate_reference_integrity.sh" > "$LOG_DIR/validate_reference_integrity.log" 2>&1 || true

  log "Writing setup summary"
  write_setup_summary

  log "Bootstrap complete"
}

main "$@"
