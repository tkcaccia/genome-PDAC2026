#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  prepare_pancreatic_shortread_workspace.sh [options]

Core path options:
  --install-root PATH           Root directory for refs/configs/docs/results. Default: /media/user/SEQ
  --work-dir PATH               Nextflow work directory. Default: <install-root>/work
  --results-dir PATH            Analysis results root. Default: <install-root>/results
  --container-cache-dir PATH    Nextflow Singularity/Apptainer cache directory. Default: <install-root>/containers/nextflow_singularity
  --apptainer-cache-dir PATH    Optional explicit Apptainer cache directory. Default: <install-root>/containers/nextflow_apptainer
  --vep-sif-path PATH           Exact path for the Ensembl VEP .sif image. Default: <install-root>/containers/ensembl-vep_release_115.2.sif

Behavior options:
  --runtime-profile NAME        Force runtime profile: apptainer, singularity, or docker
  --sync-starfusion             Force STAR-Fusion reference sync
  --skip-starfusion             Skip STAR-Fusion reference sync
  --sync-fusioncatcher          Force FusionCatcher reference sync
  --skip-fusioncatcher          Skip FusionCatcher reference sync
  --smoke-tests                 Run Nextflow smoke tests after setup
  --no-smoke-tests              Skip smoke tests
  --dry-run                     Print the resolved configuration and exit
  --help                        Show this help text

Environment variable equivalents:
  INSTALL_ROOT
  WORK_DIR_OVERRIDE
  RESULTS_DIR_OVERRIDE
  CONTAINER_CACHE_DIR
  APPTAINER_CACHE_DIR
  VEP_SIF_PATH
  NEXTFLOW_PINNED_VERSION
  RUNTIME_PROFILE
  SYNC_STARFUSION
  SYNC_FUSIONCATCHER
  RUN_SMOKE_TESTS
EOF
}

INSTALL_ROOT="${INSTALL_ROOT:-/media/user/SEQ}"
WORK_DIR_OVERRIDE="${WORK_DIR_OVERRIDE:-}"
RESULTS_DIR_OVERRIDE="${RESULTS_DIR_OVERRIDE:-}"
CONTAINER_CACHE_DIR="${CONTAINER_CACHE_DIR:-}"
APPTAINER_CACHE_DIR="${APPTAINER_CACHE_DIR:-}"
VEP_SIF_PATH="${VEP_SIF_PATH:-}"
RUNTIME_PROFILE="${RUNTIME_PROFILE:-auto}"
SYNC_STARFUSION="${SYNC_STARFUSION:-auto}"
SYNC_FUSIONCATCHER="${SYNC_FUSIONCATCHER:-true}"
RUN_SMOKE_TESTS="${RUN_SMOKE_TESTS:-true}"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install-root) INSTALL_ROOT="$2"; shift 2 ;;
    --work-dir) WORK_DIR_OVERRIDE="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR_OVERRIDE="$2"; shift 2 ;;
    --container-cache-dir) CONTAINER_CACHE_DIR="$2"; shift 2 ;;
    --apptainer-cache-dir) APPTAINER_CACHE_DIR="$2"; shift 2 ;;
    --vep-sif-path) VEP_SIF_PATH="$2"; shift 2 ;;
    --runtime-profile) RUNTIME_PROFILE="$2"; shift 2 ;;
    --sync-starfusion) SYNC_STARFUSION="true"; shift ;;
    --skip-starfusion) SYNC_STARFUSION="false"; shift ;;
    --sync-fusioncatcher) SYNC_FUSIONCATCHER="true"; shift ;;
    --skip-fusioncatcher) SYNC_FUSIONCATCHER="false"; shift ;;
    --smoke-tests) RUN_SMOKE_TESTS="true"; shift ;;
    --no-smoke-tests) RUN_SMOKE_TESTS="false"; shift ;;
    --dry-run) DRY_RUN="true"; shift ;;
    --help|-h) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n\n' "$1" >&2; usage >&2; exit 1 ;;
  esac
done

CONTAINER_CACHE_DIR="${CONTAINER_CACHE_DIR:-$INSTALL_ROOT/containers/nextflow_singularity}"
APPTAINER_CACHE_DIR="${APPTAINER_CACHE_DIR:-$INSTALL_ROOT/containers/nextflow_apptainer}"
VEP_SIF_PATH="${VEP_SIF_PATH:-$INSTALL_ROOT/containers/ensembl-vep_release_115.2.sif}"

BASE="$INSTALL_ROOT"
BIN_DIR="$BASE/bin"
LOG_DIR="$BASE/logs"
DOCS_DIR="$BASE/docs"
PIPELINES_DIR="$BASE/pipelines"
CONFIG_DIR="$BASE/configs"
REFS_DIR="$BASE/refs"
CONTAINERS_DIR="$(dirname "$CONTAINER_CACHE_DIR")"
WORK_DIR="${WORK_DIR_OVERRIDE:-$BASE/work}"
RESULTS_DIR="${RESULTS_DIR_OVERRIDE:-$BASE/results}"
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
BOOTSTRAP_LOG="$LOG_DIR/prepare_workspace_$(date +%Y%m%d_%H%M%S).log"

GENCODE_VERSION="46"
GENCODE_LABEL="gencode_v${GENCODE_VERSION}"
VEP_CACHE_VERSION="115"
VEP_IMAGE_TAG="release_115.2"
VEP_IMAGE="docker://ensemblorg/ensembl-vep:${VEP_IMAGE_TAG}"
VEP_SIF="$VEP_SIF_PATH"
NEXTFLOW_PINNED_VERSION="${NEXTFLOW_PINNED_VERSION:-25.10.4}"

NEXTFLOW_VERSION="UNKNOWN"
HOSTNAME_FQDN="$(hostname)"
USER_NAME="$(whoami)"
DATE_ISO="$(date -Iseconds)"

if [[ "$RUNTIME_PROFILE" == "auto" ]]; then
  if command -v apptainer >/dev/null 2>&1; then
    RUNTIME_PROFILE="apptainer"
  elif command -v singularity >/dev/null 2>&1; then
    RUNTIME_PROFILE="singularity"
  elif command -v docker >/dev/null 2>&1; then
    RUNTIME_PROFILE="docker"
  else
    RUNTIME_PROFILE="missing"
  fi
fi

AVAILABLE_GB="$(df -Pk "$BASE" 2>/dev/null | awk 'NR==2 {printf "%d", $4/1048576}')"
if [[ "$SYNC_STARFUSION" == "auto" ]]; then
  if [[ "${AVAILABLE_GB:-0}" =~ ^[0-9]+$ ]] && (( AVAILABLE_GB >= 220 )); then
    SYNC_STARFUSION="true"
  else
    SYNC_STARFUSION="false"
  fi
fi

if [[ "$DRY_RUN" == "true" ]]; then
  cat <<EOF
install_root=$INSTALL_ROOT
work_dir=$WORK_DIR
results_dir=$RESULTS_DIR
container_cache_dir=$CONTAINER_CACHE_DIR
apptainer_cache_dir=$APPTAINER_CACHE_DIR
vep_sif_path=$VEP_SIF_PATH
runtime_profile=$RUNTIME_PROFILE
sync_starfusion=$SYNC_STARFUSION
sync_fusioncatcher=$SYNC_FUSIONCATCHER
run_smoke_tests=$RUN_SMOKE_TESTS
available_gb=$AVAILABLE_GB
EOF
  exit 0
fi

mkdir -p \
  "$BIN_DIR" "$LOG_DIR" "$DOCS_DIR" "$PIPELINES_DIR" "$CONFIG_DIR" \
  "$REFS_DIR/genome" "$REFS_DIR/annotation" "$REFS_DIR/gatk_bundle" "$REFS_DIR/vep" \
  "$REFS_DIR/clinvar" "$REFS_DIR/fusion" "$REFS_DIR/optional" \
  "$CONTAINERS_DIR" "$CONTAINER_CACHE_DIR" "$APPTAINER_CACHE_DIR" \
  "$WORK_DIR" "$RESULTS_DIR" "$SAMPLESHEETS_DIR" "$SCRIPTS_DIR" "$TMP_DIR" "$BACKUP_DIR" "$(dirname "$VEP_SIF")"

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

ensure_nextflow() {
  local user_bin="$HOME/.local/bin"
  local target="$user_bin/nextflow"
  export PATH="$user_bin:$BIN_DIR:$PATH"
  local current_path=""
  local current_version=""
  if command -v nextflow >/dev/null 2>&1; then
    current_path="$(command -v nextflow)"
    current_version="$(nextflow -version 2>/dev/null | awk '/version/{print $3; exit}' || true)"
    current_version="${current_version:-UNKNOWN}"
  fi
  if [[ -n "$current_version" && "$current_version" != "UNKNOWN" && "$current_version" == "$NEXTFLOW_PINNED_VERSION" ]]; then
    NEXTFLOW_VERSION="$current_version"
    return 0
  fi

  [[ -n "${JAVA_HOME:-}" ]] || command -v java >/dev/null 2>&1 || die "Java 17+ is required before bootstrapping Nextflow."

  local log_file="$LOG_DIR/install_nextflow.log"
  local installer="$TMP_DIR/get_nextflow.sh"
  download_url "https://get.nextflow.io" "$installer" "$log_file"
  chmod +x "$installer"
  mkdir -p "$user_bin"

  if [[ -n "$current_path" && -x "$current_path" ]]; then
    "$current_path" self-update >>"$log_file" 2>&1 || die "Nextflow self-update failed for $current_path"
    NEXTFLOW_VERSION="$("$current_path" -version 2>/dev/null | awk '/version/{print $3; exit}' || true)"
    NEXTFLOW_VERSION="${NEXTFLOW_VERSION:-UNKNOWN}"
    export PATH="$(dirname "$current_path"):$PATH"
    return 0
  fi

  if [[ -e "$target" || -L "$target" ]]; then
    backup_if_exists "$target"
    rm -f "$target"
  fi
  (
    cd "$user_bin"
    NXF_VER="$NEXTFLOW_PINNED_VERSION" bash "$installer"
  ) >>"$log_file" 2>&1
  [[ -f "$target" && -x "$target" ]] || die "Nextflow bootstrap did not produce $target"
  NEXTFLOW_VERSION="$("$target" -version 2>/dev/null | awk '/version/{print $3; exit}' || true)"
  NEXTFLOW_VERSION="${NEXTFLOW_VERSION:-UNKNOWN}"
}

ensure_htslib_tools() {
  export PATH="$BIN_DIR:$PATH"
  if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    return 0
  fi

  command -v gcc >/dev/null 2>&1 || command -v cc >/dev/null 2>&1 || die "bgzip/tabix are missing and no C compiler is available to build htslib locally."
  command -v make >/dev/null 2>&1 || die "bgzip/tabix are missing and make is not available to build htslib locally."
  command -v bunzip2 >/dev/null 2>&1 || command -v bzip2 >/dev/null 2>&1 || die "bgzip/tabix are missing and bunzip2/bzip2 is not available to unpack htslib."

  local ver="1.20"
  local archive="$TMP_DIR/htslib-${ver}.tar.bz2"
  local src_dir="$TMP_DIR/htslib-${ver}"
  local log_file="$LOG_DIR/install_htslib_tools.log"

  download_url "https://github.com/samtools/htslib/releases/download/${ver}/htslib-${ver}.tar.bz2" "$archive" "$log_file"
  if [[ ! -d "$src_dir" ]]; then
    tar -xjf "$archive" -C "$TMP_DIR" >>"$log_file" 2>&1
  fi

  (
    cd "$src_dir"
    make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" bgzip tabix
  ) >>"$log_file" 2>&1

  install -m 0755 "$src_dir/bgzip" "$BIN_DIR/bgzip"
  install -m 0755 "$src_dir/tabix" "$BIN_DIR/tabix"

  command -v bgzip >/dev/null 2>&1 || die "Local bgzip bootstrap did not succeed."
  command -v tabix >/dev/null 2>&1 || die "Local tabix bootstrap did not succeed."
}

preflight_requirements() {
  local missing=()
  local required=(java git curl tar pigz bgzip tabix aws samtools bcftools bwa python3)
  local cmd
  for cmd in "${required[@]}"; do
    command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
  done
  if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
    missing+=("apptainer|singularity")
  fi
  if [[ "${#missing[@]}" -gt 0 ]]; then
    die "Missing required commands for this setup: ${missing[*]}"
  fi
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

This setup was prepared under the chosen install root and used the tools already present on the host, with local bootstrap fallbacks under \`$BIN_DIR\` for Nextflow and htslib tools when needed.

## Audited and configured

- Java 17+: \`$(java -version 2>&1 | head -n 1)\`
- Nextflow: \`$(nextflow -version 2>/dev/null | awk '/version/{print $2, $3; exit}')\`
- Container runtime preference: \`$RUNTIME_PROFILE\`
- Install root: \`$BASE\`
- Work directory: \`$WORK_DIR\`
- Results root: \`$RESULTS_DIR\`
- Nextflow container cache: \`$CONTAINER_CACHE_DIR\`
- Apptainer cache: \`$APPTAINER_CACHE_DIR\`
- VEP SIF path: \`$VEP_SIF\`
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

SEQ_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SEQ_BASE="$(cd -- "$SEQ_SCRIPT_DIR/.." && pwd)"
SEQ_TMP="$SEQ_BASE/tmp"
SEQ_WORK="${SEQ_WORK:-$SEQ_BASE/work}"
SEQ_RESULTS="${SEQ_RESULTS:-$SEQ_BASE/results}"
SEQ_CONFIG="${SEQ_CONFIG:-$SEQ_BASE/configs/resources.env}"
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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

LOG_FILE="$SEQ_BASE/logs/install_prereqs.log"
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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

OUT="$SEQ_BASE/refs/gatk_bundle"
LOG_FILE="$SEQ_BASE/logs/download_gatk_bundle.log"
: > "$LOG_FILE"

download_gatk() {
  local object_path="$1"
  local dest_name="$2"
  local dest="$OUT/$dest_name"
  if [[ -s "$dest" ]]; then
    seq_log "Keeping existing $dest" | tee -a "$LOG_FILE"
    return 0
  fi
  local url="https://storage.googleapis.com/gcp-public-data--broad-references/${object_path}"
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "${dest}.part" "$url" >>"$LOG_FILE" 2>&1
  mv "${dest}.part" "$dest"
}

mkdir -p "$OUT"
resources=(
  "hg38/v0/Homo_sapiens_assembly38.fasta|Homo_sapiens_assembly38.fasta"
  "hg38/v0/Homo_sapiens_assembly38.dict|Homo_sapiens_assembly38.dict"
  "hg38/v0/Homo_sapiens_assembly38.fasta.fai|Homo_sapiens_assembly38.fasta.fai"
  "hg38/v0/gdc/dbsnp_144.hg38.vcf.gz|dbsnp_144.hg38.vcf.gz"
  "hg38/v0/gdc/dbsnp_144.hg38.vcf.gz.tbi|dbsnp_144.hg38.vcf.gz.tbi"
  "hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  "hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi|Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  "hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz|Homo_sapiens_assembly38.known_indels.vcf.gz"
  "hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi|Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  "hg38/v0/1000G_omni2.5.hg38.vcf.gz|1000G_omni2.5.hg38.vcf.gz"
  "hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi|1000G_omni2.5.hg38.vcf.gz.tbi"
  "hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz|1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  "hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi|1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  "hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz|Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  "hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi|Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
  "hg38/v0/hapmap_3.3.hg38.vcf.gz|hapmap_3.3.hg38.vcf.gz"
  "hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi|hapmap_3.3.hg38.vcf.gz.tbi"
  "hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz|af-only-gnomad.hg38.vcf.gz"
  "hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi|af-only-gnomad.hg38.vcf.gz.tbi"
)

for item in "${resources[@]}"; do
  IFS='|' read -r object_path dest_name <<<"$item"
  download_gatk "$object_path" "$dest_name"
done

if [[ ! -s "$OUT/Homo_sapiens_assembly38.fasta.fai" ]]; then
  samtools faidx "$OUT/Homo_sapiens_assembly38.fasta" >>"$LOG_FILE" 2>&1
fi
if [[ ! -s "$OUT/Homo_sapiens_assembly38.dict" ]]; then
  samtools dict "$OUT/Homo_sapiens_assembly38.fasta" -o "$OUT/Homo_sapiens_assembly38.dict" >>"$LOG_FILE" 2>&1
fi
if [[ ! -s "$OUT/Homo_sapiens_assembly38.fasta.amb" ]]; then
  bwa index "$OUT/Homo_sapiens_assembly38.fasta" >>"$LOG_FILE" 2>&1
fi

awk 'BEGIN{OFS="\t"} /^@SQ/ {sn=""; ln=""; for (i=1;i<=NF;i++) {if ($i ~ /^SN:/) sn=substr($i,4); if ($i ~ /^LN:/) ln=substr($i,4)} if (sn != "" && ln != "") print sn,0,ln}' "$OUT/Homo_sapiens_assembly38.dict" > "$OUT/Homo_sapiens_assembly38.genome.bed"

cat > "$OUT/PON_CONTROLLED_REQUIRED.md" <<'PON'
# Controlled Panel Of Normals Required For Optional Standard MuTect2 PoN Usage

Public Broad hg38 references were downloaded automatically, but the GDC MuTect2 Panel of Normals files are controlled resources and were not downloaded.

Official GDC reference files page:
https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

Relevant controlled resources named there include:
- `MuTect2.PON.4136.vcf.tar`
- `MuTect2.PON.5210.vcf.tar`
- `gatk4_mutect2_4136_pon.vcf.tar`

If you later obtain authorized access, unpack the chosen PoN VCF and index it here, then update `configs/resources.env` if needed.
PON
EOF
  chmod +x "$SCRIPTS_DIR/download_gatk_bundle.sh"
}

write_sync_fusion_script() {
  local starfusion_default="$SYNC_STARFUSION"
  local fusioncatcher_default="$SYNC_FUSIONCATCHER"
  write_with_backup "$SCRIPTS_DIR/sync_rnafusion_refs.sh" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="\$(cd -- "\$(dirname -- "\${BASH_SOURCE[0]}")" && pwd)"
source "\$SCRIPT_DIR/_common.sh"
seq_source_env

BASE_OUT="\$SEQ_BASE/refs/fusion/GRCh38"
LOG_FILE="\$SEQ_BASE/logs/sync_rnafusion_refs.log"
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
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/fusion_report_db" "\$BASE_OUT/fusion_report_db"
sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/hgnc" "\$BASE_OUT/hgnc"

if [[ "\${SEQ_SYNC_FUSIONCATCHER:-$fusioncatcher_default}" == "true" ]]; then
  sync_path "s3://nf-core-awsmegatests/rnafusion/references/GRCh38/gencode_v46/fusioncatcher" "\$BASE_OUT/gencode_v46/fusioncatcher"
else
  mkdir -p "\$BASE_OUT/gencode_v46/fusioncatcher"
  printf 'FusionCatcher references intentionally not synced on this host.\\n' > "\$BASE_OUT/gencode_v46/fusioncatcher/PLACEHOLDER.txt"
fi

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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

OUT="$SEQ_BASE/refs/clinvar"
LOG_FILE="$SEQ_BASE/logs/download_clinvar.log"
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

SCRIPT_DIR="\$(cd -- "\$(dirname -- "\${BASH_SOURCE[0]}")" && pwd)"
source "\$SCRIPT_DIR/_common.sh"
seq_source_env

OUT="\$SEQ_BASE/refs/vep"
LOG_FILE="\$SEQ_BASE/logs/setup_vep.log"
VEP_SIF="\${SEQ_VEP_CONTAINER:-$VEP_SIF}"
VEP_IMAGE="$VEP_IMAGE"
: > "\$LOG_FILE"

mkdir -p "\$OUT/cache" "\$OUT/fasta"

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  seq_die "No Apptainer/Singularity runtime available for VEP setup."
fi

# Singularity on some hosts does not auto-bind arbitrary mountpoints like /media,
# so we bind the workspace explicitly for offline cache/FASTA access.
BIND_ARGS=(--bind "\$SEQ_BASE:\$SEQ_BASE")

if [[ ! -s "\$VEP_SIF" ]]; then
  "\$CONTAINER_RUNTIME" pull "\$VEP_SIF" "\$VEP_IMAGE" >>"\$LOG_FILE" 2>&1
fi

if [[ ! -d "\$OUT/cache/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" ]]; then
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

"\$CONTAINER_RUNTIME" exec "\${BIND_ARGS[@]}" "\$VEP_SIF" vep \\
  --offline \\
  --cache \\
  --dir_cache "\$OUT/cache" \\
  --species homo_sapiens \\
  --assembly GRCh38 \\
  --fasta "\$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \\
  --input_file "\$OUT/toy_input.vcf" \\
  --output_file "\$OUT/toy_output.vep.txt" \\
  --force_overwrite \\
  --no_stats \\
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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
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
check_file "$SEQ_RNASEQ_FASTA"
check_file "$SEQ_RNASEQ_GTF"
check_file "$SEQ_RNASEQ_STAR_INDEX/Genome"
check_file "$SEQ_RNASEQ_SALMON_INDEX/versionInfo.json"
check_file "$SEQ_CLINVAR_VCF"
check_file "$SEQ_VEP_CACHE_ROOT/homo_sapiens/${SEQ_VEP_CACHE_VERSION}_GRCh38/info.txt"

if [[ -n "${SEQ_SAREK_PON:-}" && -s "${SEQ_SAREK_PON:-}" ]]; then
  printf 'OK\t%s\n' "$SEQ_SAREK_PON"
else
  printf 'WARN\tNo controlled MuTect2 Panel of Normals configured; tumor-normal runs will proceed without --pon unless you add one later.\n'
fi
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
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
Rscript "$SCRIPT_DIR/deseq2_starter.R" "$@"
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

  write_with_backup "$SCRIPTS_DIR/build_pdac_autodrafts.py" <<'EOF'
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re


WES_RE = re.compile(
    r"^(?P<label>\d+)(?P<status>[NT])_\d+WESFastq_S\d+_(?P<lane>L\d{3})_R(?P<read>[12])_001\.fastq\.gz$"
)
RNA_RE = re.compile(
    r"^(?P<label>\d+)(?P<status>[NT])_\d+RNA_S\d+_(?P<lane>L\d{3})_R(?P<read>[12])_001\.fastq\.gz$"
)
WES_PLOIDY_RE = re.compile(r"^(?P<label>\d+)(?P<status>[NT])WESFastq\.ploidy_estimation_metrics\.csv$")
TWIST_BED_NAME = "Twist_ILMN_Exome_2.5_Plus_Panel.hg38.bed"


@dataclass(frozen=True)
class SampleKey:
    patient: str
    status_code: str
    lane: str


@dataclass(frozen=True)
class IntervalConsensus:
    source_path: str
    consensus_hash: str
    consensus_count: int
    total_count: int
    outlier_paths: tuple[str, ...]
    link_path: str | None = None


def backup_if_exists(path: Path, backup_root: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = backup_root / stamp
    backup_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup_dir / path.name)


def write_text(path: Path, text: str, backup_root: Path) -> None:
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]], comments: list[str], backup_root: Path) -> None:
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for line in comments:
            handle.write(f"# {line}\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def discover_fastqs(input_root: Path, regex: re.Pattern[str], skip_dirs: set[str]) -> dict[SampleKey, dict[str, str]]:
    found: dict[SampleKey, dict[str, str]] = defaultdict(dict)
    for fastq in sorted(input_root.rglob("*.fastq.gz")):
        if any(part in skip_dirs for part in fastq.parts):
            continue
        match = regex.match(fastq.name)
        if not match:
            continue
        key = SampleKey(
            patient=match.group("label"),
            status_code=match.group("status"),
            lane=match.group("lane"),
        )
        found[key][match.group("read")] = str(fastq)
    return found


def count_matching_files(input_root: Path, suffixes: tuple[str, ...], include_parts: tuple[str, ...], skip_dirs: set[str]) -> list[Path]:
    hits: list[Path] = []
    for path in sorted(input_root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in skip_dirs for part in path.parts):
            continue
        if suffixes and not path.name.endswith(suffixes):
            continue
        if include_parts and not any(part in path.parts for part in include_parts):
            continue
        hits.append(path)
    return hits


def infer_sex_from_ploidy_metrics(path: Path) -> str | None:
    values: dict[str, float] = {}
    with path.open() as handle:
        for line in handle:
            parts = [part.strip() for part in line.rstrip().split(",")]
            if len(parts) >= 4 and parts[0] == "PLOIDY ESTIMATION":
                try:
                    values[parts[2]] = float(parts[3])
                except ValueError:
                    continue
    autosomal = values.get("Autosomal median coverage")
    x_cov = values.get("X median coverage")
    y_cov = values.get("Y median coverage")
    if not autosomal or x_cov is None or y_cov is None:
        return None
    x_ratio = x_cov / autosomal
    y_ratio = y_cov / autosomal
    if y_ratio < 0.05 and x_ratio > 0.7:
        return "XX"
    if y_ratio > 0.10 and 0.25 < x_ratio < 0.7:
        return "XY"
    return None


def infer_patient_sex_map(input_root: Path, skip_dirs: set[str]) -> tuple[dict[str, str], list[str]]:
    normal_calls: dict[str, str] = {}
    tumor_calls: dict[str, str] = {}
    for path in sorted(input_root.rglob("*.ploidy_estimation_metrics.csv")):
        if any(part in skip_dirs for part in path.parts):
            continue
        match = WES_PLOIDY_RE.match(path.name)
        if not match:
            continue
        sex = infer_sex_from_ploidy_metrics(path)
        if sex is None:
            continue
        patient = match.group("label")
        status = match.group("status")
        if status == "N":
            normal_calls[patient] = sex
        else:
            tumor_calls[patient] = sex

    notes: list[str] = []
    discordant = sorted(
        patient for patient, sex in tumor_calls.items() if patient in normal_calls and normal_calls[patient] != sex
    )
    if discordant:
        notes.append(
            "Tumor ploidy was not used for patient sex inference because the following tumor calls differed from their normal: "
            + ", ".join(discordant)
        )
    return normal_calls, notes


def discover_interval_consensus(input_root: Path, skip_dirs: set[str]) -> IntervalConsensus | None:
    candidates = [
        path
        for path in sorted(input_root.rglob(TWIST_BED_NAME))
        if path.is_file() and not any(part in skip_dirs for part in path.parts)
    ]
    if not candidates:
        return None

    by_hash: dict[str, list[Path]] = defaultdict(list)
    for path in candidates:
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        by_hash[digest].append(path)

    consensus_hash, consensus_paths = max(by_hash.items(), key=lambda item: (len(item[1]), str(item[1][0])))
    outlier_paths = tuple(
        str(path) for digest, paths in by_hash.items() if digest != consensus_hash for path in sorted(paths)
    )
    return IntervalConsensus(
        source_path=str(sorted(consensus_paths)[0]),
        consensus_hash=consensus_hash,
        consensus_count=len(consensus_paths),
        total_count=len(candidates),
        outlier_paths=outlier_paths,
    )


def ensure_symlink(path: Path, source: Path, backup_root: Path) -> None:
    if path.is_symlink():
        try:
            if path.resolve() == source.resolve():
                return
        except FileNotFoundError:
            pass
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        path.unlink()
    path.symlink_to(source)


def build_wes_rows(pairs: dict[SampleKey, dict[str, str]], patient_sex_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        rows.append(
            {
                "patient": key.patient,
                "sex": patient_sex_map.get(key.patient, ""),
                "status": "0" if key.status_code == "N" else "1",
                "sample": f"{key.patient}{key.status_code}_WES",
                "lane": key.lane,
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "bam": "",
                "bai": "",
                "cram": "",
                "crai": "",
            }
        )
    return rows


def build_rna_rows(pairs: dict[SampleKey, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        rows.append(
            {
                "sample": f"{key.patient}{key.status_code}_RNA",
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "strandedness": "REVIEW_ME",
            }
        )
    return rows


def build_rnafusion_rows(pairs: dict[SampleKey, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        rows.append(
            {
                "sample": f"{key.patient}{key.status_code}_RNA",
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "bam": "",
                "bai": "",
                "cram": "",
                "crai": "",
                "junctions": "",
                "splice_junctions": "",
                "strandedness": "REVIEW_ME",
                "seq_platform": "ILLUMINA",
                "seq_center": "LOCAL",
            }
        )
    return rows


def launch_drafts_markdown(workspace_base: Path, interval_consensus: IntervalConsensus | None) -> str:
    scripts_dir = workspace_base / "scripts"
    samplesheets_dir = workspace_base / "samplesheets"
    bulk_results_root = Path("/media/user/PDAC_SEQ_analysis/results")
    if not bulk_results_root.exists():
        bulk_results_root = workspace_base / "results"
    interval_path = interval_consensus.link_path if interval_consensus and interval_consensus.link_path else (
        interval_consensus.source_path if interval_consensus else "/REVIEW_ME/EXOME_CAPTURE_INTERVALS.bed"
    )
    return "\n".join(
        [
            "# PDAC Launch Drafts",
            "",
            "- Review the generated samplesheets before launching.",
            "- Keep patient data in place. These commands read directly from the source FASTQs.",
            f"- Suggested bulk-analysis output root: `{bulk_results_root}`",
            "- Keep RNA fusion `--outdir` on `/media/user/SEQ/results` or another Linux-native filesystem.",
            "",
            "## WES Tumor-Normal",
            "",
            "```bash",
            f"{scripts_dir}/run_sarek.sh --mode tumor-normal \\",
            f"  --samplesheet {samplesheets_dir}/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \\",
            f"  --intervals '{interval_path}' \\",
            f"  --outdir {bulk_results_root}/sarek_tumor_normal \\",
            "  -resume",
            "```",
            "",
            "## WES Germline",
            "",
            "```bash",
            f"{scripts_dir}/run_sarek.sh --mode germline \\",
            f"  --samplesheet {samplesheets_dir}/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \\",
            f"  --intervals '{interval_path}' \\",
            f"  --outdir {bulk_results_root}/sarek_germline \\",
            "  -resume",
            "```",
            "",
            "## RNA-Seq Expression",
            "",
            "Edit the RNA samplesheet first and replace every `REVIEW_ME` with `unstranded`, `forward`, or `reverse`.",
            "",
            "```bash",
            f"{scripts_dir}/run_rnaseq.sh \\",
            f"  --samplesheet {samplesheets_dir}/rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv \\",
            f"  --outdir {bulk_results_root}/rnaseq_expression \\",
            "  -resume",
            "```",
            "",
            "## RNA Fusion",
            "",
            "Edit the RNA fusion samplesheet first and replace every `REVIEW_ME` with `unstranded`, `forward`, or `reverse`.",
            "",
            "```bash",
            f"{scripts_dir}/run_rnafusion.sh \\",
            f"  --samplesheet {samplesheets_dir}/rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv \\",
            f"  --outdir {workspace_base}/results/rnafusion_pdac \\",
            "  -resume",
            "```",
        ]
    )


def inventory_markdown(
    workspace_base: Path,
    input_root: Path,
    wes_rows: list[dict[str, str]],
    rna_rows: list[dict[str, str]],
    rna_bams: list[Path],
    wes_bams: list[Path],
    wes_vcfs: list[Path],
    patient_sex_map: dict[str, str],
    sex_notes: list[str],
    interval_consensus: IntervalConsensus | None,
) -> str:
    tumor_ids = sorted({row["patient"] for row in wes_rows if row["status"] == "1"}, key=int)
    normal_ids = sorted({row["patient"] for row in wes_rows if row["status"] == "0"}, key=int)
    generated_files = [
        workspace_base / "samplesheets" / "sarek_samplesheet.PDAC_WES_fastq_autodraft.csv",
        workspace_base / "samplesheets" / "rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        workspace_base / "samplesheets" / "rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        workspace_base / "docs" / "PDAC_LAUNCH_DRAFTS.md",
    ]
    lines = [
        "# PDAC Input Inventory",
        "",
        f"- Input root: `{input_root}`",
        f"- WES FASTQ rows drafted: {len(wes_rows)}",
        f"- RNA FASTQ rows drafted: {len(rna_rows)}",
        f"- WES BAMs discovered: {len(wes_bams)}",
        f"- RNA BAMs discovered: {len(rna_bams)}",
        f"- WES hard-filtered VCFs discovered: {len(wes_vcfs)}",
        f"- WES patient sex inferred from normal ploidy metrics: {len(patient_sex_map)}/{len(normal_ids)}",
        "",
        "## Cohort IDs",
        "",
        f"- Normal labels: {', '.join(normal_ids)}",
        f"- Tumor labels: {', '.join(tumor_ids)}",
        "",
        "## Generated draft files",
        "",
    ]
    lines.extend(f"- `{path}`" for path in generated_files)
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- These drafts do not move or modify source data.",
            "- WES `sex` is filled only when a matching normal ploidy metric supports an `XX` or `XY` inference.",
            "- RNA `strandedness` is left as `REVIEW_ME`. DRAGEN replay files show `rna-library-type=A` (auto-detect), which is not a direct nf-core samplesheet value.",
            "- Existing DRAGEN BAM/VCF outputs were inventoried but not used as primary workflow inputs in these drafts.",
        ]
    )
    if interval_consensus is not None:
        interval_label = interval_consensus.link_path or interval_consensus.source_path
        lines.append(
            f"- Consensus exome interval candidate ({interval_consensus.consensus_count}/{interval_consensus.total_count} copies): `{interval_label}`"
        )
        if interval_consensus.outlier_paths:
            lines.append(f"- Interval outlier kept for review: `{interval_consensus.outlier_paths[0]}`")
    lines.extend(f"- {note}" for note in sex_notes)
    lines.extend(
        [
            "",
            "## Inferred WES Sex By Patient",
            "",
        ]
    )
    for patient in sorted(patient_sex_map, key=int):
        lines.append(f"- {patient}: {patient_sex_map[patient]}")
    lines.extend(
        [
            "",
            "## Example existing processed files",
            "",
        ]
    )
    for heading, files in (
        ("WES BAM examples", wes_bams[:5]),
        ("RNA BAM examples", rna_bams[:5]),
        ("WES VCF examples", wes_vcfs[:5]),
    ):
        lines.append(f"### {heading}")
        lines.append("")
        if files:
            lines.extend(f"- `{path}`" for path in files)
        else:
            lines.append("- None discovered")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--samplesheets-dir", required=True)
    parser.add_argument("--docs-dir", required=True)
    parser.add_argument("--backup-dir", required=True)
    parser.add_argument("--link-intervals-dir")
    args = parser.parse_args()

    input_root = Path(args.input_root)
    samplesheets_dir = Path(args.samplesheets_dir)
    docs_dir = Path(args.docs_dir)
    backup_dir = Path(args.backup_dir)
    link_intervals_dir = Path(args.link_intervals_dir) if args.link_intervals_dir else None
    skip_dirs = {"SEQ_analysis", "work", "results", "smoke_tests"}
    workspace_base = samplesheets_dir.parent

    wes_pairs = discover_fastqs(input_root, WES_RE, skip_dirs)
    rna_pairs = discover_fastqs(input_root, RNA_RE, skip_dirs)
    patient_sex_map, sex_notes = infer_patient_sex_map(input_root, skip_dirs)
    interval_consensus = discover_interval_consensus(input_root, skip_dirs)
    if interval_consensus is not None and link_intervals_dir is not None:
        link_path = link_intervals_dir / "PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed"
        ensure_symlink(link_path, Path(interval_consensus.source_path), backup_dir)
        interval_consensus = IntervalConsensus(
            source_path=interval_consensus.source_path,
            consensus_hash=interval_consensus.consensus_hash,
            consensus_count=interval_consensus.consensus_count,
            total_count=interval_consensus.total_count,
            outlier_paths=interval_consensus.outlier_paths,
            link_path=str(link_path),
        )

    wes_rows = build_wes_rows(wes_pairs, patient_sex_map)
    rna_rows = build_rna_rows(rna_pairs)
    rnafusion_rows = build_rnafusion_rows(rna_pairs)

    wes_bams = count_matching_files(input_root, (".bam",), ("1900NGS_WES_Enrichment_stats-10282273",), skip_dirs)
    rna_bams = count_matching_files(input_root, (".bam",), ("DRAGEN_RNA_03_17_2025_22_02_16-59174132",), skip_dirs)
    wes_vcfs = count_matching_files(input_root, (".hard-filtered.vcf.gz",), ("1900NGS_WES_Enrichment_stats-10282273",), skip_dirs)

    write_csv(
        samplesheets_dir / "sarek_samplesheet.PDAC_WES_fastq_autodraft.csv",
        ["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2", "bam", "bai", "cram", "crai"],
        wes_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "sex is inferred from normal-sample ploidy metrics when available; confirm before use.",
            "status uses 0 for normal and 1 for tumor.",
        ],
        backup_dir,
    )
    write_csv(
        samplesheets_dir / "rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        rna_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "Replace REVIEW_ME with unstranded, forward, or reverse before launching.",
        ],
        backup_dir,
    )
    write_csv(
        samplesheets_dir / "rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        ["sample", "fastq_1", "fastq_2", "bam", "bai", "cram", "crai", "junctions", "splice_junctions", "strandedness", "seq_platform", "seq_center"],
        rnafusion_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "Replace REVIEW_ME with the correct strandedness before launching.",
        ],
        backup_dir,
    )
    write_text(
        docs_dir / "PDAC_INPUT_INVENTORY.md",
        inventory_markdown(
            workspace_base,
            input_root,
            wes_rows,
            rna_rows,
            rna_bams,
            wes_bams,
            wes_vcfs,
            patient_sex_map,
            sex_notes,
            interval_consensus,
        ),
        backup_dir,
    )
    write_text(
        docs_dir / "PDAC_LAUNCH_DRAFTS.md",
        launch_drafts_markdown(workspace_base, interval_consensus),
        backup_dir,
    )

    print(f"Wrote {len(wes_rows)} WES rows, {len(rna_rows)} RNA rows, and inventory docs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
EOF
  chmod +x "$SCRIPTS_DIR/build_pdac_autodrafts.py"

  write_with_backup "$SCRIPTS_DIR/draft_pdac_inputs.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

INPUT_ROOT="${1:-}"
[[ -n "$INPUT_ROOT" ]] || seq_die "Usage: $0 <pdac_input_root>"
[[ -d "$INPUT_ROOT" ]] || seq_die "Input root does not exist: $INPUT_ROOT"

python3 "$SCRIPT_DIR/build_pdac_autodrafts.py" \
  --input-root "$INPUT_ROOT" \
  --samplesheets-dir "$SEQ_BASE/samplesheets" \
  --docs-dir "$SEQ_BASE/docs" \
  --backup-dir "$SEQ_BASE/tmp/manual_backups/pdac_autodrafts" \
  --link-intervals-dir "$SEQ_BASE/refs/optional"
EOF
  chmod +x "$SCRIPTS_DIR/draft_pdac_inputs.sh"
}

write_wrapper_scripts() {
  write_with_backup "$SCRIPTS_DIR/run_sarek.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

MODE=""
SAMPLESHEET=""
INTERVALS=""
OUTDIR=""
WORKDIR=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --intervals) INTERVALS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$MODE" ]] || seq_die "Missing --mode (tumor-normal or germline)"
[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
[[ -n "$INTERVALS" ]] || seq_die "Missing --intervals with the exome capture BED / interval list"

seq_require_file "$SAMPLESHEET"
seq_require_file "$INTERVALS"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline sarek --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/sarek_samplesheet.XXXXXX.csv")"
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
  SORTED_INTERVALS="$(mktemp "$SEQ_TMP/sarek_intervals.XXXXXX.bed")"
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

PON_ARGS=()
if [[ -n "${SEQ_SAREK_PON:-}" && -s "${SEQ_SAREK_PON:-}" ]]; then
  PON_ARGS+=(--pon "$SEQ_SAREK_PON")
fi
if [[ -n "${SEQ_SAREK_PON_TBI:-}" && -s "${SEQ_SAREK_PON_TBI:-}" ]]; then
  PON_ARGS+=(--pon_tbi "$SEQ_SAREK_PON_TBI")
fi
if [[ "$MODE" == "tumor-normal" && "${#PON_ARGS[@]}" -eq 0 ]]; then
  seq_log "WARNING: No controlled MuTect2 Panel of Normals is configured; continuing without --pon."
fi

OUTDIR="${OUTDIR:-${SEQ_SAREK_OUTDIR:-$SEQ_RESULTS/sarek}}"
WORKDIR="${WORKDIR:-${SEQ_SAREK_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
LOCAL_IGENOMES_BASE="${SEQ_SAREK_IGENOMES_BASE:-$SEQ_BASE/refs/igenomes_stub}"
LOCAL_SNPEFF_CACHE="${SEQ_SAREK_SNPEFF_CACHE:-$SEQ_BASE/refs/annotation/snpeff_cache}"
mkdir -p "$OUTDIR" "$WORKDIR" "$LOCAL_IGENOMES_BASE" "$LOCAL_SNPEFF_CACHE"
nextflow run "$SEQ_BASE/pipelines/sarek-3.8.1/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/sarek.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
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
  "${PON_ARGS[@]}" \
  --igenomes_base "$LOCAL_IGENOMES_BASE" \
  --snpeff_cache "$LOCAL_SNPEFF_CACHE" \
  --vep_cache "$SEQ_VEP_CACHE_ROOT" \
  --vep_cache_version "$SEQ_VEP_CACHE_VERSION" \
  --vep_genome GRCh38 \
  --vep_species homo_sapiens \
  --intervals "$INTERVALS" \
  --igenomes_ignore \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_sarek.sh"

  write_with_backup "$SCRIPTS_DIR/run_rnaseq.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

SAMPLESHEET=""
OUTDIR=""
WORKDIR=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline rnaseq --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/rnaseq_samplesheet.XXXXXX.csv")"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

OUTDIR="${OUTDIR:-${SEQ_RNASEQ_OUTDIR:-$SEQ_RESULTS/rnaseq}}"
WORKDIR="${WORKDIR:-${SEQ_RNASEQ_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
mkdir -p "$OUTDIR" "$WORKDIR"
nextflow run "$SEQ_BASE/pipelines/rnaseq-3.24.0/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/rnaseq.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
  --fasta "$SEQ_RNASEQ_FASTA" \
  --gtf "$SEQ_RNASEQ_GTF" \
  --star_index "$SEQ_RNASEQ_STAR_INDEX" \
  --salmon_index "$SEQ_RNASEQ_SALMON_INDEX" \
  --igenomes_ignore \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_rnaseq.sh"

  write_with_backup "$SCRIPTS_DIR/run_rnafusion.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

SAMPLESHEET=""
OUTDIR=""
WORKDIR=""
EXTRA_ARGS=()
RESUME=""
TOOLS="${SEQ_RNAFUSION_DEFAULT_TOOLS:-arriba,fusioncatcher,salmon}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --tools) TOOLS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline rnafusion --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/rnafusion_samplesheet.XXXXXX.csv")"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

if [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/star" ]]; then
  seq_die "RNA fusion STAR references are missing. Re-run sync_rnafusion_refs.sh first."
fi
if [[ ",$TOOLS," == *",fusioncatcher,"* ]] && [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/fusioncatcher" ]]; then
  seq_die "FusionCatcher was requested but its references are missing. Re-run sync_rnafusion_refs.sh with FusionCatcher enabled first."
fi
if [[ ",$TOOLS," == *",starfusion,"* ]] && [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/starfusion" ]]; then
  seq_die "STAR-Fusion was requested but its references are missing. Re-run sync_rnafusion_refs.sh with STAR-Fusion enabled first."
fi

OUTDIR="${OUTDIR:-${SEQ_RNAFUSION_OUTDIR:-$SEQ_RESULTS/rnafusion}}"
WORKDIR="${WORKDIR:-${SEQ_RNAFUSION_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
mkdir -p "$OUTDIR" "$WORKDIR"
nextflow run "$SEQ_BASE/pipelines/rnafusion-4.1.0/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/rnafusion.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
  --genomes_base "$SEQ_RNAFUSION_GENOMES_BASE" \
  --genome GRCh38 \
  --genome_gencode_version 46 \
  --tools "$TOOLS" \
  --no_cosmic \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
EOF
  chmod +x "$SCRIPTS_DIR/run_rnafusion.sh"
}

write_config_files() {
  local singularity_enabled="false"
  local apptainer_enabled="false"
  local rnafusion_default_tools="arriba,fusioncatcher,salmon"
  [[ "$SYNC_FUSIONCATCHER" == "true" ]] && rnafusion_default_tools="${rnafusion_default_tools},fusioncatcher"
  [[ "$SYNC_STARFUSION" == "true" ]] && rnafusion_default_tools="${rnafusion_default_tools},starfusion"
  case "$RUNTIME_PROFILE" in
    singularity) singularity_enabled="true" ;;
    apptainer) apptainer_enabled="true" ;;
  esac

  write_with_backup "$CONFIG_DIR/resources.env" <<EOF
export SEQ_BASE="$BASE"
export SEQ_WORK="$WORK_DIR"
export SEQ_RESULTS="$RESULTS_DIR"
export NXF_HOME="$TMP_DIR/nextflow"
export NXF_WORK="$WORK_DIR"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE_DIR"
export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHE_DIR"
export SINGULARITY_CACHEDIR="$CONTAINER_CACHE_DIR"
export APPTAINER_CACHEDIR="$APPTAINER_CACHE_DIR"
export SEQ_RUNTIME_PROFILE="$RUNTIME_PROFILE"
export SEQ_SYNC_STARFUSION="$SYNC_STARFUSION"
export SEQ_SYNC_FUSIONCATCHER="$SYNC_FUSIONCATCHER"

export SEQ_SAREK_FASTA="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta"
export SEQ_SAREK_FASTA_FAI="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.fasta.fai"
export SEQ_SAREK_DICT="$REFS_DIR/gatk_bundle/Homo_sapiens_assembly38.dict"
export SEQ_SAREK_DBSNP="$REFS_DIR/gatk_bundle/dbsnp_144.hg38.vcf.gz"
export SEQ_SAREK_DBSNP_TBI="$REFS_DIR/gatk_bundle/dbsnp_144.hg38.vcf.gz.tbi"
export SEQ_SAREK_KNOWN_INDELS="$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
export SEQ_SAREK_KNOWN_INDELS_TBI="$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
export SEQ_SAREK_KNOWN_SNPS="$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz"
export SEQ_SAREK_KNOWN_SNPS_TBI="$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz.tbi"
export SEQ_SAREK_GERMLINE_RESOURCE="$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz"
export SEQ_SAREK_GERMLINE_RESOURCE_TBI="$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz.tbi"
export SEQ_SAREK_PON=""
export SEQ_SAREK_PON_TBI=""

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
export SEQ_RNAFUSION_DEFAULT_TOOLS="$rnafusion_default_tools"

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
  dbsnp: "$REFS_DIR/gatk_bundle/dbsnp_144.hg38.vcf.gz"
  known_indels: "$REFS_DIR/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  known_snps: "$REFS_DIR/gatk_bundle/1000G_omni2.5.hg38.vcf.gz"
  germline_resource: "$REFS_DIR/gatk_bundle/af-only-gnomad.hg38.vcf.gz"
  pon_controlled_manual: "$REFS_DIR/gatk_bundle/PON_CONTROLLED_REQUIRED.md"
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
  default_tools: "$rnafusion_default_tools"
  starfusion_synced: "$SYNC_STARFUSION"
  fusioncatcher_synced: "$SYNC_FUSIONCATCHER"

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
  enabled = $singularity_enabled
  autoMounts = true
  cacheDir = '$CONTAINER_CACHE_DIR'
}

apptainer {
  enabled = $apptainer_enabled
  autoMounts = true
  cacheDir = '$APPTAINER_CACHE_DIR'
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
  enabled = $singularity_enabled
  autoMounts = true
  cacheDir = '$CONTAINER_CACHE_DIR'
}

apptainer {
  enabled = $apptainer_enabled
  autoMounts = true
  cacheDir = '$APPTAINER_CACHE_DIR'
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
  enabled = $singularity_enabled
  autoMounts = true
  cacheDir = '$CONTAINER_CACHE_DIR'
}

apptainer {
  enabled = $apptainer_enabled
  autoMounts = true
  cacheDir = '$APPTAINER_CACHE_DIR'
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
  local rnafusion_default_tools="arriba,fusioncatcher,salmon"
  [[ "$SYNC_FUSIONCATCHER" == "true" ]] && rnafusion_default_tools="${rnafusion_default_tools},fusioncatcher"
  [[ "$SYNC_STARFUSION" == "true" ]] && rnafusion_default_tools="${rnafusion_default_tools},starfusion"

  write_with_backup "$SAMPLESHEETS_DIR/sarek_samplesheet.csv" <<EOF
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
# VEP cache root: $REFS_DIR/vep/cache
# ClinVar VCF: $REFS_DIR/clinvar/clinvar.vcf.gz
patient,sex,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai
PANCREAS001,XY,0,PANCREAS001_N,L001,/ABSOLUTE/PATH/TO/NORMAL_R1.fastq.gz,/ABSOLUTE/PATH/TO/NORMAL_R2.fastq.gz,,,,
PANCREAS001,XY,1,PANCREAS001_T,L001,/ABSOLUTE/PATH/TO/TUMOR_R1.fastq.gz,/ABSOLUTE/PATH/TO/TUMOR_R2.fastq.gz,,,,
EOF

  write_with_backup "$SAMPLESHEETS_DIR/rnaseq_samplesheet.csv" <<EOF
# Remove comment lines before using outside the provided wrappers.
# strandedness: unstranded, forward, or reverse
# fasta: $REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa
# gtf: $REFS_DIR/fusion/GRCh38/gencode_v46/gencode/Homo_sapiens_GRCh38_46.gtf
# salmon index: $REFS_DIR/fusion/GRCh38/gencode_v46/salmon
sample,fastq_1,fastq_2,strandedness
PANCREAS001_TUMOR_RNA,/ABSOLUTE/PATH/TO/RNA_R1.fastq.gz,/ABSOLUTE/PATH/TO/RNA_R2.fastq.gz,forward
EOF

  write_with_backup "$SAMPLESHEETS_DIR/rnafusion_samplesheet.csv" <<EOF
# Remove comment lines before using outside the provided wrappers.
# One of fastq_1, bam, cram, junctions, or splice_junctions must be populated.
# genomes_base: $REFS_DIR/fusion
# available default tools on this host: $rnafusion_default_tools
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
  local fusioncatcher_note="not synced on this host"
  if [[ "$SYNC_FUSIONCATCHER" == "true" ]]; then
    fusion_tools="${fusion_tools},fusioncatcher"
    fusioncatcher_note="synced from nf-core public S3 bucket"
  fi
  if [[ "$SYNC_STARFUSION" == "true" ]]; then
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
  - \`dbsnp_144.hg38.vcf.gz\`
  - \`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\`
  - \`Homo_sapiens_assembly38.known_indels.vcf.gz\`
  - \`1000G_omni2.5.hg38.vcf.gz\`
  - \`1000G_phase1.snps.high_confidence.hg38.vcf.gz\`
  - \`Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz\`
  - \`hapmap_3.3.hg38.vcf.gz\`
- Mutect2 resources:
  - \`af-only-gnomad.hg38.vcf.gz\`
  - Controlled PoN not downloaded automatically; see \`$REFS_DIR/gatk_bundle/PON_CONTROLLED_REQUIRED.md\`

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
- FusionCatcher status: $fusioncatcher_note
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
4. Optional controlled, licensed, or token-gated resources only if you need them later:
   - GDC MuTect2 Panel of Normals if you want standard PoN-backed MuTect2 filtering
   - COSMIC
   - OncoKB
   - dbNSFP
   - CADD
5. If you plan to run the DESeq2 starter, install \`DESeq2\` for \`Rscript\` first or treat the provided script as a template.

Notes:

- Comment lines in the samplesheet templates are safe when using the provided wrappers; the wrappers strip them automatically before launch.
- If your patient data already live in a PDAC-style directory tree, you can draft cohort-specific samplesheets and launch notes with:

```bash
$SCRIPTS_DIR/draft_pdac_inputs.sh /ABSOLUTE/PATH/TO/PDAC
```

- For nf-core/rnafusion, prefer an output directory on a Linux-native filesystem. Some external \`fuseblk\`/NTFS/exFAT mounts reject the pipeline's published \`stringtie/[:]\` path, so use an ext4-style \`--outdir\` when needed.
- The current host has about ${AVAILABLE_GB}G free at setup time, so STAR-Fusion references were set to \`$SYNC_STARFUSION\`.
- The public GDC/Broad bundle no longer exposes the old \`dbsnp_146.hg38.vcf.gz\` path used by some older workflows, so this setup uses the public \`dbsnp_144.hg38.vcf.gz\` object that is currently available.
- The official GDC PoN files are controlled and were intentionally not auto-downloaded.
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
- Public references for GRCh38 DNA, RNA, VEP, and ClinVar, plus a documented placeholder for the controlled GDC PoN
- Idempotent helper scripts under \`$SCRIPTS_DIR\`
- Wrapper launchers for each pipeline
- A cohort autodraft helper for PDAC-style input trees that can pre-build draft samplesheets and launch notes

## Safety

- Patient data were not moved.
- Existing files are backed up before this setup overwrites them.
- Nextflow work is kept under \`$WORK_DIR\`.
- Results are written under \`$RESULTS_DIR\`.
- References and containers remain under \`$BASE\` unless you override them.
- For nf-core/rnafusion, prefer a Linux-native output path. Some external \`fuseblk\`/NTFS/exFAT mounts reject the pipeline's published \`stringtie/[:]\` path, so use \`--outdir\` to point at an ext4-style location when needed.
EOF

  write_with_backup "$LAUNCH_EXAMPLES" <<EOF
# Launch Examples

## Draft PDAC cohort samplesheets from an existing patient-data tree

\`\`\`bash
$SCRIPTS_DIR/draft_pdac_inputs.sh /ABSOLUTE/PATH/TO/PDAC
\`\`\`

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

If your preferred storage volume is \`fuseblk\`/NTFS/exFAT or otherwise rejects colon-containing paths, override the output root explicitly:

\`\`\`bash
$SCRIPTS_DIR/run_rnafusion.sh --samplesheet $SAMPLESHEETS_DIR/rnafusion_samplesheet.csv --outdir $BASE/results/rnafusion -resume
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
  if [[ "$RUNTIME_PROFILE" == "missing" ]]; then
    printf 'No supported Nextflow runtime detected; smoke tests skipped.\n' > "$LOG_DIR/test_sarek.log"
    cp "$LOG_DIR/test_sarek.log" "$LOG_DIR/test_rnaseq.log"
    cp "$LOG_DIR/test_sarek.log" "$LOG_DIR/test_rnafusion.log"
    return 0
  fi

  mkdir -p "$TMP_DIR/nextflow"
  export NXF_HOME="$TMP_DIR/nextflow"
  export NXF_WORK="$WORK_DIR"
  export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE_DIR"
  export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHE_DIR"
  export SINGULARITY_CACHEDIR="$CONTAINER_CACHE_DIR"
  export APPTAINER_CACHEDIR="$APPTAINER_CACHE_DIR"

  local profile="$RUNTIME_PROFILE"
  local nf
  nf="$(command -v nextflow)"
  local smoke_out="$BASE/results/smoke_tests"
  local smoke_snpeff_cache="$REFS_DIR/optional/smoke_snpeff_cache"
  local smoke_sentieon_model="$REFS_DIR/optional/SentieonDNAscopeModel1.1.model"
  mkdir -p "$smoke_out/sarek" "$smoke_out/rnaseq" "$smoke_out/rnafusion"
  mkdir -p "$smoke_snpeff_cache/WBcel235.99"
  : > "$smoke_sentieon_model"
  : > "$smoke_snpeff_cache/WBcel235.99/.keep"

  "$nf" run "$PIPELINES_DIR/sarek-3.8.1/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" --outdir "$smoke_out/sarek" --snpeff_cache "$smoke_snpeff_cache" --vep_cache "$REFS_DIR/vep/cache" --sentieon_dnascope_model "$smoke_sentieon_model" > "$LOG_DIR/test_sarek.log" 2>&1 || true
  "$nf" run "$PIPELINES_DIR/rnaseq-3.24.0/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" --outdir "$smoke_out/rnaseq" > "$LOG_DIR/test_rnaseq.log" 2>&1 || true
  "$nf" run "$PIPELINES_DIR/rnafusion-4.1.0/main.nf" -profile test,"$profile" -stub-run -ansi-log false -w "$WORK_DIR" --outdir "$smoke_out/rnafusion" --genomes_base "$REFS_DIR/fusion" > "$LOG_DIR/test_rnafusion.log" 2>&1 || true
}

write_setup_summary() {
  local starfusion_status="not synced"
  local fusioncatcher_status="not synced"
  if [[ "$SYNC_STARFUSION" == "true" ]]; then
    starfusion_status="synced"
  fi
  if [[ "$SYNC_FUSIONCATCHER" == "true" ]]; then
    fusioncatcher_status="synced"
  fi

  write_with_backup "$SETUP_SUMMARY" <<EOF
# Setup Summary

Generated: $DATE_ISO

## System audit summary

- Host: \`$HOSTNAME_FQDN\`
- User: \`$USER_NAME\`
- Runtime profile selected: \`$RUNTIME_PROFILE\`
- Work directory: \`$WORK_DIR\`
- Results root: \`$RESULTS_DIR\`
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
- Broad/GDC public hg38 bundle components for BQSR, VQSR, and the Mutect2 germline resource
- Controlled MuTect2 PoN not auto-downloaded; placeholder documentation added instead
- Ensembl VEP offline cache + FASTA + local container
- ClinVar GRCh38 VCF
- RNA fusion published references: Arriba, fusion-report, HGNC, STAR, Salmon
- FusionCatcher references: $fusioncatcher_status
- STAR-Fusion references: $starfusion_status

## Downloaded automatically

- Public nf-core RNA fusion references from \`s3://nf-core-awsmegatests/rnafusion/references/\`
- Public Broad hg38 bundle files from Google public storage
- Public ClinVar GRCh38 VCF
- Public Ensembl VEP cache and FASTA
- Public GENCODE v46 transcript FASTA
- Controlled GDC PoN intentionally skipped because it requires authorized access

## Still needs manual input

- Patient FASTQ/BAM/CRAM files
- Exome capture BED / interval list
- Optional controlled GDC MuTect2 PoN if you want PoN-backed MuTect2 filtering later
- Optional licensed resources if desired later
- DESeq2 package installation if you want to run the starter directly
- RNA strandedness review if you auto-draft samplesheets from an existing PDAC-style input tree

## Exact launch commands

### Draft PDAC cohort inputs from an existing local tree

\`\`\`bash
$SCRIPTS_DIR/draft_pdac_inputs.sh /ABSOLUTE/PATH/TO/PDAC
\`\`\`

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
  ensure_nextflow
  ensure_htslib_tools
  preflight_requirements

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

  if [[ "$RUN_SMOKE_TESTS" == "true" ]]; then
    log "Running smoke tests"
    run_smoke_tests
  else
    log "Skipping smoke tests by request"
  fi

  log "Validating references"
  bash "$SCRIPTS_DIR/validate_reference_integrity.sh" > "$LOG_DIR/validate_reference_integrity.log" 2>&1 || true

  log "Writing setup summary"
  write_setup_summary

  log "Bootstrap complete"
}

main "$@"
