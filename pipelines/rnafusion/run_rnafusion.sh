#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
pipeline_root=""
references=""
samplesheet=""
outdir=""
workdir=""
tools=arriba,fusioncatcher,starfusion,salmon
resume=()
extra_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pipeline-root) pipeline_root=$2; shift 2 ;;
    --references) references=$2; shift 2 ;;
    --samplesheet) samplesheet=$2; shift 2 ;;
    --outdir) outdir=$2; shift 2 ;;
    --work-dir) workdir=$2; shift 2 ;;
    --tools) tools=$2; shift 2 ;;
    -resume|--resume) resume=(-resume); shift ;;
    *) extra_args+=("$1"); shift ;;
  esac
done

for value in pipeline_root references samplesheet outdir workdir; do
  [[ -n ${!value} ]] || { echo "Missing required option for $value" >&2; exit 2; }
done
[[ -s $pipeline_root/main.nf ]] || { echo "Missing pinned pipeline: $pipeline_root/main.nf" >&2; exit 1; }
[[ -d $references ]] || { echo "Missing references: $references" >&2; exit 1; }
[[ -s $samplesheet ]] || { echo "Missing samplesheet: $samplesheet" >&2; exit 1; }
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline rnafusion --input "$samplesheet"

mkdir -p "$outdir/pipeline_info" "$workdir"
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-$HOME/.cache/nextflow/singularity}
export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
export NXF_OPTS=${NXF_OPTS:--Xms256m -Xmx2g}

nextflow run "$pipeline_root/main.nf" \
  -profile singularity \
  -c "$SCRIPT_DIR/rnafusion.config" \
  --input "$samplesheet" \
  --outdir "$outdir" \
  --genomes_base "$references" \
  --genome GRCh38 \
  --genome_gencode_version 46 \
  --tools "$tools" \
  --no_cosmic \
  --cram \
  --star_limit_bam_sort_ram 2000000000 \
  -w "$workdir" \
  "${resume[@]}" \
  -with-report "$outdir/pipeline_info/execution_report.html" \
  -with-timeline "$outdir/pipeline_info/execution_timeline.html" \
  -with-trace "$outdir/pipeline_info/execution_trace.tsv" \
  -with-dag "$outdir/pipeline_info/execution_dag.html" \
  "${extra_args[@]}"
