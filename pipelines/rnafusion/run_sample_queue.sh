#!/usr/bin/env bash
set -Eeuo pipefail

# Required paths are passed through environment variables so this reusable code
# never embeds private project locations or sample identifiers.
: "${RNAFUSION_PIPELINE_ROOT:?Set RNAFUSION_PIPELINE_ROOT to the pinned 4.1.0 clone}"
: "${RNAFUSION_REFERENCES:?Set RNAFUSION_REFERENCES to the verified reference root}"
: "${RNAFUSION_REFERENCE_MANIFEST:?Set RNAFUSION_REFERENCE_MANIFEST}"
: "${RNAFUSION_SAMPLESHEET:?Set RNAFUSION_SAMPLESHEET}"
: "${RNAFUSION_RESULTS_ROOT:?Set RNAFUSION_RESULTS_ROOT}"
: "${RNAFUSION_WORK_ROOT:?Set RNAFUSION_WORK_ROOT}"
: "${RNAFUSION_STATE_ROOT:?Set RNAFUSION_STATE_ROOT}"
: "${RNAFUSION_LOG_ROOT:?Set RNAFUSION_LOG_ROOT}"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
tools=${RNAFUSION_TOOLS:-arriba,fusioncatcher,starfusion,salmon}
min_available_kib=$((${RNAFUSION_MIN_AVAILABLE_GIB:-24} * 1024 * 1024))
min_disk_kib=$((${RNAFUSION_MIN_DISK_GIB:-150} * 1024 * 1024))
gate_sleep=${RNAFUSION_GATE_SLEEP_SECONDS:-300}
status_file="$RNAFUSION_LOG_ROOT/run_queue_status.tsv"

mkdir -p "$RNAFUSION_RESULTS_ROOT" "$RNAFUSION_WORK_ROOT" \
  "$RNAFUSION_STATE_ROOT" "$RNAFUSION_LOG_ROOT"
exec 9>"$RNAFUSION_LOG_ROOT/run_queue.lock"
flock -n 9 || { echo "Another RNA-fusion queue is active."; exit 0; }

[[ -s $RNAFUSION_PIPELINE_ROOT/main.nf ]] || { echo "Pinned pipeline is missing." >&2; exit 1; }
[[ -s $RNAFUSION_SAMPLESHEET ]] || { echo "Samplesheet is missing." >&2; exit 1; }
[[ -s $RNAFUSION_REFERENCE_MANIFEST ]] || { echo "Reference manifest is missing." >&2; exit 1; }
awk -F '\t' 'NR == 2 && $9 == "verified" {ok=1} END {exit !ok}' "$RNAFUSION_REFERENCE_MANIFEST" || {
  echo "Reference manifest does not report verified." >&2
  exit 1
}

if [[ ! -e $status_file ]]; then
  printf 'sample\tstatus\tstarted_at\tfinished_at\texit_code\twork_cleanup\toutdir\tlog\n' > "$status_file"
fi

wait_for_capacity() {
  local sample=$1
  while true; do
    local available disk_free
    available=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
    disk_free=$(df -Pk "$RNAFUSION_RESULTS_ROOT" | awk 'NR == 2 {print $4}')
    if (( available >= min_available_kib && disk_free >= min_disk_kib )); then
      echo "Capacity gate passed for $sample: MemAvailable=$available KiB, disk=$disk_free KiB"
      return
    fi
    echo "Waiting before $sample: MemAvailable=$available KiB, disk=$disk_free KiB"
    sleep "$gate_sleep"
  done
}

consecutive_failures=0
while IFS=, read -r sample fastq_1 fastq_2 strandedness; do
  [[ $sample == sample ]] && continue
  [[ -n $sample ]] || continue
  outdir="$RNAFUSION_RESULTS_ROOT/$sample"
  sample_work="$RNAFUSION_WORK_ROOT/$sample"
  sample_state="$RNAFUSION_STATE_ROOT/$sample"
  sample_sheet="$sample_state/samplesheet.csv"
  run_name="rnafusion_${sample}"
  run_log="$RNAFUSION_LOG_ROOT/${sample}.nextflow.log"

  if [[ -e $outdir/.completed ]]; then
    echo "Skipping completed sample $sample"
    continue
  fi
  wait_for_capacity "$sample"
  mkdir -p "$outdir/pipeline_info" "$sample_work" "$sample_state"
  printf 'sample,fastq_1,fastq_2,strandedness\n%s,%s,%s,%s\n' \
    "$sample" "$fastq_1" "$fastq_2" "$strandedness" > "$sample_sheet"
  chmod 600 "$sample_sheet"

  started_at=$(date --iso-8601=seconds)
  set +e
  (
    cd "$sample_state"
    "$SCRIPT_DIR/run_rnafusion.sh" \
      --pipeline-root "$RNAFUSION_PIPELINE_ROOT" \
      --references "$RNAFUSION_REFERENCES" \
      --samplesheet "$sample_sheet" \
      --outdir "$outdir" \
      --work-dir "$sample_work" \
      --tools "$tools" \
      -resume \
      -name "$run_name"
  ) > >(tee -a "$run_log") 2>&1
  exit_code=$?
  set -e
  finished_at=$(date --iso-8601=seconds)
  cleanup=retained

  if (( exit_code == 0 )) && [[ -s $outdir/pipeline_info/execution_report.html ]]; then
    touch "$outdir/.completed"
    if (cd "$sample_state" && nextflow clean "$run_name" -f); then
      cleanup=nextflow_cleaned
    fi
    status=success
    consecutive_failures=0
  else
    status=failed
    consecutive_failures=$((consecutive_failures + 1))
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$status" "$started_at" "$finished_at" "$exit_code" \
    "$cleanup" "$outdir" "$run_log" >> "$status_file"
  if (( consecutive_failures >= 2 )); then
    echo "Stopping after two consecutive failures; failed work is retained for diagnosis." >&2
    exit 1
  fi
done < "$RNAFUSION_SAMPLESHEET"

expected=$(($(wc -l < "$RNAFUSION_SAMPLESHEET") - 1))
completed=$(find "$RNAFUSION_RESULTS_ROOT" -mindepth 2 -maxdepth 2 -name .completed | wc -l)
echo "RNA-fusion queue finished: $completed/$expected samples complete."
(( completed == expected ))
