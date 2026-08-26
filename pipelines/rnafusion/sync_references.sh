#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <destination> <log_directory> <manifest_tsv>" >&2
  exit 2
fi

destination=$1
log_dir=$2
manifest=$3
source_uri=s3://nf-core-awsmegatests/rnafusion/references/
mkdir -p "$destination" "$log_dir" "$(dirname "$manifest")"

exec 9>"$log_dir/reference_sync.lock"
flock -n 9 || { echo "Another reference sync is active."; exit 0; }

timestamp=$(date +%Y%m%d_%H%M%S)
log_file="$log_dir/reference_sync_${timestamp}.log"
inventory_file="$log_dir/reference_s3_inventory_${timestamp}.txt"
exec > >(tee -a "$log_file") 2>&1

command -v aws >/dev/null || { echo "AWS CLI is required." >&2; exit 1; }
echo "Started: $(date --iso-8601=seconds)"
echo "Source: $source_uri"
echo "Destination: $destination"
df -h "$destination"

aws --no-sign-request s3 ls "$source_uri" --recursive --summarize > "$inventory_file"
read -r expected_objects expected_bytes < <(
  awk 'NF >= 4 && $1 ~ /^[0-9][0-9][0-9][0-9]-/ && $4 !~ /\/$/ {n++; b += $3} END {printf "%d %.0f\n", n, b}' "$inventory_file"
)
echo "Expected files: $expected_objects"
echo "Expected bytes: $expected_bytes"

# aws s3 sync is idempotent and resumes by skipping matching files.
aws --no-sign-request s3 sync "$source_uri" "$destination" --only-show-errors
local_objects=$(find "$destination" -type f | wc -l)
local_bytes=$(find "$destination" -type f -printf '%s\n' | awk '{b += $1} END {printf "%.0f", b}')

status=verified
if [[ $local_objects -ne $expected_objects || $local_bytes -ne $expected_bytes ]]; then
  status=mismatch
fi

printf 'reference\tpipeline_version\tsource\tdestination\tremote_objects\tremote_bytes\tlocal_files\tlocal_bytes\tstatus\tchecked_at\n' > "${manifest}.tmp"
printf 'nf-core/rnafusion references\t4.1.0\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$source_uri" "$destination" "$expected_objects" "$expected_bytes" \
  "$local_objects" "$local_bytes" "$status" "$(date --iso-8601=seconds)" >> "${manifest}.tmp"
mv "${manifest}.tmp" "$manifest"
echo "Verification: $status ($local_objects files; $local_bytes bytes)"
[[ $status == verified ]]
