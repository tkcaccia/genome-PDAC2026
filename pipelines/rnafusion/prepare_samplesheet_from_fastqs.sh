#!/usr/bin/env bash
set -Eeuo pipefail

# Build an nf-core/rnafusion samplesheet from one directory per sample. Existing
# STAR ReadsPerGene output is used only to determine library strandedness.
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <fastq_root> <star_count_root> <output_csv>" >&2
  exit 2
fi

fastq_root=$1
star_count_root=$2
output_csv=$3
check_file="${output_csv%.csv}.input_check.tsv"
tmp_csv="${output_csv}.tmp"
tmp_check="${check_file}.tmp"

[[ -d $fastq_root ]] || { echo "Missing FASTQ root: $fastq_root" >&2; exit 1; }
[[ -d $star_count_root ]] || { echo "Missing STAR-count root: $star_count_root" >&2; exit 1; }
mkdir -p "$(dirname "$output_csv")"

printf 'sample,fastq_1,fastq_2,strandedness\n' > "$tmp_csv"
printf 'sample\traw_r1_count\traw_r2_count\tunstranded_gene_counts\tforward_gene_counts\treverse_gene_counts\tinferred_strandedness\tstatus\n' > "$tmp_check"

sample_count=0
for sample_dir in "$fastq_root"/*; do
  [[ -d $sample_dir ]] || continue
  sample=$(basename "$sample_dir")
  if [[ $sample == *','* ]]; then
    echo "Sample names containing commas are not supported: $sample" >&2
    exit 1
  fi

  # The exact suffix excludes separately generated trimmed copies.
  mapfile -t r1_files < <(find "$sample_dir" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' -print | sort)
  mapfile -t r2_files < <(find "$sample_dir" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' -print | sort)
  count_file=$(find "$star_count_root/$sample" -maxdepth 1 -type f -name '*ReadsPerGene.out.tab' -print -quit 2>/dev/null || true)

  status=PASS
  inferred=unknown
  unstranded=0
  forward=0
  reverse=0
  if [[ -n $count_file ]]; then
    read -r unstranded forward reverse < <(
      awk 'NR > 4 {u += $2; f += $3; r += $4} END {printf "%.0f %.0f %.0f\n", u, f, r}' "$count_file"
    )
    if (( reverse > forward )); then
      inferred=reverse
    else
      inferred=forward
    fi
  else
    status=MISSING_STAR_COUNTS
  fi

  if (( ${#r1_files[@]} != 1 || ${#r2_files[@]} != 1 )); then
    status=INVALID_FASTQ_PAIR_COUNT
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "${#r1_files[@]}" "${#r2_files[@]}" "$unstranded" \
    "$forward" "$reverse" "$inferred" "$status" >> "$tmp_check"

  if [[ $status == PASS ]]; then
    printf '%s,%s,%s,%s\n' "$sample" "${r1_files[0]}" "${r2_files[0]}" "$inferred" >> "$tmp_csv"
    sample_count=$((sample_count + 1))
  fi
done

mv "$tmp_check" "$check_file"
if awk -F '\t' 'NR > 1 && $8 != "PASS" {bad=1} END {exit !bad}' "$check_file"; then
  rm -f "$tmp_csv"
  echo "Input validation failed; see $check_file" >&2
  exit 1
fi
if (( sample_count == 0 )); then
  rm -f "$tmp_csv"
  echo "No valid samples were found; see $check_file" >&2
  exit 1
fi

mv "$tmp_csv" "$output_csv"
chmod 600 "$output_csv" "$check_file"
echo "Prepared $sample_count samples: $output_csv"
echo "Validation evidence: $check_file"
