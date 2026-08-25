#!/usr/bin/env bash

# Run independent downstream RNA steps without stopping the complete workflow
# when a non-critical step fails. Patient-derived outputs must point to a
# restricted directory outside this public repository.

set +e

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${1:-$repo_root/config/config.example.yaml}"
config_env="$(mktemp)"
trap 'rm -f "$config_env"' EXIT

python3 "$repo_root/scripts/render_config_env.py" "$config_file" >"$config_env"
if [[ $? -ne 0 ]]; then
  echo "Configuration validation failed: $config_file" >&2
  exit 2
fi
# shellcheck disable=SC1090
source "$config_env"

if [[ "$PROGRAMME_GMT" != /* ]]; then
  PROGRAMME_GMT="$repo_root/$PROGRAMME_GMT"
fi
if [[ -n "$MEMORY_GATE" && "$MEMORY_GATE" != /* ]]; then
  MEMORY_GATE="$repo_root/$MEMORY_GATE"
fi

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"
status_file="$LOG_DIR/run_all_status.tsv"
printf 'step\tstatus\texit_code\tlog\n' >"$status_file"

run_step() {
  local step="$1"
  shift
  local log="$LOG_DIR/${step}.log"
  echo "[$(date --iso-8601=seconds)] starting $step" | tee "$log"
  "$@" >>"$log" 2>&1
  local exit_code=$?
  if [[ $exit_code -eq 0 ]]; then
    printf '%s\tsuccess\t0\t%s\n' "$step" "$log" >>"$status_file"
  else
    printf '%s\tfailure\t%d\t%s\n' "$step" "$exit_code" "$log" >>"$status_file"
  fi
  return 0
}

run_memory_gate() {
  if [[ -z "$MEMORY_GATE" ]]; then
    return 0
  fi
  MINIMUM_AVAILABLE_GB="$MINIMUM_AVAILABLE_GB" \
    MAXIMUM_SWAP_PERCENT="$MAXIMUM_SWAP_PERCENT" \
    "$MEMORY_GATE"
}

if [[ "$RUN_STANDARD_DE" == "true" ]]; then
  if run_memory_gate; then
    if [[ -n "$PHENOTYPE_ASSIGNMENT" ]]; then
      run_step standard_de Rscript \
        "$repo_root/pipelines/rnaseq/run_standard_de_from_star_readspergene.R" \
        "$STAR_ANALYSIS_DIR" "$OUTPUT_ROOT/standard_de" "$PHENOTYPE_ASSIGNMENT" "$GTF"
    else
      run_step standard_de Rscript \
        "$repo_root/pipelines/rnaseq/run_standard_de_from_star_readspergene.R" \
        "$STAR_ANALYSIS_DIR" "$OUTPUT_ROOT/standard_de" - "$GTF"
    fi
  else
    printf 'standard_de\tblocked_by_memory_gate\t99\t%s\n' "$LOG_DIR/memory_gate.log" >>"$status_file"
  fi
fi

counts="$OUTPUT_ROOT/standard_de/star_unstranded_gene_counts_matrix.tsv"
normalized="$OUTPUT_ROOT/standard_de/DESeq2_normalized_counts.tsv"
tpm_prefix="$OUTPUT_ROOT/gene_tpm/rna_expression"
programme_prefix="$OUTPUT_ROOT/pathway_scores/pdac_programmes"
metadata_for_scores="$METADATA"
if [[ -z "$metadata_for_scores" ]]; then
  metadata_for_scores="$OUTPUT_ROOT/standard_de/rnaseq_metadata_for_standard_DE.tsv"
fi

if [[ "$RUN_GENE_TPM" == "true" ]]; then
  if run_memory_gate; then
    run_step gene_tpm Rscript \
      "$repo_root/pipelines/rnaseq/make_gene_tpm_from_counts.R" \
      "$counts" "$GTF" "$tpm_prefix"
  else
    printf 'gene_tpm\tblocked_by_memory_gate\t99\t%s\n' "$LOG_DIR/memory_gate.log" >>"$status_file"
  fi
fi

if [[ "$RUN_PATHWAY_SCORING" == "true" ]]; then
  pathway_command=(
    Rscript "$repo_root/pipelines/pathway_scoring/run_gsva_ssgsea_programme_scores.R"
    --expression "$normalized"
    --gmt "$PROGRAMME_GMT"
    --gtf "$GTF"
    --gene-column gene_id
    --transform auto
    --out-prefix "$programme_prefix"
  )
  if [[ -n "$metadata_for_scores" ]]; then
    pathway_command+=(--metadata "$metadata_for_scores")
  fi
  if run_memory_gate; then
    run_step pathway_scoring "${pathway_command[@]}"
  else
    printf 'pathway_scoring\tblocked_by_memory_gate\t99\t%s\n' "$LOG_DIR/memory_gate.log" >>"$status_file"
  fi
fi

if [[ "$RUN_IMMUNE_DECONVOLUTION" == "true" ]]; then
  immune_command=(
    Rscript "$repo_root/pipelines/immune_infiltration/run_immune_stromal_scores_from_expression.R"
    --expression "$tpm_prefix.gene_tpm.tsv"
    --gene-column gene_symbol
    --input-scale linear
    --methods "$IMMUNE_METHODS"
    --out-dir "$OUTPUT_ROOT/immune_deconvolution"
  )
  if [[ -n "$CIBERSORT_SCRIPT" && -n "$CIBERSORT_LM22" ]]; then
    immune_command+=(
      --cibersort-script "$CIBERSORT_SCRIPT"
      --cibersort-lm22 "$CIBERSORT_LM22"
      --cibersort-permutations "$CIBERSORT_PERMUTATIONS"
    )
  fi
  if run_memory_gate; then
    run_step immune_deconvolution "${immune_command[@]}"
  else
    printf 'immune_deconvolution\tblocked_by_memory_gate\t99\t%s\n' "$LOG_DIR/memory_gate.log" >>"$status_file"
  fi
fi

if [[ "$RUN_PHENOTYPE_ASSIGNMENT" == "true" ]]; then
  if [[ -z "$TME_SCORE_TABLE" || -z "$PHENOTYPE_IMMUNE_COLUMNS" || -z "$PHENOTYPE_STROMAL_COLUMNS" || -z "$PHENOTYPE_EMT_COLUMNS" ]]; then
    printf 'phenotype_assignment\tmissing_explicit_feature_configuration\t98\t%s\n' \
      "$repo_root/pipelines/phenotype_assignment/README.md" >>"$status_file"
  elif run_memory_gate; then
    run_step phenotype_assignment python3 \
      "$repo_root/pipelines/phenotype_assignment/assign_tme_phenotype_groups.py" \
      --scores "$TME_SCORE_TABLE" \
      --sample-column sample_id \
      --patient-column patient_id \
      --immune-columns "$PHENOTYPE_IMMUNE_COLUMNS" \
      --stromal-columns "$PHENOTYPE_STROMAL_COLUMNS" \
      --emt-columns "$PHENOTYPE_EMT_COLUMNS" \
      --target-per-extreme "$PHENOTYPE_TARGET_PER_EXTREME" \
      --out-prefix "$OUTPUT_ROOT/phenotype_assignment/tme_phenotype"
  else
    printf 'phenotype_assignment\tblocked_by_memory_gate\t99\t%s\n' "$LOG_DIR/memory_gate.log" >>"$status_file"
  fi
fi

echo "Workflow status: $status_file"
cat "$status_file"
