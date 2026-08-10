#!/usr/bin/env bash
# Failure-routing A1 production replay: fixed 48 cardiac papers plus the eight
# BRCA2 diagnostic papers. This intentionally changes only the risk-to-verifier
# opening predicate in pipeline/extraction.py.
set -euo pipefail

cd "$(dirname "$0")/../../../.."

RUN_ROOT="validation_runs/20260810_failure_routing_a1_56"
OUT_ROOT="$RUN_ROOT/runs"
SEED_ROOT="$OUT_ROOT/_source_seed"
LOG_ROOT="$RUN_ROOT/logs"
FIXED_SELECTION="validation_runs/20260808_fixed48_snapshot_replay_eval/selection.json"
FIXED_PMIDS="benchmarks/codex_paper_eval/runs/20260726_fixed48_production/pmids"
BRCA2_PMIDS="benchmarks/curated_extraction_eval/pmids/BRCA2.txt"
LOCK_FILE="$RUN_ROOT/source_lock.tsv"

mkdir -p "$SEED_ROOT" "$LOG_ROOT"

copy_source() {
  local gene="$1"
  local pmid="$2"
  local source_path="$3"
  local source_dir
  local destination="$SEED_ROOT/$gene/pmc_fulltext"

  source_dir="$(dirname "$source_path")"
  mkdir -p "$destination"
  cp -a "$source_path" "$destination/${pmid}_FULL_CONTEXT.md"
  for suffix in _figures _supplements; do
    if [[ -d "$source_dir/${pmid}${suffix}" ]]; then
      cp -a "$source_dir/${pmid}${suffix}" "$destination/"
    fi
  done
}

for gene in KCNH2 KCNQ1 RYR2 SCN5A; do
  while IFS= read -r pmid; do
    [[ -n "$pmid" ]] || continue
    source_path="$(
      jq -r --arg gene "$gene" --arg pmid "$pmid" \
        '.papers[] | select(.gene == $gene and .pmid == $pmid) | .source' \
        "$FIXED_SELECTION"
    )"
    if [[ ! -f "$source_path" ]]; then
      echo "locked source missing for $gene PMID $pmid: $source_path" >&2
      exit 2
    fi
    copy_source "$gene" "$pmid" "$source_path"
  done < "$FIXED_PMIDS/$gene.txt"
done

while IFS= read -r pmid; do
  [[ -n "$pmid" ]] || continue
  copy_source \
    "BRCA2" \
    "$pmid" \
    "corpus/BRCA2/$pmid/${pmid}_FULL_CONTEXT.md"
done < "$BRCA2_PMIDS"

printf 'gene\tpmid\tsha256\ttier\tsource_path\n' > "$LOCK_FILE"
for gene in KCNH2 KCNQ1 RYR2 SCN5A; do
  while IFS= read -r pmid; do
    [[ -n "$pmid" ]] || continue
    source_path="$SEED_ROOT/$gene/pmc_fulltext/${pmid}_FULL_CONTEXT.md"
    actual_sha="$(shasum -a 256 "$source_path" | awk '{print $1}')"
    expected_sha="$(
      jq -r --arg gene "$gene" --arg pmid "$pmid" \
        '.papers[] | select(.gene == $gene and .pmid == $pmid) | .source_sha256' \
        "$FIXED_SELECTION"
    )"
    if [[ -z "$expected_sha" || "$expected_sha" == "null" || "$actual_sha" != "$expected_sha" ]]; then
      echo "source lock mismatch for $gene PMID $pmid" >&2
      exit 2
    fi
    printf '%s\t%s\t%s\tmanual_cardiac_gold\t%s\n' \
      "$gene" "$pmid" "$actual_sha" "$source_path" >> "$LOCK_FILE"
  done < "$FIXED_PMIDS/$gene.txt"
done

while IFS= read -r pmid; do
  [[ -n "$pmid" ]] || continue
  source_path="$SEED_ROOT/BRCA2/pmc_fulltext/${pmid}_FULL_CONTEXT.md"
  source_sha="$(shasum -a 256 "$source_path" | awk '{print $1}')"
  printf 'BRCA2\t%s\t%s\tdiagnostic_derived_override\t%s\n' \
    "$pmid" "$source_sha" "$source_path" >> "$LOCK_FILE"
done < "$BRCA2_PMIDS"

if [[ "$(($(wc -l < "$LOCK_FILE") - 1))" -ne 56 ]]; then
  echo "source lock does not contain exactly 56 gene-paper entries" >&2
  exit 2
fi

set -a
source ./.env
set +a

run_gene() {
  local gene="$1"
  local pmid_file="$2"
  ENABLE_TIER3_REASON_CLASS_ROUTING=1 GVF_DISABLE_LOCAL_DATA=1 MAX_WORKERS=1 \
    .venv/bin/python -m cli gvf-run "$gene" \
      --email brett.kroncke@gmail.com \
      --pmid-file "$pmid_file" \
      --output "$OUT_ROOT" \
      --no-source-recovery \
      --no-corpus-sync \
      --skip doctor \
      --skip report \
      > "$LOG_ROOT/${gene}.log" 2>&1
}

pids=()
genes=()
for gene in KCNH2 KCNQ1 RYR2 SCN5A; do
  run_gene "$gene" "$FIXED_PMIDS/$gene.txt" &
  pids+=("$!")
  genes+=("$gene")
done
run_gene "BRCA2" "$BRCA2_PMIDS" &
pids+=("$!")
genes+=("BRCA2")

status=0
for index in "${!pids[@]}"; do
  if ! wait "${pids[$index]}"; then
    echo "${genes[$index]} failed; see $LOG_ROOT/${genes[$index]}.log" >&2
    status=1
  fi
done
exit "$status"
