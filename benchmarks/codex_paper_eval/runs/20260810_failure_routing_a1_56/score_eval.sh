#!/usr/bin/env bash
# Lock and score the completed A1 replay. The cardiac arm uses the production
# projection harness; BRCA2 is scored separately because its overrides are a
# diagnostic derived reference, not the manual cardiac gold standard.
set -euo pipefail

cd "$(dirname "$0")/../../../.."

RUN_ROOT="validation_runs/20260810_failure_routing_a1_56"
PRODUCTION_ROOT="$RUN_ROOT/runs"
BASE_SELECTION="validation_runs/20260808_fixed48_snapshot_replay_eval/selection.json"
PROJECTOR="benchmarks/codex_paper_eval/runs/20260726_fixed48_production/db_to_predictions.py"
EVALUATOR="benchmarks/codex_paper_eval/run_eval.py"

project_and_score() {
  local label="$1"
  local excluded_layers="$2"
  local eval_dir="$RUN_ROOT/$label"

  if [[ -e "$eval_dir/LOCK.json" ]]; then
    echo "projection already locked: $eval_dir" >&2
    return 2
  fi
  mkdir -p "$eval_dir"
  jq --arg run_id "20260810_failure_routing_a1_56_${label}" \
    '.run_id = $run_id' "$BASE_SELECTION" > "$eval_dir/selection.json"

  projection_args=(
    --run-dir "$eval_dir"
    --production-root "$PRODUCTION_ROOT"
    --out "$eval_dir/predictions.json"
  )
  if [[ -n "$excluded_layers" ]]; then
    projection_args+=(--exclude-layers "$excluded_layers")
  fi
  GVF_DISABLE_LOCAL_DATA=1 \
    .venv/bin/python "$PROJECTOR" "${projection_args[@]}"
  GVF_DISABLE_LOCAL_DATA=1 \
    .venv/bin/python "$EVALUATOR" lock --run-dir "$eval_dir"
  GVF_DISABLE_LOCAL_DATA=1 \
    .venv/bin/python "$EVALUATOR" score --run-dir "$eval_dir"
}

for gene in KCNH2 KCNQ1 RYR2 SCN5A BRCA2; do
  if ! find "$PRODUCTION_ROOT/$gene" -maxdepth 2 -type f -name "$gene.db" | grep -q .; then
    echo "missing completed $gene database under $PRODUCTION_ROOT" >&2
    exit 2
  fi
done

project_and_score "fixed48_all_layers" ""
project_and_score "fixed48_paper_only" "clinvar,pubtator"

brca2_db="$(
  find "$PRODUCTION_ROOT/BRCA2" -maxdepth 2 -type f -name 'BRCA2.db' \
    -print | sort | tail -1
)"
GVF_DISABLE_LOCAL_DATA=1 .venv/bin/python \
  benchmarks/curated_extraction_eval/run_benchmark.py \
  --mode score \
  --genes BRCA2 \
  --db "BRCA2=$brca2_db" \
  --outdir "$RUN_ROOT/brca2_diagnostic_score"
