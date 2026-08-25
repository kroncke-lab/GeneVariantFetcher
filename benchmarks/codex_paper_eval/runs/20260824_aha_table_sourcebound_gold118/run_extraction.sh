#!/bin/zsh
set -euo pipefail
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
PYTHON=/Users/kronckbm/GitRepos/GeneVariantFetcher/.venv/bin/python
RUN_DIR=/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118
EMAIL=brett.kroncke@gmail.com
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/setup_production_eval.py check --run-dir "$RUN_DIR"
mkdir -p "$RUN_DIR/operator_logs"
typeset -a pids

echo "Starting calibrated KCNH2 extraction"
(
  "$PYTHON" -m cli gvf-run KCNH2 --email "$EMAIL" --output "$RUN_DIR/production_runs" --pmid-file "$RUN_DIR/pmids/KCNH2.txt" --no-source-recovery --no-corpus-sync --no-publish-review --gold-free-run 2>&1 | tee "$RUN_DIR/operator_logs/KCNH2.log"
) &
pids+=($!)

echo "Starting calibrated KCNQ1 extraction"
(
  "$PYTHON" -m cli gvf-run KCNQ1 --email "$EMAIL" --output "$RUN_DIR/production_runs" --pmid-file "$RUN_DIR/pmids/KCNQ1.txt" --no-source-recovery --no-corpus-sync --no-publish-review --gold-free-run 2>&1 | tee "$RUN_DIR/operator_logs/KCNQ1.log"
) &
pids+=($!)

echo "Starting calibrated RYR2 extraction"
(
  "$PYTHON" -m cli gvf-run RYR2 --email "$EMAIL" --output "$RUN_DIR/production_runs" --pmid-file "$RUN_DIR/pmids/RYR2.txt" --no-source-recovery --no-corpus-sync --no-publish-review --gold-free-run 2>&1 | tee "$RUN_DIR/operator_logs/RYR2.log"
) &
pids+=($!)

echo "Starting calibrated SCN5A extraction"
(
  "$PYTHON" -m cli gvf-run SCN5A --email "$EMAIL" --output "$RUN_DIR/production_runs" --pmid-file "$RUN_DIR/pmids/SCN5A.txt" --no-source-recovery --no-corpus-sync --no-publish-review --gold-free-run 2>&1 | tee "$RUN_DIR/operator_logs/SCN5A.log"
) &
pids+=($!)

failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    failed=1
  fi
done
if (( failed )); then
  echo "At least one gene extraction failed; do not lock or score." >&2
  exit 1
fi
