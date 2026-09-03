#!/bin/zsh
set -euo pipefail
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
PYTHON=/Users/kronckbm/GitRepos/GeneVariantFetcher/.venv/bin/python
RUN_DIR=/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260903_protocol_cont120_01_candidate
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/setup_production_eval.py check --run-dir "$RUN_DIR"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/rebind_production_sources.py \
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/db_to_predictions.py \
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs" \
  --trust-mode trusted --identity-mode trusted --paper-primary \
  --out "$RUN_DIR/predictions.json"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/run_eval.py lock --run-dir "$RUN_DIR"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/run_eval.py score --run-dir "$RUN_DIR" \
  --gold-root /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/evaluation_tiers/mixed_gold_continuation_120/answer_key
