#!/bin/zsh
set -euo pipefail
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
PYTHON=/Users/kronckbm/GitRepos/GeneVariantFetcher/.venv/bin/python
RUN_DIR=/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260905_mechanism4_final
"$PYTHON" -c 'from pathlib import Path; import json,sys; from benchmarks.codex_paper_eval.setup_production_eval import runtime_fingerprint; assert runtime_fingerprint()==json.loads((Path(sys.argv[1])/"analysis_setup.json").read_text())["runtime"], "runtime drift"' "$RUN_DIR"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/rebind_production_sources.py \
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/db_to_predictions.py \
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs" \
  --trust-mode trusted --identity-mode trusted --paper-primary \
  --out "$RUN_DIR/predictions.json"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/run_eval.py lock --run-dir "$RUN_DIR"
"$PYTHON" /Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/run_eval.py score --run-dir "$RUN_DIR" \
  --gold-root /Users/kronckbm/GitRepos/GeneVariantFetcher/gene_variant_fetcher_gold_standard/normalized
