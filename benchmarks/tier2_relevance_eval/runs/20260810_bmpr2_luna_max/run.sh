#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

GVF_LUNA_SHADOW_OUTDIR=${GVF_LUNA_SHADOW_OUTDIR:-validation_runs/20260810_bmpr2_tier2_luna_max_replay}

.venv/bin/python benchmarks/tier2_relevance_eval/run_shadow.py \
  --positive-count 50 \
  --negative-count 50 \
  --boundary-count 50 \
  --workers 12 \
  --outdir "$GVF_LUNA_SHADOW_OUTDIR"
