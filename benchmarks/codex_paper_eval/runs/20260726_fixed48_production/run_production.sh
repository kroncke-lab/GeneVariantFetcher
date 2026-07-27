#!/usr/bin/env bash
# Production-strategy run over the locked 48-paper cardiac eval set.
# Source parity with the sol run: reads corpus/ only, no acquisition, no corpus mutation.
set -uo pipefail
cd "$(dirname "$0")/../../../.."          # repo root
RUN="benchmarks/codex_paper_eval/runs/20260726_fixed48_production"
set -a; . ./.env; set +a
GENES="${GENES:-KCNH2 KCNQ1 RYR2 SCN5A}"
WORKERS="${WORKERS:-6}"
PMDIR="${PMDIR:-$RUN/pmids}"
OUT="${OUT:-$RUN/production_runs}"
mkdir -p "$OUT"
for G in $GENES; do
  echo "=== $(date -u +%FT%TZ) starting $G ($(wc -l < "$PMDIR/$G.txt") pmids)"
  .venv/bin/python -m cli gvf-run "$G" \
    --email brett.kroncke@gmail.com \
    --pmid-file "$PMDIR/$G.txt" \
    --output "$OUT" \
    --no-source-recovery \
    --no-corpus-sync \
    --skip doctor --skip report \
    --extraction-workers "$WORKERS" \
    2>&1 | tee "$OUT/${G}_run.log"
  echo "=== $(date -u +%FT%TZ) finished $G (exit ${PIPESTATUS[0]})"
done
echo "=== ALL DONE $(date -u +%FT%TZ)"
