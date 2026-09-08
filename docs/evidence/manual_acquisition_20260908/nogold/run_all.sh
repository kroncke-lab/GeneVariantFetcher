#!/bin/zsh
set -u
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
set -a; . ./.env; set +a
for g in KCNH2 KCNQ1 SCN5A RYR2; do
  echo "== $g =="
  GVF_DISABLE_LOCAL_DATA=1 .venv/bin/python scripts/recall_audit/rank_manual_acquisition.py --pmid-file docs/evidence/manual_acquisition_20260908/nogold/${g}_nongold_stubs.txt --gene $g --out-dir docs/evidence/manual_acquisition_20260908/nogold/$g --max-papers 60 --email brett.kroncke@gmail.com
done
echo "ALL DONE"
