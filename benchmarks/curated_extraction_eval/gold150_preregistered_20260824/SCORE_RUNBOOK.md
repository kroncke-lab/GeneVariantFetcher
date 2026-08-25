# Score runbook

For a returned packet, place `curation_template_FILLED.csv` beside its
`packet_meta.json`, then run:

```bash
.venv/bin/python scripts/score_curation_packet.py \
  --filled-csv benchmarks/curated_extraction_eval/gold150_preregistered_20260824/BRCA1/calibration/curation_template_FILLED.csv \
  --db results/vb_curated150_20260824/BRCA1/20260824_181400/BRCA1.db \
  --gene BRCA1 \
  --out recall_metrics/gold150_20260824/BRCA1_calibration
```

Repeat for BRCA2 and BMPR2. Do not score a holdout packet until the candidate
commit/configuration is frozen. The scorer writes an explicit curated-PMID scope
and invokes `run_recall_suite.py` with `--gold-paper-exhaustive`,
`--review-gold-sync off`, `--review-gold-tier all`, and hermetic metric-only
artifacts.
