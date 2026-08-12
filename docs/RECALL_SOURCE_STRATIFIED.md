# Source-stratified recall metrics

Generated 2026-07-26 by `scripts/recall_audit/source_stratified_metrics.py`
over the four canonical ion-channel DBs listed in `docs/RECALL_STATUS.md`.
This is a frozen pre-correction diagnostic, not the current metric authority;
use `docs/RECALL_STATUS.md` for current caveats and regenerate this file after
the accepted explicit-zero rescore.
Regenerate with:

```bash
python scripts/run_recall_suite.py --score --review-gold-sync off \
  --genes KCNH2,KCNQ1,SCN5A,RYR2 --outdir recall_metrics/current
python scripts/recall_audit/source_stratified_metrics.py \
  --metrics-dir recall_metrics/current --out-md docs/RECALL_SOURCE_STRATIFIED.md
```

The ALL GOLD stratum reproduces `run_recall_suite.py` exactly (5546/6833 rows,
precision 77.2%, carriers MAE 0.613) — it re-reads that suite's per-row
`discrepancies.csv` rather than re-implementing the matcher, so the two cannot
drift.

**These numbers predate the linked-supplement recovery landing in the DBs.**
The 2026-07-25 sweep added 4,442 supplement files and refolded 249
FULL_CONTEXTs, but no DB has been refreshed to consume them, so this is the
clean "before" measurement.

ALL GOLD is the honest end-to-end number. SOURCE-COMPLETE restricts to
papers whose full text is on disk and whose advertised supplements were
all fetched — the ceiling the current corpus can support. The gap
between them is acquisition debt, not extraction error.


## Acquisition state of gold papers

| gene | gold papers | source-complete | no source | abstract-only | too short | missing suppl. |
|---|---:|---:|---:|---:|---:|---:|
| KCNH2 | 974 | 849 | 3 | 34 | 44 | 44 |
| KCNQ1 | 1316 | 1011 | 4 | 150 | 64 | 87 |
| RYR2 | 567 | 454 | 3 | 28 | 35 | 47 |
| SCN5A | 1250 | 893 | 8 | 139 | 131 | 79 |

## ALL GOLD — every gold row

| gene | rows | row recall | pmids | pmid recall | precision | carriers MAE / RMSE | affected MAE / RMSE | unaffected MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KCNH2 | 820/991 | 82.7% | 230/262 | 87.8% | 77.8% | 0.860 / 5.595 | 0.810 / 5.061 | 3.116 / 9.281 |
| KCNQ1 | 1499/1741 | 86.1% | 285/305 | 93.4% | 82.9% | 0.935 / 7.485 | 0.776 / 8.430 | 1.051 / 2.910 |
| RYR2 | 766/973 | 78.7% | 139/178 | 78.1% | 76.4% | 0.323 / 1.847 | 0.203 / 0.964 | 0.607 / 2.212 |
| SCN5A | 2461/3128 | 78.7% | 622/757 | 82.2% | 74.2% | 0.450 / 2.817 | 0.391 / 2.182 | 0.736 / 3.718 |

## SOURCE-COMPLETE — full text on disk, nothing missing

| gene | rows | row recall | pmids | pmid recall | precision | carriers MAE / RMSE | affected MAE / RMSE | unaffected MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KCNH2 | 732/879 | 83.3% | 185/213 | 86.9% | 78.8% | 0.703 / 4.363 | 0.588 / 2.690 | 3.447 / 9.866 |
| KCNQ1 | 1041/1156 | 90.1% | 206/220 | 93.6% | 78.2% | 0.363 / 2.299 | 0.410 / 3.480 | 0.809 / 2.523 |
| RYR2 | 484/644 | 75.2% | 119/144 | 82.6% | 71.0% | 0.457 / 2.312 | 0.269 / 1.146 | 0.708 / 2.389 |
| SCN5A | 1917/2274 | 84.3% | 467/517 | 90.3% | 83.1% | 0.488 / 3.034 | 0.401 / 2.196 | 0.773 / 3.888 |

## SOURCE-INCOMPLETE — the acquisition-limited remainder

| gene | rows | row recall | pmids | pmid recall | precision | carriers MAE / RMSE | affected MAE / RMSE | unaffected MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KCNH2 | 88/112 | 78.6% | 45/49 | 91.8% | 70.4% | 2.976 / 14.025 | 3.875 / 16.718 | 0.600 / 1.000 |
| KCNQ1 | 458/585 | 78.3% | 79/85 | 92.9% | 95.8% | 4.070 / 18.285 | 3.106 / 21.125 | 2.545 / 4.632 |
| RYR2 | 282/329 | 85.7% | 20/34 | 58.8% | 87.9% | 0.112 / 0.617 | 0.104 / 0.599 | 0.000 / 0.000 |
| SCN5A | 544/854 | 63.7% | 155/240 | 64.6% | 53.8% | 0.312 / 1.798 | 0.354 / 2.131 | 0.364 / 0.953 |

## All genes pooled

| stratum | rows | row recall | precision | carriers MAE / RMSE | affected MAE / RMSE | unaffected MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|
| ALL GOLD | 5546/6833 | 81.2% | 77.2% | 0.613 / 4.774 | 0.492 / 4.518 | 1.192 / 4.777 |
| SOURCE-COMPLETE | 4174/4953 | 84.3% | 79.5% | 0.492 / 3.083 | 0.416 / 2.510 | 1.200 / 4.972 |
| SOURCE-INCOMPLETE | 1372/1880 | 73.0% | 70.9% | 1.078 / 8.578 | 0.768 / 8.452 | 1.129 / 2.845 |
