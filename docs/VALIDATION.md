# Validation And Recall Metrics

GVF is optimized for high recall: finding as many true variants and carrier
counts as possible, then letting downstream review filter false positives.

## Current KCNH2 Baseline

Latest measured KCNH2 score: 2026-05-15, using the v12 manual-recovery SQLite
database and `gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv`.

| Metric | Matched / gold | Recall |
|---|---:|---:|
| PMIDs | 184 / 262 | 70.2% |
| Variant rows | 542 / 991 | 54.7% |
| Unique variants | 323 / 530 | 60.9% |
| Patients/carriers | 1758 / 2674 | 65.7% |
| Affected | 1095 / 1635 | 67.0% |
| Unaffected | 461 / 749 | 61.5% |

The older 59.1% KCNH2 result from the 2025-12-11 database is historical only.
Use the table above as the current project baseline.

## Gold Inputs

Normalized recall inputs live in:

```text
gene_variant_fetcher_gold_standard/normalized/
```

Current files include:

- `KCNH2_recall_input.csv`
- `KCNQ1_recall_input.csv`
- `SCN5A_recall_input.csv`
- `RYR2_recall_input.csv`

The first three are exported from Variant_Browser with
`scripts/build_gold_standard_from_varbrowser.py`. RYR2 is imported from the
lab-maintained spreadsheet with `scripts/build_ryr2_gold_standard_from_xlsx.py`.

## Metric Definitions

- PMIDs: fraction of gold PMIDs with at least one matched variant row.
- Variant rows: fraction of gold `(pmid, variant)` rows matched in SQLite.
- Unique variants: fraction of distinct gold variants matched anywhere.
- Patients/carriers: sum of gold carrier counts on matched rows divided by all
  gold carrier counts.
- Affected and unaffected: same as patient recall, but for phenotype-specific
  columns.

SQLite-only rows are precision findings and do not contribute to recall
denominators.

## Running Recall

Score an existing DB:

```bash
.venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
  --db KCNH2=results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db \
  --outdir recall_metrics/kcnh2_manual_recovery_after_matchfix_20260515
```

Score multiple genes once DBs exist:

```bash
.venv/bin/python scripts/run_recall_suite.py --score \
  --genes KCNH2,KCNQ1,SCN5A \
  --db KCNH2=path/to/KCNH2.db \
  --db KCNQ1=path/to/KCNQ1.db \
  --db SCN5A=path/to/SCN5A.db \
  --outdir recall_metrics/multigene_YYYYMMDD
```

Generated `recall_metrics/` outputs are local artifacts and are intentionally
gitignored. Preserve the headline numbers in `CLAUDE.md` and `TASKS.md`.

## Interpreting Gaps

Use the generated `missing_in_sqlite.csv` and `discrepancies.csv` files to split
misses into failure modes:

- Missing source: no usable full text or supplements.
- Extraction miss: source exists, but the LLM/table path missed the row.
- Migration/sanitizer drop: extraction JSON has the variant but SQLite does not.
- Matcher miss: SQLite has an equivalent variant, but comparison did not match
  it.
- Count mismatch: the variant matched, but carrier/phenotype counts differ.

Current KCNH2 status after the 2026-05-15 matcher patch:

- 449 missing variant rows across 107 PMIDs.
- 350 variant rows still needed to reach 90% variant-row recall.
- Top remaining source/extraction losses: 15840476, 14661677, 29650123,
  24667783, 16922724, 23098067, 23631430.
- 117 count mismatches remain; many top mismatches look like SQLite
  over-counting cohort tables and should be investigated separately from recall.

## Improving Recall

Highest-value actions:

1. Run paywall recovery from a Vanderbilt VPN or institutional network.
2. Add `ELSEVIER_INSTTOKEN` to unlock Elsevier subscription XML, especially PMID
   15840476.
3. Re-extract the highest-loss PMIDs that already have full-ish context:
   29650123, 24667783, 23098067, and 17905336.
4. Audit count mismatches for cohort-total bleed into per-variant carrier counts.
5. Build fresh KCNQ1 and SCN5A DBs so multi-gene recall can be measured.

## Related Docs

- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md)
- [ARCHITECTURE.md](ARCHITECTURE.md)
- [QUICKSTART.md](QUICKSTART.md)
