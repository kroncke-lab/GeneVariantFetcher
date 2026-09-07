# Live integration check: three compendium papers, 2026-09-07

Tranche 04 never exercised a stamped value through the production persistence
path (its papers hold no compendium-style count table), so after the post-lock
fix (commit `793a456f`) a calibrated `gvf-run` was executed on the three
compendium papers the calibration replay had used, with the projection on and
the same calibrated-run safeguards as a registry arm (`--pmid-file`,
`--no-source-recovery`, `--no-corpus-sync`, `--no-publish-review`,
`--gold-free-run`, Azure only). Output: `results/table_cohort_live_check_20260907/`
(ignored scratch, this note is the record). These are opened calibration
papers: the check proves the persistence path, not generalisation.

## What the pipeline wrote

| gene / PMID | extraction | lane result | SQLite `penetrance_data` |
| --- | --- | --- | --- |
| SCN5A 20129283 | deterministic table parser, 363 rows | Table 4 case series: 293 applied; Table 2 control alleles: left as parsed | 293 rows `affected_count` set, all `affected_role = case`, 363/363 `field_trust.affected = trusted`, 361/363 rows trusted |
| SCN5A 27566755 | fixed-width parser, 51 rows | Table S1 case series: 51 applied | 51/51 affected set, role `case`, all trusted |
| KCNH2 26496715 | router + Kimi, 54 rows | Table 2 case series: 54 applied | 54/54 affected set, role `case`, all trusted |
| KCNQ1 26496715 | router + grok-4.3, verified by gpt-5.6-sol, 46 rows | Table 1 case series: 35 applied, 11 parser copies stamped | 46/46 affected set, role `case`, all trusted |

So the stamp survives `migrate_to_sqlite` (role vocabulary `case`), the Step 3.7
trust gate (no reason code fires on an `affected = N, unaffected = NULL` row
whose carriers are a labelled per-variant cell) and the trusted projection.

## Scored against the normalized gold (paper-derived rows only)

| gene / PMID | identity TP/FP/FN | affected on positive gold: supplied / exact / wrong | carriers wrong |
| --- | --- | --- | --- |
| SCN5A 20129283 | 343 / 9 / 74 | 275 / 231 / 44 | 46 |
| SCN5A 27566755 | 51 / 0 / 0 | 51 / 51 / 0 | 0 |
| KCNH2 26496715 | 53 / 0 / 1 | 53 / 53 / 0 | 0 |
| KCNQ1 26496715 | 41 / 5 / 1 | 41 / 41 / 0 | 0 |

The 44 wrong affected values on 20129283 are the per-testing-centre gold rows
already documented in the README (the 46 wrong carriers are the same rows plus
two gold rows whose counts do not sum to the table either way). The live numbers
match the zero-LLM replay, which is the point: the replay is a faithful
instrument for this lane.

Cost: three single-paper gene runs, a few minutes each; the two deterministic
papers made no model calls for identity, and the KCNH2 supplement-only paper
went through the router with Kimi. Not separately priced; well under $1 by the
repository proxy.
