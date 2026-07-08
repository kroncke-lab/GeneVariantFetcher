# Recall Status

Last updated: 2026-06-12.

This file is the current measured recall snapshot. It intentionally does not
carry the active work plan or dated session log.

- Active forward plan and task checklist: `TASKS.md`.
- Benchmark/change trajectory: `docs/RECALL_HISTORY.md`.
- Re-run procedure: `docs/RECALL_REFRESH_RUNBOOK.md`.

No other doc should restate live recall tables. If a metric conflicts with this
file, this file is authoritative.

## Current Canonical Baseline

Fresh run of `scripts/run_recall_suite.py` against the four canonical DBs after
the 2026-06-12 PDF-linearized table reconstruction, iter-2 quality-aware
gate/selector, and targeted lands of all four genes:

- `results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db`
- `validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db`

Scored artifact for the strict post-1B run:
`recall_metrics/linearized_tables_20260612_strict/`.

## Four-Gene Aggregate

| Metric | Matched / Gold | Recall | Gap to 90% |
| --- | ---: | ---: | ---: |
| PMIDs | 1274 / 1502 | 84.8% | 78 |
| Variant rows | 5518 / 6833 | 80.8% | 632 |
| Unique variants | **2591 / 3010** | **86.1%** | **118** |
| Patients/carriers | 15896 / 18719 | 84.9% | 951 |
| Affected | 10435 / 12475 | 83.6% | 793 |
| Unaffected | 3441 / 3951 | 87.1% | 115 |

Rows-mode MAE:

| Count field | Sum abs error / N | MAE |
| --- | ---: | ---: |
| Carriers | 2287 / 3718 | **0.615** |
| Affected | 1535 / 3110 | **0.494** |
| Unaffected | 323 / 271 | **1.192** |

## Per-Gene Recall

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected | carriers MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 230/262 (87.8%) | 820/991 (82.7%) | 441/530 (83.2%) | 2256/2674 (84.4%) | 1404/1635 (85.9%) | 599/749 (80.0%) | 0.860 |
| KCNQ1 | 285/305 (93.4%) | 1499/1741 (86.1%) | 563/622 (**90.5%**) | 6995/7793 (89.8%) | 3909/4306 (90.8%) | 1319/1484 (88.9%) | 0.935 |
| SCN5A | 620/757 (81.9%) | 2433/3128 (77.8%) | 1022/1183 (86.4%) | 5020/6219 (80.7%) | 3836/4876 (78.7%) | 1184/1343 (88.2%) | 0.454 |
| RYR2 | 139/178 (78.1%) | 766/973 (78.7%) | 565/675 (83.7%) | 1625/2033 (79.9%) | 1286/1658 (77.6%) | 339/375 (90.4%) | 0.323 |

## Precision Snapshot

Headline precision is `precision_vs_counted_gold_pmids`, which restricts the
denominator to extra rows on gold PMIDs that carry at least one extracted count:
`5518 / (5518 + 1631) = 77.2%`.

The looser raw proxy remains useful only as a false-positive upper bound:
`5518 / (5518 + 13021) = 29.8%`.

Why the raw proxy is pessimistic:

- 11,390 / 13,021 current extra-on-gold-PMID rows have zero patient counts and
  are ClinVar/PubTator-style linkage attributions rather than count-bearing paper
  extractions.
- Only 1,631 extra rows carry any carrier/affected/unaffected count.
- About 97% are well-formed variants absent from the count-curated gold packet,
  not malformed output.
- The scorer now rejects 64 obvious figure/regex-table junk rows before scoring
  (gene-symbol-as-variant, <=2-character protein notation, residue prose). This
  removed 41 extra-on-gold-PMID rows, including 8 counted extras, with recall and
  MAE unchanged.
- 53 structural/CNV rows are real biology but currently unmatchable by the
  variant matcher.

Current per-layer precision proxy:

| Source layer | Matched DB rows | Extra rows | Counted extra rows | precision_vs_gold_pmids | precision_vs_counted_gold_pmids |
| --- | ---: | ---: | ---: | ---: | ---: |
| clinvar | 444 | 2484 | 71 | 15.2% | 86.2% |
| figure | 236 | 465 | 39 | 33.7% | 85.8% |
| llm_table | 892 | 477 | 275 | 65.2% | 76.4% |
| llm_text | 455 | 823 | 168 | 35.6% | 73.0% |
| mixed | 1949 | 392 | 176 | 83.3% | 91.7% |
| pubtator | 12 | 159 | 0 | 7.0% | 100.0% |
| regex_table | 1256 | 3864 | 929 | 24.5% | 57.5% |
| regex_text | 213 | 4971 | 2 | 4.1% | 99.1% |

Interpretation: recall gains are mostly adding real signal; the raw proxy
overstates true false positives by roughly 7x.

## Current Failure Split

Current failure-mode split from `paper_disagreement_report.csv`:

| Failure mode | Missing rows | What it means |
| --- | ---: | --- |
| source_missing_or_stub | 568 | paper/source never landed or only a stub landed |
| source_abstract_only | 250 | abstract was available, but mutation tables/body were missing |
| available_source_underextraction | 248 | usable source exists but extraction missed rows |
| source_missing_table_bodies | 184 | supplement/full text landed without the relevant tables |
| partial_underextraction | 82 | some rows extracted, table not exhausted |
| count_semantics | 36 | variant present but carrier/affected/unaffected semantics wrong |
| overinclusive_extraction | 8 | DB has many extra rows for the PMID; residual missing rows are not the main recall lever |

## Current Barriers

Use `TASKS.md` for the active checklist. The current blocker shape is:

1. **Source/table acquisition and binding.** Missing or incomplete source bodies,
   especially supplement mutation tables, remain the largest recall surface.
2. **Available-source underextraction.** Some usable sources are present but
   mutation-list tables are not exhausted.
3. **Count semantics and regex-table precision.** `regex_table` is the dominant
   count-bearing false-positive surface (`929` counted extras, 57.5% counted
   precision), so count-role attribution should be validated there first.
4. **Matcher notation and structural lanes.** Some already-extracted variants
   fail exact matching because indels, compound rows, splice/IVS rows, or
   structural/CNV forms have no matching lane.
5. **Access-gated publishers.** Wiley/Springer/Karger/Sage/Liebert failures are
   still relevant only where source audits show they block high-yield PMIDs.

## Scope Guard

This file should stay short: current metrics, current precision/failure split,
and current blockers only.

Historical recovery details, dated session logs, and superseded plans belong in
`docs/RECALL_HISTORY.md`. Operational procedures belong in the runbooks. Active
work items belong in `TASKS.md`.
