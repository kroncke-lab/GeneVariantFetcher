# Recall Status

Last updated: 2026-07-20.

This file is the current measured recall snapshot. It intentionally does not
carry the active work plan or dated session log.

- Active forward plan and task checklist: `TASKS.md`.
- Benchmark/change trajectory: `docs/RECALL_HISTORY.md`.
- Re-run procedure: `docs/RECALL_REFRESH_RUNBOOK.md`.

No other doc should restate live recall tables. If a metric conflicts with this
file, this file is authoritative.

Live view: the published status dashboard renders these numbers at
<https://kroncke-lab.github.io/GeneVariantFetcher/dashboard/> (built by
`scripts/build_status_dashboard.py` into `docs/dashboard/`).

## Metrics scope: cardiac four only

**Recall, precision, and MAE are computed only against the four cardiac genes —
KCNH2, KCNQ1, SCN5A, RYR2 — because only those have a fully human-curated,
manually derived gold standard** (`gene_variant_fetcher_gold_standard/`). The
non-cardiac genes (APOE, BRCA1, BRCA2, MYBPC3) have only curator/LLM-derived
`gold_overrides/` answer keys, which are useful for review but are **not** a
manual gold standard, so they are **excluded from every headline metric here and
on the dashboard**. Score them for spot checks if you like, but never fold them
into the reported recall/precision/MAE. To reproduce the headline numbers,
restrict scoring to the four cardiac genes, e.g.
`run_benchmark.py --genes KCNH2,KCNQ1,SCN5A,RYR2`.

## Current Canonical Baseline

Fresh run of `scripts/run_recall_suite.py` against the four canonical DBs after
the 2026-07-12 four-gene supplement reconciliation, fold-gap closure, and the
strictly gated SCN5A supplement-source land:

- `results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db`
- `validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db`

Scored artifact:
`recall_metrics/fulltext_supplements_20260712/`.

## Targeted pfs12 enforcement spot check (not the canonical baseline)

The enforced paper-final-check protocol was replayed on fresh single/few-paper
DBs for KCNH2 (PMIDs 15840476 and 33013630), KCNQ1 (30758498), RYR2
(28404607), and SCN5A (29325976). Metrics below are restricted to those PMIDs;
they validate trust projection behavior and must not be compared with the
four-gene aggregate as if they were a fleet refresh.

| Gene | Variant recall | Affected recall | Matched carrier MAE / RMSE | End-to-end carrier MAE / RMSE | Grounded missing groups | Applied facts |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 85/86 (98.8%) | 107/108 (99.1%) | 0.035 / 0.243 | 0.047 / 0.264 | 1 | 1 |
| KCNQ1 | 14/93 (15.1%) | 68/228 (29.8%) | 5.833 / 10.206 | 3.064 / 5.207 | 91 | 0 |
| RYR2 | 235/244 (96.3%) | 242/276 (87.7%) | 0.004 / 0.065 | 0.143 / 1.676 | 1 | 1 |
| SCN5A | 20/91 (22.0%) | 37/118 (31.4%) | n/a | 1.297 / 1.460 | 18 | 0 |

Trusted and raw count metrics were identical for KCNH2, KCNQ1, and SCN5A. For
RYR2, a source-quoted phenotype contradiction quarantined one affected field;
matched affected MAE changed from 0.00426 raw to 0.00427 trusted. The source
explicitly says that carrier had no VT/CPVT, while the current gold row expects
affected status, so this is a gold-adjudication discrepancy rather than evidence
that the source-grounded gate selected the wrong field. The grounded missing
groups make the KCNH2, KCNQ1, RYR2, and SCN5A runs fail acceptance and route to
re-extraction instead of allowing their high spot-check recall gaps to pass
silently. A full canonical four-gene rerun remains required before changing the
headline baseline.

## Four-Gene Aggregate

| Metric | Matched / Gold | Recall | Gap to 90% |
| --- | ---: | ---: | ---: |
| PMIDs | 1276 / 1502 | 85.0% | 76 |
| Variant rows | 5546 / 6833 | 81.2% | 604 |
| Unique variants | **2596 / 3010** | **86.2%** | **113** |
| Patients/carriers | 15944 / 18719 | 85.2% | 904 |
| Affected | 10483 / 12475 | 84.0% | 745 |
| Unaffected | 3441 / 3951 | 87.1% | 115 |

Rows-mode MAE:

| Count field | Sum abs error / N | MAE |
| --- | ---: | ---: |
| Carriers | 2287 / 3724 | **0.614** |
| Affected | 1535 / 3116 | **0.493** |
| Unaffected | 323 / 271 | **1.192** |

## Per-Gene Recall

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected | carriers MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 230/262 (87.8%) | 820/991 (82.7%) | 441/530 (83.2%) | 2256/2674 (84.4%) | 1404/1635 (85.9%) | 599/749 (80.0%) | 0.860 |
| KCNQ1 | 285/305 (93.4%) | 1499/1741 (86.1%) | 563/622 (**90.5%**) | 6995/7793 (89.8%) | 3909/4306 (90.8%) | 1319/1484 (88.9%) | 0.935 |
| SCN5A | 622/757 (82.2%) | 2461/3128 (78.7%) | 1027/1183 (86.8%) | 5068/6219 (81.5%) | 3884/4876 (79.7%) | 1184/1343 (88.2%) | 0.452 |
| RYR2 | 139/178 (78.1%) | 766/973 (78.7%) | 565/675 (83.7%) | 1625/2033 (79.9%) | 1286/1658 (77.6%) | 339/375 (90.4%) | 0.323 |

## Full-Text and Supplement Coverage

The consolidated corpus currently contains 1,312 KCNH2, 2,396 KCNQ1, 656
RYR2, and 1,590 SCN5A paper contexts. Among papers with an on-disk convertible
supplement, the fold audit is **1,577/1,577 folded**: KCNH2 355/355, KCNQ1
604/604, RYR2 233/233, and SCN5A 385/385. The pre-refresh gap was 289 papers.

The Elsevier per-file reconciliation recovered 64 missing `mmc` files across 49
papers without re-downloading article bodies, then updated 427 folded contexts.
Remaining source gaps are papers or referenced supplements that are not locally
available through current publisher access; they are not an on-disk fold gap.

## Precision Snapshot

Headline precision is `precision_vs_counted_gold_pmids`, which restricts the
denominator to extra rows on gold PMIDs that carry at least one extracted count:
`5546 / (5546 + 1629) = 77.3%`.

The looser raw proxy remains useful only as a false-positive upper bound:
`5546 / (5546 + 13036) = 29.8%`.

Why the raw proxy is pessimistic:

- 11,407 / 13,036 current extra-on-gold-PMID rows have zero patient counts and
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

1. **Source/table acquisition and binding.** The local fold backlog is closed;
   missing or access-gated source bodies and referenced supplements remain the
   largest recall surface.
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
