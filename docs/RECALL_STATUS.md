# Recall Status

Last updated: 2026-08-08.

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

## Current protocol note (2026-08-08)

The canonical DBs were rescored from scratch on 2026-08-08 with local-data
fallback disabled and source-heavy disagreement artifacts explicitly skipped.
All six recall numerators/denominators below reproduced exactly. The rescore did
correct slightly stale matched-count denominators and source-layer precision
totals in this file. It made no extraction or canonical-DB mutation. The current
end-to-end strategy and resolved workstation route are described in
`docs/ARCHITECTURE.md`; per-run `run_manifest.json` and `llm_traces/` remain the
authoritative execution record.

## Current Canonical Baseline

Fresh run of `scripts/run_recall_suite.py` against the four canonical DBs after
the 2026-07-12 four-gene supplement reconciliation, fold-gap closure, and the
strictly gated SCN5A supplement-source land:

- `validation_runs/canonical_baseline/KCNH2.db`
- `validation_runs/canonical_baseline/KCNQ1.db`
- `validation_runs/canonical_baseline/SCN5A.db`
- `validation_runs/canonical_baseline/RYR2.db`

Current local scored artifact: `recall_metrics/current/`.

## Targeted pfs12 enforcement spot check (not the canonical baseline)

> **Historical as of 2026-07-26.** The per-paper final check (Steps 3.8/3.9) is now **parked — default off** on cost/latency grounds; see `docs/EXTRACTION_CONTRACT.md` and `docs/PROTOCOL_CHANGELOG.md`. The numbers below record what the enforced `pfs12` protocol measured when it ran; they are not a description of the current default pipeline.

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
| Carriers | 2287 / 3729 | **0.613** |
| Affected | 1535 / 3121 | **0.492** |
| Unaffected | 323 / 271 | **1.192** |

End-to-end count MAE (all asserted gold rows, treating a missed extraction as
zero rather than dropping it from the denominator): carriers **1.744**
(10413/5971), affected **1.414** (7566/5349), unaffected **3.862** (3240/839).

## Per-Gene Recall

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected | carriers MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 230/262 (87.8%) | 820/991 (82.7%) | 441/530 (83.2%) | 2256/2674 (84.4%) | 1404/1635 (85.9%) | 599/749 (80.0%) | 0.860 |
| KCNQ1 | 285/305 (93.4%) | 1499/1741 (86.1%) | 563/622 (**90.5%**) | 6995/7793 (89.8%) | 3909/4306 (90.8%) | 1319/1484 (88.9%) | 0.935 |
| SCN5A | 622/757 (82.2%) | 2461/3128 (78.7%) | 1027/1183 (86.8%) | 5068/6219 (81.5%) | 3884/4876 (79.7%) | 1184/1343 (88.2%) | 0.450 |
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
`5546 / (5546 + 1637) = 77.2%`.

The looser raw proxy remains useful only as a false-positive upper bound:
`5546 / (5546 + 13036) = 29.8%`.

Why the raw proxy is pessimistic:

- 11,399 / 13,036 current extra-on-gold-PMID rows have zero patient counts and
  are ClinVar/PubTator-style linkage attributions rather than count-bearing paper
  extractions.
- Only 1,637 extra rows carry any carrier/affected/unaffected count.
- The scorer rejected 63 obvious figure/regex-table junk rows before this
  rescore (gene-symbol-as-variant, <=2-character protein notation, residue
  prose); recall and MAE are computed after that deterministic filter.
- 53 structural/CNV rows are real biology but currently unmatchable by the
  variant matcher.

Current per-layer precision proxy:

| Source layer | Matched DB rows | Extra rows | Counted extra rows | precision_vs_gold_pmids | precision_vs_counted_gold_pmids |
| --- | ---: | ---: | ---: | ---: | ---: |
| clinvar | 441 | 2416 | 79 | 15.4% | 84.8% |
| figure | 231 | 460 | 39 | 33.4% | 85.6% |
| llm_table | 883 | 477 | 274 | 64.9% | 76.3% |
| llm_text | 480 | 823 | 168 | 36.8% | 74.1% |
| mixed | 1899 | 393 | 176 | 82.9% | 91.5% |
| pubtator | 12 | 159 | 0 | 7.0% | 100.0% |
| regex_table | 1442 | 3734 | 899 | 27.9% | 61.6% |
| regex_text | 158 | 4574 | 2 | 3.3% | 98.8% |

Interpretation: the raw proxy counts roughly eight times as many extra rows as
the count-bearing proxy, so neither should be presented as clean precision.

## Last Source-Backed Failure Split (2026-07-20)

The 2026-08-08 run was intentionally metric-only. The last completed
source-backed split from `paper_disagreement_report.csv` is retained below and
should be regenerated before claiming that its categories moved:

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
   count-bearing false-positive surface (`899` counted extras, 61.6% counted
   precision), so count-role attribution should be validated there first.
4. **Matcher notation and structural lanes.** Some already-extracted variants
   fail exact matching because indels, compound rows, splice/IVS rows, or
   structural/CNV forms have no matching lane.
5. **Access-gated publishers.** Wiley/Springer/Karger/Sage/Liebert failures are
   still relevant only where source audits show they block high-yield PMIDs.

## Scope Guard

This file should stay short: current metrics/precision, the latest dated
source-backed failure split, and current blockers only.

Historical recovery details, dated session logs, and superseded plans belong in
`docs/RECALL_HISTORY.md`. Operational procedures belong in the runbooks. Active
work items belong in `TASKS.md`.
