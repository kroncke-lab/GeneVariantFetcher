# Recall Status

Last updated: 2026-08-15.

This file is the current measured recall snapshot. It intentionally does not
carry the active work plan or dated session log.

- Active forward plan and task checklist: `TASKS.md`.
- Benchmark/change trajectory: `docs/RECALL_HISTORY.md`.
- Re-run procedure: `docs/RECALL_REFRESH_RUNBOOK.md`.

No other doc should restate live recall tables. If a metric conflicts with this
file, this file is authoritative.

Audit note: the scorer now preserves explicit adjudicated zero counts instead
of conflating them with null. The tables below remain the last published
pre-correction baseline and must be re-scored before any changed headline is
claimed. Source-reachable strata are secondary reader diagnostics; ALL GOLD
remains the primary turnkey acceptance denominator.

The deterministic parsers also now leave affected/unaffected null when the
source supplied only a carrier total; explicit phenotype/control cells still
preserve zero. Mutating zeros in an already-extracted DB does not reproduce that
policy and changes the matched-count denominator, so promotion requires a paired
live cardiac re-extraction and rescore with assertion coverage plus conditional
and end-to-end error reported together. That gate is tracked in `TASKS.md`.

A first paired rescore ran on 2026-08-13 against the 73-paper curated cardiac
fixture and is recorded in `docs/RECALL_HISTORY.md`. It does **not** change the
headline below: its KCNQ1 arm was degraded by a source-recovery failure, and its
"before" side is an older DB rather than a re-extraction with the old code. Its
one load-bearing result is that matched-row MAE and end-to-end count error move
in **opposite** directions — conditional accuracy improves 16–57% while
end-to-end carrier/affected error rises, because the pipeline now declines to
answer more often. Any future report of one without the other is misleading.
The lead approved advancement from Gate 1 on 2026-08-13 without changing the
headline. Gate 2 (`gold_120`) then ran as a fixed, gold-value-blinded sample of
30 source-available, count-eligible papers per cardiac gene (120 attempts / 116
unique PMIDs; seed 2026081301). The patched-system revalidation **passed the
precision gate**: `precision_vs_counted_gold_pmids` is 95.70% raw and 95.87%
trusted, above the current 77.3% floor. Variant recall is 84.09%; carrier MAE is
0.308 raw / 0.299 trusted, versus the 0.614 canonical all-paper baseline. The
immediately preceding stochastic sample was lower at 0.266 / 0.243, so the
accepted revalidation improves precision while showing a small conditional-MAE
regression relative to that one run. The experimental BMPR2, BRCA1, and BRCA2
queues remain 50, 50, and 46 papers; their no-publish extraction/QC run is
complete. It used 972 calls, 4.261M tokens, 3.712 summed provider-hours, and a
$23.664 public-price proxy. The two lead-approved BRCA2 gold papers recovered
all 7/7 curated variant identities; carrier-count coverage is 3/7 with MAE
1.333. That sparse seven-record stratum is not an exhaustive paper-level
precision denominator. After this gate passed, the fixed BMPR2 50 / BRCA1 50 /
BRCA2 46 queues were refreshed in collaborator-facing Variant Browser staging
from the current-system run. Live trusted evidence counts are 482, 7,260, and
2,346, with exact manifest membership and order. All 111 historical BRCA2
adjudications remain auditable and now require re-review; none is eligible for
the default adjudication or gold export. A bounded two-paper BRCA2 count-
recovery dry-run added no inferred counts (0/26 gaps grounded). Public
annotations remain unchanged.

Historical view: the published dashboard at
<https://kroncke-lab.github.io/GeneVariantFetcher/dashboard/> is the archived
2026-07-08 pre-correction snapshot. It does **not** render the current table
below. Regenerate it with `scripts/build_status_dashboard.py` only after the
scorer rescore and parser re-extraction gate are accepted; until then this file
remains authoritative.

## Metrics scope: cardiac four only

**Recall, precision, and MAE are computed only against the four cardiac genes —
KCNH2, KCNQ1, SCN5A, RYR2 — because only those have a fully human-curated,
manually derived gold standard** (`gene_variant_fetcher_gold_standard/`). The
non-cardiac genes do not have a complete manually derived gold standard. BRCA2
now has two papers with collaborator review and lead approval in Variant
Browser, while its remaining override rows and the APOE/BRCA1/MYBPC3 keys are
internally derived. All remain **excluded from every headline metric here and
on the dashboard**. Score them as explicitly labeled strata, but never fold them
into the reported cardiac recall/precision/MAE. To reproduce the headline numbers,
restrict scoring to the four cardiac genes, e.g.
`run_benchmark.py --genes KCNH2,KCNQ1,SCN5A,RYR2`.

## Protocol context (landed 2026-07-20)

The trust/provenance/gold-integrity arc **#161–#165** has landed (see
`docs/RECALL_HISTORY.md`). These are gold-free / additive / scorer-invariant
changes — they harden trustworthiness and BRCA-readiness and **do not change the
four-gene headline** below (`#165` was verified to add 0 new cardiac quarantine
on KCNH2). The current end-to-end protocol is described in
`docs/ARCHITECTURE.md`. Dated cost-and-quality measurements live in
`docs/PROTOCOL_COST_EVAL.md`; they describe the routes used for those samples,
not current defaults. The next acceptance sequence lives only in `TASKS.md`.

## Canonical rollout tiers

The active rollout population is governed by exactly three manifests in
`benchmarks/evaluation_tiers/`: 50 gold-scored attempts, a 119-attempt cardiac
gold expansion (115 unique PMIDs; KCNH2 29 after removing 10086972), and 546
full reviewer-backlog attempts (507 unique PMIDs). These tiers govern
evaluation/review scope, not the authoritative four-gene headline cohort below. The full tier includes BMPR2 and ranked
50-paper LMNA/TTN subsets; it keeps BMPR2 and BRCA1 at 50 papers and BRCA2 at 46
rather than expanding the experimental genes.

## Gate 2 patched gold-120 result (not a headline; passed rollout gate)

The locked run is
`benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/`. Exact production
trace manifests for all four genes are bound into both the predictions and the
pre-gold lock. The primary score preserves all raw count observations; the
read-only trusted diagnostic masks only persisted quarantined fields.

| Projection | Variant precision | Variant recall | Precision vs counted extras | Count-bearing-only precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Raw locked primary | 534/1438 (37.13%) | 534/635 (84.09%) | 534/(534+24) (**95.70%**) | 148/(148+24) (86.05%) | 146/635 (22.99%) | 0.308 |
| Trusted diagnostic | unchanged | unchanged | 534/(534+23) (**95.87%**) | 139/(139+23) (85.80%) | 137/635 (21.57%) | 0.299 |

The acceptance metric counts every matched gold row as signal and restricts
only the extra-row denominator to predictions carrying a patient count. The
86.05% diagnostic instead restricts the numerator to count-bearing matches and
is intentionally stricter, but it cannot be compared with the 77.3% floor. The
trust projection removes one counted extra and masks nine matched carrier
assertions. The generic parser fix retains vertical-table variant identity but
declines to invent one carrier per row unless the table proves patient/subject
row semantics; the fresh run confirms all 42 PMID 26746457 classification-table
identities survive with null carrier counts. Relative to the prior locked run,
counted extras fell 58→24 while carrier absolute error moved 41/154→45/146 raw
(41/137 trusted). This sample does not replace the canonical all-paper
four-gene headline below.

A 2026-08-15 **diagnostic** rescore of the same locked predictions against
the live gold snapshot and current matcher lives in
`benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/`
(545/633 = 86.10% recall). It is not a new lock and does not replace the
table above.

## Active 50-paper collaborator-grounded count-semantics cohort

As of 2026-08-11, the active scored count evaluation and new strategy comparisons use the fixed
cardiac 48 plus only the two BRCA2 papers with lead-approved Variant Browser
adjudications by Nate (PMIDs 26833046 and 26848529). The six internally derived
BRCA2 papers have been removed from active membership without changing the
historical run.

The live BRCA2 reviewer queue is now 46 papers: four of the six were removed
from the prior 50-paper snapshot and two were already absent. Nate's two papers
and their 87 current gold records remain intact.

Projecting the already locked predictions onto this 50-paper cohort required no
new model calls. Carrier MAE is **0.0608 (20/329)** after count-scope repair
versus 0.9058 (298/329) under the legacy answer key; count recall is 32.67%
(329/1007). Report the strata alongside the combined result: cardiac 48 is
0.0491 (16/326), while the two-paper BRCA2 stratum is 1.3333 (4/3) and retains
known scope limitations. Details:
`benchmarks/count_semantics_eval/runs/20260811_collaborator_gold_50/`.

## Historical locked 56-paper count-semantics audit (not the active cohort)

The preserved 2026-08-10 count-semantics study used a deliberately difficult, fixed set
of 48 cardiac papers plus 8 BRCA2 papers. It is an error-analysis benchmark,
not a replacement for the four-gene canonical baseline below. The cardiac
papers use the manual cardiac gold standard; the BRCA2 arm used a
mixed-provenance internal answer key and remains a historical diagnostic only.

Most importantly, the **predictions were locked and did not change**. Compact
source cards and a blind independent audit showed that the largest apparent
carrier errors were stale or inconsistent answer-key definitions of the
current study cohort, plus one scorer that ignored an existing `gold_v2`
adjudication. After correcting those reference/scoring defects, carrier MAE
among the same 378 supplied predictions changed from **0.8148 (308/378)** to
**0.0794 (30/378)**. This is a 90.3% reduction in measured absolute error, but
it is **not** a 90.3% extraction improvement.

| Fixed 56-paper metric | Before | After |
| --- | ---: | ---: |
| Carrier MAE | 0.8148 | **0.0794** |
| Carrier count recall | 34.02% | **34.05%** |
| Affected MAE | **0.7869** | 0.7902 |
| Unaffected MAE | 0.1802 | **0.1307** |

The slight affected-MAE regression is source-supported and is retained as a
negative control against adjudicating toward the predictions. Cardiac-only
carrier MAE on this selected set is 0.0491 (16/326); BRCA2 is 0.2692 (14/52).
Neither value is directly comparable to the canonical all-paper cardiac MAE of
0.614 because the cohorts and answer-key maturity differ.

The risk-ranked claim-verification implementation is merged and enabled inside
Tier-3 extraction, but the seven-call **Luna** experiment was a shadow analysis.
The measured A1 production trace used GPT-5.6 Sol for claim verification. The
planned post-recovery failure router has not yet been implemented or measured.
Reproducible metrics and locked digests live in
`benchmarks/count_semantics_eval/runs/20260810_luna_xhigh_56/`; the publication-
oriented design record is `benchmarks/count_semantics_eval/METHODS_20260810.md`.

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
