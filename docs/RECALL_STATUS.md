# Recall Status

Last updated: 2026-08-21.

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
headline. Gate 2 (`gold_120`) is a fixed, extraction-blinded sample of 30
source-available, count-eligible papers per cardiac gene, now 118 attempts / 114
unique PMIDs after the non-genetics KCNH2 PMIDs 10086972 and 14642689 were
quarantined. The 2026-08-20 no-inference
revalidation **passed the precision gate**:
`precision_vs_counted_gold_pmids` is 97.51%, above the current 77.3% floor, and
variant recall is 86.57%. Carrier coverage is 30.49%, but conditional carrier
MAE regressed to 0.425 from 0.308 on the accepted 2026-08-13 run. Affected and
unaffected coverage fell to 7.58% and 5.06% because unsupported partitions are
now deliberately null. This is an identity recall/precision improvement, not a
new penetrance or affected/unaffected baseline. The registered Tier-3 queues
remain BMPR2 50, BRCA1 50, and BRCA2 45; a separate no-publish collaborator
candidate now completes the requested full process for exactly **50 BRCA1 + 50
BRCA2 + 50 BMPR2 papers**. For every gene, pinned manifest membership, staged
JSON membership, and final SQLite paper membership agree exactly. PMID 19944633
(the canine BRCA2 paper that exposed the prior harness defect) is absent from
all three. Mandatory VariantFeatures enrichment and high-confidence
wrong-gene/out-of-range quarantine leave 3,582 BRCA1, 722 BRCA2, and 260 BMPR2
live variants. The final DBs contain zero placeholder titles, wrong-gene live
rows, nameless identities, negative counts, impossible unquarantined
partitions, or duplicate penetrance strata. Typed family observations remain
as raw evidence and the trust gate masks them from carrier-facing totals.

The source-evidence replay gate is gold-free and independent of the ordinary
row-count gate. It requires a unique compatible identity, never aliases BIC
digits to cDNA, never matches by codon position, fills only null or identical
untyped values, and restores the backup on absent, ambiguous, or conflicting
untyped evidence. A shared ref+position+alt rule folds abbreviated
substitution-looking frameshift labels only when the explicit frameshift proves
the relationship. Protein-only evidence remains separate from cDNA-rich rows
because a later cDNA can make streaming uniqueness false. Grok 4.6 `xhigh`
adversarial review found the ambiguity waiver, insertion-order split, and a
live regex/source-notation gap; all are now pinned by regression tests. Full
evidence and residual conflicting-protein groups are recorded in
`docs/evidence/brca_bmpr2_full150_audit_20260821.md`. None of these three genes
has a new source-adjudicated precision estimate, so public publication remains
on hold.
Staging also exposed and closed a Variant Browser fail-open path: a reviewed
subject whose evidence vanished was detached but its gold revision could remain
current. The importer now disputes that revision and a regression test pins the
behavior. The BRCA2 replay has zero active canine evidence, zero current canine
gold, and zero current gold for all 111 detached/re-review subjects. Historical
revisions remain auditable. None of this is evidence that the new extraction is
publication-ready; public annotations remain unchanged.

The first 2026-08-21 current-code cardiac attempt is also not evidence. It was
terminated before lock when the operator observed that `gvf-run` had
auto-discovered the local KCNH2 gold CSV for intermediate recovery-layer
scoring. No `RUN_STATUS.json`, prediction lock, or score was accepted. The
replacement harness disabled gold discovery for the entire production run via
`--gold-free-run` and completed from a clean scaffold. The locked trusted result
is 489 TP / 438 FP / 144 FN (77.25% recall), counted-extra precision 97.80%, and
carrier MAE 0.247. A real SCN5A Q1077/Q1077del lookup mismatch over-held true
identities; the constrained post-lock diagnostic is 546 TP / 452 FP / 87 FN
(86.26% recall). Because that correction was informed by the locked residuals,
the diagnostic is not a new blind lock or public headline. The subsequent
candidate-local run at
`benchmarks/codex_paper_eval/runs/20260821_candidate_local_gold119/` is a fresh
gold-free production run locked before scoring. It records 546 TP / 290 FP /
87 FN (86.26% recall, 65.31% raw precision), 98.03% counted-extra precision,
94.44% count-bearing-only precision, and carrier MAE 0.255. All four production
statuses completed with zero stage failures and 1,245 trace records verified
against their write-time hashes. Its explicit empty for unrelated PMID
14642689 was correct; that record was quarantined from the live tier only after
the immutable run was scored.

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
`benchmarks/evaluation_tiers/`: 50 gold-scored attempts, a cardiac gold
expansion now containing 118 attempts (114 unique PMIDs; KCNH2 28 after
quarantining 10086972 and 14642689), and 545
full reviewer-backlog attempts (506 unique PMIDs). These tiers govern
evaluation/review scope, not the authoritative four-gene headline cohort below. The full tier includes BMPR2 and ranked
50-paper LMNA/TTN subsets; it keeps BMPR2 and BRCA1 at 50 papers and BRCA2 at 45
rather than expanding the experimental genes.

## Current-code blind candidate-local lock

The accepted clean blind lock is
`benchmarks/codex_paper_eval/runs/20260821_candidate_local_gold119/`. It used
the then-live 119-attempt cohort and remains immutable after the subsequent
14642689 cohort quarantine.

| Projection | TP | FP | FN | Variant recall | Variant precision | Counted-extra precision | Count-bearing-only precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Locked trusted production | 546 | 290 | 87 | 86.26% | 65.31% | 98.03% | 94.44% | 184/633 | 0.255 |

The candidate-local table fix was specified from a production-path failure and
reviewed before this run; no locked residual was fed back into its predictions.
The source binding, prediction digest, four production trace-manifest digests,
and post-lock score are recorded in the run directory.

## Prior current-code lock and post-lock identity diagnostic

The immutable current-code lock is
`benchmarks/codex_paper_eval/runs/20260821_current_changes_gold119/`; its
provenance erratum and the SCN5A analysis are linked from that run and
`docs/evidence/scn5a_q1077_identity_diagnostic_20260821.md`.

| Projection | TP | FP | FN | Variant recall | Variant precision | Counted-extra precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Locked trusted primary | 489 | 438 | 144 | 77.25% | 52.75% | 97.80% | 178/633 | 0.247 |
| Post-lock trusted diagnostic | 546 | 452 | 87 | 86.26% | 54.71% | 98.03% | 188/633 | 0.250 |
| Post-lock raw diagnostic | 551 | 578 | 82 | 87.05% | 48.80% | 98.04% | 202/633 | 0.535 |

The lock-to-diagnostic delta is mixed: +57 TP and +14 FP, so the precision of
net added identity matches is 80.3%, not the 98.03% counted-extra metric. The
lookup is SCN5A-only, proves the local Q1077del reference, applies only after
exact protein/cDNA matching fails, trusts simple missense/stop calls, and never
rewrites the literature identity. Structural and generic offsets remain held.
Grok 4.6 xhigh gave a GO for merge/trusted staging and a HOLD for public
publication. A disjoint pre-registered validation is required.

## Prior Gate 2 no-inference gold-119 result (not current-code acceptance)

The locked run is
`benchmarks/codex_paper_eval/runs/20260820_gold119_noinference/`. Exact
production trace manifests for all four genes are bound into both the
predictions and the pre-gold lock. The extraction code and prompt hashes are
recorded in the run manifests; the worktree was intentionally dirty with the
reviewed no-inference patch.

| Projection | Variant precision | Variant recall | Precision vs counted extras | Count-bearing-only precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026-08-20 raw locked primary | 548/1216 (45.07%) | 548/633 (86.57%) | 548/(548+14) (**97.51%**) | 195/(195+14) (93.30%) | 193/633 (30.49%) | 0.425 |
| 2026-08-13 accepted raw comparison | 534/1438 (37.13%) | 534/635 (84.09%) | 534/(534+24) (**95.70%**) | 148/(148+24) (86.05%) | 146/635 (22.99%) | 0.308 |

The acceptance metric counts every matched gold row as signal and restricts
only the extra-row denominator to predictions carrying a patient count. The
93.30% diagnostic instead restricts the numerator to count-bearing matches and
is intentionally stricter, but it cannot be compared with the 77.3% floor. The
fresh run confirms all 42 PMID 26746457 classification-table identities survive
with null carrier counts. It also records KCNH2 PMID 14642689 and SCN5A PMID
22685113 as explicit empty outcomes after circuit-breaker/source failure rather
than silently omitting them. Affected coverage is 48/633 (7.58%, MAE 0.688) and
unaffected coverage is 32/632 (5.06%, MAE 1.000); these are abstention-driven
coverage changes and must accompany any comparison. This sample does not
replace the canonical all-paper four-gene headline below. It predates the
current active-DB, trace-presence, and trusted-identity hardening. The
2026-08-21 lock above supersedes that outstanding rerun task but does not
replace this historical comparison.

## Active 50-paper collaborator-grounded count-semantics cohort

As of 2026-08-11, the active scored count evaluation and new strategy comparisons use the fixed
cardiac 48 plus only the two BRCA2 papers with lead-approved Variant Browser
adjudications by Nate (PMIDs 26833046 and 26848529). The six internally derived
BRCA2 papers have been removed from active membership without changing the
historical run.

The active BRCA2 reviewer queue is now 45 papers: four of the six provenance
exclusions were removed from the prior 50-paper snapshot, two were already
absent, and the explicitly canine BRCA2 ortholog paper PMID 19944633 was removed
on 2026-08-20. Nate's two papers and their 87 current gold records remain intact.

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
