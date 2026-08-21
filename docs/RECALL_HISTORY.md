# GVF Recall History

**Append-only.** This file is the durable, chronological record of recall
benchmarks and the changes that moved them. Add new dated entries at the **top of
the Timeline**; never delete or rewrite past entries (correct a past number only
by adding a new entry that supersedes it, with a note). Live/working status lives
in [`RECALL_STATUS.md`](RECALL_STATUS.md) (the single source of truth for the
*current* numbers) and the re-run procedure in
[`RECALL_REFRESH_RUNBOOK.md`](RECALL_REFRESH_RUNBOOK.md); this file is the *history*.

Grant target: **90% unique-variant recall**, submission **October 2026**.
Validation surface: four cardiac genes **KCNH2 / KCNQ1 / SCN5A / RYR2**
(gold standards in `gene_variant_fetcher_gold_standard/normalized/`). All
recall numbers below are **figures-skipped, DB-observed** scoring via
`scripts/run_recall_suite.py` unless noted, and depend on the Vanderbilt
Elsevier insttoken.

## Headline benchmark trajectory (four-gene aggregate)

| Date | Milestone | Unique-variant recall | Variant-row recall | carriers MAE |
|---|---|---:|---:|---:|
| 2025-11-18 | Project start (initial commit) | — | — | — |
| 2026-05-25 | First multi-gene scored artifact | ~1859/3013 (61.7%) | — | — |
| 2026-05-26 | Pre-SUA-sweep active DBs | 1884/3010 (62.6%) | 3714/6833 (54.4%) | — |
| 2026-05-26 | Post-SUA source-replay sweep | 2206/3010 (73.3%) | 4173/6833 (61.1%) | — |
| 2026-05-26 | Post-rollback + ClinVar/PubTator recovery | 2450/3010 (81.4%) | 5065/6833 (74.1%) | 1.055 |
| 2026-05-29 | cDNA↔protein matcher bridge + count-outlier guard | 2473/3010 (82.2%) | 5160/6833 (75.4%) | ~0.90 |
| **2026-06-05** | **Elsevier mmc supplement acquisition (landed)** | **2523/3010 (83.8%)** | **5350/6833 (78.3%)** | **0.882** |
| 2026-06-06 | refresh_recall re-extract (recall up, MAE regressed) | 2572/3010 (85.4%) | 5423/6833 (79.4%) | 1.274 |
| **2026-06-06** | **Duplicate-penetrance idempotency fix** | **2572/3010 (85.4%)** | **5423/6833 (79.4%)** | **0.614** |
| **2026-06-12** | **PDF-linearized table reconstruction + iter-2 quality gate + targeted KCNQ1 land** | **2590/3010 (86.0%)** | **5514/6833 (80.7%)** | **0.634** |
| **2026-06-12** | **+ targeted lands KCNH2/SCN5A/RYR2 (all four genes)** | **2591/3010 (86.1%)** | **5518/6833 (80.8%)** | **0.615** |
| **2026-07-12** | **Four-gene supplement reconciliation + gated SCN5A land** | **2596/3010 (86.2%)** | **5546/6833 (81.2%)** | **0.614** |
| 2026-07-20 | Trust/provenance/gold-integrity arc (#161–#165) — *no recall change by design* | 2596/3010 (86.2%) | 5546/6833 (81.2%) | 0.614 |
| — | **Target** | **2709/3010 (90.0%)** | — | → 0 |

Gap to the 90% unique-variant target: **113** variants (was 1126 at the 62.6%
starting point; ~89% of the original gap has been closed).

Per-gene unique-variant recall, latest canonical (2026-07-12):
KCNH2 **83.2%**, KCNQ1 **90.5%** (crossed the 90% target), SCN5A **86.8%**,
RYR2 **83.7%**.

---

## Timeline (newest first)

### 2026-08-21 — Current-code gold-119 lock exposed and constrained an SCN5A isoform mismatch

The clean `--gold-free-run` replacement completed all 119 attempts and locked
the trusted count/identity projection at 489 TP / 438 FP / 144 FN (77.25%
recall), 97.80% counted-extra precision, and 0.247 carrier MAE. The unexpectedly
large trusted-only recall loss localized to VariantFeatures identity projection:
the local SCN5A reference is the 2,015-aa Q1077del isoform while clinical papers
commonly number the 2,016-aa Q1077 isoform.

The post-lock correction preserves paper identity and uses a rigorously proven
SCN5A `N-1` coordinate only as a VariantFeatures lookup key for simple
missense/stop calls. It also fixes bare uppercase proline parsing. On the same
locked predictions the trusted diagnostic is 546 TP / 452 FP / 87 FN (86.26%
recall); raw is 551 / 578 / 82 (87.05%). This is in-sample calibration, not a
new lock: the mixed delta is +57 TP / +14 FP (80.3% extra-match precision), and
the original lock remains unchanged. Grok 4.6 xhigh returned GO for trusted
staging and HOLD for public publication. See
`docs/evidence/scn5a_q1077_identity_diagnostic_20260821.md`.

The BRCA1 50 / BRCA2 45 / BMPR2 50 exact manifests were imported privately via
the trusted projection. A Variant Browser fail-open gold bug found during the
BRCA2 replay was fixed: vanished reviewed subjects now make their current gold
revision disputed/non-current. Live BRCA2 has zero canine evidence, zero current
canine gold, and zero current gold for all 111 detached/re-review subjects.
Public annotations were not changed.

### 2026-08-20 — No-publish BRCA1/BRCA2/BMPR2 refresh completed; BRCA2 failed structural scope

After the gold-119 acceptance lock, the fixed collaborator manifests ran with
publication disabled and mandatory VariantFeatures enrichment plus quarantine.
BRCA1 completed 50/50 papers with 50 full texts, BRCA2 46/46 with 45 full texts,
and BMPR2 50/50 with 45 full texts. All three exited 0 with no stage failures or
warnings and write-time-verified traces (582, 374, and 192 LLM calls).

| Gene | Final variants | Nameless | Penetrance percentages | VF quarantined | Structural gate |
| --- | ---: | ---: | ---: | ---: | --- |
| BRCA1 | 5,348 | 0 | 0 | 324 | PASS |
| BRCA2 | 1,947 | 0 | 0 | 156 | **FAIL** |
| BMPR2 | 554 | 0 | 0 | 32 | PASS |

BRCA2 failed because the primary extractor correctly abstained on an explicitly
canine BRCA2 paper (PMID 19944633), but downstream PubTator recovery attached
four cDNA links. VariantFeatures did not quarantine them because in-range canine
coordinates can look valid against the human gene. The reusable readiness audit
now detects the strong species-scoped title pattern and fails the run, with a
tested fallback to staged abstract metadata when the database title is blank.
This establishes the need for paper/link-level relevance quarantine rather than
variant-global residue quarantine alone.

The two lead-approved BRCA2 papers retained 7/7 approved identities; carrier
coverage remains 3/7 with MAE 1.333. That seven-record positive subset is not
exhaustive paper-level gold, so unmatched predictions cannot be labeled false
positives. No new precision claim is made for any of the three genes. Evidence:
`docs/evidence/collaborator_readiness_20260820.md` and
`benchmarks/codex_paper_eval/runs/20260820_brca2_gold2_noinference/`.

### 2026-08-20 — No-inference gold-119 acceptance lock passed

The live 119-attempt / 115-unique cardiac manifest was re-extracted after the
penetrance and phenotype-partition no-inference changes. Production sources and
all four write-time-verified trace manifests were rebound and SHA-256 locked
before the scorer read gold. The primary projection contains 1,216 predicted
variants over 119 papers; four papers are explicit empty predictions. KCNH2
10086972 remains excluded from the live manifest.

| Metric | 2026-08-13 accepted lock | 2026-08-20 no-inference lock |
| --- | ---: | ---: |
| Variant recall | 534/635 (84.09%) | **548/633 (86.57%)** |
| Raw variant precision | 534/1438 (37.13%) | **548/1216 (45.07%)** |
| Counted-extra precision | 534/(534+24) (95.70%) | **548/(548+14) (97.51%)** |
| Count-bearing-only precision | 148/(148+24) (86.05%) | **195/(195+14) (93.30%)** |
| Carrier supplied / MAE | 146/635 / **0.308** | **193/633** / 0.425 |
| Affected supplied | 70/635 (11.02%) | 48/633 (7.58%) |
| Unaffected supplied | 53/634 (8.36%) | 32/632 (5.06%) |

The 97.51% counted-extra result passes the 77.3% rollout floor and variant
recall improves. Carrier conditional MAE regresses, so this is not a universal
metric win. Affected/unaffected supply falls because diagnosis, enrollment, and
arithmetic completion no longer manufacture phenotype partitions; that is a
deliberate abstention policy and must not be reported as a new penetrance
baseline. Exact artifacts:
`benchmarks/codex_paper_eval/runs/20260820_gold119_noinference/`.

### 2026-08-13 — Fixed experimental queues refreshed in collaborator staging

After the accepted gold-120 gate and paper/source QC, the current-system
BMPR2 50, BRCA1 50, and BRCA2 46 outputs were imported into the existing
Variant Browser review snapshots under dataset label
`collaborator_reextract_current_system_20260813`. Every live PMID and review
position exactly matches its pinned manifest; no paper was added or removed.
Trusted live evidence changed from 528 to 482 for BMPR2, 6,299 to 7,260 for
BRCA1, and 2,735 to 2,346 for BRCA2. The refresh exposes 470, 3,663, and 591
individual records plus 3,871, 27,172, and 4,920 exact provenance facts.

BRCA2's 111 prior adjudications and 111 gold-record keys were preserved exactly
against before/after audit exports. Every adjudication is detached from the new
evidence and marked `needs_re_review`: 98 subjects have a changed evidence
fingerprint and 13 no longer have a current evidence subject. The 98 matching
gold records are disputed/non-current; the 13 vanished-subject revisions remain
historical but detached. Default export returns zero BRCA2 adjudications and
zero gold records until collaborators review the new extraction. One reviewed
frameshift identity was retained while five empty legacy aliases were moved to
a reversible archive namespace before import. A source-grounded recovery dry-
run on PMIDs 12942367 and 22382806 used two GPT-5.6 Sol `high` calls, found 26
count gaps, grounded none, and wrote nothing. No public annotations were
published. Full hashes and reconciliation details are in
`runs/20260813_experimental_146/VARIANT_BROWSER_PUBLICATION.md`.

### 2026-08-13 — Fixed BMPR2/BRCA1/BRCA2 experimental run completed

The approved no-publish experiment stayed fixed at BMPR2 50, BRCA1 50, and
BRCA2 46. All 146 attempts completed with write-time-verified trace manifests:
972 calls (970 successful), 4,261,341 tokens, 3.712 summed provider-hours, and a
$23.664 public-list-price proxy. The three jobs ran concurrently and finished
in 70.7--97.8 minutes. Source integrity passed at 45/50, 50/50, and 45/46 full
text. Trust gate tg5 retained 231/252, 2062/2527, and 416/479 count rows;
family-count carrier fields remain raw-but-masked instead of being imported as
patient counts.

Only the two BRCA2 papers with collaborator-reviewed, lead-approved gold were
scored. All 7/7 curated variant identities were recovered; carrier counts were
supplied for 3/7 assertions with MAE 1.333. The seven selected records are not
an exhaustive inventory of the papers, so the raw extra-variant denominator is
not presented as true paper-level precision. Report-only somatic/germline QC
found 12.94% somatic+ambiguous BRCA1 records and 18.31% ambiguous BRCA2 records;
the initial run completed with publication disabled pending paper-level review.
The subsequent collaborator-staging refresh is recorded above; public
annotations remain unchanged. Detailed extraction evidence is in
`runs/20260813_experimental_146/EXPERIMENTAL_COST_AND_QC.md`.

### 2026-08-15 — Diagnostic gold-120 rescore (locked predictions, live gold)

Free re-score of `runs/20260813_gold120_verticalfix` predictions against the
live gold snapshot and current matcher. No new extraction. Locked parent
reports were not rewritten; artifacts are in
`diagnostics/current_gold_matcher_20260815/`.

| Metric | Locked 2026-08-13 | Diagnostic 2026-08-15 |
| --- | ---: | ---: |
| Variant recall | 534/635 (84.09%) | **545/633 (86.10%)** |
| Raw precision | 534/1438 (37.13%) | 545/1335 (40.82%) |
| Counted-extra precision | 534/(534+24) (95.70%) | 545/(545+10) (98.20%) |
| Carrier supplied | 146/635 (22.99%) | 160/633 (25.28%) |
| Carrier MAE | 0.308 | 0.2875 |

KCNH2 10086972 is empty on the live key, so the 120-paper and live 119-paper
totals match. This is not a new Gate 2 lock and does not change the four-gene
headline. Remaining 88 FNs are classified in that diagnostic `NOTES.md`.

### 2026-08-13 — Patched gold-120 revalidation accepted; experimental launch opened

The fresh current-system run in `runs/20260813_gold120_verticalfix` locked 120
attempts, exact run-local source digests, and four write-time-verified production
trace manifests before scoring. It reproduced 534/635 variant recall (84.09%)
with 1,438 predicted rows. The vertical-table fix reduced count-bearing extras
to 24 raw / 23 trusted, yielding **95.70% / 95.87%** counted-extra precision;
the separate count-bearing-only diagnostic is **86.05% / 85.80%**.

Carrier MAE is 0.308 raw / 0.299 trusted. That remains roughly half the 0.614
canonical all-paper baseline but is worse than the immediately preceding
stochastic sample's 0.266 / 0.243, so this is recorded as a precision gain with
a small conditional-MAE regression, not a universal quality improvement. Raw
carrier absolute error is 45 across 146 supplied matched counts (prior: 41/154).
The accepted lock opened the fixed BMPR2 50 / BRCA1 50 / BRCA2 46 experimental
launch with publication disabled pending QC.

### 2026-08-13 — Gate 2 precision-definition correction: passed

The initial Gate 2 decision below compared 158 count-bearing matches / 216
count-bearing predictions (73.15%) against the canonical 77.3%
`precision_vs_counted_gold_pmids` floor. That comparison was invalid: the floor's
numerator is **all matched gold rows**, while only its extra-row denominator is
restricted to rows carrying a patient count. On the correct like-for-like
definition, gold-120 is **534 / (534 + 58) = 90.20% raw** and
**534 / (534 + 54) = 90.82% trusted**, so Gate 2 passed. The scorer now emits
both metrics with explicit numerators to prevent recurrence.

The audit also found a gold-independent deterministic defect: a vertical-table
parser assigned one carrier to each row of a laboratory classification table.
The patched parser retains all 42 variant identities but emits no patient count
unless the table proves patient/subject row semantics. On the locked predictions
the exact counterfactual reaches **95.36% raw / 96.04% trusted** counted-extra
precision; raw/trusted carrier MAE are 0.278/0.254. A fresh patched-system
gold-120 run is the final revalidation before the approved BMPR2 50, BRCA1 50,
and BRCA2 46 launch.

### 2026-08-13 — Gate 2 gold-120 initial decision (superseded above)

The current full `gvf-run` route processed the fixed, gold-value-blinded
`tier2_gold_120.tsv`: 30 count-eligible manual-gold papers per cardiac gene,
120 gene-paper attempts / 116 unique PMIDs. Predictions plus the four exact
production LLM trace manifests were locked before the scorer read gold. The
run-local source snapshot was bound before scoring without changing cohort
membership or reading gold values.

| Metric | Raw primary | Trusted count projection |
| --- | ---: | ---: |
| Variant precision | 534/1440 (37.08%) | unchanged |
| Variant recall | 534/635 (84.09%) | unchanged |
| Count-bearing-only precision | 158/216 (**73.15%**) | 144/198 (**72.73%**) |
| Carrier assertion coverage | 154/635 (24.25%) | 140/635 (22.05%) |
| Carrier MAE | 0.266 | **0.243** |

This initial interpretation is retained for audit but is **not current**. It
mistook the count-bearing-only diagnostic for the repository's differently
denominated counted-extra precision metric. The correction above supersedes its
gate decision. Exact score, telemetry, and gate record:
`benchmarks/codex_paper_eval/runs/20260813_gold120_current/`.

### 2026-08-13 — First paired rescore of the phenotype-null arc (not a headline)

The re-extraction that `RECALL_STATUS.md` and `TASKS.md` have been demanding
since the 2026-08-12 phenotype-null commit finally ran, together with five
defect fixes in that same arc (see `PROTOCOL_CHANGELOG.md` for the defects).

**Instrument.** The 73-paper cardiac arm of `benchmarks/curated_extraction_eval`
through full default `gvf-run` (79 min), scored against the frozen gold subset.
Both sides scored by the **current** scorer: the committed
`expected_baseline.json` was written by the *old* scorer, so comparing to it
would have conflated the scorer rewrite with the extraction rewrite. The "before"
column is the four canonical baseline DBs re-scored on the same fixture today.

| Metric | Before | After |
| --- | ---: | ---: |
| Unique variants | 1584/1721 (92.0%) | 1559/1721 (**90.6%**) |
| Variant rows | 2435/2714 (89.7%) | 2391/2714 (**88.1%**) |
| Affected | 5233/5615 (93.2%) | 5181/5615 (**92.3%**) |
| Matched-row MAE — carriers | 0.393 | **0.168** |
| Matched-row MAE — affected | 0.322 | **0.270** |
| Matched-row MAE — unaffected | 0.377 | **0.201** |
| End-to-end error — carriers | 1.289 | **1.654** |
| End-to-end error — affected | 1.196 | **1.749** |
| End-to-end error — unaffected | 0.285 | **0.137** |
| Counted extras on gold PMIDs | 850 | **552** |
| Counted precision | 74.1% | **81.2%** |

**The two error families move in opposite directions, and that is the finding.**
Matched-row MAE falls by 16–57% while end-to-end carrier/affected error *rises*.
The pipeline stopped guessing: when it asserts a count it is far more often
right, and it fabricates a third fewer counted rows, but it declines to answer
more often, so total error against gold rows grows. Reporting only matched-row
MAE would have shown a clean win that is not the whole truth. Whether
null-instead-of-a-guess is worth that coverage cost is a deliberate policy
choice; it currently exists as a side effect rather than a decision.

**The recall delta is mostly acquisition, not parsing.** Of 96 lost rows, 57 come
from three papers with recorded source failures — RYR2 33606749
(`source_missing_or_stub`, zero bytes on disk), SCN5A 26746457 and 19251209
(`source_missing_table_bodies`) — and only 2 rows classify as `count_semantics`.
Real gains landed too: RYR2 19398665 went 2→26 matched rows as the figure path
began contributing.

**Why this is not a headline and not a gate pass.** The KCNQ1 arm was degraded
(`fetch_paywalled.py` exited 1), and the before side is an older DB rather than a
re-extraction with the old code, so the recall delta is not a controlled A/B.
`expected_baseline.json` was deliberately not rewritten and the four-gene
canonical headline in `RECALL_STATUS.md` is unchanged. Gate 1 (`gold_50`) in
`TASKS.md` remains the next step.

### 2026-08-13 — Gate 1 approved; 120-paper gold comparison retained

The lead approved advancing the current extraction system after reviewing the
paired-rescore MAE and counted-precision evidence. The requested 100--150-paper
scale refers to comparison with the already manually curated cardiac gold. The
prior `cardiac_120` manifest was a reviewer queue and only 52/120 attempts had
manual-gold rows, so it was replaced before measurement by `gold_120`: a fixed,
gold-value-blinded sample of 30 source-available, count-eligible papers per
cardiac gene (120 attempts / 116 unique PMIDs; seed 2026081301). The experimental
BMPR2, BRCA1, and BRCA2 queues remain 50, 50, and 46 papers, including all six
BRCA2 provenance exclusions. A mistakenly widened BMPR2 100-paper run and an
initial two-paper probe of the old cardiac reviewer manifest were stopped during
extraction as soon as their scope mismatches were identified; neither was scored
or published. Their partial artifacts remain under
`results/bmpr2_brca_300_20260813/` and `results/cardiac_120_20260813/` as aborted
historical runs. No headline metric changed.

### 2026-08-12 — Phenotype-null and figure-link hardening (rescore required)

Deterministic table paths now preserve NULL for affected/unaffected partitions
the source did not assert, while explicit source zero remains an observation.
The scorer carries assertion presence through aggregation and adjudication so
zero and missing are no longer conflated. Figure observations now enrich a
compatible existing variant-paper link, fail closed on ambiguous multi-cohort
parents, and retain role/locator provenance pending trust review. No live
extraction or headline rescore was performed. The 2026-07-12 headline therefore
remains the pre-correction baseline; acceptance requires a paired live cardiac
re-extraction and rescore with assertion coverage, matched-row error, and
end-to-end error reported together.

### 2026-08-11 — Full reviewer tier expanded to BMPR2, LMNA, and TTN

The full private-review tier now includes every populated reviewer workspace:
BMPR2 retains its existing 50-paper cohort, while LMNA and TTN are narrowed
from their 99-paper temporal snapshots to 50 each by descending count-bearing
evidence, descending total variant evidence, then PMID. Neither narrowed
workspace had reviewer adjudications or current gold records. The full tier is
now 546 gene–paper attempts / 507 unique PMIDs across 11 gene–disease pairs.
The live queues match the committed manifests, every paper has stable review
order and a non-empty paper-specific reviewer summary, and the LMNA/TTN source
runs pin the selection policy and manifest checksums.
This supersedes the earlier 396-attempt cohort definition without changing the
50-paper scored gate, the 120-attempt cardiac gate, or any recall metric.

### 2026-08-11 — Evaluation rollout consolidated into three canonical tiers

The active cohort surface was consolidated without changing predictions,
answer keys, or headline metrics: 50 gold-scored gene–paper attempts, then 120
cardiac reviewer attempts (98 unique PMIDs), then 396 attempts in the eight
established private-review workspaces (357 unique PMIDs). The remembered
“about 340 papers” referred to unique articles; repeated papers under different
gene–disease pairs remain distinct extraction/review attempts. BMPR2 and the
99-paper LMNA/TTN temporal runs remain historical experiments rather than
silently widening the production cohort. Exact manifests and derivation tests
live in `benchmarks/evaluation_tiers/`. No model calls or rescoring occurred.

### 2026-08-11 — Active cohort narrowed to the two collaborator-reviewed BRCA2 papers

A provenance audit found that only BRCA2 PMIDs 26833046 and 26848529 in the
eight-paper override were traceable to lead-approved Variant Browser
adjudications by Nate. PMIDs 10398279, 15365993, 18489799, 21356067, 22655046,
and 25802882 originated in an internal candidate re-derivation/review pass.
Those six were excluded from active Variant Browser review publishing and new
strategy-comparison membership. Four were present and removed from the live
BRCA2 snapshot; two were already absent. Nate's two papers and all 87 associated
current gold records were retained. The dated 56-paper locked predictions and
machine-readable metrics remain unchanged for reproducibility.

The replacement active manifest is the cardiac 48 plus the two approved BRCA2
papers (50 total). A PMID-filtered projection of the same locked predictions
gives carrier MAE 0.9058→0.0608 (298→20 absolute error over 329 supplied
counts), with count recall 32.64%→32.67%. Strata are essential: cardiac is
0.0491 (16/326), whereas BRCA2 collaborator-2 is 1.3333 (4/3). No model calls
were made for this scope correction. Full methods and arithmetic:
`benchmarks/count_semantics_eval/METHODS_20260811_COLLABORATOR_GOLD.md` and
`runs/20260811_collaborator_gold_50/metrics.json`.

### 2026-08-10 — Count-scope adjudication repaired the historical 56-paper error estimate; extraction predictions unchanged

This was a prospective-policy, locked-prediction error analysis, not a fresh
extraction experiment. The fixed cohort was the 48 cardiac papers in the
standard paper-eval set plus 8 BRCA2 curated-benchmark papers. Cardiac rows used
manual gold; BRCA2 used a mixed-provenance internal answer key. The 2026-08-11
audit subsequently removed six BRCA2 papers from active membership.
Prediction files were SHA-256 locked before adjudication, and all 378 supplied
carrier predictions remained byte-identical.

The intervention had six parts: define the count scope prospectively; rank
count-bearing multi-cohort/large-count claims; build compact source cards;
review the largest conflicts with GPT-5.6 Luna at xhigh; centralize the
fail-closed `gold_v2` loader used by both scorers; and audit 15 randomized,
prediction-blind source cards independently with Grok 4.5 High, Gemini 3.1 Pro
High, and Claude Fable 5 Max. The three reviewers agreed on all five headline
carrier totals and six corrected controls. Two additional controls deliberately
retained gold=3 against prediction=2, and one corrected affected count worsened
affected error by one, guarding against prediction-following adjudication.

Measured carrier MAE on the same 378 supplied counts changed from **0.8148 to
0.0794** (absolute error **308 to 30**, −90.3%). Count recall was essentially
flat (**34.02% to 34.05%**) because predictions were not added; the small change
comes from excluding one duplicate gold row. Affected MAE changed
**0.7869 to 0.7902** and unaffected MAE **0.1802 to 0.1307**. The selected
cardiac-48 slice ended at 0.0491 carrier MAE (16/326); BRCA2 ended at 0.2692
(14/52) and remains provisional.

What worked: compact ambiguity cards, explicit current-cohort rules, exact
closed-vocabulary v2 status handling, notation-aware evidence routing, and
blind negative controls. What did not work: broad missing-count recovery
grounded **0/162** completed gap checks while consuming 153,010 tokens, and two
blind cards lacked decisive source excerpts. Seven compact Luna calls consumed
27,682 tokens and 80.9 seconds, but the workloads are not an apples-to-apples
efficiency comparison.

The scorer/claim-verification changes are in `4a23b42` and are contained by
`origin/main`. The Luna route itself remains shadow-only, and compact
verification still runs before downstream figure/recovery layers; therefore
this result does not replace the canonical four-gene production headline.
Design, equations, locks, decisions, negative results, and reproducibility
pointers are consolidated in
`benchmarks/count_semantics_eval/METHODS_20260810.md`.

### 2026-08-10 — Strategy-comparison standard extended to 56 papers (+BRCA2 arm); first BRCA2 production baseline
Benchmark-surface change; the four-gene canonical headline is untouched. The
blinded paper-eval standard (`benchmarks/codex_paper_eval/`) grew from the
fixed cardiac 48 to a **56-paper manifest**
(`highcarrier48_plus_brca2_20260810.tsv`): the frozen 48 plus the 8 BRCA2 gold
papers from the curated extraction benchmark. The harness gained per-gene gold
resolution — cardiac genes score against the manual gold standard, BRCA2 falls
back to the curator-adjudicated `gold_overrides` answer key, and the resolved
path is recorded per gene in `selection.json`/`report.json` (`gold_sources`) so
the provenance difference is never silent. Seeded random sampling stays
cardiac-only. Schema-1 (production-import) predictions now lock without
per-paper token telemetry, which gvf-run structurally does not aggregate;
schema-2 harness-native runs still require it.

First BRCA2 production baseline (`runs/20260810_brca2_8_production`, same
calibrated replay as the cardiac `20260726_fixed48_production`): **variant
recall 80.9%** (89/110; cardiac production is 78.8%), paper-only view 73.6%.
Precision 18.2% is dominated by the two known subset-gold papers (26833046,
26848529 — 198 of 400 FPs; curator follow-ups open), so treat BRCA2
precision/count-MAE as provisional. Counts: carriers 41% supplied / MAE 0.31,
affected 34% / 1.38, unaffected 34% / 0.32 — same commit-rarely-but-accurately
shape as cardiac. Detail: `runs/20260810_brca2_8_production/COMPARISON.md`.

### 2026-07-20 — Trust, provenance, and gold-integrity hardening (#161–#165); recall headline unchanged
A run of protocol changes that harden **trustworthiness and BRCA-readiness**
rather than four-gene recall. All were designed gold-free / additive /
scorer-invariant, so the canonical four-gene headline is **unchanged**
(86.2% unique-variant, 81.2% variant-row, 0.614 carriers MAE — see the trajectory
table). The point of this arc is *honesty and generalization*, not a recall
number: it makes the pipeline safe to trust on non-cardiac genes (BRCA1/BRCA2)
where the old cardiac-tuned heuristics silently fabricated penetrance.

- **#161 — Live gold sync from Variant_Browser.** GVF now pulls the human-
  adjudicated approved gold from the Azure-backed review machine API into a
  versioned SQLite cache (immutable snapshots, per-sync change log,
  reviewer/approver identity, reversible exclusions; native JSON fields
  normalized on ingest). Raw/disputed/withheld/stale/checksum-invalid inputs
  fail closed.
- **#162 — Generalized table-role validation.** The table router now rejects
  row-ID, population-frequency, cohort-denominator, and clinical-measure columns
  by class (not per-gene) and records the selected column + count type — directly
  attacks "table numbers pulled from the wrong column."
- **#163 — Versioned Azure gold snapshots + scoring tiers.** Cardiac / all /
  noncardiac scoring tiers; required-sync scoring reads the selected tier for
  recall/precision/MAE/RMSE; hardened tier filtering + bulk exclusions.
- **#164 — Source-grounded per-paper final checks.** The default-on Step 3.8
  final check must quote the source; returned quotes are programmatically verified
  against the paper (same-table validation for header/row fragments) before a
  finding can affect trust. Step 3.9 deterministically composes only source-
  verified objective count/phenotype contradictions; weak "unsupported count"
  findings stay advisory and never mutate raw counts.
- **#165 — Study-design-aware counts + provenance honesty** (the coworker BRCA
  critique). `trust_gate` **tg3**: `negative_count`, full-partition arithmetic
  (a fully-specified affected+unaffected+uncertain that ≠ total is quarantined),
  and `implied_unaffected_zero` — a *derived* 100%-penetrance claim
  (unaffected=0 with affected=total, unsourced) in an affirmatively
  cohort/biobank/case-control/cascade study is soft-quarantined on the
  **unaffected field only** (dormant on proband/unknown-design papers, so 0 new
  cardiac quarantine, verified on KCNH2). Plus `variant_papers.source_notation`
  (verbatim as-reported variant string beside the normalized IDs); a prompt that
  decouples count from phenotype (label AFFECTED only when disease is stated;
  never assume an unaffected count); a penetrance/segregation/cascade discovery
  lane + priority signal; and the vertical supplement-catalogue parser no longer
  asserting "pathogenic" per row. Full offline suite 1214 passed; CI green.

Deferred to a **measured** pass (flip source-free `unaffected=0`→null and
evidence-gate `affected=patient_count` in the deterministic parsers; flip the
count-classifier/guard defaults off→flag) — believed gold-neutral but they touch
always-on, *scored* count defaults, so they need a live cardiac re-extraction +
rescore first. Tracked in `docs/AUTONOMY_ROADMAP.md`.

### 2026-07-12 — Four-gene idempotent supplement reconciliation and SCN5A gain
Closed the local supplement-fold gap across KCNH2/KCNQ1/RYR2/SCN5A from 289
papers to **0** (1,577/1,577 papers with convertible supplements folded).
Elsevier acquisition now reconciles individual `mmc` files rather than treating
"any supplement exists" as paper completion, reuses matching supplements across
sibling gene corpora, folds changed sources immediately, and keeps supplement
recovery separate from article-body recovery. The live pass recovered 64 missing
files across 49 papers and updated 427 contexts without re-downloading bodies.

Acceptance-gated re-extraction promoted only SCN5A (PMIDs 12193783, 20541041,
23973953, 29017927): SCN5A unique variants **1020→1027**, rows **2432→2461**,
and PMID coverage **619→622**, with carriers/affected MAE slightly improving and
unaffected MAE holding. KCNH2, KCNQ1, and RYR2 candidates regressed recall or
count coverage and were not promoted. Fresh aggregate: unique variants
**2596/3010 (86.2%)**, rows **5546/6833 (81.2%)**, carriers MAE **0.614**.
Grok and Claude CLI reviews informed the per-file completion model and the
supplement-only acquisition route.

### 2026-06-22 — Annotation-table guard: blank-cell wrong-gene + gnomAD-count-as-carrier (PMID 33013630)

**Garbage-source warning.** PMID **33013630** ("Are Variants Causing Cardiac
Arrhythmia Risk Factors in SUDEP?") is a **hypothesis/review**. Its Table 1 has
columns `Gene | Variant | gnomAD allele count | SIFT | PolyPhen-2 | References`
and **zero patient/carrier data** — it is variant *annotation*, never primary
carrier evidence (extraction rule 4). Treat any carrier counts sourced from this
PMID as invalid.

Two deterministic `pipeline/table_router.py` bugs were turning this annotation
table into carrier evidence (both hit *every* gene, including no-gold targets):
1. **Blank Gene cells (markdown rowspan) were not inherited.** Gene-grouped
   tables name the gene once and leave continuation rows blank; those rows were
   stamped with the *target* gene instead of the real one. So HCN4 `Val759Ile`,
   HCN2 `Pro802Ser`, SCN5A `Ala572Asp`, AKAP9 `Arg2607Gly` … all leaked into
   KCNH2. The flagged symptom was **KCNH2 `Val759Ile` = 870 affected / 0 / 870**
   — the 870 is HCN4's *gnomAD allele count*.
2. **`_looks_like_row_level_clinical_list` fired on the caption's clinical cue
   alone**, minting 1 inferred carrier per annotated row even with no subject
   column.

**Fix (landed in code, gold-free, generalized to the whole class):**
forward-fill blank Gene cells before the gene-filter; and an **annotation
table** — one whose only quantitative columns are variant annotation
(population-frequency allele counts like gnomAD/ExAC/TOPMed **or** in-silico
predictors like SIFT/PolyPhen/CADD/REVEL) with **no subject column** — no longer
infers a carrier per row, regardless of a clinical caption or "score" cue. The
current extractor returns **0** variants from this table (and from a
prediction-score-only table such as `REVEL score | CADD score`). 4 regressions
added in `tests/unit/test_table_router.py` (incl. a guard rail that a genuine
proband list with no annotation columns still infers); full offline suite green
(958 passed / 53 skipped). The literal 870 was already nulled in canonical
`02_strict` by the earlier count-outlier guard; the code fix stops the whole
class at the root.

**Did a good source do a good job? Yes — dropping 33013630 loses zero signal.**
The only genuine KCNH2 variants this review name-drops (`Arg744*`, `Gly924Ala`,
"identified in a SUDEP patient") are redundantly and far better captured by real
primary sources GVF already extracts well, e.g.:
`Arg176Trp` 25+ papers (incl. 19160088=128, 16754261=32/60, 21244686=88),
`Arg744*` (11802537=7, 22338672=1/3), `Gly924Ala` (**38370760**=65),
`Gly749Ala` (38370760=41), `Arg1047Leu` (38370760=67, +10 more). PMID
**38370760** is a standout high-quality cohort table GVF extracted cleanly with
affected/unaffected splits.

**Residual (not yet cleaned — see `TASKS.md`):** the May-29 canonical KCNH2
(21 vars from this PMID: 9 sole-source FP incl. `Val759Ile`, 12 with bogus
counts) and KCNQ1 (14 vars) DBs still carry this contamination; SCN5A/RYR2 are
clean. The live `Val759Ile=870/0/870` row also persists in the stale
`validation_runs/20260518` copies. Un-ingesting PMID 33013630 (the fixed
extractor yields 0 for it) is a precision/FP cleanup deferred to the next refresh.

### 2026-06-22 — Docs source-of-truth cleanup
Consolidated the recall docs so there is no competing "next-run plan" surface:
`TASKS.md` now owns the active Exact-Match Recovery Plan, `RECALL_STATUS.md` is a
short current metrics/failure-split snapshot, and this file owns the dated recall
trajectory. The embedded `RECALL_STATUS.md` session log and stale 2026-05-29
plan were removed from the live status file; their substantive benchmark history
is represented by the timeline entries below.

### 2026-06-12 — Targeted lands extended to the remaining three genes
After the KCNQ1 headline land, `scripts/targeted_land.py` was run on the other
three canonical genes. SCN5A PMID 19716085 promoted a cleaner re-extraction
(+1 unique variant, +4 variant rows; SCN5A **86.3%→86.4%**, and its carriers MAE
improved **0.489→0.454** as over-counted rows were replaced). KCNH2 PMID 32681117
and RYR2 PMID 24136861 each found a candidate whose re-extraction *held* recall —
promoted as cleaner, non-regressive (gold-gated, no recall/row/MAE regression).
Four-gene aggregate **86.0%→86.1%** unique (2591/3010), rows **80.7%→80.8%**
(5518/6833), carriers MAE **0.634→0.615**. Each gene backed up as
`<GENE>.db.before_targeted_land.db` before promotion.

### 2026-06-12 — PDF-linearized table reconstruction → iter-2 quality gate/selector → fast targeted land
PDF supplement tables that fold in linearized (one cell per line) are now
reconstructed into delimited rows before extraction (`ExpertExtractor.
_augment_pdf_linearized_tables`), with a post-extraction cDNA↔protein backfill.
The replay gate/selector were generalized to be **no-gold-safe**: the gene-scoped
deterministic table parse is the structural baseline — a fewer-row re-extraction
is accepted only if it does not lose paired variants and covers ≥85% of the
deterministic parse's positions (over-extracted prior → drop the excess; lossy
under-extraction → reject). Root insight: the blocker was never that extraction
dropped protein — the `refresh_run_db` count gate preferred a stale over-counted
cDNA-only extraction (KCNQ1 30758498: 182 rows, recovers 34/183 gold) over a
clean paired one (100 rows, recovers 85/183). Validated cross-gene (fires on
KCNQ1 + an untuned RYR2 paper; declines SCN5A 0/120). `scripts/targeted_land.py`
lands a single paper's win in **minutes not ~1h** (scan → bridge only candidates
→ re-extract+gate → surgical layer-preserving inject → gold-gated promote).
Landed KCNQ1 PMID 30758498 → **KCNQ1 87.9%→90.5% unique** (crosses target),
aggregate **85.5%→86.0% / rows 79.9%→80.7%**, gated (recall+row-recall+MAE),
backup `KCNQ1.db.before_targeted_land.db`. 787 unit tests incl. a non-cardiac
generalization test. Codex (gpt-5.5, read-only) reviewed the gate/selector design.

### 2026-06-12 — Source-layer backfill + cheap notation junk gate
Added a shared source-layer classifier, explicit `variant_papers.source_layer`
backfill for the four canonical DBs, and a scorer/migration reject for obvious
figure/regex-table junk: gene-symbol-as-variant, <=2-character protein strings,
and residue prose. Backups were written as
`.before_source_layer_20260612_093534` before any DB mutation. Re-scoring the
canonical DBs and the backups (fallback inference, no `source_layer` column)
produced identical aggregate recall/MAE/precision blocks, confirming the scorer
prefers the explicit column without changing semantics.

Recall and MAE were unchanged: uniqV **2572/3010 (85.4%)**, rows
**5423/6833 (79.4%)**, carriers MAE **0.614**, affected MAE **0.523**,
unaffected MAE **1.219**. The junk gate dropped 64 raw rows before scoring
(41 extra-on-gold-PMID rows, 8 counted extras). Headline counted precision is now
`5423/(5423+1660)` = **76.6%**; the loose raw gold-PMID upper bound is
`5423/(5423+13642)` = **28.4%**. The old `manual_or_legacy` bucket is now
`llm_text`.

### 2026-06-11 — Canonical rescore + precision decomposition
Fresh `scripts/run_recall_suite.py` run against the four canonical DBs confirmed
the current aggregate: uniqV **2572/3010 (85.4%)**, rows **5423/6833 (79.4%)**,
PMIDs **1274/1502 (84.8%)**, carriers MAE **0.614**, affected MAE **0.523**,
unaffected MAE **1.219**. Precision-vs-gold-PMIDs remains a false-positive
upper bound at **28.4%**, but decomposition shows only **1668/13683** extra
rows on gold PMIDs carry counts, giving `precision_vs_counted_gold_pmids`
**76.5%**. Added `variant_papers.source_layer` plus per-layer precision output
so ClinVar/PubTator/figure/linker rows are auditable without a manual sample.
Confirmed the MAE non-regression land gate already existed; fixed the stale
history note that claimed it was still missing.

### 2026-06-06 — Duplicate-penetrance idempotency fix (carriers MAE 1.274→0.614, recall preserved)
The 2026-06-05 `refresh_recall` re-extraction lifted recall (uniqV 83.8→85.4%,
rows 78.3→79.4%, patients 80.6→82.1%) but **spiked aggregate carriers MAE
0.882→1.274** — 94% of it KCNQ1, ~98% of that a single paper (PMID 32893267, the
Lahrouchi 2020 LQTS exome cohort). Root cause was **not** junk variants from
no-patient-data papers: the migration wrote one `penetrance_data` row per
supplement-table appearance of a variant (that paper lists each variant across
Table S4 + Table S14 × two ACMG schemes), and the scorer **sums** the rows linked
to each (pmid, variant) — so V254M scored 4× gold (100 vs 25). The count-outlier
guard can't catch it (uniform inflation lifts the per-paper median too, so nothing
trips the >10× rule) and the land gate checks only unique-variant recall, never
MAE, so it promoted silently.
- Fix is idempotency at the write path: exact-duplicate insert guards in
  `migrate_to_sqlite.insert_variant_data` for penetrance/individual/functional/
  phenotypes/variant_metadata, so re-migration or a variant repeated across table
  cells never writes a second identical row. Back-fill for DBs built before the
  guards: `dedup_existing_rows` (`scripts/dedup_db.py`).
- Canonical DBs deduped (backups `<GENE>.db.before_dedup.db`): aggregate carriers
  MAE **1.274→0.614** (below the 0.882 baseline), recall unchanged. Per-gene
  carriers MAE: KCNH2 0.860, SCN5A 0.489, KCNQ1 3.61→**0.897**, RYR2 0.323.
  Removed exact dups: penetrance 6883, phenotypes 7431 across the four DBs.
- Tests: `tests/unit/test_migrate_idempotent.py` (4) + the migration/scoring suite
  (125) green.
- Follow-up landed 2026-06-06: `refresh_recall` now has an MAE non-regression
  land gate (`_promotion_decision`) requiring recall hold plus no carriers,
  affected, or unaffected MAE regression before promotion.

### 2026-06-05 — Supplement acquisition: the real recall lever (+1.6pp uniqV)
Recon (`docs/SUPPLEMENT_ACQUISITION_PLAN.md`) re-framed the problem: the
Cloudflare-blocked publishers (Karger 0.3% / Sage 0.0% of the gap) are
near-irrelevant; the leak was that the **Elsevier full-text API fetches body
only**, dropping the `mmc` supplement mutation tables (~31% of the gap), and the
re-fold layer that makes on-disk supplements visible to Tier-3 was never wired
into the refresh path.
- Wired the supplement re-fold into `refresh_run_db` (Phase 0, no-network) + a
  nested-zip fix + a fold-gap QC counter.
- Added `ElsevierAPIClient.download_supplements` (mmc refs from the authenticated
  XML → open `ars.els-cdn.com` CDN) + `scripts/fetch_elsevier_supplements.py`.
- **Solved the #1 single blocker, SCN5A 29325976** — flagged 2026-06-01 as 87
  missing variants behind a "Cloudflare-blocked `mmc1.docx`". The supplement was
  in fact openly served on the ScienceDirect CDN; fetching + folding +
  re-extracting recovered **64 of 87** variants.
- Batch over the addressable Elsevier papers: 16 gained supplements; re-extract
  recovered 103/126 of their missing variants (per-paper proxy: +175 total).
- **Landed in the canonical DBs** via surgical injection (preserving
  clinvar/pubtator/figure layer rows): aggregate uniqV **82.2%→83.8%** (+50),
  rows **75.5%→78.3%** (+190), patients 77.8%→80.6%, carriers MAE 0.910→0.882.
  Nothing regressed. Backups: `{gene}.db.before_supplements_20260605.db`.
- Figures (the prior day's experiment) gave **0 recall lift** — the gaps are
  splice/intronic variants in supplement tables, not figures. Honest negative.
- T6 (Wiley/Springer supplements) investigated → access-blocked: Springer API key
  is off in `.env`; addressable Wiley DOIs return TDM 403. Reactivates when keys/
  access are restored.

### 2026-06-04 — Corpus consolidation, dashboard, figure fetch (0 recall)
- Consolidated all fetched source into `corpus/<GENE>/<PMID>/` + INDEX
  (6,382 papers), idempotent never-downgrade builder, corpus-as-cache for
  `gvf-run`.
- Static provenance/coverage/adjudication dashboard (`gvf dashboard`).
- Caption-triage à-la-carte figure fetch: 313 OA figures fetched; **0 gold
  recall lift** — confirmed the recall gap is supplements, not figures.

### 2026-06-01 — Acquisition replay & no-gold source QC (per-gene diagnostics)
Per-gene acquisition experiments on separate run dirs measured KCNH2 ~84.15% and
SCN5A ~84.19% uniqV (figures-skipped) after fetched-source replay + ClinVar/
PubTator. Diagnosis flagged **SCN5A 29325976** as the single largest gap (87
variants in a blocked supplement) — solved 2026-06-05. These were diagnostic runs
on `01_off`/turnkey dirs, not all consolidated into the canonical DBs.

### 2026-05-29 — Matcher bridge + count-outlier guard (+0.8pp uniqV, MAE 1.055→~0.90)
- `cli/compare_variants.py`: cDNA↔protein substitution bridge (matches `Y51X`
  gold to `c.153C>A` extracted by implied codon).
- `pipeline/extraction.py`: scope gene-column-less tables by caption (stops
  multi-gene leakage; PMID 26669661 161→24 SCN5A).
- Count-outlier guard cleared study-wide-N values (96 across 28 PMIDs): MAE win,
  no recall change. KCNQ1 was the dominant offender (carriers MAE 2.916).

### 2026-05-26 — SUA source-replay sweep + rollback + recovery (62.6%→81.4% uniqV)
The single biggest jump. Acceptance-gated source replay over the
`source_unbound_available` pool (147 PMIDs), then per-PMID rollback (restore the
higher-variant extraction), DB rebuild, and DB-observed ClinVar/PubTator
recovery. KCNQ1 moved 43.1%→84.2% (LQT supplement parser + native table shapes);
KCNH2 68.3%→82.6%.

### 2026-05-21 — Elsevier insttoken unblock (foundational acquisition win)
Vanderbilt institutional `X-ELS-Insttoken` installed; 242 of 246 previously-
paywalled Elsevier articles across the cardiac genes returned full text. This is
the credential the current cardiac-gene baselines depend on.

### 2025-11-18 → 2026-05 — Pipeline build-out (pre-systematic-benchmark)
Discovery → harvest → Tier1/2 filter → Tier3 extraction → migrate → recall
scoring pipeline; publisher routes (PMC, Elsevier/Wiley/Springer APIs, browser
recovery); gold-standard builders; the recall suite + MAE foundation. Systematic
four-gene recall benchmarking began ~2026-05-25.
