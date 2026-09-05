# Recall Status

Last updated: 2026-09-05.

This file is the current measured recall snapshot. It intentionally does not
carry the active work plan or dated session log.

- Active forward plan and task checklist: `TASKS.md`.
- Benchmark/change trajectory: `docs/RECALL_HISTORY.md`.
- Re-run procedure: `docs/RECALL_REFRESH_RUNBOOK.md`.

No other doc should restate live recall tables. If a metric conflicts with this
file, this file is authoritative.

## Accepted locked snapshot (legacy linkage-assisted headline)

The accepted headline is the immutable, gold-blind
`benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118` lock.
It covers 118
gene-paper attempts / 114 unique PMIDs over KCNH2, KCNQ1, RYR2, and SCN5A and
was prediction-locked before the 632-row human cardiac gold values were read.
This historical table is retained for audit, but it is not the extraction-task
headline going forward: its 554 TP include ClinVar/PubTator citation linkage.
The paper-derived replay below (512/125/120) is the primary interpretation of
whether the protocol found variants in the papers.

| TP / FP / FN | Recall | Raw precision | F1 | Counted-extra P | Count-bearing-only P |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 554 / 283 / 78 | **87.658%** | **66.189%** | **75.425%** | **97.880%** | **94.595%** |

| Count | Supplied / gold | Coverage | Conditional MAE | RMSE |
| --- | ---: | ---: | ---: | ---: |
| Carrier | 207 / 632 | 32.753% | **0.193** | 0.742 |
| Affected | 56 / 632 | 8.861% | **0.321** | 1.102 |
| Unaffected | 28 / 631 | 4.437% | **0.321** | 0.779 |

### Latest phenotype-zero/provenance validation (not the headline)

The fresh gold-blind lock
`benchmarks/codex_paper_eval/runs/20260902_false_zero_recovery_gold118` applies
the same source-closed phenotype protocol to all 118 paper attempts. It closes
the recurring cross-stage failure in which a deterministic recovery could be
lost because later verification, trust, or refresh code did not recognize its
provenance. Code-owned stamps are now scrubbed from model JSON and may be added
only by the corresponding deterministic audit. The phenotype denominator also
comes from dedicated penetrance data before the ambiguous legacy patient-count
mirror.

| TP / FP / FN | Recall | Raw precision | F1 | Counted-extra P | Count-bearing-only P |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 543 / 275 / 89 | **85.918%** | **66.381%** | **74.897%** | **98.192%** | **95.413%** |

| Count | Supplied / gold | Coverage | Conditional MAE |
| --- | ---: | ---: | ---: |
| Carrier | 207 / 632 | 32.753% | **0.179** |
| Affected | 120 / 632 | 18.987% | **0.508** |
| Unaffected | 37 / 631 | 5.864% | **2.081** |

This validates the repaired data path but does not replace the accepted
headline: identity recall is 1.740 percentage points lower, and fresh stochastic
extraction supplied fewer carrier and affected values overall. The acceptance
cases all survive extraction, verification, trust, database migration, lock,
and score: RYR2 W4645R **4/2/2**, C2277R **8/7/1**, G357S
**185/97/62** (+26 uncertain), P2328S **62/17/42** (+3 uncertain), and KCNH2
Y652X **6/6/NULL**. This is evidence for the protocol path, not evidence that
the overall coverage problem is solved.

The locked figure has 1,085 matched phenotype rows. Raw null is evaluated as
zero only in that visualization: **423 affected and 505 unaffected** values are
null, not literal model zeros. There are just **11** stored zeros (1 affected,
10 unaffected). Eight agree with zero-valued gold. The three positive-gold
disagreements are source-definition conflicts: KCNQ1 V141M PMID 28491547 is a
three-child AF/SQTS series (**3/3/0**); KCNH2 C566R PMID 30246897 is one carrier
with marked QT prolongation but no syncope/palpitations/arrest (**1/1/0** under
the paper-target LQTS definition); and KCNH2 P926AfsX14 PMID 20181576 is a
four-carrier family whose ECG and symptom subsets do not match the legacy gold
partition. None is a source-positive count silently rewritten to zero. Grok 4.6
`xhigh` independently audited these cases and the high-value null rows.

### Current two-cohort failure diagnostic (2026-09-04)

The current diagnostic reads the latest completed locks near the requested
120-attempt scale. The identity columns are **not a head-to-head comparison**:
Gold 118 is scored from its legacy trusted projection, which can include
ClinVar/PubTator linkage, while Mixed 120 uses the paper-derived primary lane.
They remain useful within-lock snapshots and failure inventories.

| Cohort | Scored lane | TP / FP / FN | Recall | Carrier supply | Conditional MAE | End-to-end MAE |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Gold 118 | Legacy trusted projection | 543 / 275 / 89 | 85.918% | 207 / 632 (32.753%) | 0.179 | 1.112 |
| Mixed 120 | Paper-derived primary | 261 / 141 / 123 | 67.969% | 130 / 384 (33.854%) | 0.231 | 1.359 |

Recall is principally source-bound: acquisition plus figure/unsearchable-
notation failures account for **66/89 (74.2%)** Gold 118 false negatives and
**112/123 (91.1%)** Mixed 120 false negatives. Only 23 and 11 misses,
respectively, are reachable by downstream parser/model/projection changes on
the bytes already acquired.

Carrier error is also an omission problem. Gold 118 has 703 absolute-error
units: 540 (76.8%) come from matched variants with no supplied count, 126
(17.9%) from missed identities, and 37 (5.3%) from wrong supplied values. Mixed
120 has 522 units: 290 (55.6%) from missed identities, 202 (38.7%) from omitted
counts, and 30 (5.7%) from wrong supplied values.

The largest false-negative attempts are KCNH2 PMID 29650123 (21; acquired bytes
omit 15 substitutions and six notations are not searchable), SCN5A 15898185
(10; nine behind a stub), and RYR2 27114410 (9; no source) in Gold 118; and
SCN5A 25163546 (20; stub), SCN5A 20031634 (13; no source), and BRCA1 10528853
(13; figure-only) in Mixed 120. The ranked worklists and row-level traces are in
[`evidence/current_protocol_diagnostic_20260904/`](evidence/current_protocol_diagnostic_20260904/).

Three locked paper-agnostic reading/count comparisons were tested. V1 improved
Mixed 01 carrier end-to-end MAE by 0.731 but raised recall only 0.826 percentage
points, below the identity gate; V2 changed recall by -0.330 points on Mixed 02
and -0.521 points on Mixed 120, so it did not replicate. A temporary-copy replay
of the `tg6` source-backed table-outlier exception found zero eligible rows and
zero aggregate changes in these two locks. No candidate is promoted. The formal
`gold_120b` replication tier is registered but still unrun and has no metric.

### Affected/unaffected source diagnostic (2026-09-05; no new score)

The [22-paper failure panel](evidence/phenotype_failure_panel_20260905/README.md)
examines 30 attempts selected from four opened locks: 14 papers have some A/U
capture, including two exact controls, and eight have none. It distinguishes
missing components, acquired but unconsumed clinical tables, lost count facts,
and disagreement over symptoms, diagnosis, ECG or follow-up endpoints. Its
deliberate enrichment is not a population estimate. The earlier identity-source
ceiling does not measure completeness of phenotype evidence.

Full bodies for SCN5A 20031634 and 25163546 are now recovered from author
repositories and synced to the external corpus, with incomplete supplement
status retained. PMID 20031634 Table 1 was visually checked: 115 carriers split
into 54 ECG-positive, 55 not ECG-positive and six undetermined. The source
partition differs from this lock's gold and requires endpoint adjudication.
These are source acquisitions, not replay results; all locked scores and the
accepted headline above remain unchanged. The panel report retains before/after
source provenance, row evidence and both CLI reviews; `TASKS.md` owns the next
implementation order.

The subsequent [general repository-fallback validation](evidence/repository_fallback_20260905/README.md)
reproduced both requested downloads through the normal harvester. A fixed
16-paper network panel (those two plus 14 other ranked gold papers) downloaded
3 bodies: 20031634, 25163546 and 29650123. This is **1/14** additional-paper
download successes, not 3 new corpus recoveries; 29650123 already had a body
and still lacks its mutation roster. The corpus builder upgraded 25163546's
body and retained the other existing sources. That first validation did not score extraction. The subsequent targeted
calibration below measures the utility of both recovered bodies.

### Targeted source and reader calibration (2026-09-05)

The [implementation campaign](evidence/recall_campaign_20260905/README.md)
uses opened papers and preserves the accepted headline. Three fresh locked
runs cover 15 attempts over nine unique papers; overlapping arms are not pooled.

| Calibration | Attempts | TP / FP / FN | Test cost proxy |
| --- | ---: | ---: | ---: |
| Initial prototype | 10 | 348 / 41 / 12 | $2.6282 |
| Prototype restricted to final ablation attempts | 4 | 55 / 11 / 6 | Included above |
| Final code, identical source and asset hashes | 4 | 55 / 6 / 6 | $2.0381 |
| Second recovered body, 25163546 | 1 | 0 / 0 / 20 | $0.0254 |

SCN5A **20031634** changes from its old paper-derived **0 / 0 / 13** to
**11 / 2 / 2**. Eight matched variants supply A/U, with zero conditional MAE
for both fields; five +1 carrier disagreements reflect unknown phenotype
included by the source but excluded by gold. Replaying the identical model
response isolates the insertion validator's one recovered row; the fresh run
retains that row through trusted projection. Three matched identities still
lack A/U, and two source notation ambiguities remain unresolved.

The final four-attempt run supplies affected/carrier values for 49/61 gold
rows and unaffected values for 30/61. Conditional MAE is 0.224 / 0.122 / 0.000
(affected / carriers / unaffected); end-to-end MAE is 0.967 / 1.098 / 0.197.
Count changes outside the insertion example also reflect fresh-model variation
and are not credited to this code change. All three speculative parser rules
were withdrawn after their precision or incremental-yield checks failed.

25163546's recovered body still lacks useful variant/count material: all 20
gold identities remain missed, despite successful PDF acquisition. Its all-miss
figure exposed and now verifies an empty-agreement reporting repair. Missing
supplements and unconsumed clinical tables remain the main unresolved work.
The [final difference figure](../benchmarks/codex_paper_eval/runs/20260905_mechanism4_final/figures/gold_difference.png)
and [all-miss figure](../benchmarks/codex_paper_eval/runs/20260905_repository1_final/figures/gold_difference.png)
retain these limitations. Test proxy total **$4.69165**, separate from
**$3.07193** reported Claude consultation and unpriced Grok/Agy CLI reviews.

### How to report these numbers (2026-08-25)

Two independent maximum-effort reviews (Grok 4.6 `xhigh`, GPT-5.6 Sol `max`)
audited the headline. Full writeup:
`docs/evidence/generalization_consult_20260825.md`. Three corrections bind here.

**1. Counted-extra 97.880% is not a conventional precision.** It is
`554/(554+12)`: every matched row in the numerator, only count-bearing extras in
the denominator. Report it as the Gate-2 score, or transparently as **2.17
counted extras per 100 matched rows**. The conventional count-bearing precision
is `210/(210+12)` = **94.595%**. Raw 66.189% is precision against an incomplete
fixture under unmatched=false; because some gold rows are themselves out of
scope, it is not a clean lower bound either. Never report 97.9% to a reader who
will hear "the system is 98% precise".

**2. The headline mixes paper extraction with database linkage.** The locked
projection includes ClinVar/PubTator rows. Locked-DB replay:

| Lane | TP / FP / FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: |
| Linkage-assisted (the table above) | 554 / 283 / 78 | 66.189% | 87.658% | 75.425% |
| Paper-derived only | 512 / 125 / 120 | 80.38% | 81.01% | 80.69% |

Linkage FPs: KCNQ1 79, SCN5A 72, KCNH2 7, **RYR2 0** — the reason RYR2's
per-gene row looks clean. Do not delete linkage rows; report the lanes
separately, paper-source primary for the extraction task.

**3. Conditional MAE is not the decision metric.** It is computed only over rows
the pipeline chose to answer, so abstaining improves it. End-to-end error
(already implemented at `cli/compare_variants.py:2828`) tells the other half:

| Field | Conditional MAE | End-to-end MAE |
| --- | ---: | ---: |
| carriers | 0.193 | **1.541** |
| affected | 0.321 | **1.528** |
| unaffected | 0.321 | **0.512** |

Never quote one without the other, and never quote either without coverage.
Carrier MAE 0.193 comes from only 23 errors across 207 rows and is near its
practical floor; the real deficit is coverage. The real non-zero carrier
coverage gap is **268 rows**, not the 425 implied by pooled coverage — SCN5A
32533946 alone contributes 83 gold assertions that are all explicit zeros, which
is the gold-zero/NULL convention gap rather than an attribution miss.

**Standing caveat on gene-generality.** `gvf_data/kcnh2_variant_aliases.json`
declares `"source": "Gold standard Excel + generated forms"` and is read at
runtime, so KCNH2 normalizes with access to the answer key. Re-matching all 554
locked pairs with the lookup stubbed gives 554/554 — no score-time result
depends on it — but extraction-time influence is unrecoverable. Do not cite
KCNH2 as evidence of gene-general behaviour until the map is rebuilt from public
resources.

The grouped/structural and source-recovery clean locks are retained as failed
validations, not candidate headlines.
`20260824_grouped_structural_gold118` improved raw/count precision
and all count errors but lost four TP (542/278/90, 85.759% recall).
`20260824_source_recovery_gold118` restored 546 TP and improved carrier MAE to
0.198, but carrier coverage fell 206→197 and count-bearing-only precision fell
93.694%→93.396%. No best-replicate selection is allowed.

The promoted lock repairs a separately verified AHA acquisition gap. PMID
15466642's complete publisher DOM contains a collapsed patient table that the
prior renderer discarded. Its source-bound production result is 8 TP / 0 FP /
0 FN, and the full fresh lock improves TP, FP, FN, raw/count precision, all
three count coverages, and all three conditional MAEs versus the postfix lock.
Full per-gene metrics, integrity hashes, reviewer disposition, and the active
`$100` ledger are in
`docs/evidence/gold118_aha_table_lock_20260824.md`.

The exact BRCA1 50 + BRCA2 50 + BMPR2 50 extraction remains a staged candidate,
not a metric cohort: only 3/150 papers overlap approved human gold and BMPR2 has
zero. Its precision, recall, and MAE are **undefined**, not zero.

**Its human curation was descoped by Brett on 2026-08-25** — the cardiac
gold-120 is the human-curated standard being worked against. Exact-150 metrics
therefore stay undefined indefinitely, and no human curation time should be
spent there without Brett reopening it.

The packets are parked, not deleted. Note that the original
`gold150_preregistered_20260824/` split is **defective** and marked superseded:
because assignment ranked `sha256("<seed>|<GENE>|<PMID>")`, six PMIDs were
calibration for one gene and holdout for another, putting five of the twenty
BRCA1 holdout papers into BRCA2 calibration, and seven of the eight cross-gene
PMIDs were bound to different source bytes per gene. If that cohort is ever
reopened, use `gold150_preregistered_20260825_amended/`, whose firewall is
verified by `scripts/audit_split_firewall.py`. Run that audit on any future
multi-gene packet tree.

### Acquisition ceiling of the named-variant gold (2026-09-03)

`scripts/recall_audit/gold_source_presence_sweep.py` classifies every gold row
against everything on disk for its paper (body, converted binary supplements,
article PDFs), blind to predictions. On the 1,422 runnable mixed-gold attempts
(7,107 rows) **15.8% of gold rows sit behind a hard acquisition ceiling**
(source absent, abstract-length stub, glyph-garbled PDF text, or a
substitution absent from every searchable byte) and **28.7% behind the wide
ceiling** (adds figure-only and non-searchable notation). Only 9 rows live in a
binary supplement the body lacks. The rows are concentrated in table-body
capture failures (KCNQ1 14678125 37/41, 17192539 51/57; SCN5A 21273195,
24631775, 11901046). Report tranche recall on all gold first and the
reachable-gold figure as a labeled diagnostic; never drop the two unknown
classes from a denominator. Details and the acquisition worklist:
`docs/evidence/gold_source_presence_sweep_20260903.md`.

Mixed-gold tranche 01 (paper-derived, frozen `506a949c` baseline, 49
attempts): **155 / 61 / 87**, recall 64.0%, precision 71.8%; 70 of the 87
misses are acquisition, 12 notation-unknown, 5 reachable
(`docs/evidence/mixed_gold_tranche01_20260903.md`). Tranche 01 is unusually
acquisition-bound (31% hard ceiling vs 15.8% suite-wide).

## Historical audit context

The following material preserves older gates and all-paper baselines for audit;
it does not override the authoritative lock above. The scorer now preserves
explicit adjudicated zero counts instead of conflating them with null.
Source-reachable strata are secondary reader diagnostics; ALL GOLD remains the
primary turnkey acceptance denominator.

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

## Shared reliability change set (2026-09-02) — headline unchanged

A repo-wide parser/data-quality repair landed on 2026-09-02. **It does not move
any number in this file.** The cardiac curated gold gate is *bit-identical*
before and after (`diff` of the two reports is empty): pmids 72/73, variant rows
2436/2714, unique variants 1585/1720, carriers MAE 0.393, affected 0.321,
unaffected 0.377. It adds 157 tests and the offline suite passes with zero
failures.

A separate session was editing the phenotype-count and prompt path concurrently
(`pipeline/count_classifier.py`, `phenotype_count_guard.py`, `claim_verifier.py`,
`prompts.py`, `benchmarks/codex_paper_eval/run_eval.py` and their tests). That
work is preserved and untouched here, but it means the working-tree suite total
covers both change sets, and the fresh BMPR2 run below executed while those
files were being edited. Python binds modules at process start, so the run used
the code as of its own start time and is internally consistent; its recorded
fingerprint should not be re-derived from the current tree.

Seven shared, gene-agnostic defects were repaired: nested archives that reached
no converter, four mutually inconsistent text normalizers, cross-route identity
keyed on the raw notation string, denominators reported from step intent rather
than the finalized paper list, unpersisted count roles and zero provenance,
incomplete protocol provenance, and a re-fold that could replace real converted
text with a converter placeholder. Full evidence, the measured before/after, and
three adversarial reviews:
[`evidence/shared_pipeline_reliability_20260902.md`](evidence/shared_pipeline_reliability_20260902.md).

Three consequences bind how earlier artifacts may be read:

1. **Source-completeness artifacts written before this date can be wrong about
   which papers had full text.** The completeness report was written from the
   download step's return value, so a run served from the corpus cache reported
   every paper as abstract-only while its own per-paper records said full text.
   One audited run shipped four artifacts disagreeing about the same 50 papers.
   Re-derive that field with `pipeline/source_ledger.py` rather than quoting it.
2. **`single_carrier_papers` in every pre-2026-09-02 artifact is meaningless.**
   The predicate read two keys that do not exist in the extraction schema, so it
   was constant and flagged every paper. On the audited run the true figure is
   26, not the reported 50.
3. **Duplicate identities in existing databases are largely a linkage artifact.**
   The audited 50-paper database holds 81 pure-spelling duplicate pairs, and
   database linkage contributed the largest share because each ingest carried
   its own raw-string identity helper. 51 further pairs are genuine source
   contradictions that need human adjudication, not merging.

### Fresh 50-paper BMPR2 validation run (2026-09-02)

A gold-free 50-paper BMPR2 run was completed on the repaired protocol
(`results/bmpr2_shared_fixes_20260902/BMPR2/20260902_074142`, fingerprint
`f6787618`, 85.8 min, 264 LLM calls, trace integrity `write_time_verified`).
**No precision, recall, or MAE is defined for it** — BMPR2 has no gold standard,
and nothing in this work was tuned to an answer key. What the run establishes is
internal consistency, measured against the prior run of the same 50-paper
roster:

| Property | Prior run (2026-08-26) | Fresh run (2026-09-02) |
| --- | ---: | ---: |
| Source-ledger checks passed | n/a (artifact self-contradictory) | 10 / 10 |
| Papers full text / abstract-only / unverified | 0 / 50 / — | 50 / 0 / 0 |
| Source-class discrepancies | not reported | 0 |
| Single-carrier papers | 50 (list of 30) | 25 (list of 25) |
| Variants | 522 | 436 |
| Foldable spelling duplicates | 87 | 5 |
| Identity conflicts held for adjudication | 51 | 35 |
| Sparse partials deliberately not completed | 167 | 91 |
| Chimeric identities | not measured | 0 |
| Cross-gene rows in `variants` | not measured | 0 |
| Archives on disk expanded | nested never expanded | 31 / 31 |

The five residual foldable pairs are one characterized class: a linkage row
carrying a genomic coordinate beside a paper row without one, with byte-identical
cDNA and protein. That was an over-applied sparse rule in the fold predicate,
since a derived coordinate is an annotation rather than an identity axis. It is
fixed and tested, so current code no longer creates them, but the fix postdates
the run artifact. **Those five are not curation items.**

One stage failure occurred: the paywall fetch exited 1 on the DOI
`10.3390/ijms25052734Submission`, where a rendered landing page printed the
identifier with the next word glued on and four separate DOI cleaners each
lacked a word boundary. No data was affected — all 50 papers already had full
text, the fetch queue held one row, and the override file was empty. Fixed with
one canonical DOI path.

One paper's nested supplement archive had never been read by any route. PMC
lists it as that paper's only supplement (392.7 KB, byte-exact against the
on-disk archive), so its entire declared supplementary material was previously
invisible. Its four extracted variants match the published abstract, with no
rows from the second gene's table inside the same file.

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

## 2026-08-24 gold-118 locks and source-backed validation

The accepted blind result is the immutable
`benchmarks/codex_paper_eval/runs/20260824_postfix_gold118` lock: 546 TP,
284 FP, and 86 FN over 632 gold rows (86.39% recall, 65.78% raw precision,
74.69% F1). Counted-extra precision is 97.50% and count-bearing-only precision
is 93.69%. Conditional count results are 68/206 carrier absolute-error units
(MAE 0.330), 27/49 affected (0.551), and 9/18 unaffected (0.500). This is the
headline current measurement; its trace-derived public-list-price proxy is
$11.32 for 574 calls and 2.760M tokens.

Two general repairs were initially staged and were **not a blind headline**.
First,
the trusted DB projection retains source-grounded structural-only identities,
and migration preserves exact structural cDNA breakpoints rather than
collapsing them into a broad exon label. Second, the deterministic table parser
recognizes a parent count header split into opposing phenotype subcolumns and
sets carrier total to their sum while preserving the two supplied partitions.
The latter was replayed from source for SCN5A PMID 29709101 in a staged DB and
recovered all 12 table totals exactly; no accepted extraction or score lock was
mutated.

Applied as a source-backed candidate projection before the fresh locks,
the repairs give 548 TP / 284 FP / 84 FN (86.71% recall, 65.87% raw precision,
74.86% F1), while carrier MAE falls to 47/208 = 0.226, affected MAE to 27/59 =
0.458, and unaffected MAE to 9/28 = 0.321. Counted-extra precision remains
97.51% and count-bearing-only precision is 93.75%. These numbers are diagnostic:
promotion requires a fresh gold-free extraction, hash lock, and post-lock score.
The exact unblinded overlay and machine-scored report are stored under the
lock's `diagnostics/` directory and are bound to the baseline predictions,
selection, overlay, and four gold-file SHA-256 digests. The exact evidence,
reproduction command, and budget accounting are in
`docs/evidence/gold118_grouped_header_candidate_20260824.md`.

Two fresh runs tested those repairs. The immutable grouped/structural lock is
542 TP / 278 FP / 90 FN (85.76% recall, 66.10% raw precision), with 98.19%
counted-extra precision, 95.41% count-bearing-only precision, and carrier MAE
0.222 over 207 supplied rows. It failed recall/identity non-regression.

The subsequent source-recovery lock is 546 TP / 285 FP / 86 FN (86.39% recall,
65.70% raw precision), with 97.50% counted-extra precision, 93.40%
count-bearing-only precision, and carrier MAE 0.198 over 197 rows. It failed the
accepted count-bearing precision and carrier-coverage floors, so it also was
not promoted. Selecting a different stochastic run for each metric is
prohibited.

The later source-bound AHA lock supersedes the postfix headline after a fresh
gold-free extraction and post-prediction lock: 554 TP / 283 FP / 78 FN (87.66%
recall, 66.19% raw precision), 97.88% counted-extra precision, 94.59%
count-bearing-only precision, and carrier MAE 0.193 over 207 rows. Affected and
unaffected coverage also rise to 56 and 28 while both MAEs fall to 0.321. It
passes every preregistered floor and is the authoritative result at the top of
this file.

KCNH2 PMID 29650123 alone contributes 20 remaining FNs. Live publisher
inspection found only the already-acquired `mmc1.docx`, containing cohort and
beta-blocker tables plus outcome figures but no mutation roster. This is a
source/provenance blocker, not permission to infer twenty singleton variants.

The exact 50 BRCA1 + 50 BRCA2 + 50 BMPR2 candidate remains a structural review
set, not a metric cohort. Only three of its 150 papers overlap an approved
curated fixture and BMPR2 has none. Precision, recall, and MAE for that set are
therefore undefined until the preregistered 90-paper calibration and 60-paper
holdout receive exhaustive human answer keys; they must never be reported as
zero or inferred from the structural audit.

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

## Historical all-paper baseline (pre-lock; not current headline)

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
