# Mixed-gold protocol campaign, 2026-09-03: what moved, what could not, and why

Brett's goal for the day: run three or four mixed-gold tranches to raise
paper-derived variant recall and improve count MAE using only the papers,
quantify how much gold is missing because the data are absent from the
acquired text/supplements/figures, root-cause the largest failures per
tranche, and take advice from Grok 4.6 `xhigh` and the agy CLI throughout.
This document is the standalone summary; the per-tranche evidence documents
hold the tables.

## 1. The acquisition ceiling (answer to "how many are we missing?")

`scripts/recall_audit/gold_source_presence_sweep.py` classifies every gold
row against every byte on disk for its paper (article body, converted binary
supplements, article PDFs), blind to predictions
(`gold_source_presence_sweep_20260903.md`).

| stratum | gold rows | hard ceiling | wide ceiling | present in body |
|---|---:|---:|---:|---:|
| all 1,533 attempts | 7,316 | 1,321 (18.1%) | 2,236 (30.6%) | 5,071 (69.3%) |
| 1,422 runnable attempts | 7,107 | 1,122 (15.8%) | 2,037 (28.7%) | 5,061 (71.2%) |

Hard ceiling = no source, abstract-length stub, glyph-garbled PDF text, or a
substitution absent from every searchable byte. Wide ceiling adds figure-only
and non-searchable notation (unknown, not proven absent). Only 9 rows live in
a binary supplement the body lacks: **supplements are not where the missing
gold is**; the rows are table-body capture failures on identifiable papers
(KCNQ1 14678125 37/41, 17192539 51/57; SCN5A 21273195, 24631775, 11901046),
which are an acquisition worklist, not extraction work.

## 2. Tranche results

| tranche | attempts / gold rows | baseline (frozen `506a949c`) | candidate | identity rule | count effect |
|---|---:|---|---|---|---|
| 01 | 49 / 242 | 155/61/87, R 64.0%, P 71.8% | v1 `b56f469f`: 157/54/85, R 64.9% (+0.83 pp), P 74.4% | reject (bar +1.0 pp) | carriers supplied 48→125, cond. MAE 0.81→0.10, e2e MAE 2.68→1.95, counted-extra P 82→91% |
| 02 | 49 / 303 | 268/54/35, R 88.4%, P 83.2% | v2 `32ec857c`: 267/55/36, R 88.1% (−0.33 pp) | reject | flat (carriers 232→227); all differences model variance |
| cont120-01 | 120 / 384 | pending | v2 (identical runtime) | pending | pending |

Root cause of the misses (`scripts/recall_audit/fn_root_cause.py`, new): on
tranche 01, 70 of 87 misses were acquisition, 12 notation the probe cannot
search, and 5 reachable; on tranche 02, 13 acquisition, 20–21 notation, 2
reachable. The +1.0 pp discovery bar is 2–3 rows on these tranches, i.e.
essentially every reachable row.

## 3. What was changed in the protocol (all gene-agnostic, all $0 per paper)

Traced from the five reachable tranche 01 misses and its count/extra audits
(commit `b56f469f`), hardened after the tranche 01 result and two reviewer
counterexamples (commit `32ec857c`):

1. Publisher access shells whose folded supplements carry the tables are
   usable/reusable corpus sources (22 of 70 refused corpus full texts).
2. Corpus reuse prefers a much richer `CLEANED` sibling; preprocessing keeps
   a staged richer rendering.
3. Source-only legacy protein spellings (`L860fsx89`, `1795insD`) become
   identities; the projection falls back to `legacy_notation`.
4. `HGVSp`/`HGVSc` headers and `p.(Glu101Lys)` values are recognised.
5. A closed case-series count column sets carriers = counted cases, typed
   `case` so the phenotype guard keeps the affected value; a bare case-control
   `Cases` column never fires.
6. ClinVar/PubTator are a row's origin only via the ingest markers; notes
   mentioning ClinVar no longer move a paper row into the linkage lane.
7. Scorer: same-codon bridges do not fire between conflicting concrete cDNA
   alleles.
8. Trust gate: genotype-frequency labels are population counts.

Offline suite: 2,791 passed after v2.

## 4. What the campaign learned about the measurement itself

- **Identity recall on the gold is acquisition-capped, not reading-capped.**
  80% of tranche 01's misses and 37% of tranche 02's are absent from every
  acquired byte; most of the rest are figure-only or unsearchable notation.
- **Single-extraction paired arms carry provider nondeterminism of the same
  size as the effects being sought.** On tranche 02 every arm difference was
  the model emitting or omitting a protein notation or a carrier value at
  temperature 0. A ±1–2 TP / ±5 count swing per 49 attempts is noise.
- **Counts are where the protocol moved.** Tranche 01's candidate raised
  carrier supply by 77 rows and cut the end-to-end carrier error by 27%
  while the identity score moved two rows. A secondary count endpoint was
  preregistered before the tranche 02 candidate lock
  (`mixed_gold_count_endpoint_preregistration_20260903.md`) with identity
  non-inferiority as a hard guard; tranche 01 would have passed it (locked
  diagnostic only), tranche 02 does not.
- **Grok's cancelled-run trap** and the concurrent-session fingerprint hazard
  are recorded in the reviewer-CLI memory notes.

## 5. Reviewer input (Grok 4.6 `xhigh`, agy = Gemini 3.1 Pro `high`)

Both reviewers endorsed the arm design (frozen baseline vs cumulative
candidate; confirmation only for the identical runtime on the next unopened
tranche) and both ranked the case-series rule and the fold-marker gate as the
top regression risks of v1; their concrete counterexamples became v2 tests.
Both judged the secondary count endpoint legitimate if registered before the
next score, motivated by the blinded ceiling, guarded by identity
non-inferiority, and never applied retroactively. Grok located the registry
digest binding that forces the rule into a sibling file, and asked for the
observed-recall ≥ 0 guard that was adopted.

## 6. Cost

Public-price proxy: tranche 01 $8.575 (both arms), tranche 02 $8.961; ledger
$62.055 used / $37.945 remaining before the continuation tranche, whose paired
estimate is $25.11.

## 7. Continuation tranche 01

Pending; see `mixed_gold_cont120_01_20260903.md`.
