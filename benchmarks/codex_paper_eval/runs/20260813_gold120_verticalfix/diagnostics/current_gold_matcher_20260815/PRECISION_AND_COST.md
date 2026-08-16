# Gold-120: how to raise precision and cut cost

Analysis of the 2026-08-15 diagnostic rescore of the locked
`20260813_gold120_verticalfix` predictions. No new extraction. Locked parent
`report.json` is not a new Gate 2 lock and was not rewritten.

## 1. Two precision numbers — say which one you are moving

Live diagnostic (`report.json` `overall`):

| Denominator | Formula | Value |
| --- | --- | ---: |
| Raw identity precision | 545 / (545+790) = 545/1335 | **40.82%** |
| Counted-extra precision (Gate 2) | 545 / (545+10) | **98.20%** |

Gate 2 already uses **counted-extra** precision
(`precision_vs_counted_gold_pmids`): every matched gold row is signal; only
extra predictions that assert a patient count enter the false-positive
denominator. The locked 2026-08-13 Gate 2 pass was 534/(534+24) = 95.70%
against a 77.3% floor. This diagnostic is 98.20% on the same definition.

Those two numbers can move in opposite directions. Collapsing them into one
"better precision" target is not a plan:

- Dropping the 780 identity-only extras would raise raw precision and leave
  Gate 2 unchanged — and most of those extras are gold-incomplete, not
  invented.
- Adding a counted extra that is a real patient count would *lower* Gate 2
  even if the identity is correct.

**Any precision action below is labeled raw or counted-extra.**

## 2. Remaining extras: 10 counted vs 790 raw

`overall.fp` = 790. `counted_precision.counted_extra_rows` = 10.
Identity-only extras = 780. Twin-merge already collapsed 103 notation-twin
rows before scoring.

### 2a. The 10 counted extras (the only Gate 2 FP surface)

Every paper with `counted_precision.counted_extra_rows > 0`:

| Gene | PMID | Counted extra | Asserted counts | Class |
| --- | --- | --- | --- | --- |
| KCNH2 | 18808722 | L187P | carriers=24 | **Other-gene gold.** KCNQ1 gold on this same PMID is `L187P`. KCNH2 gold is the eight KCNH2 alleles (all TP). |
| SCN5A | 23538678 | c.693delCA | carriers=2 | **Near-twin of the TP.** Gold/match is `c.692_693delCA`. Matcher left a second counted row. |
| KCNH2 | 21779290 | p.Leu564Leu | 1/1/0 | Synonymous; gold on this paper is only K897T. |
| KCNH2 | 21779290 | p.Tyr652Tyr | 1/1/0 | Synonymous. |
| KCNH2 | 20181576 | P926A | affected=2 | Gold has `P926fsX` (matched as TP). Missense vs frameshift — not the same allele. |
| KCNH2 | 10086971 | c.3108insG | 4/1/3 | Not in the four-row gold. Sibling extras are other frameshifts (D1037fs, G965fs, L953fs) without counts. |
| RYR2 | 22677073 | c.2398+5G>T | carriers=1 | Not in the 25-row RYR2 gold (24 TP, leftover FN is Q2958R). |
| SCN5A | 17675083 | p.Arg34Cys | carriers=1 | Gold is G400A + H558R (both TP). R34C is a common SCN5A polymorphism shape. |
| SCN5A | 19406494 | p.Ala1949Pro c.5845G>C | 3/1/2 | Gold is only R808C (TP). Eight other extras on this paper are count-null. |
| SCN5A | 29709101 | I1660V | carriers=2 | Gold has 11 rows, all TP after the `c.4813+3_4813+6dup` rename. Sibling extra `p.Ile1659Val` is count-null. |

`docs/GOLD_CURATION_QUEUE_2026-08-14.md` §6 records **three**
`GOLD_GAP_REAL_VARIANT` counted extras: source-quoted, correctly attributed,
absent from gold. They sit in this list of 10. They are curator decisions, not
extractor defects.

If Gate 2 is the target, the whole remaining FP budget is these 10 rows.
Removing the one other-gene count (L187P) and merging `c.693delCA` would make
counted-extra **545/(545+8) = 98.55%** with no new extraction. Curating the
three GOLD_GAP rows as gold would move them from extra to TP and raise both
recall and counted-extra.

### 2b. The 790 raw FPs, exclusive class split

Each extra is in exactly one row. Other-gene matching used the live cardiac
gold on the same PMID (`matches` / `twin_identical`).

| Class | n | What it is |
| --- | ---: | --- |
| Other-gene gold spill | 69 | Extra matches another cardiac gene's gold on the same PMID: 52 KCNH2 + 8 KCNQ1 + 7 SCN5A + 2 RYR2. Bulk is KCNQ1 **26496715** (multi-gene supplement; gold correctly holds 54 KCNH2 + 4 SCN5A under this PMID). Includes counted L187P. |
| Named gold-completeness (KCNH2 **26746457** eTable) | 34 | Gold registers 10 / ~43 eTable-5 cohort variants; all 10 match. Queue: 24 of these extras are correct extractions. All 34 are count-null (vertical-table fix). **Raw-only.** |
| Synonymous / polymorphism-shaped identity | 10 | `A29A`-style or `p.Xxx=` extras that did not hit other-gene gold. |
| Unmerged notation twin | 1 | Twin of a TP that `merge_notation_twins` refused (count conflict or ambiguity). |
| Counted extras not in the rows above | 9 | The Gate 2 list minus L187P (already in other-gene). |
| **Unallocated remainder** | **667** | Stated remainder. Not inventoried allele-by-allele. |
| **Total** | **790** | |

The 667 remainder is not mystery invention. **584 / 790 extras sit on 69
papers that already matched every gold row** (`tp == gold_n`, `fp > 0`). Gold
registered a subset; the pipeline kept reading. Largest gold-exhausted
clusters:

| Gene | PMID | FP | TP = gold_n | Note |
| --- | --- | ---: | ---: | --- |
| KCNQ1 | 19632626 | 67 | 1 | Catalogue vs 1-row gold |
| SCN5A | 32533946 | 52 | 83 | Extra identities after a complete gold hit |
| RYR2 | 29350269 | 49 | 3 | 20 twins already merged; 49 remain |
| RYR2 | 33536282 | 44 | 1 | |
| KCNH2 | 26746457 | 34 | 10 | Counted in the named eTable class above |
| SCN5A | 29544605 | 28 | 6 | |

The other 206 extras sit on 14 mixed papers (`tp < gold_n`). The two large
ones are RYR2 **19926015** (78 FP / 2 FN, caption-stub) and KCNQ1 **26496715**
(72 FP / 2 FN, multi-gene; most of the 69 other-gene hits live here).

**Do not raise raw precision by deleting gold-exhausted extras.** That hides
true extractions (26746457 is the worked example). Raw precision moves by
*adding those identities to gold* (FP→TP) or by teaching the extractor not to
emit off-gene / polymorphism rows. Counted-extra barely moves either way
because these clusters are count-null.

## 3. Cost is three measured shares, not a new model list

Locked gold-120 public-price proxy from `docs/PROTOCOL_COST_EVAL.md`
(same table as `COST_AND_GATE.md`):

| Model role | Calls (successful) | Cost proxy | Share of $9.774 |
| --- | ---: | ---: | ---: |
| GPT-5.6 Sol verification/figure vision | 393 (388) | **$7.019** | 71.8% |
| Grok 4.3 primary extraction | 115 (114) | **$2.587** | 26.5% |
| Kimi K2.6 table routing | 19 (19) | **$0.168** | 1.7% |
| **Total** | **527 (521)** | **$9.774** | 100% |

Locked trace manifests split the 393 Sol calls without inventing a second
price table: **244 `figure_variant_reader` + 148 `VariantClaimVerifier` + 1
extraction-adjudication**. Figure vision is 62% of Sol *calls*; claim
verification is 38%. Token-weighted dollars inside Sol are not published
separately; levers below attach to those call shares of the **$7.019** Sol
line, not to invented sub-prices.

### Cheaper / faster levers, each tied to a share

1. **Sol figure vision (244 of 393 Sol calls → the $7.019 line).** Largest
   cost share. A skip/reuse rule (no figures on disk, text already holds the
   paper's variants, caption-only stubs) cuts this call count. Measuring
   which of the 244 calls added a unique identity is $0 on the locked traces.
   Implementing a skip is a later extract; it is not required to rank the
   lever.
2. **Sol claim verification (148 of 393 Sol calls → the $7.019 line).**
   Already risk-ranked inside `ExpertExtractor`. Two $0-to-specify, later-to-
   measure cuts: (a) run the same cards *after* merge so twins are not
   verified twice (TASKS already lists this move); (b) tighten the rank so
   identity-only or gold-exhausted catalogue rows do not get a Sol card.
   This is fewer Sol calls on the existing verifier, not a new model.
3. **Grok extraction (115 calls → $2.587).** Already ~one primary extract per
   attempt. A successful deterministic short-circuit *saves* a Grok call; the
   2026-08-15 short-circuit *refuse* (census saw count columns, parse emitted
   none) *adds* a Grok call on the 26496715 class. Do not advertise that
   refuse as a cost win. Skipping Grok on abstract-only / zero-source papers
   is the only other cut on this share, and it is small (115 is already the
   floor for "read every paper once").
4. **Kimi routing (19 calls → $0.168).** 1.7%. Not a project.

What does **not** belong on this list: shopping a cheaper extract model,
un-parking the per-paper final check, or turning `COUNT_RECOVERY_ENABLED` on.
Those are either new-model or already-rejected cost.

A full 119-paper re-extract is another ~$10 of the same mix. It is the wrong
tool for a precision ranking: Gate 2 is 98.20% on these locked predictions,
and the dominant raw-FP clusters are gold-complete-subset catalogues that a
re-extract will emit again.

## 4. Ordered next actions

### $0 — do these; they answer the precision/cost question from the lock

1. **Curate the 10 counted extras** (counted-extra / Gate 2). Decide other-gene
   L187P, matcher `c.693delCA`, the two synonymous counted rows, and the three
   GOLD_GAP identities. This is the only list that can still move 98.20%.
2. **Do not delete the 780 identity-only extras to "fix" 40.82%.** 584 of 790
   FPs are on papers that already hit every gold row. 26746457's 34 extras
   are the documented eTable debt.
3. **If raw precision is the actual target, expand gold** on the named
   completeness papers, starting with 26746457 (24 already source-confirmed).
   That is FP→TP on the raw denominator and does not touch Gate 2 (count-null).
4. **Keep / slightly extend the scorer** ($0): `c.693delCA` ↔ `c.692_693delCA`,
   and RYR2 `c.169-198_273+820del` ↔ EXON 3. Do not add per-variant aliases.
5. **Trace-only Sol audit** ($0): which of the 244 figure-vision calls added an
   identity the text extract already had. That sizes lever 3.1 before anyone
   spends another $7 on vision.

### Needs new extraction or acquisition — do not start the paid 119

6. **Paid gold-120 re-extract (~$10 / ~40 min)** answers *recall* (KCNQ1
   21956039 11/11 on the 4-paper check) and the short-circuit refuse. It is
   not required for any precision or cost claim in this note. Do not run it
   for this ranking.
7. Figure-OCR on 29650123 TIFFs and caption-stub reconvert are recall/source
   work, not precision work.
8. The regex_table column binder is a count-quality / MAE lever on 26496715,
   not a counted-extra precision lever (that paper currently has **0**
   counted extras).

## What "better and cheaper" means on this sample

- **Counted-extra / Gate 2:** already 98.20%. The remaining 10 rows are a
  curator + one matcher gap. Do not spend $10 to chase this number.
- **Raw identity 40.82%:** mostly gold-subset catalogues and one multi-gene
  supplement. The honest raise is gold completeness (26746457 first) and
  not emitting other-gene rows, not another extract.
- **Cost:** 72% is Sol, and 244/393 Sol calls are figure vision. The cheap
  next measurement is "which vision calls were redundant on the lock." Grok
  is already one call per paper. Kimi is noise.
