# Cardiac gold-standard curation queue — 2026-08-14

Label errors and completeness debt in the cardiac gold standard, found while
taxonomizing every false negative (101) and false positive (904) of the locked
`benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix` run. Every item
below was adversarially re-verified against PubMed metadata or the locked
on-disk source before listing. **These are curator decisions, not code.** The
source-verified subset listed under **Applied 2026-08-15** is now in the live
gold snapshot; the rest is still unapplied. Gate 2 counted extras (10 rows)
and the precision/cost ranking are in §8 and
[`PRECISION_COST_LEVERS.md`](PRECISION_COST_LEVERS.md).

Scale: the items in §1–§3 account for ~15 of the run's 101 false negatives and
2 scored false positives. Correcting them moves measured variant recall by
roughly +1.2 points with zero pipeline change.

## Applied 2026-08-15 (source-verified subset only)

Re-verified against PubMed, gold raw rows, and on-disk source, then applied to
the live snapshot (`gene_variant_fetcher_gold_standard/` plus the curated-eval
RYR2 recall copy). Locked benchmark runs were not rewritten.

| Item | Action |
|---|---|
| KCNH2 `c.2592+1G->A` PMID 10086972 | remapped to **10086971** (Berthet Circulation 1999; 10086972 is a blood-pressure paper) |
| RYR2 19398665 N4101I | renamed **N4104I** (Figure 2) |
| RYR2 30403697 G4722S | renamed **G4772S** (source text; sibling-row notes updated) |
| RYR2 32866913 R85I | renamed **T85I** (title/body `p.Thr85Ile`); counts left at 2/1/1 |
| SCN5A 29709101 `c.4813+3_4813_6dup` | renamed **`c.4813+3_4813+6dup`** |
| RYR2 16517285 R169Q ×2 | dropped the second identical `1/1/0` row |
| RYR2 22677073 V2113M ×2 | dropped the second identical `1/1/0` row |

Not applied to the underlying gold snapshot: 19216760 EXON 3 DELETION ×2 (two
families, 4/4/0 and 2/2/0); 14642689; editorial PMIDs; compound row; any count
edits. KCNH2 14642689 was quarantined from the live Tier 2 manifest on
2026-08-21 after a second blind run correctly returned no KCNH2 variant; its
historical gold row remains unchanged until the intended source is identified.

`tier2_gold_120.tsv` no longer includes KCNH2 10086972 (10086971 remains) or
KCNH2 14642689. Live tier membership is 118 attempts / 114 unique PMIDs. Locked
gold-120 runs still contain the historical attempts.

## 1. PMID errors (2 rows)

| gene | gold PMID | variant | evidence | proposed fix |
|---|---|---|---|---|
| KCNH2 | 10086972 | c.2592+1G->A (5 carriers / 4 affected / 0) | PubMed: 10086972 is "Relation of weight... Minneapolis Children's Blood Pressure Study" (Circulation 1999). The row's variant and counts exactly match the run's extraction from PMID **10086971** (Berthet et al., same issue, previous pages) — an off-by-one typo. The pipeline extracted c.2592+1G>A from 10086971 and was charged both an FP (10086971) and an FN (10086972). | change PMID to 10086971 |
| KCNH2 | 14642689 | A561T | PubMed: 14642689 is "Expression of angiotensin II receptors... atrial fibrillation" (JACC 2003) — not a genetics paper. Correct target PMID unknown; A561T is a common KCNH2 LQT2 variant with many candidate sources. | quarantined from live Tier 2; curator to locate intended paper before any restoration |

## 2. Editorial mis-attributions (11 rows)

| gene | gold PMID | rows | evidence |
|---|---|---|---|
| SCN5A | 15898185 | 9 rows (H558R, R34C, S524Y, ...) | PMID is R. Brugada's one-page editorial "Genetics, ethics and ethnicity" (Heart Rhythm 2004;1:608-9, Editorial/Comment per PubMed). Only S1103Y and R1193Q appear in its prose (both already score TP). The other 9 rows belong to the Ackerman ethnicity paper it comments on. |
| SCN5A | 22966897 | 2 rows (P1177L, R190W) | Mazzanti/Priori editorial ("Molecular autopsy... Time to discuss pros and cons"); neither token appears in any on-disk file for this PMID. |

## 3. Variant transcription errors (4 rows, source-verified)

| gene | PMID | gold has | source has | evidence |
|---|---|---|---|---|
| RYR2 | 19398665 | N4101I | **N4104I** | the paper's mutation-map figure plainly prints N4104I; extraction produced N4104I and scored FP+FN |
| RYR2 | 30403697 | G4722S | **G4772S** | source text has G4772S 9x, G4722S 0x |
| RYR2 | 32866913 | R85I | **T85I** | the paper's title is "Rare RYR2 p.Thr85Ile variant..." (R85I appears nowhere) |
| SCN5A | 29709101 | c.4813+3_4813_6dup | **c.4813+3_4813+6dup** | gold dropped the second '+'; paper text carries the intronic form |

## 4. Duplicate gold rows (scorer consumes gold 1:1, so each duplicate is a structural FN)

- RYR2 16517285: `R169Q` appears twice (RYR2_recall_input.csv lines 51 and 474).
- RYR2 22677073: `V2113M` appears twice (lines 62 and 73).
- RYR2 19216760: `EXON 3 DELETION` appears twice. Note the extractor DOES find
  this variant as `p.Asn57_Gly91del`; crediting it needs the exon-coordinate
  matcher extension (tracked in TASKS §5), but the duplicate row is gold-side.

## 5. Compound rows

- KCNQ1 33082985 carries `W248F`, `L347R`, **and** a compound `W248F + L347R`
  row for the same double-heterozygous patient. Under 1:1 consumption the
  compound row can never match once the components have.

## 6. Completeness debt (informational — these depress measured precision, not recall)

- KCNH2 26746457: eTable 5 lists 43 KCNH2 cohort variants; gold registers 10.
  24 correct extractions score FP. The eTable's cohort-MAF column implies ~1
  carrier for most rows if the curator wants counts.
- 3 additional count-bearing, source-quoted, correctly-attributed extractions
  score FP because gold does not register them (GOLD_GAP_REAL_VARIANT class in
  the 2026-08-14 FP taxonomy).
- KCNQ1 26496715 is a multi-gene supplement: gold correctly registers 54 KCNH2
  and 4 SCN5A rows under this PMID — no change needed, listed here because any
  per-gene gold export must keep those rows out of the KCNQ1 arm.

## 7. Source-availability flags

- KCNQ1 31293497 `A590T`: after the circuit-breaker fix the paper extracts (5
  rows), but A590T appears nowhere in the acquirable source text. Verify the
  gold row's provenance (supplement-only? different notation?).
- RYR2 28798025 `G1885E`: present in source (2 mentions) but missed by the
  extractor even after the fix — genuine residual extraction FN.

## 8. Counted extras on the 2026-08-15 diagnostic (Gate 2 FP surface)

The diagnostic rescore left **10 counted extras** (extras that assert
carriers/affected/unaffected). That is the entire remaining Gate 2
false-positive budget. Paper-level write-up:
`benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/PRECISION_AND_COST.md`.
Do not treat the 780 identity-only extras as this list.

| gene | PMID | extra | asserted | proposed reading |
|---|---|---|---|---|
| KCNH2 | 18808722 | L187P | carriers=24 | Other-gene: KCNQ1 gold on this PMID is L187P. Not a KCNH2 gold miss. |
| SCN5A | 23538678 | c.693delCA | carriers=2 | Near-twin of gold/TP `c.692_693delCA`. Scorer bridge, not a new gold row. |
| KCNH2 | 21779290 | p.Leu564Leu | 1/1/0 | Synonymous; gold on this paper is only K897T. |
| KCNH2 | 21779290 | p.Tyr652Tyr | 1/1/0 | Synonymous. |
| KCNH2 | 20181576 | P926A | affected=2 | Gold/TP is `P926fsX`. Distinct allele unless source says otherwise. |
| KCNH2 | 10086971 | c.3108insG | 4/1/3 | Not in the four-row gold. |
| RYR2 | 22677073 | c.2398+5G>T | carriers=1 | Not in the 25-row gold (leftover FN is Q2958R). |
| SCN5A | 17675083 | p.Arg34Cys | carriers=1 | Gold is G400A + H558R. Common polymorphism shape. |
| SCN5A | 19406494 | p.Ala1949Pro | 3/1/2 | Gold is only R808C. |
| SCN5A | 29709101 | I1660V | carriers=2 | All 11 gold rows already TP. |

Three of these ten are the §6 `GOLD_GAP_REAL_VARIANT` identities (source-quoted,
correctly attributed, absent from gold). Identify which three against source
before adding gold rows. Ranking: [`PRECISION_COST_LEVERS.md`](PRECISION_COST_LEVERS.md).

## Provenance

FN/FP taxonomies and adversarial verifications: 2026-08-14 nine-agent analysis
of the locked verticalfix run (session artifacts; per-FN classes and per-paper
dossier retained in the session scratchpad). PubMed identities verified against
live metadata on 2026-08-14.
