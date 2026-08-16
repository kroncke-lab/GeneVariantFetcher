# Diagnostic rescore — current gold + current matcher (2026-08-15)

Free re-score of the locked `20260813_gold120_verticalfix` predictions against
the live gold snapshot and the current `run_eval.py` matcher (including
`utils/structural_alleles.py`). No new extraction, no LLM calls. Locked
`report.json` / `evidence.csv` / `paper_metrics.csv` in the parent run dir were
not rewritten; this directory holds a copy of the locked `selection.json`,
`predictions.json`, and `LOCK.json` plus a fresh score.

Live `gold_120` membership is 119 attempts after dropping KCNH2 10086972. That
paper now has 0 gold rows and 0 predictions, so the 120-paper and 119-paper
totals are identical on this rescore.

## Headline (not a replacement for the locked Gate 2 table)

| Projection | Variant recall | Raw precision | Counted-extra precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: |
| Locked 2026-08-13 (then-gold, then-matcher) | 534/635 (84.09%) | 534/1438 (37.13%) | 534/(534+24) (95.70%) | 146/635 (22.99%) | 0.308 |
| This diagnostic | 545/633 (86.10%) | 545/1335 (40.82%) | 545/(545+10) (98.20%) | 160/633 (25.28%) | 0.2875 |

Net: **+11 TP**, **−2 gold rows** (the two source-verified duplicate drops),
**−13 FN**, **−114 FP**, counted extras **24→10**. Twin-merge collapsed 103
notation-twin prediction rows.

## What moved (21 papers)

Gold edits: remapped `c.2592+1G->A` onto 10086971; renamed N4104I / G4772S /
T85I / `c.4813+3_4813+6dup`; dropped duplicate R169Q and V2113M.

Matcher: KCNQ1 25471708 recovered M159X / Q530X / R192CFS91X / R518X (splice and
stop bridges); KCNQ1 26496715 recovered A344A / R562S and absorbed 59 cDNA
twins; several other papers absorbed cDNA/protein twins.

Structural alleles did **not** add TPs on this prediction set. The matcher
equates `EXON 3 DELETION` ↔ `p.Asn57_Gly91del` and `P.K1505_Q1507DEL` ↔ `ΔKPQ`,
but those identities were not present in the locked predictions for the FN
papers (19216760 extracted only N3308S; 28087622 never emitted ΔKPQ). RYR2
19926015 did extract `c.169-198_273+820del`, which is the genomic span of
coding exon 3 (cds 169–273); that bridge is still missing.

## Remaining 88 FNs (22 papers)

Classified against `docs/GOLD_CURATION_QUEUE_2026-08-14.md` and the locked
source, not re-extracted.

| Class | FN | Papers | Can a re-extract close it? |
| --- | ---: | --- | --- |
| Gold-side (wrong/editorial PMID, compound row) | 13 | 14642689, 15898185 (9), 22966897 (2), 33082985 compound | No. Drop from the answer key. |
| Source not readable as text (TIFF tables) | 20 | KCNH2 29650123 | No, not from the same XML/markdown. Needs figure OCR. |
| Caption-stub / thin table body (TASKS acquisition list) | 17 | 31520628 (7), 14678125 (6), 24667783 (2), 19926015 (2) | Only after reconverting the real table bodies. |
| Extraction-code already landed, not in these predictions | 11 | KCNQ1 21956039 (11/11 on the 2026-08-14 4-paper check) | Yes — this is the paid gold-120 arm. |
| Residual extraction / zero-output / hard notation | 27 | 28798025 (6, 0 pred), 27114410 RYR2 (5), 19216760 (4), 15466642 (3), 26496715 (2 DUP/fs), 31293497, 17971661, 21288276, 28087622 ΔKPQ, 22677073 Q2958R, 27114410 KCNQ1, 10086971 C1117fsX | Maybe some. Not promised by current code. |

The two 19216760 `EXON 3 DELETION` rows are real families (4/4/0 and 2/2/0), not
duplicates. They stay FN until extraction emits the exon event.

## How much closer can extraction get?

These are ceilings on the **current gold denominator (633)**, not new
measurements.

- **Already realized, $0:** 86.10%.
- **If the remaining 13 gold-side rows are dropped:** 545/620 = **87.90%**, still
  no new extraction.
- **If the paid gold-120 arm repeats the 21956039 recovery:** +11 → 556/633 =
  **87.84%** (or 556/620 = **89.68%** after gold-side drops). Other
  circuit-breaker / gene-less-table papers may add a few more; the earlier
  4-paper check recovered 19/22 previously-missed rows, 11 of them here.
- **If 29650123 TIFF tables become readable:** +20 more, to ~91% before residual
  hard cases.
- **Realistic residual after source + gold cleanup:** the 27 hard FNs above.
  That is about **95.7%** (606/633) as an optimistic identity ceiling, not a
  forecast.

This diagnostic cannot measure the large-table short-circuit refuse
(26496715-class). That only appears on a new extraction.

Do not treat these figures as a new Gate 2 lock or a four-gene headline.
