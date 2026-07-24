# Validation notes — Codex high-carrier 48-paper run

## Outcome

The primary, post-lock notation-audited score is **702 TP / 52 FP / 299 FN**:
**93.1% precision, 70.1% recall, and 80.0% F1** over 1,001 gold variant
rows. The locked extraction contains 754 predictions, each with non-empty
evidence and a source location.

This is an extraction-blinded evaluation, not a selection-blinded one. The user
supplied the same 48 gold-ranked high-carrier PMIDs used by Claude. Gold values
and gold row counts were absent from routing and extraction. The prediction
file was made read-only and SHA-256 locked before scoring.

The supplied set contains **1,001 gold variant rows** and the gold carrier
column sums to **5,949**. Thus “~1,000” describes variant/count assertions, not
the summed number of carriers.

## Count coverage and conditional accuracy

| field | supplied / 1,001 | count recall | exact among supplied | MAE | RMSE |
|---|---:|---:|---:|---:|---:|
| carriers | 597 | 59.6% | 93.3% | 0.508 | 5.027 |
| affected | 339 | 33.9% | 95.0% | 0.375 | 2.894 |
| unaffected | 215 | 21.5% | 87.0% | 0.549 | 2.769 |

MAE and RMSE are conditional on a matched variant for which both gold and
Codex supplied the field. Count recall must therefore be read alongside the
small error magnitudes; abstentions and missed variants do not enter MAE/RMSE.

## Scorer sensitivity

The preserved pre-adjudication matcher scored **611 TP / 143 FP / 390 FN**
(81.0% precision, 61.0% recall, 69.6% F1). Review after the prediction lock
identified 91 notation-equivalent pairs, primarily:

- verbose parenthetical three-letter protein labels versus compact one-letter
  labels;
- cDNA frameshift/deletion rows versus the corresponding protein event;
- legacy `splice`, `sp`, `stop`, deletion, and insertion shorthand.

The audited matcher raised the semantic score to 702 TP without editing any
variant, count, evidence, or source location. Every adjudicated pair is listed
in `matcher_adjudication.csv`; both raw matcher reports remain preserved.

## Comparison with Claude's completed runs

| system | precision | recall | F1 | extracted | total tokens | tokens/paper | wall minutes |
|---|---:|---:|---:|---:|---:|---:|---:|
| Codex, audited | 93.1% | 70.1% | 80.0% | 754 | 856,768 | 17,849 | 58.2 |
| Codex, pre-adjudication | 81.0% | 61.0% | 69.6% | 754 | 856,768 | 17,849 | 58.2 |
| Claude Fable 5 | 81.0% | 74.0% | 77.3% | 915 | 3,129,963 | 65,208 | 14 |
| Claude Sonnet | 82.0% | 72.8% | 77.1% | 889 | 3,590,702 | 74,806 | 28 |
| Claude Haiku | 78.5% | 64.9% | 71.1% | 828 | 2,612,908 | 54,436 | 11 |

Codex audited recall is 3.9 percentage points below Fable, 2.7 below Sonnet,
and 5.2 above Haiku. Codex used 72.6% fewer tokens than Fable, 76.1% fewer than
Sonnet, and 67.2% fewer than Haiku. The wall-clock comparison is not
like-for-like: this Codex harness processed papers serially, while Claude used
parallel subagents. Token definitions also differ (Codex exact API telemetry
versus Claude subagent accounting), and Claude used its production matcher
rather than the Codex notation-adjudicated matcher.

By gene, Codex recall was SCN5A 83.6%, KCNH2 69.2%, KCNQ1 43.2%, and RYR2
91.3%. Relative to Claude, Codex's distinctive result is much stronger RYR2
recall (Claude: 75.7% Fable, 64.7% Sonnet, 74.3% Haiku) but weaker KCNQ1 recall
(Claude: 60.1%, 59.8%, 48.0%).

## Error concentration

Seven papers account for **263 of 299 false negatives (88.0%)**:

| gene | PMID | tool | FN | dominant issue |
|---|---:|---|---:|---|
| KCNQ1 | 17192539 | text | 56 | narrative names only one mutation; remaining mutation identities are absent |
| KCNQ1 | 17470695 | table | 56 | downloaded structured preview lacks the variant-level rows implied by the paper |
| KCNQ1 | 14678125 | text | 41 | only aggregate mutation-region counts are preserved |
| SCN5A | 26746457 | text | 30 | route chose text; the known Figure 1 carrier matrix required figure/OCR fallback |
| KCNH2 | 14661677 | OCR | 29 | partial article and unusable publisher artifacts |
| KCNH2 | 19038855 | text | 28 | variant identities reside in a table image not transcribed into the text |
| KCNH2 | 24667783 | table | 23 | selected table representation exposes aggregate phenotype data but no usable identifiers |

The first three KCNQ1 papers alone contribute 153 false negatives, 51.2% of
all misses. Nine papers contribute 276 misses, 92.3% of the total. This is a
whole-paper/source-routing tail, not a uniform per-variant degradation.

The largest count disagreements also expose scope ambiguity:

- KCNH2 PMID 19160088: Codex used 16 R176W carriers in the present population
  cohort; gold uses 112 LQTS patients mentioned as prior/unpublished clinical
  experience.
- SCN5A PMID 20470418: Codex used the 26 genotyped infant hearts (17 SIDS,
  9 other-cause deaths); gold uses 85 carriers and 39 affected.
- SCN5A PMID 18451998: Codex used the pedigree/table split of 41 carriers,
  9 symptomatic and 32 asymptomatic; gold uses 50/47/3.
- KCNQ1 PMID 33141630: Codex used all 124 identified T224M carriers and 34
  affected among 88 deeply phenotyped carriers; gold carrier count is 88.
- RYR2 PMID 25814417: Codex combined 179 living mutation-positive relatives
  with six mutation-positive SCD cases (185 total); gold carrier count is 179.

These remain scored as errors, but they should be adjudicated as
population/denominator-definition disagreements before changing extraction
logic.

## Validation checks

- Selection, prediction, and report each contain 48 unique papers, 12 per gene.
- The manifest exactly matches the user-supplied 48 PMID/gene pairs.
- `selection.json` exposes only source metadata, representations, and hashes;
  it contains no gold values or row counts.
- All 48 selected source files still match their recorded SHA-256 hashes and
  byte sizes.
- This historical run predates all-representation hashing: its selection records
  hashes for the primary markdown sources, not artifact JSON, PDFs, or figures.
  New runs hash and verify every selected representation before extraction and
  again before locking.
- `predictions.json` and `selection.json` match the hashes in `LOCK.json`.
- Per-paper tokens sum exactly to 856,768; per-paper elapsed time sums to
  3,489.436 seconds; wall time is 3,490.596 seconds.
- TP/FP/FN sums, micro precision/recall/F1, per-gene totals, count recall,
  MAE, and RMSE were independently recomputed from the detailed output.
- `evidence.csv` has exactly 754 data rows and no blank evidence or
  source-location fields.
- The report's DuckDB queries executed successfully and returned the expected
  1/4/12/3/48/5 row shapes for overview, gene, chart, count, paper, and model
  datasets.
- The bounded seven-dataset MCP report payload passed artifact validation and
  rendered successfully.

## Recommended next change

Add a data-rich-paper fallback: when the selected representation yields zero
or implausibly few variants, retry table, figure/OCR, and PDF representations
before accepting the paper result. The failure concentration suggests this
will improve recall more than spending additional tokens on already successful
variant rows.
