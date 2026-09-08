# Five-item pass: gold adjudication, title tier, paper tier, manual acquisition, tranche 05 — 2026-09-08

Brett's instruction (2026-09-08), after the review of the 2026-09-07 Astra
work: take the five next steps in order, give the recorded data and a proposed
answer for the two curator questions, use best judgment on the phenotype lever
but keep it general, treat non-cardiac tables as out of scope, replace the
EZproxy route with a ranked manual download list of one to two hundred papers,
and run the confirmation tranche. Starting point: `main` at `9954b23d`.

## 1. Gold adjudication: what is recorded, and a proposed answer

Two different papers, two different questions. The 46 wrong affected values
that keep the projection under its 90% gate all come from the first; the
definition of "affected" that any future phenotype reader must target comes
from the second. Nothing in gold was changed; the tables below are proposals.

### SCN5A 20129283 (Kapplinger 2010, Brugada compendium, nine testing centres)

Table 4 prints one row per nucleotide change with **No. of unrelated
individuals** and the list of **Testing center** numbers (293 rows; N sums to
451 = 438 genotype-positive patients + 13 double-mutation carriers). Gold
(`SCN5A_recall_input.csv`) holds 417 rows for 347 variants: 300 variants with
one row, 47 with two to five rows. The extra rows follow the centre list: R878H
(5 individuals, centres 1, 2, 4, 5, 7) is five gold rows of `1 / 1 / 0`; G1661R
is two rows because the paper prints two nucleotide changes. The paper never
reports per-centre counts, so when N exceeds the number of centres the split
was invented: D1243N is `1 / 1 / 0 | 1 / 1 / 0 | 5 / 5 / 0` against a printed
5, I1660V is `4 / 4 / 0 | 4 / 4 / 0` against a printed 5, G1408R is `1 | 1 | 2`
against 7. **14 of the 47 multi-row variants do not sum to Table 4** (three
more could not be joined on notation); K1493X is `2 / 2 / 0` where the table
prints 1 (the 2 belongs to K1493del, which gold also holds).

The pipeline emits one row per printed nucleotide change with the table's N.
The scorer's one-to-one assignment pairs it with one of the per-centre rows,
so every multi-row variant scores as wrong although the value is the paper's.

**Proposed answer.** Collapse 20129283 to one gold row per printed nucleotide
change: `carriers = affected = Table 4 N`, `unaffected = 0` (Brugada-ascertained
cases), keeping distinct nucleotide changes with the same protein effect (G1661R
G>A 2 and G>C 1; K1493del 2 and K1493X 1) as separate rows. The per-centre
split adds no information the paper contains and is internally inconsistent
for 14 variants. Proposal per variant, with the current rows, the Table 4 value
and the status: [`gold_20129283_multirow_proposal.csv`](item1_gold_adjudication/gold_20129283_multirow_proposal.csv)
(47 variants); the full join for all 347 gold variants is
[`gold_20129283_vs_table4.csv`](item1_gold_adjudication/gold_20129283_vs_table4.csv).
Under this collapse the projection's 46 wrongs become 45 exact plus one K1493X
correction, and its new-value precision on the opened locks is 100%.

### SCN5A 30059973 (Baruteau 2018, 442 children with SCN5A mutations)

Gold holds 185 rows, one per distinct variant in Table 5, each recorded as
`N / N / 0` with every carrier in `affected_lqt3` (E1784K `69 / 69 / 0`). The
paper's own tables say otherwise. Table 1 stratifies the 442 by baseline ECG
phenotype: **196 negative ECG phenotype**, 47 isolated LQT3, 8 isolated BrS,
113 isolated PCCD, 6 SSS, 3 DCM, 69 overlap; two-thirds were diagnosed through
family screening. Table 14 gives the split for the 12 most common variants
(150 of the 445 occurrences); for E1784K: 29 ECG-negative, 13 isolated LQT3, 17
isolated PCCD, 10 overlap.

| variant | total | gold C/A/U | ECG-negative | LQT3 | BrS | PCCD | SSS | DCM | overlap |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| E1784K | 69 | 69/69/0 | 29 | 13 | 0 | 17 | 0 | 0 | 10 |
| V411M | 10 | 10/10/0 | 2 | 5 | 0 | 2 | 0 | 0 | 1 |
| I1768V | 9 | 9/9/0 | 6 | 1 | 0 | 1 | 0 | 0 | 1 |
| Q1507_P1509del | 9 | 9/9/0 | 4 | 5 | 0 | 0 | 0 | 0 | 0 |
| G1743E | 8 | 8/8/0 | 6 | 0 | 0 | 1 | 0 | 0 | 1 |
| G1319V | 8 | 8/8/0 | 6 | 0 | 0 | 2 | 0 | 0 | 0 |

Full table with source line numbers and the proposed columns:
[`gold_30059973_table14_vs_gold.csv`](item1_gold_adjudication/gold_30059973_table14_vs_gold.csv).

**Proposed answer.** For these 12 variants record `carriers = row total`,
`unaffected = ECG-negative`, `affected_lqt3 = isolated LQT3`, `affected_brs1 =
isolated BrS`, `affected_other = PCCD + SSS + DCM`, and the overlap count as
ambiguous. For the other 173 variants the paper gives no per-variant phenotype:
keep `carriers = N` and leave affected/unaffected unknown rather than `N / 0`.
This matters beyond the benchmark: 29 of 69 E1784K carriers with a negative
ECG counted as affected inflates penetrance. It also settles what a phenotype
reader must produce for this table shape (variant x phenotype-category counts):
the negative-phenotype column is `unaffected`, disease columns are `affected`
by disease, overlap is ambiguous. The pipeline's current behaviour, 185
carrier rows and no affected value, is the correct abstention until then. Gold
already follows this convention where a paper prints it (RYR2 28237968 has
seven rows with unaffected relatives), so the change is consistency, not a new
rule.

## 2. Phenotype lever: where the gap is, what was built, what it measures

Decomposition of the 255 positive-gold affected values still null after the
projection on the two opened cardiac locks
([`gap_analysis_locks_before_title_tier.txt`](item2_measurement/gap_analysis_locks_before_title_tier.txt)):

| still null | rows | shape |
| --- | ---: | --- |
| model-authored rows without any carrier count | 95 | prose or model-read tables |
| model-authored rows with carriers, no affected | 51 | model rows; the lane refuses them by design |
| SCN5A 28341781 Table 2 (fixed-width, one proband per row) | 52 | count column refused as "not a people count" |
| KCNH2 27871843 and 11854117 rows without carrier counts | 22 | tables without a recognised count column |
| KCNH2 11854117 Table T1 "No. of Subjects" | 18 | caption names no cohort; paper tier would rescue |
| RYR2 28237968 per-person proband/relatives table | 13 | symptom column, relatives counts |

The one general shape worth a deterministic rule was the 28341781 block: a
**mutation catalogue with one proband per row in a paper whose title defines
the cohort** ("... Characteristics of Probands With Brugada Syndrome: A Japanese
Multicenter Registry"). Three things were wrong with it: the fixed-width
clinical mutation parser labelled the implicit one-proband count with the text
column "Coding Effect", so the classifier saw a non-people column; the second
page of the table carried the caption "Table 2. Continued", which names no
cohort; and no tier could use the title.

Built (all general, no paper-specific rule):

- **Title-ascertainment tier** (`title_ascertains`, `TIER_TITLE`, on by default
  via `GVF_TABLE_COHORT_TITLE_ASCERTAINMENT`): a title that binds a people noun
  to a disease ("probands with Brugada syndrome", "patients referred for long QT
  syndrome", "LQTS patients") certifies a count table whose caption names no
  other cohort (as the sentence tier already did) and, only from the title, a
  one-mutation-per-proband catalogue: implicit one carrier per row, a
  mutation/variant noun in the caption, no roster wording (clinical,
  characteristics, subjects, carriers, relatives, screened, follow-up ...), no
  phenotype/status/symptom/control column, exactly one carrier on the row. A
  caption that itself names the disease case series certifies such a catalogue
  without the title. Titles naming exome, genome, population, autopsy or
  relatives designs never qualify. The sentence-level tier still cannot
  certify per-person rows.
- **Continuation headings** ("Table 2. Continued", "(cont.)") resolve to the
  table's first printed caption and are ignored as competing captions.
- **Parser provenance**: the fixed-width clinical mutation parser labels an
  implicit one-proband row as such; archived rows with the old label are read
  the same way so replay and refresh classify them as a fresh extraction would.

Measurement, zero LLM calls ([`item2_measurement/`](item2_measurement/)):

| | positive-gold affected exact, cardiac four | new values |
| --- | ---: | --- |
| projection off | 203 / 1,118 (18.2%) | |
| projection on, before today | 567 / 1,118 (50.7%) | 364 exact / 46 wrong / 13 null-reference |
| projection on, title tier on | **615 / 1,118 (55.0%)** | **412 exact / 46 wrong / 13** |

The title tier's own contribution is 48 exact and 0 wrong, all SCN5A 28341781.
Identity, carriers, unaffected and counted extras are unchanged; 0 new values
fall on gold real-split rows. New-value precision is 412 / 458 = 89.96%, still
under the 90% gate by the 46 Kapplinger rows of section 1 and 100% outside
them. On tranche 04 the tier changes nothing. Across 497 archived cardiac
papers ([`smoke_cardiac_title_tier.json`](item2_measurement/smoke_cardiac_title_tier.json))
it stamps exactly one table, the same one. Unit suite: 3,020 passed.

What remains is not deterministic: 146 of the 255 null values are
model-authored rows. That is the derived-count contract gate already at the top
of `TASKS.md`. Per-person clinical rosters with a "Clinical diagnosis" column
(RYR2 29925740, 51 rows) need a reader of the row's own status cell; the
projection refuses them on purpose.

## 3. Sentence-level paper tier: stays off

Brett: non-cardiac tables are out of scope. On the cardiac calibration locks the
sentence tier adds 18 exact and 0 wrong (KCNH2 11854117). On the cardiac
archive ([`smoke_cardiac_sentence_tier.json`](item3_paper_tier/smoke_cardiac_sentence_tier.json))
it stamps two tables: 11854117 through a **footnote** ("... per unit time for
patients with the factor present ...") that happens to contain "patients with"
and "LQTS", and SCN5A 29709244 (23 rows, "Number of persons", a sodium-channel-
blocker challenge cohort) where gold records a real `2 / 1 / 1` split the
projection would overwrite. One accidental right and one wrong is not a
mechanism. The 11854117 rows remain a gap.

## 4. Manual acquisition worklist

[`docs/evidence/manual_acquisition_20260908/`](../manual_acquisition_20260908/README.md):
363 cardiac and BRCA2 gold papers have rows behind the acquisition ceiling
(1,314 rows). The top 150 hold 83.5%; 57 of them need no institutional access
(open-access copy exists, nothing usable on disk, or supplements missing from
an open-access body), 80 need a library browser session, 11 have no DOI, 2 need
a clean PDF. Built by `scripts/recall_audit/rank_manual_acquisition.py`, which
also ranks a no-gold PMID list by the abstract-only yield predictor.

## 5. Tranche 05: one live arm of current main

Design preregistered in [`PLAN.md`](PLAN.md) before opening, including the
registry state: tranche 04's never-run candidate slot is closed by an
append-only `abandon_arm` ledger event (new `setup_production_eval.py abandon`
subcommand, tested) so that 05 opens next in order. Run
`20260908_protocol_cont120_05_baseline` (121 attempts, 111 PMIDs). Results are
appended below once locked and scored.

### Result

Run `20260908_protocol_cont120_05_baseline`: 121 attempts / 112 PMIDs, seven
gene processes, all `completed`, 48.2 minutes of production, locked and scored
through the registry's `lock_and_score.sh`; the ledger records the arm. API
cost by the list-price proxy **$7.07** ([`budget.json`](item5_tranche05/budget.json):
gpt-5.6-sol $4.74 on 276 calls, grok-4.3 $2.13 on 140 calls of which 19 failed
without usage, Kimi $0.21).

Identity, paper-derived lane (all genes): **297 TP / 145 FP / 158 FN**, recall
65.3%, precision 67.2%, F1 66.2%, counted-extra precision 97.7% (7 extra rows
with counts). Cardiac four (118 attempts): 270 / 141 / 157, recall 63.2%,
precision 65.7%.

| gene | TP | FP | FN | recall | precision |
| --- | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 81 | 62 | 15 | 84.4% | 56.6% |
| KCNQ1 | 42 | 13 | 28 | 60.0% | 76.4% |
| RYR2 | 24 | 14 | 38 | 38.7% | 63.2% |
| SCN5A | 123 | 52 | 76 | 61.8% | 70.3% |
| BRCA1 / BRCA2 / MYBPC3 | 20 / 6 / 1 | 4 / 0 / 0 | 1 / 0 / 0 | | |

Four papers hold 73 of the 158 misses with few or no hits: RYR2 27452199 (34),
SCN5A 24721456 (14, rank 16 on the manual-acquisition worklist), SCN5A 22360817
(13, rank 18) and KCNQ1 32470535 (12). The per-gene RUN_STATUS files flag 20
attempts (SCN5A 13, KCNQ1 5, KCNH2 2) whose on-disk source scan is abstract-only
while the finalized extraction records say full text; the scorer reports every
attempt as `corpus_as_locked`. Recall is bounded by acquisition here as on every
prior tranche.

Counts, cardiac four ([`ablation_stdout.txt`](item5_tranche05/ablation_stdout.txt)):
carriers supplied on 66 of 374 positive-gold rows (47 exact), affected on 43 of
344 (35 exact, **10.2%** exact recovery), unaffected on 18 of 50 (3 exact).

**The projection lane, including the title tier, was inert on this tranche**, as
it was on tranche 04. Its metadata block is present on all 121 extractions with
the title tier enabled: 0 tables classified, 0 rows stamped. The rows it saw
were 575 deterministic rows without any carrier count (identity-only tables)
and 286 model-authored rows; not one attempt short-circuited on a deterministic
or router table (every paper went to the grok-4.3 primary, 50 of them through
verification), so the fast-path coverage audit had nothing to record either.
The strip ablation ([`ablation/`](item5_tranche05/ablation/)) is therefore
identical on every metric, off and on.

Preregistered rules: (1) newly supplied exactness, **uninformative** (nothing
supplied); (2) no manufactured split, pass (0); (3) nothing else moves, pass;
(4) supply on the addressable pool, **uninformative** by the rule's own clause,
the tranche contains no compendium-style count table and no per-proband
catalogue; (5) no wrong values to list; (6) title tier, no stamps to list. This
is the second consecutive null with a reason, not a failure and not a
confirmation.

What it says about the confirmation path: the count-table shapes the projection
targets sit in a minority of gold papers (the opened locks 02/03 held several;
04 and 05, 241 attempts, held none). Drawing whole mixed tranches will keep
producing nulls. A confirmation cohort has to be selected by a **gold-free
source property**, "the frozen source contains a deterministic count table or
a one-proband-per-row catalogue", applied to the still-unopened tranches under a
preregistered rule before any score is read. That design is recorded as the
next step in `TASKS.md`; it was not run today.

The canonical stratified figure regenerates only after a candidate-labelled
arm, so `run_eval.py score` left it unchanged; the run's own figures are
[`phenotype_count_recovery.png`](../../../benchmarks/codex_paper_eval/runs/20260908_protocol_cont120_05_baseline/figures/phenotype_count_recovery.png)
and [`gold_difference.png`](../../../benchmarks/codex_paper_eval/runs/20260908_protocol_cont120_05_baseline/figures/gold_difference.png).

## Reproduction

```bash
.venv/bin/python -m pytest tests/unit -q
.venv/bin/python scripts/replay_table_cohort_phenotype.py --title-tier \
  --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate \
  --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_03_candidate \
  --out-dir tmp/title_tier_replay
.venv/bin/python scripts/table_cohort_smoke.py --root results \
  --genes KCNH2 KCNQ1 SCN5A RYR2 --title-tier --out tmp/title_tier_smoke.json
.venv/bin/python scripts/recall_audit/rank_manual_acquisition.py \
  --out-dir docs/evidence/manual_acquisition_20260908 --max-papers 250
```
