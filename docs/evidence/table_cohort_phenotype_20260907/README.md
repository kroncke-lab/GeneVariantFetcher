# Table-cohort phenotype projection: affected-count supply, 2026-09-07

Brett's request: get affected-count exactness "way up", adversarially, using agy
and grok to think it through, with a $100 API budget for testing. GPT's 2026-09-07
assessment framed the target as exact clinical-fact recovery: on the two locked
`cont120_02/03` candidate runs the pipeline supplied an affected value on
221 of 1,118 positive-gold rows (18.2% exact), although 92% of what it supplied
was exact.

## 1. Where the missing affected values were

Row-level decomposition of the 1,118 positive-affected gold rows (four cardiac
genes, both locks; `results/performance_assessment_20260907/` is the source of
the aggregate numbers):

| rows | state |
| ---: | --- |
| 477 | identity matched, `affected` NULL, pipeline **carriers == gold affected** |
| 232 | identity missed |
| 203 | affected exact |
| 126 | matched, carriers and affected both NULL |
| 57 | affected NULL, carriers wrong (45 are 20129283's per-centre gold rows) |
| 18 | affected supplied but wrong |
| 5 | affected NULL, carriers == gold carriers != gold affected (a real split) |

Gold convention: affected == carriers on 1,018 / 1,118 positive rows and
unaffected == 0 on 1,023: a variant found in N disease-ascertained probands is
recorded as N carriers, N affected, 0 unaffected.

The 477 rows sit almost entirely in deterministic-parser tables whose caption or
count column names the case series, with the phenotype fields left NULL because
the tables have no phenotype column:

| paper | rows | table evidence |
| --- | ---: | --- |
| SCN5A 20129283 (Kapplinger 2010) | 239 | Table 4 "Compendium of Brugada syndrome-associated SCN5A mutations", column "No. of unrelated individuals"; the column sums to 451 = 438 genotype-positive cases + 13 double-mutation patients stated in the results |
| KCNH2 + KCNQ1 26496715 | 83 | "Summary of putative LQT1/LQT2-associated mutations", column "No. of patients" (supplement-only source) |
| SCN5A 28341781 | 52 | fixed-width Table 2 "Included SCN5A Mutations and Variants"; count is a parenthetical, column label garbled; cohort named only in the title |
| SCN5A 27566755 | 47 | "Table S1: List of Mutations by Coding Effect, Location, and Frequency in 406 LQT3 Patients", column "COUNT" |
| KCNH2 11854117 | 18 | body is an empty shell; column "No. of Subjects" under "HERG mutations" |
| RYR2 28237968 (adversarial) | 13 | one proband per row with relatives columns; gold counts the family, so copying the row count would be wrong on 6 of 13 |

The parser already recognised 20129283's Table 2 controls ("Rare control" /
"Polymorphism" status cells) and emitted `unaffected = N`, `affected = 0`; the
always-on guard then cleared the zero as unsourced. The LQTS-compendium parser
(19841300) already emits `affected = count` for its case rows with
`count_type = case`, which the guard accepts. The projection generalises that
existing, gold-consistent behaviour to any table whose printed caption or header
names the cohort.

## 2. Reviewer consultation before implementation

Both CLIs received the same self-contained brief (`reviews/brief_*.txt`);
grok additionally read up to twelve repository files. Full outputs are in
`reviews/`.

| point | agy | grok | disposition |
| --- | --- | --- | --- |
| Is a case-series count column a source-backed affected count? | Yes when the table structurally names patients/cases; a direct read of the authors' claim | Yes; draw the line at the count column's own noun (and the caption), not at paper-level enrollment | Implemented as tier 1 only |
| Paper-level (title/abstract) ascertainment | Crosses the contract line; do not ship | Do not ship; it is paper-level NLP that fires on every registry and autopsy series | Implemented behind `GVF_TABLE_COHORT_PAPER_ASCERTAINMENT`, default off; adds 18 rows on the locks (`replay_tier2_diagnostic/`) |
| Control tables | affected = 0 is sound when the source defines the cohort as healthy; never treat allele counts as people | Leave affected NULL unless the count is people; do not turn "2,600 reference alleles" into 408 unaffected humans | People-unit control tables get the closed N/0/N partition; allele-unit tables are left exactly as the row parser wrote them |
| Exclusions | add in-vitro/cells/assays/HEK/Xenopus, SCD/SUD, relatives, sub-clinical; "phenotype" and "genotype-positive" are too aggressive | add autopsy/SIDS/post-mortem, literature/published, occurrences, hom/het, incidental, "characteristics of", case AND control columns; "frequency"/"allele"-before-control/"carriers" too aggressive | All additions taken except "cohort" and "putative" (both kill legitimate compendia); "phenotype" dropped from caption exclusions, kept for headers; "frequency" only excluded as an allele-frequency column; control classification runs before exclusions |
| Measurement | within-arm on/off ablation is the right instrument; a paid baseline arm is wasteful; require ≥ 90% conditional exactness | same; forecast from the $0 replay, gate the newly supplied set at ≥ 90%, zero new affected on gold real-split rows, counted extras non-decreasing | Adopted in `PLAN.md`; one paid arm on tranche 04 with a strip-mode ablation |
| Next lever | fixed-width identity fast path that skips phenotype tables (30059973) | same, plus caption reconstruction and role rebinding | Queued in TASKS.md |
| Unproven claim | the 55–60% forecast assumes every one of the 477 converts | the exclusion list as drafted would eat 20129283 and 27566755 | Checked in the replay below: 50.7%, and neither paper was eaten |

Both reviewers rejected "referred for genetic testing means affected" as a
biological claim; the implementation follows the lab's gold convention for
disease-ascertained probands and says so in the provenance rather than hiding it.

## 3. What was built

`pipeline/table_cohort_phenotype.py`, hooked at the extraction success boundary
after the audited patient-row lane and before the always-on guard. Rules and
provenance are specified in `docs/EXTRACTION_CONTRACT.md` ("Deterministic
table-cohort projection") and `docs/ARCHITECTURE.md`. Supporting changes:
`TABLE_COHORT_PHENOTYPE_SOURCE` joins the code-owned stamp set that is scrubbed
from model output; the guard accepts the code-owned control zero;
`refresh_run_db.py` preserves stamped observations; two settings flags. Tests:
`tests/unit/test_table_cohort_phenotype.py` (21, including the extractor
success-boundary hook). Offline suite before the final exclusion tweak: 2,989
passed.

## 4. Calibration: zero-LLM replay on the two opened locks

`scripts/replay_table_cohort_phenotype.py` re-runs the module over each archived
per-PMID extraction JSON with the run's own frozen source text, re-applies the
guard, fills the derived values into a copy of the locked `predictions.json`
(fill-null-only), and rescores both arms with the harness matcher. Nothing in the
locked run directories was written.

Cardiac four, pooled (`replay_tier1/replay_summary.json`):

| field | positive gold | supplied off/on | positive exact off/on | exact recovery off/on | wrong off/on | conditional exact off/on |
| --- | ---: | --- | --- | --- | --- | --- |
| carriers | 1,280 | 867 / 867 | 778 / 778 | 60.8% / 60.8% | 89 / 89 | 89.7% / 89.7% |
| affected | 1,118 | 248 / 658 | 203 / 567 | **18.2% / 50.7%** | 44 / 90 | 82.3% / 86.3% |
| unaffected | 220 | 176 / 176 | 101 / 101 | 45.9% / 45.9% | 16 / 16 | 90.9% / 90.9% |

Identity TP/FP/FN 1189/241/289 and 42 counted-extra rows in both arms. Newly
supplied affected values: 366 exact, 46 wrong, 11 on unmatched or gold-null
rows; 0 on rows where gold reports `affected != carriers`. Per paper: KCNH2
26496715 0→53 supplied / 0 wrong; KCNQ1 26496715 11→41 / 0 wrong; SCN5A
27566755 0→44 / 0 wrong; SCN5A 20129283 0→283 / 46 wrong.

**All 46 wrong values are one gold convention.** Kapplinger's Table 4 prints the
pooled "No. of unrelated individuals" with a "Testing center" column, and gold
stores one row of `1` per centre (A735V ×4, E1225K ×4, G1408R 1+1+2 ...). The
pipeline's carriers are already scored wrong on those rows for the same reason,
and the extra gold rows are identity misses. Excluding that convention the newly
supplied set is 366 / 367 exact (99.7%); including it, 88.8%. This is a curator
question about row granularity, not something to patch in code
(`replay_tier1/new_values_*.csv` lists every row).
`gold_20129283_per_centre_rows.csv` quantifies it: 47 variants carry 117 gold
rows; for 31 of the 46 with a Table 4 entry the pooled "No. of unrelated
individuals" equals the sum of the per-centre rows (a pure granularity
convention), while for 15 the gold rows do not sum to the table count either
way (D1243N 1+1+5 against 5; D1275N 1+1 against 3), so those rows need a
curator, not a parser.

`guard_clear_audit.md` closes a tempting alternative: restoring every value the
always-on guard cleared on these locks would recover 10 exact model-row
affected values and add 6 wrong ones, and its zero restorations land only on
zero-gold rows. Relaxing the guard is not a lever; the projection recovers the
deterministic-table rows without touching model rows.

## 5. Gold-free transfer smoke

`scripts/table_cohort_smoke.py` runs the module over every archived extraction
JSON under `results/` and the evaluation runs and lists every table it would
touch (`transfer_smoke_noncardiac.json`, `transfer_smoke_cardiac_results.json`).

- Non-cardiac (BMPR2, BRCA1, BRCA2, MYBPC3, APOE, TTN, LMNA; 1,909 papers):
  every classification came from a column-labelled case count ("Number of
  patients", "Number of Index Patients", "Patientsnumber", "No. of patients",
  "Case count"). One shape was wrong in spirit: BMPR2 26820968 "Gene mutations in
  Chinese CTEPH patients and PE without PH patients" mixes two patient groups;
  captions containing "without", "versus" or "compared with" now refuse.
- Cardiac (493 historical papers): ten distinct tables. Nine are defensible
  (Tester 2005 "No. of patients" compendia, Kapplinger Table 4, Ware 2013
  "Cases (n=2111)", a 2025 RYR2 "Number of Patients" sheet, 26496715, 27566755).
  **One is a bug:** SCN5A 25904541's fixed-width header "BrS (2111) / LQT
  (2888) / Control (8975)" was joined by the parser into "BrS + LQT + Control"
  and classified as a control column (342 rows stamped). A count label mixing
  case/disease and control nouns must refuse; the fix landed in the first commit
  after the tranche-04 lock (runtime files cannot change while an arm's
  fingerprint is pinned), and the re-run smoke (`transfer_smoke_cardiac_results.json`)
  now refuses those 342 rows under `count_column_mixes_cases_and_controls` and
  classifies nothing new.
- Refusal tallies are dominated by `per_person_clinical_row`, `no_carrier_count`,
  `model_authored_row` and `caption_mixes_cases_and_controls`, which is the
  conservative shape intended.

## 6. Tranche 04: one live arm, within-arm ablation

Design fixed in `PLAN.md` before opening. Run
`20260907_protocol_cont120_04_baseline` (registry arm label `baseline`; current
`main` at `e85dcff6` with the projection on; 120 attempts, 110 PMIDs, six gene
processes, all `completed`/`ok`, 50 minutes). Locked and scored through the
registry's own `lock_and_score.sh`; the consumption log records the arm. The
`strip` replay (`tranche04_ablation/`) nulls every field stamped
`table_cohort_phenotype_v1` in the archived extraction output and rescores.

Arm headline (all genes): identity TP 491 / FP 145 / FN 231, recall 68.0%,
precision 77.2%; KCNQ1 36/4/81 and SCN5A 103/64/84 carry most of the misses.
Cardiac four: affected exact recovery 51 / 545 positive-gold rows (9.4%),
carriers 349 / 665 (52.5%), unaffected 6 / 71.

**The projection was inert on this tranche.** Its metadata block is present on
all 120 extractions: 0 tables classified, 0 rows applied or stamped. The
deterministic rows it saw were 321 per-person clinical rows (refused by design)
and 4 rows under a "Symptomatic (Yes/No)" column; the other rows were 247
model-authored and 81 without a carrier count. Off and on arms are therefore
identical on every metric.

Preregistered rules: (1) newly supplied exactness, **uninformative** (nothing
supplied); (2) zero new affected on gold real-split rows, pass (0); (3)
identity, carriers, unaffected and counted extras unchanged, pass; (4) supply on
the addressable pool, **uninformative** by the rule's own clause: the tranche
contains no compendium-style count table. This is a null with a reason, not a
failure, and not a confirmation. Per `PLAN.md` no second arm is spent on
tranche 05.

## 7. What the tranche-04 gap looks like instead

Per-person rows dominate: 270 matched rows where gold records the row's patient
as `1 / 1 / 0` and the pipeline holds carriers = 1 with affected null, plus 8
rows where gold affected > 1.

- **RYR2 28404607 (229 rows): a gold-quality flag.** "Supplemental Table 1:
  Compendium of variants identified by WES testing" lists RYR2/CASQ2 variants
  found in 6,517 individuals undergoing clinical whole-exome sequencing for any
  indication; the paper's point is the background frequency of CPVT-associated
  variants (8.8% of the WES cohort, 97.7% VUS). Gold stores 263 of its 269 rows
  as `1 / 1 / 0`. Deriving affected there would match gold and be wrong; the
  projection's `exome` / `sequencing` exclusions refuse it. This is exactly the
  "enrollment is not diagnosis" case both reviewers warned about, and it needs a
  curator decision before anything counts those rows as affected.
- **RYR2 29925740 (41 + 8 rows).** "Table S1. Subject Clinical and Genetic
  Characteristics": one CPVT/LQT1 patient per row with "Most severe symptom" and
  "Clinical diagnosis" columns. Gold treats each diagnosed subject as affected.
  A per-person extension gated by the same caption/header classifier would take
  these, but the header exclusion refuses tables with clinical columns on
  purpose (28237968's relatives table is the counter-example), so this stays a
  documented candidate, not a shipped rule.

Cost (`budget.json`): 412 calls, $7.75 by the repository list-price proxy
(gpt-5.6-sol $5.62, grok-4.3 $2.04, Kimi $0.08; 12 failed calls without usage
records). $92.25 of the $100 testing budget remains unspent.

## 8. Post-lock fix and live integration check

The cardiac smoke's one bug (§5) is fixed in the commit after the lock: a
count label that names both a case/disease group and controls
("BrS + LQT + Control", SCN5A 25904541) now refuses with
`count_column_mixes_cases_and_controls`; the smoke and the two-lock replay were
re-run unchanged elsewhere. Because tranche 04 never exercised a stamped value
through migrate, trust gate and projection live, a three-paper calibrated
`gvf-run` (SCN5A 20129283 + 27566755, KCNH2 26496715, KCNQ1 26496715) was run
after the fix into `results/table_cohort_live_check_20260907/`; its outcome is
recorded in `live_check.md`. These are opened calibration papers: the check
proves the persistence path, not generalisation.

## 9. Pre-push verification and corrections

Independent local verification reproduced the aggregate scores in both opened
replays and the tranche-04 null, and checked all 15 pinned tranche-04
prediction/selection/setup, run-status and trace-manifest hashes. No paid
extraction or new tranche was used. The recorded tranche-04 figures were
visually inspected; the canonical candidate membership remains unchanged.

Four review findings were corrected before pushing:

- Control classification returned before the caption, clinical-header and
  people-count checks. A synthetic "Functional assays in control cells" table
  therefore received a trusted 0/N phenotype partition. Controls now pass the
  same exclusions as cases, including mixed disease/control captions and
  adjacent case/control columns. Allele-unit controls still certify nothing.
- `Table 2` did not resolve `Table 2. ...` on one line, hiding exclusion words
  such as "patients and controls". Inline captions now resolve without matching
  `Table 20` or a body sentence such as "Table 2 shows ...".
- Strip-mode ablation erased pre-existing counts when the projection merely
  certified them. Each new audit records the phenotype values the ordinary
  guard would keep without the projection; strip restores those values. Legacy
  stamped records lacking this audit now refuse reconstruction. Tranche 04 has
  no stamps and its original null remains reproducible.
- The new-value CSV performed a fresh fuzzy gold lookup instead of using the
  scorer's one-to-one assignment. That reused a reference row for K1493del and
  counted two unscored notations as exact. The audit now uses scored reference
  assignments and includes the actual gold notation.

**Corrected row audit:** 364 new exact affected values, 46 wrong, and 13
unmatched/merged/gold-null log entries. The 410 scored additions reconcile with
the aggregate supplied-count increase (248 to 658), and 364 exact additions
reconcile with exact recovery (203 to 567 of 1,118 positive-gold rows). Identity,
carriers, unaffected, counted extras and the zero real-split hits are unchanged.
The earlier 366-exact row-audit total in this historical report is superseded;
the aggregate recovery table was already correct. Corrected artifacts:
[`review_replay/`](review_replay/).

The original assertion that *all 46* disagreements were per-centre duplicates
was too broad. Forty-five match repeated reference variants; K1493X has one
reference row asserting 2 carriers/affected while Table 4 prints 1 unrelated
individual. The source-frozen table also prints K1493del = 2, while the scorer
assigns it a per-centre reference row of 1. These counts remain source/reference
curation questions, without a paper-specific runtime patch. Conditional
exactness is 364/410 (88.8%) unexcluded, or 363/364 (99.7%) after the plan's
duplicate-reference exclusion removes 45 wrong and one exact value.
This is calibration, and tranche 04 supplies no
independent confirmation.

Regression tests cover the unsafe control shapes, inline caption resolution,
preservation of existing case/control counts during ablation, refusal of legacy
stamps without prior-count evidence, and reference assignment with duplicate
gold rows. Verification commands and scratch outputs are under the ignored
`tmp/verify_table_cohort_20260907/`; historical locked inputs were not rewritten.

Final validation: **2,996 unit tests passed**, plus nine bounded end-to-end and
negative benchmark tests; repository-wide Ruff lint/format and pre-commit hooks
passed. The 497-record cardiac archive scan leaves all shared classifications
unchanged and still refuses the 342 pooled case/control rows. Its five additional
table entries are the already-stamped live-check records and supply no new values.
The 1,911-record non-cardiac scan adds no classification and removes the three
older BMPR2 mixed-caption entries already discussed in §5. Machine-readable
checks, source hashes and scan differences: [`review_verification.json`](review_verification.json).
