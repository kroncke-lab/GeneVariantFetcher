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
  case/disease and control nouns must refuse; fix queued for the first commit
  after the tranche-04 lock, because runtime files cannot change while the arm's
  fingerprint is pinned.
- Refusal tallies are dominated by `per_person_clinical_row`, `no_carrier_count`,
  `model_authored_row` and `caption_mixes_cases_and_controls`, which is the
  conservative shape intended.

## 6. Tranche 04 (pending)

See `PLAN.md` for the design fixed before opening. Results, the strip-mode
ablation, the preregistered rule verdicts and the cost receipt are appended in
§7 once the arm is locked and scored.

## 7. Results

_Pending._
