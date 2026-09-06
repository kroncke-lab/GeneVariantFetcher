# Source checks after tranche 02 locked

These checks explain observed errors; they did not change the reader, predictions,
reference values, or registered decisions. The second tranche uses the identical
candidate. All locations below refer to the candidate's actual extraction input,
which was SHA-identical to the paired baseline input for all 120 attempts.

## SCN5A, PMID 27566755: carrier gain is not phenotype capture

Carrier absolute error fell 198 → 70, accounting for 128 of the tranche's net
131-unit reduction. Affected error stayed 406 and neither arm supplied any
A/U values. This paper was previously scored and is excluded from the
predeclared previously-unscored subset.

The recovered supplement starts at line 171. Table S1, lines 177 onward, lists
mutations and a COUNT column in 406 LQT3 patients. Its caption establishes the
study cohort, but this alone does not identify the per-variant clinical endpoint
needed to interpret affected/unaffected. The body also describes asymptomatic
subjects, cardiac-event outcomes, and a 391-person follow-up subset. A generic
carrier-to-affected copy would silently choose an endpoint. Next work should
retain the table's cohort and endpoint context and reconcile it with the stated
extraction contract before filling phenotype counts.

Source: `benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate/production_runs/SCN5A/20260905_175554/pmc_fulltext/27566755_FULL_CONTEXT.md`.

## RYR2, PMID 25814417: phenotype definition and denominator

Both arms output G357S with N/A/U = 185/97/62; the locked v2 reference is
185/73/106. The carrier value is exact. The reference's phenotype split sums to
179 living relatives, while its carrier total includes six genotyped sudden-death
cases; see `benchmarks/count_semantics_eval/ADJUDICATIONS_20260810.md`, line 33.
The six-person difference is deliberate and must not be “corrected” by forcing
A + U = all carriers.

The actual input includes Supplementary Table 3 (line 226), with 179 distinct
patient rows. Its separate columns are previous symptoms, VA in a basal test,
and CVA in a basal test. Direct tabulation gives:

| Field | Source values |
| --- | --- |
| Previous symptoms | 133 No, 45 with a named symptom, 1 missing |
| VA in basal test | 69 Yes, 81 No, 29 missing |
| CVA in basal test | 42 Yes, 108 No, 29 missing |

Neither simple VA nor CVA tabulation reproduces the reference 73/106. The
existing `pipeline/patient_row_phenotype.py` already applies an audited rule:
any selected symptom/objective field positive, all selected fields negative,
or missing without a positive finding. This produces 91/62/26 among living
relatives, with six confirmed SCD cases added to affected for 97/62/26 overall.
The 26 uncertain people stay unknown. This was already diagnosed in the
22-paper audit; it is not a new parser opportunity. The exported prediction's
short citation does not expose that arithmetic, but the implementation does.
Keep this case out of forecasts for reference-matching parser gains. Reconcile
the phenotype endpoint and cohort with the reference through a separate
adjudication; do not force the 26 unknown people into unaffected. The locked
reference and score remain unchanged.

Source: `benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate/production_runs/RYR2/20260905_175554/pmc_fulltext/25814417_FULL_CONTEXT.md`.
SHA-256: `a131c2488e8250c5383aba2ff98bd8ae2b681ae90b1e32c50ec1707b23c090ef`.

## BRCA2, PMID 26848529: incomplete reference and source conflicts

The reference provenance is `collaborator_approved_nonexhaustive`: it contains
three identities, while the article describes many more. Identity extras fell
82 → 74, but carrier-bearing extras rose 0 → 74. These rows account for the
apparent collapse in any-count-bearing precision on the full tranche. A/U
values for this paper remained blank, and phenotype-bearing extra rows fell
across the tranche overall (affected 11 → 6; unaffected 6 → 4).

Supplementary Table 4 begins at line 279; line 281 explicitly distinguishes
**Total cases**, **Carrier**, and **Population frequency**. Rows below it contain
variant-specific carrier integers: e.g. c.262_263delCT has Carrier 2 versus
Total cases 518. Thus many extra identities/counts can be supported by source
without appearing in this non-exhaustive reference. We retain all official
false-positive accounting and do not claim that all 74 supplied values are
validated. One specific conflict remains: c.2808_2811delACAA is predicted as 2
from main Table 3, while Supplementary Table 4 lists 3 (line 289). Evidence must
retain table/cohort provenance and surface such conflicts, rather than choosing
the maximum or summing duplicate appearances.

Source: `benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate/production_runs/BRCA2/20260905_175554/pmc_fulltext/26848529_FULL_CONTEXT.md`.
SHA-256: `4888932616b5cef48066009f20cec9705fcce83297bfeb90cc6ed8cc80e391b8`.

## Previously unscored examples: actionable reading and representation gaps

The following additional checks were selected after scoring from the largest
remaining A/U errors in the predeclared previously-unscored subset. They are
root-cause examples, not another validation sample or new measured improvement.

**RYR2 18929323:** actual source line 21 explicitly associates P2328S with
13 carriers and V4653F with 6. Baseline emits both integers. Candidate emits
both identities but leaves all count fields blank; its extraction JSON still
quotes the exact per-variant carrier numbers in notes and key quotes, while
saying the clinical split is not available by variant. This adds 19 carrier
error units. Recovering source is unnecessary. A field-specific evidence rule
should preserve explicit carrier counts while abstaining separately on A/U.
Do not distribute the aggregate 13/19 symptom history across the two variants.
This is an observed extraction-record regression; it does not establish that
a particular changed code line caused the stochastic model output.

**MYBPC3 20433692:** source contains Additional file 1 twice, including the
folded supplement at lines 347 onward. It has mutation, family, patient ID,
mutation-positive flag, and age-at-diagnosis columns; the legend explicitly
defines “No Dx” as unaffected/healthy. Continued rows have lost leading merged
cells and shifted patient IDs/values left (e.g. after D75N), and some wrapped
ECG cells create spurious rows. Current A/U error is 41 versus baseline 35,
with affected supply 3 → 0 and unaffected 4 → 2. This supports restoring the
original DOC table grid, carrying forward only valid mutation/family cells,
and deduplicating by family + patient + mutation before interpreting diagnosis.
Do not infer genotype from symptom status: the table includes affected
non-carriers. Main Table 3 also pools historical families and must not be
substituted for the current-study patient roster just because it has tidy
HCM/healthy counts.

**KCNH2 15364333:** both arms miss the intronic identity entirely despite a
28.8-kB body with the old T1945+6C notation and explicit carrier discussion.
Phenotype capture cannot improve until the identity is represented. The body
also distinguishes genotyped/obligate carriers and baseline/expanded ECG
phenotyping, so a recovered identity still needs count-endpoint context.

**SCN5A 22885917 and KCNH2 9693036:** actual sources are 33.2 kB and 13.7 kB,
respectively; neither size nor a FULL_CONTEXT filename proves that the
count-bearing tables are present. These remain examples of why acquisition
success should be measured by required table/roster evidence, not body bytes.
No new download success is claimed from these frozen-source tests.
