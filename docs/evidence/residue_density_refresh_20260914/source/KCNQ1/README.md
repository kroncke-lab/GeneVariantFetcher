# KCNQ1 source and endpoint correction — 2026-09-14

The G589D count from PMID 21244686 was incorrectly treated as **492 affected
carriers**. The paper reports 492 genotype-positive carriers and continuous
QTc measurements, without the variant-specific binary clinical/ECG phenotype
partition needed here. The archived database explicitly describes its affected
count as inferred from the patient count. Both that partition and the single
compound Y171X observation are quarantined, with their total observed carriers
preserved in the decision ledger.

This review is bounded to the **three papers contributing the largest frozen
affected counts**. It does not establish that the remaining KCNQ1 literature
uses one consistent endpoint or contains independent people.

## Inputs for the new population union

| Input | Canonical missense clinical A | Canonical missense clinical U |
|---|---:|---:|
| Original frozen input | 6,514 | 1,757 |
| [clinical_input.csv.gz](clinical_input.csv.gz) | 5,848 | 1,595 |
| [clinical_input_cardiac_events.csv.gz](clinical_input_cardiac_events.csv.gz) | 6,022 | 1,757 |

Both ledgers preserve all **699 original clinical keys** and identity columns.
They include original counts, exact removals, corrected counts, review status,
and endpoint labels. The primary quarantines the reviewed cardiac-event splits;
the sensitivity retains them. **The latter is a mixed clinical-LQTS plus
cardiac-event sensitivity, not a pure cardiac-event cohort or LQTS-unaffected
denominator.**

Read corrected `affected_literature` and `unaffected_literature`, and filter
`n_literature > 0` before rebuilding the population union. Release population
alleles from unsupported clinical aggregates through the usual identity rules;
no population allele or gnomAD count is removed by this ledger.

## Exact decisions

| PMID | Frozen evidence | Primary decision | Sensitivity decision |
|---|---|---|---|
| 21244686 | G589D A492/U0 | Quarantine unsupported phenotype partition | Same |
| 21244686 | Y171X A1/U0 | Quarantine compound carrier without a resolved phenotype partition | Same |
| 23856471 | A341V A111/U41 | Preserve separately as cardiac-event evidence | Retain with explicit event endpoint |
| 23856471 | G589D A63/U121 | Preserve separately as cardiac-event evidence | Retain with explicit event endpoint |
| 23856471 | c.386+18089C>T A560/U0 | Quarantine whole-study size assigned to a modifier allele | Same |
| 32893267 | Additional A344= A1/U0 | Remove duplicate; retain verified A29/U0 | Same |
| 32893267 | Additional c.477+5G>A A1/U0 | Remove duplicate; retain verified A7/U0 | Same |

The seven decisions remove primary **A1,229/U162 across all types**, including
missense A666/U162. Some removed rows are real observations of a different
endpoint, so this is not a count of fabricated people. No unknown phenotype is
relabelled unaffected, and no clinically/ECG-affected LQTS case is removed merely
for being symptom-free.

[observation_corrections.csv](observation_corrections.csv) preserves each original
CSV line, database variant ID, count, source location, source URL, source hash,
reason and sensitivity disposition. The original source DB rows are available
in [selected_source_db_rows.csv.gz](selected_source_db_rows.csv.gz).

## Why these distinctions matter

**PMID 21244686:** the Methods select all available Finnish founder-mutation
carriers, including 492 G589D carriers and a person also carrying Y171X. The
paper studies QTc as a quantitative trait and reports separate ascertainment for
syncope history, medication and devices. Syncope information is available for
488 G589D carriers, not all 492. No specified clinical/ECG threshold count can
be recovered from the reported means or plots. Neither lack of symptoms nor
genotype positivity supplies the missing clinical LQTS A/U split. See
[Methods, Patients and Statistical analyses](https://pmc.ncbi.nlm.nih.gov/articles/PMC3032654/).

**PMID 23856471:** the founder-population splits explicitly concern cardiac
events: symptomatic carriers experienced events before age 35, whereas the
event-free group was older and untreated or treated only after age 35. The
symptom-free participants are still described as LQTS patients. Therefore the
111/41 and 63/121 partitions are valid event-risk evidence but cannot be silently
used as LQTS present/absent. The separate
[event observation ledger](preserved_cardiac_event_observations.csv) retains them.
The unrelated A560 for a protective modifier equals the total study enrollment
(224+152+184); that is not evidence that every participant carries its T allele.
See [Methods and replication populations](https://pmc.ncbi.nlm.nih.gov/articles/PMC3864834/).

**PMID 32893267:** this source contains genuine clinically diagnosed LQTS cases,
whose ECGs were centrally assessed. Its Supplementary Table S4 identifies
KCNQ1 cases separately from the SCN5A Brugada study. All 270 frozen records map
unambiguously to a KCNQ1 S4 row by nucleotide or protein identity. The 268
count-matching records equal the sum of the explicit European and Japanese
case columns. The other two records are duplicate one-count emissions for
variants already represented by larger correct counts. Thus the high affected
counts in this source are retained rather than removed for clinical
ascertainment alone. See [cohort methods and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC7790744/)
and [the complete 270-row count check](32893267_case_count_checks.csv).

No atrial-fibrillation or short-QT endpoint was identified among these three
reviewed KCNQ1 sources. That is a bounded result, not a screen of every remaining
paper. A variant associated with another phenotype elsewhere is not excluded
from an actual LQTS case count in this review.

## Reproduction and limits

[build_source_ledger.py](build_source_ledger.py) verifies the exact archived
database SHA-256, reproduces every original clinical-key count from the frozen
observations, checks the paper ranking, matches S4 counts, and writes
deterministic GZip/CSV with LF newlines. Run it with the
BayesianPenetranceEstimator Python environment. All evidence files are smaller
than 1.2 MB. [source_checks.json](source_checks.json) pins the source database,
observations, identity ledger, primary source texts, XLSX and script.

This does not change source databases or old frozen results. Repeated founder
families across publications, remaining symptom-versus-ECG partitions,
multivariant genotypes, age and treatment still need fuller adjudication. The
resulting neighborhood remains an empirical-posterior feature under the stated
gnomAD-unaffected assumption, not a calibrated patient probability.
