# Source adjudication for the incremental Astra test

Prepared from already-opened frozen sources before this test's API outputs.
This file is investigator reference material, not model input. No unopened
tranche was inspected. The corpus symlink guard passed. Paths below are relative
to `validation_runs/model_routing_20260906/frozen_sources/`; line numbers refer
to the frozen `GENE/PMID/PMID_FULL_CONTEXT.md` files, not reformatted packets.

The reference concerns source-supported counts for the stated population and
endpoint. It is separate from legacy gold agreement and acceptance by the
production literal-evidence validator. A count derived from people identified
across sentences may be correct while failing that validator. A null does not
mean zero.

## SCN5A 20031634: table columns and unknown phenotype

Recommended packet: lines 109–175 and 184–216. The table header at 109–112
distinguishes mutation carriers, mutation-positive BrS-ECG positive,
mutation-positive not BrS-ECG positive, undetermined phenotype, and
mutation-negative BrS-ECG positive. The methods at 164–175 define the target
BrS phenotype at baseline or after provocation; the cohort is not a symptom-only
cohort.

| Target | Carriers | Affected | Unaffected | Undetermined | Row |
| --- | ---: | ---: | ---: | ---: | --- |
| p.Gly1408Arg | 14 | 4 | 9 | 1 | 114, family B |
| c.3963+2T>C | 10 | 2 | 8 | 0 | 117, family E |
| p.Ala665GlyfsX16 | 9 | 6 | 2 | 1 | 121, family I |

These are explicit columns, not subtraction. Source OCR renders the family-E
nucleotide string as `c.396312T.C`; preserve that raw source alongside the
normalized target. The header and adjacent rows establish its column ownership.
Table 2 at 184–199 combines six undetermined carriers into a 61-person
not-BrS group; Table 1 separately reports 55 negatives and six undetermined.
Do not overwrite variant-specific unknown phenotype using the aggregate table.
Mutation-negative relatives in the rightmost column are not carriers.

## RYR2 30403697: relatives and per-variant inheritance

Recommended packet: lines 93–124. Carrier references count people explicitly
documented in the supplied table and family-history cells, including relatives.
They are derived counts; none is an explicit scalar total printed beside the
variant. Do not imply an exhaustive pedigree census beyond the documented
people. Do not score affected/unaffected in this packet because diagnostic
status, symptoms, and registry enrollment are different endpoints.

| Target | Documented carriers | Supporting rows | Derivation |
| --- | ---: | --- | --- |
| p.R417L | 3 | 107–108 | Subjects 1 and 2 plus their shared carrier father; paternal cis inheritance |
| p.R2028H | 2 | 113 | Subject 7 plus mother, from whom variant 1 was inherited |
| p.Y4721C | 2 | 113 | Subject 7 plus father, from whom variant 2 was inherited |
| p.G4772S | 2 | 120 | Subject 14 plus the explicitly variant-positive symptomatic cousin |

The father of subjects 1/2 is described only as asymptomatic; that does not
independently prove a negative diagnostic test. Subject 7's parents are
phenotypically silent heterozygous carriers, but each parent carries the
specified one of two variants, not both variants. The cousin's G4772S does not
establish carriage of subject 14's additional intronic variants at line 124.
Family-history relatives must not also be counted as new numbered subjects.

## RYR2 25435091: initial diagnosis versus treatment response

The frozen file has three lines, with source metadata, table, and letter
flattened into line 3. Include the table beginning `Clinical Characteristics of
the Cohort of Eight Family Members Carrying the Mutation in RYR2`, its footnotes,
and the letter from `To the Editor,` to immediately before `FUNDING`.

The source explicitly reports eight carriers of C2277R, seven with the CPVT
phenotype by exercise stress testing. The eight people are II:1, II:3, II:6,
II:8, II:9, II:15, III:4, and III:9. II:9 underwent an additional adrenaline test
but did not meet diagnostic criteria. Reference for **initial diagnostic EST
classification**: carriers 8, affected 7, unaffected 1. The source calls this
87.5% penetrance. The negative is a closed source-defined partition and also
identified by the explicit non-diagnostic finding, not a silent remainder.

Only three carriers had prior arrhythmic symptoms. Later treatment suppresses
arrhythmia in three, flecainide is added for five, and all are asymptomatic at
final follow-up. None of those timepoints replaces initial diagnostic status.
The source's four sudden deaths lack confirmed genotypes; do not add them to
the eight carriers. II:1 is already in the eight; do not add the proband twice.

## RYR2 18929323: per-variant aggregates and repeated measurements

Recommended packet: lines 19–25. Line 21 explicitly identifies P2328S in
13 gene-positive patients and V4653F in six. Line 23 describes a repeated Holter
assessment in six carriers of each variant; those are not additional people.
The 19 control participants are a separate selected cohort, not unaffected
carriers of either target variant.

| Target | Carriers | Affected | Unaffected |
| --- | ---: | --- | --- |
| p.P2328S | 13 | null | null |
| p.V4653F | 6 | null | null |

For the planned **clinical manifestation** endpoint, 13 historical symptomatic
people and two asymptomatic people with negative exercise tests are reported
across the combined cohort without a variant assignment. They cannot supply
either variant's exact split. The paper calls the overall group CPVT patients;
that cohort label must not silently change the planned manifestation endpoint.
Table 1 is cited but absent from the frozen source.

## RYR2 19398417: abstract-only clinical evidence

The frozen `.md` file is actually JSON. Read its complete `abstract` value; do
not mistake the file extension for a full paper. The results say that CPVT was
diagnosed in two of nine evaluated relatives, and W4645R was found in four
relatives, including two without symptoms. The variant is said to segregate
with disease with incomplete penetrance. The conclusion describes silent
transmission across two generations before symptoms in the next two.

Reference carriers = 4 and symptom-negative carriers = 2 are explicit. A
symptom-positive reference of 2 requires the cross-sentence disease-segregation
linkage to the two diagnosed relatives; it is less direct than the other two
fields. If the frozen task requires an exact literal per-variant affected
statement, affected must be null instead. Keep the raw source-adjudicated
clinical interpretation separate from literal-validator acceptance, and do not
penalize an explained abstention as a numeric false positive. Do not use nine
screened relatives as carriers or infer seven unaffected carriers.

## SCN5A 32533946: literature and population counts in a functional study

Recommended packet: lines 222–255, 345–350, and 1099–1117. Suitable targets are
p.Thr220Ile, p.Gly752Arg, and p.Glu1784Lys. The table contains BrS, LQT,
unaffected, and gnomAD counts, but the caption calls them literature/gnomAD
counts and the footnote at 348–350 explicitly says the individuals' disease
counts come from a literature curation. Methods describe automated patch
clamping, not enrollment of these clinical cohorts.

All three targets have **null current-study clinical carrier, affected, and
unaffected counts**. Their variant identities and functional data remain valid.
Do not turn the assay's 83 variants into people, add BrS/LQT populations whose
overlap is unspecified, or interpret gnomAD allele counts as unaffected people.

## RYR2 25814417: living cohort, unknown tests, and timepoints

Recommended packet for the living enrollment cohort: lines 31–33 and 49–61,
plus the complete Supplementary Table 3 at 225–408. Carriers = **179 living
screened relatives**, explicitly at line 33; 1,404 is screened N. Six confirmed
historical SCD carriers are a distinct group. The 53 people who initially
refused testing and the three subsequently genotyped SCD/aborted-SCD cases
belong to a later follow-up narrative. Do not sum these groups into the requested
living enrollment count or assume that every historical SCD was genotyped.

The table contains 179 distinct, contiguous patient IDs. Independently counted
source values:

| Field | Positive | Negative | Missing |
| --- | ---: | ---: | ---: |
| Previous symptoms | 45 named symptoms | 133 No | 1 |
| VA in basal test | 69 Yes | 81 No | 29 |
| CVA in basal test | 42 Yes | 108 No | 29 |

The current deterministic production rule uses symptoms OR ventricular
arrhythmia: **91 positive / 62 negative / 26 uncertain among living people**.
That rule already exists; recovering it is not unique Astra capability or new
parser upside. Adding six separate genotyped fatal cases yields the existing
broader protocol observation 185 / 97 / 62 with 26 uncertain. It is a different
cohort from the living-only request. Neither a CVA-only test nor symptom-only
count equals that combined endpoint. If the new model task does not freeze an
endpoint/predicate, score only the 179 carrier count and retain phenotype
outputs as unscored endpoint sensitivity. Do not force the legacy 73/106 split.

## Reserve source: SCN5A 25163546

The recovered 53-page supplement is not an individually identified SCN5A
patient roster. Supplemental Table 6 lists 20 SCN5A variants at lines 1378–1397,
under a gene/transcript/exon/cDNA/protein header at 1180–1183. It gives no
per-variant person frequency or patient IDs. The 639 enrolled DCM patients at
96–103 cannot be copied onto each variant. This is a useful source-availability
abstention control, not a source-complete carrier count opportunity.
