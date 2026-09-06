# Independent review of the frozen source references

Reviewed the eight prepared packets, query wording, original source files,
`reference_values.json`, and the frozen plan. No model responses or scores were
inspected for this review. The packet and reference files were not changed.

**Finding:** no outright numeric reference error was found for the frozen
population and endpoint queries. One affected count is interpretive and deserves
the separately reported sensitivity described below. Two representation limits
also affect how model abstentions should be interpreted.

## Integrity and scope

All eight packet SHA-256 values match `plan.json`. The reference SHA-256 matches
the plan (`8c41c3289d3d897703c9d134d9ee1fd7592f71c8975fc9ee0fb4c14571f8f3a8`).
Every packet's original source hash matches the current frozen source file, and
every `L` source unit matches Python `splitlines()` at the indicated location.
Query IDs match the reference keys: **43 queries, 24 numeric references, and
19 null references**. No null reference means zero.

The manually selected variants and specified endpoints make this a bounded
source-reading component test. They do not establish whole-paper discovery,
automatic endpoint selection, production acceptance, or population accuracy.

## Reference findings by packet

| Packet | Review finding |
| --- | --- |
| `scn200_table` | All nine values agree with Table 1: Gly1408Arg 14/4/9; c.3963+2T>C 10/2/8; Ala665GlyfsX16 9/6/2. The separate unknown values are 1/0/1; unaffected is not carrier minus affected. |
| `ryr304_relatives` | Documented-person carrier counts 3/2/2/2 are supported by subjects plus explicitly linked relatives. Each parent of subject 7 contributes only to the stated inherited variant. The father shared by subjects 1/2 is counted once. The cousin's G4772S does not imply the cousin has subject 14's other variants. No phenotype query is scored. |
| `ryr254_timepoint` | 8/7/1 matches the initial diagnostic phenotype, supported by S11, S21, and S22. The one negative is an explicitly described non-diagnostic carrier, not an unexamined remainder. Later treatment response does not replace the frozen endpoint. |
| `ryr189_aggregate` | Carrier counts 13 and 6 are explicit at L21. All four phenotype values should be null for the requested variant-specific manifestation endpoint: the symptomatic and asymptomatic subsets are not assigned to variants. Repeated Holters at L23 are not new people. |
| `ryr193_abstract` | Carrier 4 and symptom-negative 2 are explicit at A26–A27. Symptom-positive 2 is a defensible cross-sentence interpretation, but weaker than the other references; retain primary and show the sensitivity below. |
| `scn325_prior_counts` | All nine nulls are correct for newly studied human clinical participants. Literature and gnomAD counts at L222–L255 remain excluded by the explicit footnote at L348–L350 and functional-study methods. |
| `ryr258_roster` | 179/45/133 is correct for living baseline people and the specifically frozen **Previous symptoms** column. It is not the clinical diagnosis or combined symptoms/VA endpoint. One missing symptom record remains unknown. |
| `scn251_variant_list` | All six nulls are correct. The two target variants are listed, but the supplied list supplies no per-variant count or person identity. One variant-list row cannot establish one carrier. |

## RYR2 25814417 recount

The packet includes all 179 unique, contiguous IDs from 1 to 179. Directly
counting its `Previous symptoms` cells gives 133 `No`, 27 `Syncope`, 13
`Dizziness`, three `Syncope/Dizziness`, one `Syncope/Dizzeness`, one
`Syncope/ Dizziness`, and one missing `-`. Thus named symptoms total 45;
patient 37 at L266 is the sole missing record. The misspelling does not make an
otherwise explicit symptom disappear. The reference correctly avoids treating
VA or CVA findings as symptoms in this query.

The existing production rule's 91/62/26 living partition and 97/62/26 broader
partition answer different questions. Comparing the new 45/133 result with
those values as an improvement or regression would be an endpoint error.

## Interpretive sensitivity: RYR2 19398417 affected count

The abstract directly says four relatives carry W4645R, including two without
symptoms. It separately describes two CPVT diagnoses, segregation with disease,
and symptom emergence in later generations. Reading those statements together
supports the frozen affected value 2, given that the task permits derived
counts. However, the abstract does not print an explicit sentence assigning two
symptomatic carriers to W4645R, and `including two without symptoms` alone is
not a formal exhaustive partition. A cautious model can explain why it withholds
the affected integer.

Keep the 43-query frozen primary result unchanged. Report an additional
**42-query sensitivity excluding only `ryr193_abstract/q1_affected`**. Do not
change the reference after observing model results, turn a conservative null
into a numeric false positive, or use this interpretive field alone to justify
routine Astra inclusion. A literal-validator rejection can be correct even when
the source-adjudicated reading is reasonable.

## Representation limits

1. SCN5A 20031634's target c.3963+2T>C appears in the OCR source as
   `c.396312T.C`. The expected row values are unambiguous after identifying the
   splice row, but its identity requires OCR interpretation. A model's refusal
   to equate the malformed string is an identity/source-quality abstention, not
   evidence that it cannot read the row's count columns. Preserve this detail
   when classifying errors and validation failures.
2. SCN5A 25163546 uses Python `splitlines()` source IDs, whereas the earlier
   investigator note used newline-based `nl` numbering. Two embedded line
   separators create an offset by the relevant table. The packet includes the
   Table 6 title and both target rows, but its introductory span ends at L1184,
   just before the column header at L1185. This omission does not create a
   supported count or make any null reference incorrect; it does mean the
   packet is not a complete caption-plus-header representation. Leave the
   frozen packet intact and disclose this small completeness limitation.
