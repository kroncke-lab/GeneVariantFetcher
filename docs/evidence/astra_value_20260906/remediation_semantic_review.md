# Independent semantic review of the completed Grok remediation

Reviewed all eight locked remediation responses against the frozen packets,
with detailed checks of the nine formerly Astra-only numeric answers, the six
RYR2 18929323 questions, and both negative-control packets. All eight receipt
digests match `remediation_locked.json`; all eight frozen original source
digests still match their packets. No requests were made, and no primary or
remediation outputs, references, or scores were changed.

**Finding: no new semantic count error was found behind the matching numbers.**
The review supports 23 numeric answers and 19 abstentions for the stated packet
questions. The remaining numeric answer, RYR2 19398417 symptom-positive 2,
retains the previously documented interpretive uncertainty. Thus the unchanged
43-query primary and 42-query sensitivity remain the appropriate separate
readouts. This does not establish whole-paper extraction or production
acceptance.

## Nine answers formerly returned only by Astra

| Packet | Grok values, C/A/U | Independent finding |
| --- | --- | --- |
| `ryr254_timepoint` | 8/7/1 | Strong support. S11 and S22 explicitly state eight C2277R carriers and seven with the EST phenotype. S21 and the S5–S6 table identify II:9 as the one non-diagnostic carrier. The answer uses the initial diagnostic assessment, not prior symptoms, treatment response, or final asymptomatic follow-up. |
| `ryr193_abstract` | 4/2/2 | Carrier 4 and symptom-negative 2 are explicit at A26–A27. Symptom-positive 2 remains interpretive; see below. |
| `ryr258_roster` | 179/45/133 | Strong support for the frozen living baseline **Previous symptoms** endpoint. The cited rows and missing-person exclusion are correct; detailed recount below. |

Grok recovered all nine frozen numeric answers under the changed package.
Eight have strong source support; the same ninth remains uncertain for both
models. This supports a cheaper-reader recovery result on this selected panel,
not an inference that any single changed transport or prompt setting caused
the recovery.

### RYR2 25814417: people, endpoint, and source ownership

An independent parse of the packet found **179 unique, contiguous patient IDs
1–179**. The `Previous symptoms` column contains:

- 133 `No` cells;
- 27 `Syncope`, 13 `Dizziness`, three `Syncope/Dizziness`, one
  `Syncope/Dizzeness`, and one `Syncope/ Dizziness`: 45 symptom-positive people;
- one missing `-`, patient 37 at L266.

The affected claim cites exactly all 45 positive patient rows, plus the table
title/header at L225/L229. It cites no negative or missing row. The unaffected
claim cites all 133 explicit-negative rows **and L266**, the missing record
that its reason explicitly subtracts. Including that extra source is justified
by the stated exclusion; it does not count patient 37 as unaffected. The
reason's arithmetic is 179 minus 45 minus 1, consistent with direct counting.

L33 binds the 179 living mutation-positive cohort to p.G357S, and L225 identifies
the table as basal features of living mutation-positive subjects. L229 labels
the endpoint columns. The carrier claim cites L33, while the two derived
claims cite the table title/header and rows. The complete response therefore
has the required variant/cohort ownership. For portable, independently
validated claims, repeating L33 in each derived claim's provenance would make
that ownership more self-contained; the omission is not evidence of an
incorrect count in this packet.

Neither VA nor CVA cells contribute to this symptom tally. The response keeps
the six deceased genotyped carriers, other historical sudden deaths, and later
follow-up cohorts outside the living baseline denominator. The 45/133 split
must not be presented as a replacement for the existing combined symptoms/VA
91/62 living partition or 97/62 broader partition: those are different questions.

### RYR2 19398417: retain the uncertainty and quote failure

The affected claim cites A23 as well as A26–A27 and says four carriers minus two
without symptoms matches two CPVT diagnoses. This is better citation coverage
than a subtraction supported only by the carrier sentence, but it does not
fully resolve the semantics. Two diagnoses among nine evaluated relatives are
not a literal statement that exactly two of the four carriers had symptoms;
`including two without symptoms` does not formally close the partition. The
abstract's broader segregation and symptom narrative supports the clinical
reading, with the uncertainty already recorded before remediation.

The quote combines two nonadjacent passages with `...`. It is not one
verbatim contiguous quote and correctly fails the mechanical quote rule. This
is a representation failure plus an independently existing interpretive
limitation, not a new fabricated count. Keep the frozen primary score and the
42-query sensitivity excluding only this affected question.

## Newly returned RYR2 18929323 and negative controls

`ryr189_aggregate` correctly returns 13 P2328S and 6 V4653F carriers from L21.
It does not count repeat Holter recordings at L23 as additional people. The
four phenotype abstentions remain conservative and defensible: L21 gives
pooled historical symptoms; L23 discusses treatment, repeat-test participation,
and two asymptomatic people, rather than a complete per-variant manifestation
partition. Treatment use and repeat-test membership must not silently stand in
for the queried phenotype. The 19 unaffected controls at L25 are a separate
control cohort, not unaffected carriers of these two variants.

`scn325_prior_counts` correctly abstains on all nine current-study clinical
questions. The response cites the Table S1 literature/gnomAD title, the three
target rows, and the L348–L349 footnote assigning the disease counts to prior
literature curation. Functional patch-clamp measurements do not convert those
historical people into newly enrolled clinical participants.

`scn251_variant_list` correctly abstains on all six questions. L1380/L1381
identify the two requested variants in the category list, but the packet
contains no per-variant number of people or person identifiers. Its nulls mean
the supplied evidence cannot establish the counts, not that these variants
have no carriers. The previously documented omitted list header does not
create a numeric patient count.

## Other returned claims and acceptance limits

The nine SCN5A 20031634 table answers preserve the correct genotype-positive
columns, including the separate undetermined counts. The four RYR2 30403697
carrier totals correctly use the two siblings plus one shared father, one
variant-specific parent for each inherited allele, and the explicitly
genotyped cousin. They count documented people within the packet, without
claiming an exhaustive pedigree census. The remediated R417L quote is now a
contiguous source substring and retains both sibling rows as sources.

This review does not promote these answers through the current production
validator. Its replay still accepts only the explicit C2277R carrier count 8,
which was already retained in the saved Grok 4.3 baseline. See
[baseline_comparison.md](baseline_comparison.md) for the separate marginal-value
audit. Agreement on investigator-selected questions cannot by itself establish
new regular-protocol yield, safe automatic endpoint selection, or performance
on unseen papers.
