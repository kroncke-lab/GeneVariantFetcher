# Frozen methods for results consultation

This summary was written from the frozen plan, source packets, investigator
reference review, and scheduling amendment. No model response bodies or scores
were read while preparing it. It is a methods brief, not a result or decision.

## Question and design

The operational question is whether Astra adds enough correct accepted counts
beyond a cheaper reader to justify regular use. The frozen experiment is an
**opened-source component diagnostic**, not a held-out or whole-paper benchmark.
It contains eight distinct papers and 43 requested count fields: 24 numeric
source references and 19 null references. A null means unavailable or
unsupported, never zero. Multiple count fields from one paper are correlated;
43 fields are not 43 independent paper observations.

Every packet has three initial independent calls: Astra low, Grok 4.6 low first
pass, and Grok 4.6 low repeat. Thus the fixed panel contains 24 initial calls.
The Grok repeat receives the same source, prompt and strict JSON schema as the
first Grok call. It does not receive the first answer, Astra's answer, reference
values or a corrective hint. It is a fixed inexpensive repeat comparator,
not a gold-triggered retry or an answer-adjudicating workflow. The two Grok
labels are experimental labels; the repeat may be dispatched before the first
pass under the seeded within-packet order.

All calls use low reasoning effort, a 4,096-token completion cap, strict JSON
schema, a 180-second deadline and no automatic retries. Each requested claim
contains its ID, value or null, explicit/derived/unknown basis, exact allowed
source IDs, quotation and concise reason. Quote and provenance validity must be
scored separately from numeric agreement. Both models receive the same generic
extraction rules; the selected target variants, population and clinical
endpoint are provided. Consequently this does not test autonomous discovery of
every variant, retrieval of missing source, or automatic endpoint selection.

The frozen plan requires completion of the fixed panel before outcome
inspection, subject to the conservative $2.17 envelope. All failures and
undispatched cells remain reportable. After lock and scoring, any remaining
allowance may support same-prompt replication of Astra-only wins. Such
replications must be labeled separately from the initial panel.

## Task coverage

| Paper and packet | Queries | What the task tests |
| --- | ---: | --- |
| SCN5A 20031634, `scn200_table` | 9 | Variant-specific carrier/BrS-ECG columns; separates explicitly negative from undetermined phenotype and excludes mutation-negative relatives. |
| RYR2 30403697, `ryr304_relatives` | 4 | Carrier counts include documented relatives outside numbered subjects, shared parents once, and variant-specific parental inheritance. No phenotype split is queried. |
| RYR2 25435091, `ryr254_timepoint` | 3 | Initial CPVT diagnostic testing, distinct from prior symptoms, final treatment response and ungenotyped deaths. |
| RYR2 18929323, `ryr189_aggregate` | 6 | Explicit per-variant aggregate carriers, repeated Holter tests and unassigned phenotype subsets when a cited table is absent. |
| RYR2 19398417, `ryr193_abstract` | 3 | Abstract-only carrier/symptom claims, screened relatives versus carriers, and cross-sentence linkage. |
| SCN5A 32533946, `scn325_prior_counts` | 9 | Abstention on clinical counts originating in literature/gnomAD within a functional study. |
| RYR2 25814417, `ryr258_roster` | 3 | A complete 179-person living baseline table; previous symptoms, explicit No and missing cells remain distinct from VA/CVA and historical deaths. |
| SCN5A 25163546, `scn251_variant_list` | 6 | Abstention when a genuine article-specific variant list lacks per-variant patient frequencies; a variant row is not one carrier. |

## Reference and scope safeguards

The source reference was frozen before live answers, and an independent
source-only review checked packet hashes, source lines, query correspondence
and the numeric/null references without inspecting model output. This is
investigator reference material, not API input. Agreement with it must remain
distinct from agreement with legacy gold and acceptance by current production
validation. A source-correct derived sum need not be admissible under the
unchanged literal-count path.

The independent reference review found one interpretive field:
`ryr193_abstract/q1_affected`. Its value links separate abstract statements;
the abstract does not print one explicit variant-specific symptomatic total.
Keep the frozen 43-query primary and separately show a 42-query sensitivity
excluding this field. An explained null on it is not a numeric false positive,
and this field alone cannot establish Astra's value.

Two representation limits are recorded before result interpretation. The
SCN5A 20031634 splice-variant target has malformed OCR in the packet, so refusal
to equate the string can be an identity/source-quality abstention. The SCN5A
25163546 introductory span misses the column-header line because embedded
Unicode line separators changed numbering; the target variant rows and table
title remain present, and none establishes a patient count. The frozen packet
and reference remain unchanged.

RYR2 25814417 explicitly asks for **previous symptoms** among living baseline
participants. That answer cannot be compared as an improvement or regression
against the existing combined symptoms/VA production endpoint or a broader
cohort that adds historical fatal cases. The relatives task is a count of
documented people, not proof of an exhaustive whole-paper pedigree census.

## Scheduling amendment

The original serial runner stopped between completed calls after three calls
had settled. The recorded amendment reports no pending reservation and no
canceled or retried API call. Transport latency, observed without inspecting
clinical answers, motivated resuming the remaining unchanged requests with
three worker processes. Prompt, model, cap, deadline, reference and scoring
were unchanged; per-call budget locking remained in force.

This changes scheduling, so concurrent service conditions and shared endpoint
limits are relevant to latency interpretation. It does not transform the
diagnostic into a clean latency benchmark. The amendment's plan hash is
`fb378be9745daa46085e2853c622125153f01e047da07ca94097b3b53f943131`;
the reference hash is
`8c41c3289d3d897703c9d134d9ee1fd7592f71c8975fc9ee0fb4c14571f8f3a8`.

The root experiment owns final scoring, production acceptance replay, cost
accounting, replication and the inclusion decision. No result is asserted here.
