# Semantic review of Astra's nine unique numeric responses

Reviewed after all 24 initial outputs were locked, using the frozen packets,
responses, `source_scores.json`, and the earlier source/reference adjudication.
No API calls were made and no primary artifacts were changed. This review
assesses source meaning beyond agreement with the numeric reference and the
mechanical quote check.

**Eight of the nine claims have strong source support for their precisely
specified endpoint. One is interpretive and is not a secure independent gain.**
All nine arose where both experimental Grok 4.6 attempts failed to return a
claim; they therefore demonstrate observed response availability on these
requests, not that Grok returned the wrong clinical interpretation.

## RYR2 25435091: three supported claims

| Claim | Assessment |
| --- | --- |
| Carriers 8 | Strong, explicit. S11 directly names C2277R and eight carriers; the short quote preserves the variant and count role. |
| Initial diagnostic affected 7 | Strong, explicit across the cited source. S22 states seven of eight carrier subjects have the EST CPVT phenotype; S11 supplies the named variant and S24 anchors initial assessment. The count is not inferred from symptoms or follow-up response. |
| Initial diagnostic unaffected 1 | Strong, explicit across the cited source. S21 identifies II:9 as not meeting diagnostic criteria, S5 identifies that person in the carrier table, S6 explains the non-diagnostic extrasystoles, S17 defines the phenotype, and S20 establishes the mutation. Prior palpitations do not invalidate diagnostic-phenotype negativity. |

All three quotes are verbatim contiguous excerpts in the cited material. There
is no missing bibliography or substantive citation needed for these three
packet claims. The latter two short quotes omit the literal C2277R string, but
their cited supporting units supply it. The current literal validator rejects
them because it requires local quote-to-variant binding. That is a difference
between the supplied multi-unit evidence and the validator's accepted evidence
form, not a false clinical count.

The original production comparison already contains these counts. In
`benchmarks/codex_paper_eval/runs/20260906_model12_grok43/predictions.json`,
the `papers` entry for RYR2 25435091 has `p.Cys2277Arg` with carriers 8,
affected 7, and unaffected 1 in its primary `variants` array, with `llm_table`
source evidence. Consequently, Astra's sole accepted literal carrier fill in
the new experiment is **not a new count or a null fill relative to that retained
Grok 4.3 production baseline**. It is an experimental fallback success relative
to two failed Grok 4.6 requests.

## RYR2 19398417: two supported claims and one uncertain inference

Carriers 4 and symptom-negative 2 are strong explicit claims. Both cite A26–A27,
which name W4645R and the four relatives, including two without symptoms.
The newline inside each quote preserves the original line wrapping. The output
correctly distinguishes symptom silence from negative diagnostic testing. The
validator's requested-role rejection does not erase the human-readable
variant/count/symptom relationship.

The affected 2 claim is **semantically uncertain**, even though it matches the
frozen reference. Astra's stated reason is four carriers minus two without
symptoms, supported by a general statement about symptom emergence across
generations. `Including two without symptoms` does not by itself prove that
every remaining carrier was assessed as symptomatic. The reason and citations
omit A23, the separate statement that two evaluated relatives were diagnosed
with CPVT. That omission matters to the strongest available cross-sentence
interpretation, rather than being decorative citation formatting. Even A23
would still require disease-segregation linkage and care about the symptom
endpoint.

The answer acknowledges the lack of individual symptom assessments and
timepoints but still emits the complement. Preserve its numeric primary score
and the previously specified 42-query sensitivity excluding only
`ryr193_abstract/q1_affected`. Do not count this claim among the eight strongly
source-supported unique values, and do not use it to justify automatic protocol
inclusion. No bibliography lookup can replace the absent individual linkage in
this abstract-only packet.

## RYR2 25814417: all three claims supported, with a scoped symptom endpoint

Carriers 179 is explicit at L33. L51 links the living cohort to Supplementary
Table 3, while L225/L229 identify the table population and columns. The quoted
screening sentence distinguishes 1,404 screened people from 179 living carriers.
Later or deceased cohorts are correctly excluded.

The symptom-positive 45 claim is supported beyond the scalar total:

- The cited patient rows are exactly the 45 rows with named previous symptoms;
  no positive row is omitted and no negative or missing row is included.
- Expanding the patient-ID ranges in Astra's reason gives exactly 45 unique
  IDs, identical to those source rows.
- The caption, variant/cohort statement, and `Previous symptoms` header are
  cited. The example quote preserves the source's `Syncope/Dizzeness` spelling.

The symptom-negative 133 claim cites all 179 patient rows, plus cohort and
header evidence. Directly recounting those cited rows gives 133 explicit `No`
entries, 45 named-symptom entries, and the missing `-` for patient 37 at L266.
Although its reason also displays subtraction, the complete cited roster
independently proves the negative count. This is not an open-cohort complement
that silently turns missing assessments into negatives.

The answer also notices that its row-derived fraction differs from the prose
percentage and obeys the specifically requested column. It does not confuse
VA/CVA results with previous symptoms or symptom-negative with disease-negative.
These claims are therefore strong for this test's specified endpoint. They are
not new clinical-diagnosis totals, and the corresponding production endpoint
must not be replaced by 45/133. The existing patient-row implementation already
computes the symptom-positive tally within its separate combined-phenotype
audit; a stronger model is not required to perform that deterministic tally.

For claim-level downstream use, retain the complete roster audit, not merely
the single example quote. The positive claim cites all positive rows; the
negative claim supplies the complete roster needed to independently verify
exhaustiveness across the returned evidence set.

## Grok mechanical quote failure: RYR2 30403697, not SCN5A 20031634

The sole mechanically invalid quote among Grok repeat's numeric-exact claims is
`ryr304_relatives/q1_carriers`, whose value is 3. It concatenates truncated
excerpts from the subject-1 and subject-2 table rows into a single quote, omitting
the remaining columns between them. The combined string is not a contiguous
verbatim source span. Its cited L105/L107/L108, shared-carrier-father reasoning,
and carrier count are substantively correct. This is a quotation-format failure,
not an invented relative or a duplicate-person error. It must remain failed in
the frozen mechanical score unless separately analyzed as a sensitivity.

All nine SCN5A 20031634 Grok-repeat quotes are verbatim and cite the relevant
headers and target row. There is no mechanical quote failure in that packet.
Its OCR mapping of `c.396312T.C` to c.3963+2T>C is openly acknowledged and remains
the previously disclosed identity-quality limitation. Current literal-validator
rejections of count roles or target binding are separate from mechanical quote
validity and from numeric correctness. The aggregate 12/13 mechanically valid
numeric claims combines the nine SCN5A values with four RYR2 relative counts;
the missing pass is the spliced RYR2 quote.

## Implication for inclusion

This semantic review strengthens the claim that Astra returned useful,
reviewable evidence during Grok 4.6 failures. It does not establish nine new
accepted production gains: one inference remains uncertain, only one claim
passes this experiment's current literal-validator lane, and that one already
exists in the retained production baseline. A paid reliability fallback and
routine improvement over the existing protocol are different decisions.
