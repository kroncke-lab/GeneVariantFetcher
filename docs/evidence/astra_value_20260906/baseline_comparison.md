# Marginal benefit against the retained Grok 4.3 baseline

**The new test demonstrates no additional comparable, accepted count over the
saved regular-protocol baseline.** Of 24 positive frozen source questions, 17
already have the same stored baseline number: 16 are strongly comparable and
one is the previously flagged interpretive affected count. Four carrier counts
are missing in that baseline, but both experimental Grok 4.6 arms and Astra
recover the same documented-relative candidates. The remaining three questions
use a different cohort or phenotype endpoint and cannot be scored as baseline
errors or Astra gains.

The machine-readable companion is `baseline_comparison.json`. It records each
question, unique strict allele match, original primary prediction, scope and
endpoint qualification, source evidence, initial experimental result, and input
hashes. The audit reads only
`benchmarks/codex_paper_eval/runs/20260906_model12_grok43/predictions.json`'s
primary `papers[].variants` array; it does not substitute linkage-assisted or
comparison-lane predictions. No paid calls or remediation outputs were used.

## All positive reference questions

Values below are carriers/affected/unaffected where three fields are shown.

| Frozen packet and targets | Saved Grok 4.3 baseline | Frozen reference | Comparability and marginal value |
| --- | --- | --- | --- |
| SCN5A 20031634: Gly1408Arg | 14/4/9 | 14/4/9 | All three already present, Table 1 family B, same phenotype columns. |
| SCN5A 20031634: c.3963+2T>C | 10/2/8 | 10/2/8 | All three already present, Table 1 family E. The OCR identity limitation remains, but the stored identity is already the normalized target. |
| SCN5A 20031634: Ala665GlyfsX16 | 9/6/2 | 9/6/2 | All three already present, Table 1 family I. Unknown phenotype is not counted as unaffected. |
| RYR2 30403697: R417L, R2028H, Y4721C, G4772S | All four carrier fields null | 3, 2, 2, 2 carriers | New source-supported documented-relative candidates, shared by both experimental Grok arms and Astra. The packet's explicitly documented set is not an independently proved exhaustive whole-paper count; derived-count acceptance also remains outstanding. |
| RYR2 25435091: C2277R | 8/7/1 | 8/7/1 | All three already present with initial EST evidence. Astra's accepted experimental carrier fill is therefore not new over the retained protocol. |
| RYR2 18929323: P2328S, V4653F | 13 and 6 carriers | 13 and 6 carriers | Both already present with the same current-study Holter cohort evidence. Astra's initial request failed on this packet. |
| RYR2 19398417: W4645R | 4/2/2 | 4/2/2 | Carriers 4 and symptom-negative 2 already present. Affected 2 has the same interpretive ambiguity in both outputs; it is not counted as a strong safe exact reference. |
| RYR2 25814417: G357S | 185/97/62 | 179/45/133 | Not comparable: baseline includes six genotyped SCD cases and a combined symptoms/VA phenotype; packet asks only living baseline people and the Previous symptoms column. |

The SCN5A 20031634 baseline export often says `no quote captured`, but it
explicitly records each Table 1 family location. The source table independently
supports those exact values and roles. This audit establishes that the counts
already exist; it does not claim their exported provenance is ideal or that
every stored count would pass the newly tested literal recovery validator.

## Astra's nine unique values relative to the failed Grok 4.6 attempts

Those nine divide as follows:

- RYR2 25435091: all three already exist in the retained baseline.
- RYR2 19398417: all three already exist; two are directly supported and one
  remains interpretive.
- RYR2 25814417: all three answer a differently scoped question.

Thus six of nine have the same baseline values, including five strongly
comparable values and one interpretive value. Three have an incompatible scope
or endpoint. None demonstrates a new accepted regular-protocol count. This
does not negate Astra's measured availability advantage against these Grok 4.6
requests; it limits what that advantage means for the decision to pay for Astra
in the existing protocol.

## Counts and limits

| Classification | Positive questions |
| --- | ---: |
| Same stored value, strong endpoint/scope comparability | 16 |
| Same stored value, interpretive reference | 1 |
| Baseline null, documented-subset carrier candidate shared by cheaper Grok | 4 |
| Different cohort or phenotype endpoint | 3 |
| Total | 24 |

This is a read-only audit of a saved production run, not a new same-day Grok 4.3
rerun or a population-level benchmark. The selected targets and detailed query
endpoints were supplied in advance. Do not report 17/24 as a regular-protocol
accuracy rate: the denominator includes four qualified subset candidates and
three incompatible questions. Retain the separate 42-query sensitivity for
RYR2 19398417's affected count and the original locked primary scores.
