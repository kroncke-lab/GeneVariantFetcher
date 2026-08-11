# Count-semantics adjudications — 2026-08-10

These decisions repair answer-key scope errors exposed by the fixed 56-paper
evaluation. They do not change locked predictions. Each decision was checked
against the local paper source, reviewed with `azure_ai/gpt-5.6-luna` at
`xhigh`, and then independently audited using blinded source cards by Grok 4.5,
AGY Gemini 3.1 Pro, and Claude Fable 5 at their strongest available effort.

## Prospective count policy

- Count the focal variant carriers in the current paper's study or analysis
  cohort, not recruitment totals, prior cohorts, unpublished observations,
  relatives without the variant, or allele copies.
- For treatment/effect papers, use the enrolled analysis cohort. For population
  genotyping papers, use all focal genotyped carriers.
- Record phenotype partitions only for the explicitly assessed subset. They may
  sum to less than the carrier total.
- If phenotype was not systematically ascertained, use null. Never infer zero or
  `carriers - affected` without direct support.
- Exclude duplicate rows for the same current cohort.

## Decisions

| Gene | PMID | Variant | Legacy | Adjudicated | Source-grounded rationale |
|---|---:|---|---|---|---|
| KCNH2 | 19160088 | R176W | 112/112/0 plus duplicate 16/0/16 | 16/null/null; duplicate excluded | The present Health 2000 survey reports 16 carriers. The 112 LQTS patients are prior/unpublished clinical observations. Symptoms were not systematically ascertained in the population cohort. |
| SCN5A | 28339995 | D1790G | 85/85/0 | 30/30/0 | The current flecainide study population consisted of 30 D1790G carriers described as LQT3 patients. The 85 genotype-confirmed people are the broader source-family pool. |
| KCNQ1 | 33141630 | T224M | 88/34/54 | 124/34/54 | The paper reports 124 genotyped carriers. Detailed phenotype follow-up covered 88, of whom 34 met the clinical LQTS definition and 54 did not; the other 36 remain unpartitioned. |
| SCN5A | 26746457 | R1193Q | 7/7/0 | 19/7/12 | The variant occurred in 19 participants and 7/19 had arrhythmia phenotypes. |
| SCN5A | 20470418 | S1103Y | 85/39/46 | 26/17/9 | The existing v2 adjudication represents the paper-defined S1103Y case/control cohort and is now honored by the scorer. |
| KCNH2 | 10841244 | L552S | 44/12/25 | 42/12/25 | The paper reports 40 heterozygous people plus two homozygous people. The legacy 44 is an allele-copy total. Phenotype was reported for a subset. |
| SCN5A | 18451998 | E1784K | 50/47/3 | 41/38/3 | The current multi-center cohort contains 41 carriers; 38 had long QTc and three had normal QTc. |
| RYR2 | 25814417 | G357S | 179/73/106 | 185/73/106 | Six genotyped sudden-death cases plus 179 living mutation-positive relatives yield 185 carriers. The phenotype partition describes the living cohort. |
| RYR2 | 33606749 | R420Q | 1/1/0 | 2/2/0 | The genotyped symptomatic proband inherited the variant from his symptomatic mother. |
| RYR2 | 33606749 | T1223A | 1/1/0 | 2/1/1 | The genotyped symptomatic proband inherited the variant from his asymptomatic mother. |
| RYR2 | 33606749 | T2390I | 1/1/0 | 2/1/1 | The genotyped symptomatic proband inherited the variant from his asymptomatic father. |

## Audit boundaries

The blind audit used 15 randomized evidence cards and withheld both locked
predictions and existing gold values. All three reviewers independently agreed
on the carrier totals for the five headline rows and the six newly corrected
controls. Two further controls lacked enough evidence in the generated cards
and were left unchanged. Two controls (SCN5A PMID 26921764 R1193Q and M369K)
unanimously retained gold=3 even though the locked prediction was 2, providing a
direct check against prediction-following adjudication.

Reviewer disagreement was confined mainly to phenotype axes whose definitions
or denominators differed. The authoritative partitions above therefore follow
the paper's explicit phenotype definition, not a model-completed arithmetic
partition. Full reviewer outcomes are in
`benchmarks/count_semantics_eval/MULTIMODEL_REVIEW_20260810.md`.

## Framework findings

- Broad missing-count recovery was rejected after 0 grounded additions in the
  first 162 completed gap checks across two KCNH2 papers. False-positive variant
  gaps made that route both low-yield and costly.
- GPT-5.6 xhigh verification requires a 64k output budget because hidden
  reasoning consumes the same budget before JSON is emitted. The larger budget
  is restricted to GPT-5.6; other xhigh models retain their normal cap.
- Three-letter prediction notation must route to one-letter source notation,
  including prefixless and stop-codon aliases. Long converted paragraphs must
  be excerpted around the variant mention rather than truncated from character
  zero.
- Gold v2 status is an exact, closed vocabulary consumed through one shared
  fail-closed helper by both scorers and the claim-verification pilot.
- Arithmetic completion is unsafe for nested cohorts. Unsupported phenotype
  partitions remain null.

Runtime LLM traces are retained under
`validation_runs/20260810_failure_routing_a1_56/luna_count_semantics_shadow_traces*`.
