# Failure-routing A1 — locked 56-paper replay (2026-08-10)

## Decision

**Implemented behind `ENABLE_TIER3_REASON_CLASS_ROUTING`, default OFF.** A1
materially reduced the intended Sol verifier lane and its changed-decision slice
did not lose a variant, but the full stochastic replay failed the predeclared
aggregate carrier-MAE gate. Keep the flag off in production until a replicate or
paired same-primary-output ablation resolves that failure.

A1 changes only the condition that opens claim verification. The existing risk
score, threshold, verifier model, evidence cards, ranking, 20-card cap, prompts,
trust logic, recovery layers, and primary extraction route are unchanged.

- Claim verification signals: denominator ambiguity, repeated suspicious count
  tuples, arithmetic mismatch, count-bearing regex/figure/mixed rows,
  non-per-variant provenance, unknown large-count provenance, carrier sum above
  the paper census, and study-wide-count suppression.
- Completeness-only signals route to `completeness_rescue`; missing-count signals
  route to `count_recovery`; missing-table/source blockers route to
  `source_recovery`. A1 records these routes but does not implement those future
  recovery stages.

## Run integrity

The replay covered 56 gene-paper entries: the same locked 48 cardiac papers as
the 2026-08-08 baseline plus all eight BRCA2 curated diagnostic entries. All five
gene jobs completed with zero extraction failures, stage failures, or warnings.
The jobs ran in parallel and the slowest finished in 17.4 minutes.

- Source lock: `validation_runs/20260810_failure_routing_a1_56/source_lock.tsv`
- Source-lock SHA-256: `7a88befc07e9f25d0fdabdc45566d3ba9c67d9a3ad56c83aa2c3eb7dff03b6a1`
- Every final `*_FULL_CONTEXT.md` re-hashed to the lock after completion.
- Cardiac inputs match the prior selection's source hashes exactly.
- `GVF_DISABLE_LOCAL_DATA=1` prevented implicit corpus, sibling-DB, or
  VariantFeatures discovery during extraction and final scoring.

The first local scoring attempt was stopped before it emitted predictions or a
lock because the projection command lacked `GVF_DISABLE_LOCAL_DATA=1` and began
reading the mounted VariantFeatures database. The final `score_eval.sh` applies
the guard to projection, locking, and scoring; only those guarded outputs below
are authoritative.

## Fixed-48 cost

| Route | Baseline calls / tokens | A1 calls / tokens | Token change |
|---|---:|---:|---:|
| Kimi table routing | 13 / 46,379 | 13 / 45,884 | -495 |
| Grok paper extraction | 43 / 843,407 | 49 / 850,549 | +7,142 |
| Sol claim verification | 146 / 394,420 | 106 / 285,472 | **-108,948 (-27.6%)** |
| Sol figure reading | 76 / 86,467 | 76 / 88,626 | +2,159 |
| Sol paper adjudication | 1 / 2,169 | 0 / 0 | -2,169 |
| **Total** | **279 / 1,372,842** | **244 / 1,270,531** | **-102,311 (-7.5%)** |

Summed provider-call time fell from 3,670.5s to 3,522.4s (-4.0%). The 27.6%
verifier-token reduction clears the predeclared >=25% cost gate. A1 routed 17 of
25 above-threshold cardiac papers to verification and skipped eight whose risks
were completeness/missing-count only. Those eight consumed no verifier calls;
the other 17 produced 106 evidence-card calls.

BRCA2 added 59 calls / 267,609 tokens (6 Grok extraction, 3 Kimi routing, 36 Sol
verification, and 14 Sol figure calls). None of its three above-threshold papers
was completeness-only, so A1 skipped no BRCA2 verifier work. The complete
56-paper replay used 303 calls / 1,538,140 tokens.

## Fixed-48 quality

| Projection | Baseline TP/FP/FN | A1 TP/FP/FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| all layers | 831 / 1060 / 170 | 829 / 1065 / 172 | 43.95% -> 43.77% | 83.02% -> 82.82% | 57.47% -> 57.27% |
| paper-derived only | 710 / 592 / 291 | 708 / 598 / 293 | 54.53% -> 54.21% | 70.93% -> 70.73% | 61.66% -> 61.38% |

The aggregate variant gates passed: precision -0.18pp, recall -0.20pp, all-layer
F1 -0.20pp, and paper-only F1 -0.28pp. The paper-only KCNH2 F1 moved -0.82pp;
the other per-gene paper-only F1 changes were +0.07pp (KCNQ1), 0 (RYR2), and
-0.10pp (SCN5A).

Count coverage rose, but the raw carrier-MAE gate failed:

| Projection / field | Baseline coverage, MAE | A1 coverage, MAE | MAE change |
|---|---:|---:|---:|
| all carriers | 307/1001, 0.723 | 326/1001, 0.902 | **+0.179 (FAIL)** |
| all affected | 260/1001, 0.673 | 267/1001, 0.708 | +0.035 |
| all unaffected | 225/1001, 0.373 | 246/1001, 0.159 | -0.215 |
| paper-only carriers | 306/1001, 0.686 | 325/1001, 0.868 | +0.181 |

The carrier regression is not concentrated in papers A1 rerouted. On the eight
changed-decision cardiac papers, baseline and A1 have exactly the same 58 TP / 203
FP / 30 FN. Carrier coverage increased 7 -> 10 while total absolute carrier
error improved 43 -> 42 (MAE 6.14 -> 4.20). The non-rerouted 40-paper slice
accounts for all aggregate degradation. In particular, SCN5A PMID 20470418
newly extracted 26 p.Ser1103Tyr carriers versus gold 85 (59 absolute error); it
still triggered semantic verification, so A1 did not bypass its verifier. RYR2
PMID 25814417 is route-sensitive: A1 retained the source-supported total of 185
(179 living carriers + six SCD) while gold records 179; the old verifier removed
the count entirely because its compact evidence packet omitted the count line.
That improves coverage but adds six matched-row error under the current gold
semantics.

This causal slice is favorable to A1, but it does not override the predeclared
raw aggregate gate. The rollout therefore remains default off.

## BRCA2 diagnostic arm

The eight BRCA2 answer keys are curator/derived `gold_overrides`, **not a fully
manual gold standard** and not eligible for headline program metrics. Against
that diagnostic reference, A1's fresh run scored:

- PMID recall: 8/8 (100%)
- variant rows: 99/110 (90.0%); unique variants: 93/103 (90.3%)
- patients: 922/1088 (84.7%); affected: 652/809 (80.6%); unaffected: 188/188
  (100%)
- matched carrier MAE: 0.192 (10 absolute error over 52 matched values)
- matched affected MAE: 0.0 (15 matched values)

The older canonical BRCA2 DB scored 97/110 rows, 91/103 unique variants,
carrier MAE 4.04, and affected MAE 1.97 on the same diagnostic subset. This is
directional context, not a paired A0/A1 claim: the canonical DB differs in age,
primary extraction, and recovery history. The gold-PMID precision proxy moved
18.8% -> 19.3%; its counted-extra proxy moved 45.3% -> 36.9%, largely from
count-bearing regex-table extras and subject to the benchmark's known
non-exhaustive-gold limitation.

## Locked artifacts

- All-layer report:
  `validation_runs/20260810_failure_routing_a1_56/fixed48_all_layers/report.json`
- All-layer predictions SHA-256:
  `9c61cd0263131a822a1d616dfd01d10d47310d39012dbf459fb8c5afe85022a5`
- Paper-only report:
  `validation_runs/20260810_failure_routing_a1_56/fixed48_paper_only/report.json`
- Paper-only predictions SHA-256:
  `161c5e5d7a52937ad601343ea0a57a214ef03f20a6e9146c90995dd8456bf1f9`
- BRCA2 diagnostic report:
  `validation_runs/20260810_failure_routing_a1_56/brca2_diagnostic_score/summary.json`
- Exact route telemetry:
  `validation_runs/20260810_failure_routing_a1_56/route_telemetry.tsv`

Offline validation after adding the default-off rollout flag: focused routing,
settings, and provenance suites passed 43/43; the full unit suite passed 1,688
with one skipped.
