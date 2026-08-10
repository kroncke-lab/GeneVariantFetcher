# Fixed-48 source-snapshot replay — 2026-08-08

This is the durable compact record for the fresh production replay of
`highcarrier_48_papers_20260723.tsv`: 12 papers each for KCNH2, KCNQ1, RYR2,
and SCN5A. BRCA2 is not part of this benchmark. The extraction used the exact
preserved source selection and `GVF_DISABLE_LOCAL_DATA=1`, preventing the run
from discovering a different local corpus or sibling database.

## Run identity and integrity

Code base at launch: `5bfb41cdc3d25db341f3170baab598401b0549b9` plus the
dirty-worktree changes captured by each `run_manifest.json`.

| Gene | Run directory | Status | Duration | Calls | Trace integrity |
|---|---|---|---:|---:|---|
| KCNH2 | `KCNH2/20260808_121502` | completed / ok | 965s | 66 | write-time verified |
| KCNQ1 | `KCNQ1/20260808_121506` | completed / ok | 812s | 66 | write-time verified |
| RYR2 | `RYR2/20260808_121509` | completed / ok | 1014s | 70 | write-time verified |
| SCN5A | `SCN5A/20260808_121513` | completed / ok | 1004s | 77 | write-time verified |

Local production root:
`validation_runs/20260808_fixed48_snapshot_replay_hermetic/`.
Every gene reported 12/12 full-text sources, zero stage failures, and zero stage
warnings. Raw call records are intact. The trace coverage audit has one known
linkage gap: `clinical_table_routing` decision records do not carry an
`accepted_response_trace_id`; this is an observability defect, not a missing
request or response.

## Locked score projections

| Projection | Lock directory | Predictions SHA-256 | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| all layers | `validation_runs/20260808_fixed48_snapshot_replay_eval` | `803786badfe0aa1be806926ce9cb41f76339b86ac5305f6ad531eb9396c1027f` | 831 | 1060 | 170 | 43.9% | 83.0% | 57.5% |
| paper-derived only | `validation_runs/20260808_fixed48_snapshot_replay_paper_only` | `bb546581d477e26906b4581f4a39bb4b966825ff68dc41b200b62194d0f61988` | 710 | 592 | 291 | 54.5% | 70.9% | 61.7% |

The all-layer selection SHA-256 is
`f2192713bd0ce20806ef2c98b906af735578a39125392adce6a6092dd456600f`;
the paper-only selection SHA-256 is
`62ddf729cf85fdd697f72b1623c593f81240fdcaa277366bdc7833bfe4c3c04c`.
The selections differ only because each projection has its own run identity
and lock material; both cover the same 48 gene/PMID sources.

Count coverage and matched-value error:

| Projection | Carriers | Affected | Unaffected |
|---|---|---|---|
| all layers | 307/1001 (30.7%), MAE 0.723 | 260/1001 (26.0%), MAE 0.673 | 225/1001 (22.5%), MAE 0.373 |
| paper-derived only | 306/1001 (30.6%), MAE 0.686 | 259/1001 (25.9%), MAE 0.676 | 224/1001 (22.4%), MAE 0.321 |

## Change from the 2026-07-26 production lock

The July all-layer lock was 789 TP / 985 FP / 212 FN (44.5% precision, 78.8%
recall, 56.9% F1). The fresh replay is +42 TP / +75 FP / -42 FN: +4.2 recall
points, -0.5 precision points, and +0.6 F1 points.

The July paper-only view was 670 / 588 / 331 (53.3%, 66.9%, 59.3%). The fresh
view is +40 TP / +4 FP / -40 FN: +4.0 recall points, +1.2 precision points, and
+2.4 F1 points.

Almost all of the gain is RYR2 PMID 19926015. The July markup circuit breaker
discarded the paper, producing 0 TP / 0 FP / 40 FN. The current circuit breaker
accepts it: paper-derived scoring is 37 / 8 / 3, while all-layer scoring is
38 / 78 / 2. The approximately 70 extra false positives in the latter are
ClinVar linkage rows and should be treated as a paper-identity-linkage precision
problem rather than a paper-reading regression.

Carrier MAE improved from 1.424 to 0.723, partly because two old large errors
disappeared. Affected MAE worsened from 0.137 to 0.673. The clearest actionable
case is KCNQ1 PMID 18713323: carrier-only table totals were also written as
affected counts without an explicit phenotype split. `TASKS.md` tracks the
deterministic parser change and trust-gate backstop; neither was changed before
this measurement.

## Exact model telemetry

| Route | Calls | Input tokens | Output tokens | Total tokens | Summed duration |
|---|---:|---:|---:|---:|---:|
| `Kimi-K2.6-1` table routing | 13 | 15,916 | 30,463 | 46,379 | 121.6s |
| `grok-4.3` paper extraction | 43 | 567,021 | 196,948 | 843,407 | 2,030.9s |
| `gpt-5.6-sol` claim verification | 146 | 325,510 | 68,910 | 394,420 | 1,029.5s |
| `gpt-5.6-sol` figure reading | 76 | 68,519 | 17,948 | 86,467 | 485.6s |
| `gpt-5.6-sol` paper adjudication | 1 | 2,006 | 163 | 2,169 | 2.9s |
| **Total** | **279** | **978,972** | **314,432** | **1,372,842** | **3,670.5s** |

Explicit fixed PMIDs bypass discovery and Tier 1/2 relevance, so this run does
not test Luna. Luna and Terra connectivity both passed on the same workstation,
but they remain opt-in. The first controlled Luna A/B should be Tier 2 relevance
or extraction-priority triage on a labeled discovery set; count attribution,
claim verification, figure interpretation, and adjudication require separate
quality evidence.

## Blinding boundary

The final projection contents are DB-derived, gold-PMID enrichment was not
enabled, and the external scorer only read gold after each `predictions.json`
was finalized. However, `gvf-run` automatically creates read-only per-layer
scorecards when registered gold is present. Those scorecards did not feed back
into extraction, but they can precede the external lock. This production
projection is therefore prediction-content-independent of gold, not a strict
demonstration that no process read gold before the lock. A future `--no-score`
or blind-recovery projection is required for that stronger claim.
