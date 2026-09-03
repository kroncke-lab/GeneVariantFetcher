# Mixed-gold harness: Grok 4.6 `xhigh` design review

Date: 2026-09-03

This is a design review, not a new extraction or score. The local Grok CLI was
run as `grok-4.6` with `xhigh` reasoning, web disabled, no subagents, and no
tools. It received a compact evidence packet describing the lane, inventory,
tranche, locking, and cost design. It did not inspect repository code. No paid
GVF extraction was run; Grok consultation billing was not exposed and is not
asserted to be zero.

## Findings and dispositions

| Review concern | Disposition |
|---|---|
| ClinVar/PubTator must not count as a variant found in a paper | Accepted. `papers[].variants` is paper-origin only; linkage identities are hash-locked separately and appear only in the secondary `linkage_assisted` diagnostic. |
| `mixed` is too ambiguous to call paper-derived | Accepted. `mixed`, `manual_or_legacy`, and unknown origins are retained as `unattributed_variants` and excluded from both scored lanes. |
| A source label must record evidence origin, not merely the module that emitted a row | Existing origin/witness separation is load-bearing: `source_layer` remains the stable first origin and `observed_source_layers` records later corroboration. The projector classifies the stable origin and tests composite legacy forms. |
| Gold-derived identity resources are another leakage path | Accepted. `--gold-free-run` now disables every file-backed alias map even when its cache was warmed earlier in the process, records the state in `RUN_STATUS.json`, and paper-primary projection refuses a run without that proof. The public hard-coded ClinVar spelling map remains normalization-only. |
| Unknown origins must not silently disappear | Accepted. They remain inspectable in the locked unattributed lane, with a top-level held-row count. |
| Reused extraction databases/predictions contaminate a confirmation run | Accepted. The registered design permits frozen source-byte reuse only. Generated commands create new run roots and never pass `--resume-dir`; prior predictions, extraction JSON, databases, and traces are forbidden. |
| One-arm results on different tranches are not a clean protocol comparison | Accepted. The registered comparison is a paired frozen-baseline/candidate run on the same PMID-clustered manifest. Cost metadata now includes both one-arm and paired estimates. |
| The endpoint and pooling rule were underspecified | Accepted. The registry names paper-derived micro variant-identity recall as primary, paper-derived precision as a non-inferiority guardrail, and requires per-gene/per-gold-provenance reporting without treating the pooled number as a scientific cross-provenance claim. |
| Cheap tranches invite holdout reuse and metric shopping | Accepted. Registry order is frozen. `setup_production_eval.py` refuses out-of-order tranches, a candidate before its baseline, and a repeated arm. Scoring appends the tier, arm, runtime digest, and locked prediction digests to `consumption_log.jsonl` before exposing the report. |
| Same PMID under two genes is not independent | Accepted. Assignment is article-atomic and tests enforce that a PMID appears in exactly one tranche. The registered cluster unit is PMID. |
| Incomplete BRCA2 and mixed-provenance override gold cannot support one pooled scientific headline | Accepted. Provenance is carried per attempt and scored separately. These strata remain operational diagnostics, not exhaustive paper-coverage claims. |
| The 111 locally unavailable sources are a selected hard tail | Accepted limitation. They remain in `inventory.tsv` with reasons and are not mislabeled as extraction failures. The runnable suite is explicitly 1,422 of 1,534 attempts. |
| Mean-token cost is not an invoice or tail bound | Accepted. Estimates exclude acquisition/review, use a dated price proxy, carry 25% headroom, and must be recalibrated against the first paid tranche's actual charge. Candidate changes can alter cost. |

Grok's narrative also printed a conflicting `$0.98–$1.30` per-tranche estimate
despite the supplied `$4.98–$5.30` figures. That arithmetic was not adopted.
The registry recomputes costs from the checked source manifests and exact model
token totals: `$4.98–$5.30` per arm, `$9.96–$10.61` paired, and
`$12.45–$13.26` paired with headroom.

## Result

The original design was a no-go for paid confirmation because ambiguous origin,
runtime alias leakage, unpaired comparisons, endpoint choice, cache policy, and
tranche burn were not mechanically closed. The implementation now closes those
items. It is suitable for the first registered paired measurement, with two
limits that remain explicit: unavailable-source selection and heterogeneous,
partly non-exhaustive gold provenance. The first paid pair must validate the
cost proxy; it must not be promoted as a cross-provenance scientific headline.
