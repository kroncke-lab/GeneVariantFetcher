These results justify **keeping Grok 4.3 as the primary full-paper reader**, **not promoting the tested Astra-medium full-paper config**, and **not treating the literal overlay as an A/U upgrade**. They also justify a **bounded, source-linked roster/transcription path** as the next experiment—not as a production swap. Grok 4.6 is **not ranked**; corrected Azure probes failed, earlier SDK pilots silently dropped effort and stay unscored, and CLI success is a different service path.

## What is not a model ranking

The scored comparison is not model-weights-only. Control used Grok 4.3 with shared Sol-medium verification/adjudication; Astra-medium primary changed model, effort, and output cap (32k vs 15k). Inputs were frozen source only. Two papers (PMID 20129283, 30059973) used deterministic shortcuts: 602/788 asserted identities and 521 control TPs. Their identical TP/FP/FN does not show model equivalence. The post-hoc non-shortcut slice (178/26/8 → 87/11/99; recall 95.70% → 46.77%) is descriptive, not a replacement gate.

Overall identity (699/36/89 vs 608/21/180) is dominated by a few papers. Astra’s extra 91 FN are mostly empty/length-limit failures (32533946: 83 FN; 30403697: 8 FN; 25814417: 1 FN), partly offset where Astra recovered more identities. Astra’s higher identity precision (96.66% vs 95.10%) is partly a failure looking like precision: 32533946 produced 20 control FPs and 0 Astra predictions. Do not read empty output as conservative accuracy.

The overlay’s identity scores copy the baseline because the reader **cannot add identities or overwrite non-null counts**. The one accepted fill (RYR2 G1886S, 30403697, 2 carriers) is a null-field, validator-passing current-study literal. Zero accepted A/U fills is the contract working, not a measured “Astra cannot read phenotypes.”

## Denominator and causal errors to refuse

Exact/supplied is among **emitted values on matched rows**; abstention is not a correct zero. End-to-end error treats misses and abstentions as zero without rewriting stored nulls. That is internally consistent, but **conditional MAE and supplied-exact fractions are success-conditional**. Astra’s 37/37 exact affected values are 37 of 650 positive affected fields, not complete recovery. “Adds 15 exact affected, loses two exact unaffected, loses 91 identities, carrier error +51.9%” is a joint operational-plus-model outcome, not an Astra-medium effect.

Identity extras with counts (13/3/1 vs 15/1/0) sit **outside** the matched-value denominator and are not proven hallucinations. Rejection tallies mix units (54 identity-level vs field-level reasons) and must not be summed into one field-error rate. Provenance-bucket MAE among supplied matches is not end-to-end MAE.

Raw lane: 75 unvalidated proposals, 71 numerical gold matches, 5.9% A/U error drop (1,132 → 1,065) is an **opened diagnostic**, not an accepted improvement or forecast. Gold disagreement does not by itself prove a proposed count wrong; gold agreement does not waive derivation, join, or current-study binding. Sign of count-error on the preflagged-exclusion slice (combined A/U 628 → 585) depends on influential/disputed rows. Keep the full panel and the sensitivity; neither is a population forecast.

The compact Astra-low probe changed **task, schema, source scope, and effort together**. Transcribing 48 people with matching layout fields is not clinical scalar-count accuracy. No Dx ≠ ECG-normal ≠ asymptomatic; two blank family cells need prose; duplicate tables are one cohort.

## 75%-of-human ceiling

This does **not** measure a 75%-of-human ceiling, a corpus rate, or a promotion gate. The 12 papers were previously opened and failure-enriched. There is no random sample, blind holdout, or blinded human reread. Missing Table 1 (18929323), unfinished 32k responses, endpoint ambiguity, identity retention, and the literal validator all cap the workflow before a human-relative ceiling could be estimated.

## Next test to isolate

One factor at a time, same frozen sources, arms locked before scores:

1. Bounded person-record transcription (fixed schema, one table copy, captions/headers/footnotes/nearby scope).
2. Source completeness/retrieval, separate from reader strength.
3. Code-side dedup and aggregation under an explicit endpoint contract; do not weaken the literal validator to accept model sums.
4. Empty-output / retry / table-hint fallback, so a failed JSON body does not discard already-filtered table candidates.

Do not re-confound model, effort, cap, and prompt. Do not score Grok 4.6 until health probes return usable completions with effort actually on the wire.

## Budget and operational caveats

Mostly honest. $38.58 returned proxy + $81.46 unknown/retry reserves + $25 margin = $145.04 against ~$150, covering **all API tests including failures and pilots**; CLI excluded; no Anthropic extraction. Unknown usage was not zeroed. Returned prices are proxies, not reconciled Azure charges; cache, invoice, and unreturned SDK attempts remain open. Overlay “$6.41 incremental” is reader returned-proxy, not full-stack cost. Several arms show $0 returned with large reserves; that is correctly conservative if the envelope is not sold as an invoice cap. Claims that effort is now serialized, prediction arrays were unchanged by post-extraction accounting, and offline tests passed are **documented here, not re-verified in this review**.
