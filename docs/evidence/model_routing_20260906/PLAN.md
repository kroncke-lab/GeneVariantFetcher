# New Azure readers and selective clinical reading: opened calibration

Started 2026-09-06. The user authorizes up to approximately $150 in extraction
API spend for this test; CLI consultation costs are excluded. This is a new
envelope, separate from previous campaign ledgers. Azure is preferred; no
Anthropic extraction is planned. This document records the experiment, while
TASKS.md remains the active checklist.

The question is whether Grok 4.6 should replace Grok 4.3 as primary reader and
whether Astra should read independently on papers with unresolved clinical
counts, or should itself be the primary reader. New deployment access is not
evidence of an accuracy gain.

1. Repair model compatibility and run tiny deployment probes. Preserve older
   routes/defaults. Record actual request parameters, model identifiers and
   usage. Reserve $5 including failed-call uncertainty.
2. Freeze 12 already-opened gene-paper attempts: RYR2 18929323, SCN5A 30059973,
   MYBPC3 20433692, MYBPC3 21302287, RYR2 30403697, SCN5A 20031634,
   RYR2 19398417, RYR2 25435091, SCN5A 32533946, SCN5A 25163546,
   RYR2 25814417 and SCN5A 20129283. These include omissions, joins, exact
   controls, legitimately unavailable counts and endpoint/unit disputes.
   Bind actual previously acquired components and their hashes before calls.
   Do not tune prompts with gold values; exposure to the earlier paper audit
   makes this calibration, not a blinded population sample.
3. Run the current production extraction path with Grok 4.3, Grok 4.6 and
   Astra as the primary model, keeping source and other stages fixed. Record
   which papers actually reach the primary model; deterministic shortcuts
   are part of the current workflow, not evidence that models are equivalent.
   Grok 4.6 and Astra request high reasoning; the older deployment retains its
   current default. This compares model configurations, not weights alone.
   Initial reservation: $60.
4. Test a bounded independent Astra clinical reading on source-bearing papers
   where the Grok 4.6 result is incomplete. The independent reading does not
   receive Grok's proposed answers or the gold. Preserve raw proposals and
   separate source-validated accepted fields from rejected/unresolved fields.
   Do not accept a count because two models agree. Record the trigger and
   compare with the ordinary Grok result. Initial reservation: $35.
5. Repeat the useful contrasts on the same fixed panel, or extend to another
   opened panel if failure mechanisms justify it. Reserve $35 for this and
   $15 for retries/uncertainty. Do not spend just to exhaust the envelope.

The roster and initial design are fixed before new paper outputs. Any revision
must be dated and state which outputs were already inspected. No original
source, gold, lock or scored output may be overwritten. All arms are locked
before scoring. Use the established paper-derived trusted projection and
scoring/figure tools; report raw-reader output as a separate diagnostic.

Report identity TP/FP/FN, positive-reference count coverage, exact supplied
counts, count error including omissions, wrong supplied counts, per-paper
changes, endpoint disputes, API usage and cost per additional correct count.
Report the SCN5A 20129283 influence separately without replacing official
scores. Preserve cardiac/manual versus MYBPC3 reference provenance. A useful
small-panel result does not promote defaults, reset the existing gates or
authorize an unearned confirmation claim on tranche 04.

Pricing uses explicit, dated token-rate proxies until provider charges can be
reconciled; unknown or failed-call usage is never silently zero. Reserve a
conservative cost before dispatch and record all returned usage, including
reasoning, retries and operational failures. CLI advice is adversarial input,
not independent human validation.

Pre-paper amendment (2026-09-06, after tiny deployment probes only): Astra's
Chat Completions endpoint requires `max_completion_tokens`; compatibility now
translates GVF's legacy cap. Grok 4.6's first two Azure probes timed out; run
baseline and Astra while checking that endpoint, and omit a failed deployment
from accuracy comparisons. If Grok remains unavailable, use the current Grok
4.3 baseline to trigger the independent clinical reader. This is an operational
fallback and cannot establish a Grok 4.6 accuracy result.

Pin non-primary adjudication and verification to Sol `medium` in every arm to
prevent inherited primary reasoning from confounding them; the previous default
was unspecified, so the fresh Grok 4.3 arm is a controlled comparator, not a
byte-identical replay of historical defaults. Tier 2 Luna retains `xhigh`.
Request 32,000 primary output tokens; current per-model clamping gives Grok
15,000 and Astra 32,000. Record actual request caps and truncations. Vision stays
Sol with its existing provider-default effort. Count recovery stays off.

The additive clinical-reader experiment only fills missing count fields on
identities already retained by the baseline's trusted paper projection.
Independently proposed values for already-filled fields are retained as
contradictions for audit, never silently overwritten. The reader receives the
frozen paper text and target gene only, no earlier answer, variant list, or gold.
Require a verbatim source quote and variant/cohort/person-unit binding; preserve
rejected proposals separately. A table cell spanning variants, allele-frequency
number, assay count, or ambiguous cohort total cannot become a patient count.

Operational update before Grok 4.6 paper outputs: all three configurations passed
exact shared-client probes, including Grok 4.6 high reasoning with JSON output.
The original three-arm design proceeds. Early timeouts remain operational
failures in the smoke receipts, not paper-level false negatives.

Accounting update during the first baseline/Astra calls: Azure Grok Chat
Completions reports reasoning separately from `completion_tokens` (unlike the
Astra response). The frozen hook's returned-usage proxy is conservatively
reconciled once per second using `max(completion_tokens, total_tokens -
prompt_tokens)` for Grok, retaining all raw counters. A further $20 uncertainty
reserve stays inside the $150 envelope during these runs. The normalizer is
`reconcile_live_budget.py`; neither model input nor extraction runtime changes.

Independent-reader trigger and raw diagnostic frozen while production runs were
in progress, before scores or trusted paper projections were inspected: trigger
when the chosen baseline retains at least one identity with a null count and its
frozen source contains a clinical person/count signal (`carrier`, `patient`,
`proband`, `symptomatic`, `asymptomatic`, `affected`, `unaffected`). This permissive
trigger deliberately permits negative controls; measure its firing rate and cost.
No retained identities means no additive recovery. The independent reader can
propose explicit counts or individually enumerated patient-row derivations.
Only the existing literal-evidence validator's accepted, current-study, explicit
values enter the primary additive result. Raw proposed fills are a separately
labelled unvalidated diagnostic lane. A disagreement or ambiguous identity match
is never resolved using gold. No validator relaxation is allowed in this test.

Pre-score expansion (2026-09-06): the complete Grok 4.3 comparator is now locked;
several high-reasoning primary calls on the new resource remain pending after
10 minutes, while tiny health probes still return successfully. No gold scores
have been inspected. Test the **same predeclared clinical-reader policy on the
locked Grok 4.3 baseline as well as Grok 4.6**. This adds a useful current-reader
workflow comparison and starts independent clinical work while primary runs
finish. It does not replace the three primary arms or change any trigger/merge
rule. The two independent readings, where both parents trigger, also provide
an opened same-source repeat diagnostic. Keep both fresh call sets and actual
routing costs separate; the $150 total envelope is unchanged.

Operational revision, before any gold scores (2026-09-06): Astra high used all
32,000 completion tokens as reasoning with **no visible extraction** on
MYBPC3 20433692 and RYR2 25814417. The production empty-content retry remains
part of that arm. Grok 4.6 high also reached its 15,000 visible-output cap on
SCN5A 20031634 and MYBPC3 20433692, with partial text. These are observed
request-configuration failures, not evidence of a human-performance ceiling.
Add an Astra **medium**, otherwise identical 12-paper primary configuration
with the same 32,000 cap. Launch only when the ledger has room for its bounded
reservations. Retain the original high arm and all failed/partial attempts.
This revision responds to visible completion status and usage, not gold scores;
it does not change sources, prompts, downstream models or acceptance rules.
The original launcher is preserved as `campaign_primary_v1.py`.

Reservation refinement during execution, without changing scientific calls:
for pending/unknown Astra calls whose input **byte upper bound** is below
272,000, the live accountant removes the unnecessary 2x long-context premium
from the original reservation. It still reserves three complete SDK attempts,
including each attempt's full 32,000-token reasoning/output cap. Original
reservations are retained as `original_reserved_usd`. This prevents overly
conservative concurrent reservations from blocking work while actual spend is
low. The $150 ceiling and extra $20 uncertainty margin remain unchanged.

## Superseding transport correction before scoring

An offline inspection of LiteLLM's actual parameter transformation proved that
`reasoning_effort` was omitted for both new model names under `drop_params=True`.
The earlier labels were requested efforts, **not transmitted efforts**. Sol's
shared-stage `medium` and Grok 4.3's default were unaffected. The remaining
new-model integration pilots (including the mistaken medium pilot and partial
clinical run) were stopped, preserved, and explicitly **abandoned unscored**;
they cannot support a high-versus-medium or accuracy comparison. In-flight
usage remains reserved. See `integration_pilot_abandonment.json` and
`reasoning_parameter_transport_audit.json`.

The corrected comparisons are:

| Arm | Primary | Transmitted effort | Primary output request cap |
| --- | --- | --- | --- |
| Completed, locked control `grok43` | Grok 4.3 | provider default | 15,000 |
| Fresh `grok46_verified` | Grok 4.6 | high | 32,000 |
| Fresh `astra_medium_verified` | Astra | medium | 32,000 |
| Independent reader on `grok46_verified` | Astra | medium | 32,000 |

The SDK allow-list now explicitly preserves the supported effort parameter.
Mock-transport tests inspect the actual OpenAI SDK HTTP request for both new
models at medium/high; they verify effort, token-cap spelling and Astra's lack
of sampling parameters. Grok 4.6 gets a separate 32,000 local request cap after
observed truncation at the inherited 15,000 cap; this is a GVF policy, not a
claim about the provider's maximum. Both new configurations must also pass a
live 32,000-cap micro-probe before full-paper dispatch.

The old control remains valid because the compatibility changes affect only
new model names or malformed selected resource routes; its actual model
parameters, shared stages and source inputs are unchanged. Runtime fingerprints
remain distinct and recorded. This is a practical configuration comparison,
not a weights-only contrast. The independent clinical prompt and acceptance
rules are unchanged; only effective effort is now medium. No gold scores have
been read and all final scored arms must lock first. The original $150 total
includes every pilot, cancellation and corrected call, with CLI costs excluded.

The v2 budget hook counts separately-reported Grok reasoning, waits for temporary
concurrent reservations to clear, and captures OpenAI SDK retry events per call.
A successful final response after a retry retains an additional unknown-usage
reservation for that earlier attempt. The original hook is preserved as
`budget_guard_pilot_v1.py`; the original clinical script is preserved as
`clinical_reader_pilot_v1.py`. The live reserve normalizer is stopped before v2
runs so it cannot halve a reservation that already omits the long-context premium.

Corrected-reader fallback before any scores: Grok 4.6 failed the verified live
32k probe and subsequent bounded SDK/direct-HTTP probes, including two effort
levels, both JSON/text modes, both cap spellings and omission of the stream
field. No corrected Grok paper batch is dispatched without a passing probe.
The corrected independent Astra-medium run therefore uses the completed Grok
4.3 control, under the originally allowed unavailable-deployment fallback. Its
fresh run name ends `_astra_clinical_verified`; the earlier partial default-
effort clinical pilot remains unscored. Grok 4.6 availability is not an accuracy
or ceiling result. A final bounded Responses/low-effort health check may still
permit the corrected Grok arm; no further full batch is started speculatively.

Before this corrected independent run, identity matching was made symmetric:
explicit cDNA/protein components in a proposed combined notation are matched
against the baseline's explicit components. No sequence or position-only alias
is inferred. Multiple matching baseline identities remain a rejection. Two
synthetic tests cover the unique combined-notation and ambiguous-compound cases.
No count acceptance/derivation gate is relaxed and no gold score was consulted.

Final availability decision before scores: the last 55-second direct HTTP checks
failed for both Responses at low effort and Chat Completions at provider-default
effort, each with a 32,000 cap and no retries. The corrected Grok 4.6 paper arm
remains prepared but undispatched. Complete and score the control, corrected
Astra primary, and corrected Astra clinical overlay on the control. This is an
operationally limited comparison, with no Grok 4.6 accuracy ranking.

Bounded source-roster diagnostic before scoring: corrected Astra medium again
exhausted 32k entirely as reasoning on MYBPC3 20433692; the unchanged production
retry is retained. The clinical medium arm has returned one 8-field raw proposal
set for MYBPC3 21302287 and no accepted fills. No gold score is open. After all
three CLIs independently reviewed MYBPC3 20433692's source, add one Astra-low
8,192-output-token diagnostic: transcribe the single deduplicated patient table
into compact person records with original line IDs. Do not emit/score scalar
count predictions or change an existing arm. This tests bounded transcription
operability, not model-only accuracy or a new clinical acceptance rule. Preserve
source hashes, exact prompt, returned usage and ambiguities. CLI suggestions
are not reference truth. This call remains inside the same $150 envelope.

Pre-score telemetry amendment: SCN5A 30059973's independent call timed out at
1,200 seconds with no usage returned. Its baseline predictions remain intact
and the full possible-call reservation is retained. The native lock contract
previously required exact counters even for such a traced failure, encouraging
failure exclusion or fabricated zeros. Add a narrow exception requiring null
counters, an explicit unknown-failed-call status, and an actual same-paper,
SHA-bound failed API trace with absent provider usage. Native write-time trace
manifest verification remains required. No source, prediction, count gate or
score arithmetic changes. Attach the already-recorded failed-call reference
before locking; disclose incomplete cost telemetry in the generated report.
The final report uses the campaign ledger for known usage and unknown reserves.

Pre-score dispatch-accounting clarification: a queued request may be refused
by the campaign budget guard before any API call is issued. The pre-lock binder
may mark incremental usage as exact zero only when the same-paper write-time
decision trace records that specific refusal, no API-call trace exists, and the
campaign ledger contains no reservation for that paper/clinical arm. Genuine
timeouts retain unknown usage and their full reservations. This metadata-only
case keeps a budget-blocked paper in the comparison with its baseline intact;
it is not an evaluated model response. Source/count rows are unchanged.
