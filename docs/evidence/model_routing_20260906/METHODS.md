# Interpretation and reproduction notes

This test uses 12 deliberately selected, previously opened gene-paper attempts.
They include missed counts, readable patient tables, source-scope disputes,
correct controls, count-free sources and an abstract-only source. Models see
frozen source text and the target gene. Investigators had earlier exposure to
these papers, so this is opened calibration with no population/human-accuracy
or successful discovery/confirmation claim.

The valid primary comparison is the current Grok 4.3 configuration versus
Astra medium, using the same sources and shared stages. Grok has a 15k local
output cap; Astra requests 32k. Sol-medium verification/adjudication, Kimi table
routing and Luna-xhigh Tier 2 are fixed. Two papers take deterministic extraction
paths. This is a configuration/workflow comparison, not an isolated model-weights
comparison or a replay of every historical default. API completion, valid JSON,
extracted identities and accepted counts are separate outcomes.

The independent Astra-medium overlay reads each triggered source without the
baseline answer or gold. The trigger requires a retained baseline identity with
a null count plus a clinical signal in the source. It only fills null fields on
unambiguously matched baseline identities. Explicit current-study values must
pass the existing literal-evidence validator. Existing nonnull values stay put;
contradictions and rejected proposals are retained. Raw additive proposals form
an unvalidated diagnostic lane. The overlay cannot recover an identity the
baseline missed or automatically repair a wrong nonnull baseline count.

All final arms lock before any new scoring. Pilot calls made before the actual
SDK payload correctly transmitted reasoning effort are abandoned unscored.
Their requested-effort labels do not establish a high-versus-medium contrast.
Corrected Grok 4.6 failed bounded deployment probes and was not dispatched as a
paper arm. The Grok CLI is a different service path and its working reviews do
not establish Azure deployment availability. See PLAN.md's dated amendments.

Identity TP/FP/FN compare variant-paper rows to reference membership. For each
carrier/affected/unaffected field, report exact supplied values, wrong supplied
values, positive-reference coverage, and absolute error over every asserted
reference row. Missing identities and matched abstentions contribute zero only
in the evaluation arithmetic; their stored values remain null. Supplied counts
on identity extras are reported separately so a prettier matched-count score
cannot hide count-bearing false positives. Disagreement with the reference is
not by itself source adjudication.

Report both preflagged sensitivity sets: omit SCN5A 20129283, then also omit
RYR2 25814417. Do not replace official panel scores with these subsets. Preserve
cardiac/manual versus MYBPC3 reference provenance and per-paper deltas. No new
human reread was conducted and no population-weighted sample was drawn.

Budget accounting uses returned-usage price proxies, includes reasoning tokens,
retains unknown failed/canceled-call reservations and observed successful-call
SDK retry reservations, and adds a $25 uncertainty margin within the $150
API envelope. CLI cost is excluded. Cache-price adjustments and Azure invoice
reconciliation are not applied. Request durations include transport retries;
whole-run elapsed time also includes waits for budget reservations, so it is
not a pure model-latency measure. Native failed-call telemetry remains null
with a hash-bound failure trace and known returned subset, never fabricated zero.

The separate compact-roster diagnostic changes source selection, task, schema
and effort together. Its transcription audit is against converted source text
and the original DOC's table structure, with explicit qualifications for two
blank family cells and competing phenotype definitions. It does not emit new
accepted scalar predictions or measure a general effect of lower reasoning.

Reproduction entry points (the archived runs remain immutable):

- `campaign.py` prepares/extracts/projects/locks/scores the primary arms with
  source/runtime fingerprints and the experiment budget hook.
- `clinical_reader.py` builds the independent overlay; `bind_failed_usage.py`
  attaches the already-recorded failure metadata before its native lock.
- `analyze.py` summarizes locked scored arms and the campaign cost ledger.
- `operational_summary.py` summarizes primary/clinical call outcomes.
- `audit_compact_roster.py` and `audit_original_layout.py` validate the diagnostic.
- Original executed helper bytes are retained in `executed_code_snapshots.json`
  where formatting or a superseding launcher would otherwise obscure them.

The canonical stratified figure retains its registered historical/opened-tranche
sampling frames. This selected model panel receives its own per-run and paired
companion figures; it is not pooled into the canonical cohort denominator.

A useful way to locate the limit is to measure successive conditions: acquired
count-bearing evidence, correctly bound identities/people, the intended phenotype
definition, and a validated accepted assertion. These are conditional stages,
not independent percentages to multiply. A model upgrade can improve reading
while end-to-end recovery stays flat because an earlier source is missing or a
later derivation gate rejects the output. The experiment does not quantify an
intrinsic human-relative ceiling for any of these stages.

## Pre-score accounting amendment

The final runtime check detected a changed benchmark file. Exact reconstruction
of the prepared 251-file fingerprint proved the initial drift was confined to
`run_eval.py`'s failed-call usage contract and reporting. Before lock, inspection
also found that the production exporter silently omitted failed-call usage.
`db_to_predictions.py` now publishes null totals, failed trace IDs and the known
returned subset, including within model buckets. Both changes concern accounting.
No extraction/source-processing runtime changed while the API jobs ran.

The first locking driver was stopped before lock. Its audit receipt is preserved
in `runtime_telemetry_amendment.json`; the final two-file reconstruction and
re-export receipt are `production_usage_amendment.json`. Original setup hashes
are unchanged. The final export must have identical scientific prediction arrays
to the earlier export before lock is allowed. Prior source bytes are retained in
JSON snapshots, and all final arms must lock before scores are opened. The
clinical binder independently proves that its timeout has unknown usage while
the final budget refusal was never dispatched, with no call/ledger row. These
are pre-score measurement amendments, not an unchanged-harness claim.

A further **post-hoc descriptive subset** omits both deterministic-shortcut
papers, SCN5A 20129283 and 30059973. Their extraction path was known before
scoring; this additional reporting slice was added after the control score
showed their large share of the identity denominator. It does not replace the
12-paper result or introduce a passing gate. The two preflagged source/reference
sensitivity sets above remain distinct.
