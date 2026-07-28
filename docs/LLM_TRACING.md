# LLM trace and curation adjudication

GVF records durable, per-call LLM traces so a curator can reconstruct what the
pipeline asked, what the provider returned, and why a later routing or selection
decision was made.

## The hidden-chain-of-thought boundary

**Model APIs do not expose private hidden chain-of-thought, and GVF does not
claim to record it.** Three different things are recorded, and the recorder keeps
them apart on purpose:

| What | Where | What it means |
|---|---|---|
| Reasoning **content** | `reasoning_capture.response_paths` | The provider actually returned a reasoning summary, a `reasoning_content` string, or a thinking block. `provider_exposed_reasoning_available` is true only for these. |
| Reasoning **token counts** | `reasoning_capture.reasoning_token_usage` | Billing telemetry (`usage.completion_tokens_details.reasoning_tokens`). A count — zero or positive — is **never** evidence that any reasoning content was exposed and can never set the flag. |
| Explicit **rationales** | response envelope + `decision_event` records | Model-authored fields (`tool_rationale`, `inclusion_rationale`, `count_rationale`, `selection_policy`). Ordinary model output; the appropriate artifact for scientific adjudication. |

## What is recorded

Each `llm_call` trace contains:

- gene, PMID, pipeline stage, component, operation, model, and attempt;
- the exact textual messages/input sent to the provider;
- safe request parameters, including temperature, output budget, response
  format, and requested reasoning effort;
- the provider response envelope and visible output text;
- timing, token usage, completion state, and API/parse failure details;
- `reasoning_capture` as defined above.

Each `decision_event` trace records a pipeline decision derived from those
calls, plus the linkage a curator needs to follow retries:

- `accepted_response_trace_id` — the deterministic **primary** link: the exact
  call whose parsed content became the stored result (a JSON repair, when a
  repair produced the data);
- `accepted_response_trace_ids` — every contributing call, chronologically. A
  result assembled from a primary call plus a continuation names both;
- `discarded_trace_ids` — calls that reached the provider but whose content was
  not used;
- `failed_trace_ids` — parse failures and provider errors, kept **separate**
  from discarded so the two stay distinguishable;
- `attempt_trace_links` — every attempt as `{attempt, role, trace_id, outcome}`,
  where `role` is `primary`, `empty_content_retry`, `json_repair`,
  `continuation`, `adjudication`, `fallback`, `table_router`, `figure_ocr`,
  `figure_variant_read`, `pedigree_detect`, `pedigree_extract`,
  `count_recovery` or `claim_debate`, and `outcome` is `accepted`,
  `discarded`, `parse_failed` or `error`;
- `repaired`, `parse_failures`, `provider_attempts`;
- `decision_source` — which branch produced the outcome, including fail-open
  branches, which distinguish a provider **API failure** from an
  **unparseable response**.

### How acceptance is decided

Acceptance is **derived, not mutated in place**. A provider call whose content
parses is recorded as `parsed` — a *candidate*. Which candidate's content became
the stored result is decided by the stage:

- Stages with a single logical call (filters, triage, claim verification) have no
  selection step, so the last candidate is the accepted one.
- Stages that pick among several models call `finalize_attempt_selection(...)`
  **after** the winner is chosen, naming the calls that produced the retained
  result. Everything else that reached the provider becomes `discarded`.

This matters for `paper_variant_extraction`: the model-fallback loop can run
model B successfully and still keep model A's higher-yield result. Letting each
parsing call overwrite the accepted link reported B — the last call made — as the
accepted response, and never marked A's discarded.

Attempt numbers come from the ledger, not from a counter local to the calling
method, because `@llm_retry` re-enters `call_llm_json` wholesale on a retriable
provider error. A method-local counter restarted at 1 on each re-entry, so one
logical call emitted several records all labelled `attempt 1`.

A provider success is **not** acceptance. Count recovery marks its call `parsed`
only after `parse_response` succeeds, so a batch-failure decision cannot claim an
accepted response; the same call can be `parsed` while the decision reports zero
grounded counts, because quote/role validation is a separate gate.

### Redaction

Key names are compared with hyphens and underscores folded together, so
`api-key` (Azure), `x-api-key` (Anthropic), `X-ELS-Insttoken` (Elsevier),
`set-cookie`, `bearer`, `token`, `refresh_token` and `azure-ai-api-key` all
redact. Any mapping stored under `headers` / `extra_headers` /
`default_headers` / `request_headers` / `response_headers` is
**redact-by-default**: header *names* are retained so the request stays
reviewable, but values are dropped unless the header name is on a small safe
allow-list (`content-type`, `api-version`, `user-agent`, `anthropic-version`, …).
Unknown SDK objects serialize to `{"unserializable_object": true, "type": …}` —
never `repr()`, because an SDK client's `__repr__` routinely renders its
credential.

Inline base64 images are represented by media type, size, and SHA-256 digest;
the image bytes are not duplicated into every JSON trace. Text is retained in
full because it is the auditable curation input.

## Production runs

`gvf gvf-run` configures tracing once for the **whole run lifetime** — not just
extraction — so post-extraction stages (count recovery, the per-paper final
check, claim verification) are traced on every path, including
`--skip extract`. Traces live beside the run:

```text
<run_dir>/
  llm_trace_report.html
  llm_traces/
    <GENE>/<PMID>/*.json
    _by_pmid/<PMID>/*.json     # PMID known, gene not
    _unscoped/*.json
    trace_index.jsonl
    trace_manifest.json
```

### Run isolation

A `--skip extract` run reuses an **older** run's directory. Writing into that
run's `llm_traces/` and rebuilding its manifest under the current id would
destroy the earlier manifest and present stale traces as this run's artifact.
So:

- when this run created the run dir → `llm_traces/` + `llm_trace_report.html`;
- when the run dir is reused → `llm_traces_<gvf run id>/` +
  `<gvf run id>_llm_trace_report.html`, and the earlier artifacts are untouched;
- `build_trace_manifest` **refuses** (raises `TraceRunMismatchError`) to rebuild
  a manifest whose records carry a different `run_id`, unless a caller passes
  `allow_mixed_runs=True`, which records `mixed_run_trace_records` as an
  integrity error instead of hiding the mixture.

**`GVF_LLM_TRACE_DIR` is a storage BASE, not one run's directory.** It names
durable storage (an encrypted volume, a large disk); pointing every run straight
at it mixed records across runs and then made every later manifest rebuild raise
`TraceRunMismatchError`. So each run gets a per-run child named after its
sanitized run id:

```text
$GVF_LLM_TRACE_DIR/
  gvfrun-20260728T101500Z-a1b2c3d4/     # one gvf-run
  gvfrun-20260728T134500Z-e5f6a7b8/     # the next one
  recover-counts-…/                     # a standalone recovery
```

`GVF_LLM_TRACE_RUN_ID` is how an already-selected child is handed down. When it
is set, the configured directory is used **verbatim** — that is what stops a
nested stage (the extraction workflow inside `gvf-run`, or a subprocess) from
appending a second level to a path its parent already picked, and what keeps one
`gvf-run`'s extraction and post-extraction records in ONE tree under ONE id. Both
variables are exported for nested and subprocess stages that make model calls
(`scripts/refresh_run_db.py` re-extraction) and **restored** when the top-level
in-process run ends, so a sequential run never inherits the previous run's
identity.

`scripts/recover_counts.py` selects its own run-scoped child for the same reason;
`--trace-dir` overrides it with an exact directory.

### Integrity levels

`trace_index.jsonl` is appended the moment each record lands, so its `sha256` is
the record's **write-time** digest. That is the only thing that can detect
tampering which happened *before* a manifest existed. `build_trace_manifest`
cross-checks every record against it and reports one of three honest levels:

| Level | Meaning |
|---|---|
| `generated_now` | The manifest was built from the files on disk just now (or a digest check failed). It proves internal consistency, **not** that anything is unmodified since capture. |
| `write_time_verified` | Every record still matches the digest recorded when it was written. |
| `locked` | Write-time verified **and** covered by a lock file (`benchmarks/codex_paper_eval` sets this after `lock`). |

Errors surfaced by name: `write_time_digest_mismatch`,
`missing_write_time_digest`, `mixed_run_trace_records`, `trace file changed after
manifest`, `trace_index.jsonl changed after manifest`, `trace file set mismatch`.
The HTML banner prints the level it has; it never says "verified" for a manifest
it generated itself.

### Report size policy

Trace bodies are full prompts. Measured at 0.49 MB of HTML per paper, an
unbounded report reaches ~400 MB at 800 papers, and Tier-2 filter calls fire on
every candidate PMID (1000–3000/gene) into the same file. So:

- embedded strings are capped (`--max-field-chars`, default 24 000). A truncated
  body carries a visible marker stating how many characters were cut and naming
  the on-disk record with the full text;
- records per paper are capped (`--max-records-per-paper`, default 400);
- above `--max-papers-per-file` (default 60) trace groups, the report shards:
  an index page plus one file per paper under `<report stem>_papers/`;
- **nothing is dropped silently.** Every omission appears in the report's
  "Omitted from this page" panel with its extent and record path, and is logged.

### Reading the report

`llm_trace_report.html` is a self-contained, searchable per-paper timeline. Open
it directly in a modern browser; it does not upload data or require a server.
It shows the run id, the integrity level, per-stage **route coverage** (calls,
decisions, accepted links, failures) with an explicit gap column for stages that
made model calls but emitted no decision event, and per-card labels for
`Accepted call` / `Discarded call` / `Retry / repair`. Filters: all, model calls,
decisions, accepted, retries & repairs, failures. Dark mode and print supported.

Rebuild manifest + report after an interrupted or manual workflow:

```bash
.venv/bin/python scripts/build_llm_trace_manifest.py /path/to/run/llm_traces
.venv/bin/python scripts/build_llm_trace_html.py /path/to/run
```

## Route coverage

Every production LLM call goes through `capture_llm_call`. The table below is the
per-route status; **"linkage verified" means a test asserts the accepted /
discarded / failed distinction for that route**, not merely that the code emits
an event.

| Stage | Decision event | Linkage verified by |
|---|---|---|
| `tier2_relevance_filter` | `tier2_relevance_decision` (both fail-open branches) | `test_llm_trace.py` (concurrency) |
| `clinical_data_triage` | `clinical_data_triage_decision` | — emitted, single-call derivation |
| `extraction_priority_triage` | `extraction_priority_triage_decision` | — emitted, single-call derivation |
| `paper_variant_extraction` | `paper_extraction_selection` | `test_route_observability.py::TestExtractionSelectionLinkage` |
| `paper_extraction_continuation` | folded into the parent selection event | `TestStageScopedCoverage` (not a gap) |
| `paper_extraction_adjudication` | folded into the parent selection event | `TestStageScopedCoverage` (not a gap) |
| `clinical_table_routing` | `table_routing_decision` | `TestTableRoutingLinkage` (4 cases incl. deterministic) |
| `variant_claim_verification` | `variant_claim_verification_decision` | — emitted, single-call derivation |
| `variant_claim_debate` | `claim_debate_decision` | `TestClaimDebateRoute` |
| `paper_final_check` | `paper_final_check_decision` | — emitted, single-call derivation |
| `paper_source_grounded_summary` | `paper_source_grounded_summary_decision` | — emitted, single-call derivation |
| `count_recovery` | `count_recovery_decision` | `TestCountRecoveryAcceptTiming` (3 cases) |
| `paper_relevance` | `paper_relevance_decision` | `test_llm_trace.py` (Anthropic path) |
| `synonym_relevance` | `synonym_relevance_decision` | — emitted, single-call derivation |
| `figure_text_extraction` | `figure_text_extraction_decision` | `TestFigureTextRoute` (3 cases) |
| `figure_variant_read` | `figure_variant_read_decision` | `TestFigureVariantRoute` (3 cases) |
| `pedigree_detection` | `pedigree_detection_decision` | `TestPedigreeRoute` (3 cases) |
| `pedigree_extraction` | `pedigree_extraction_decision` | `TestPedigreeRoute` |
| `representation_route` (benchmark) | `representation_route_decision` | `test_codex_paper_eval.py` |
| `paper_curation` (benchmark) | `paper_curation_decision` | `test_codex_paper_eval.py` |

Routes marked "— emitted, single-call derivation" make exactly one logical
provider call, so the accepted link is derived from the single candidate. They are
covered by the shared derivation tests in `test_llm_trace.py` rather than by a
route-specific test.

**Still incomplete:** `scripts/smoke_azure_models.py` is a connectivity probe with
no scope and no decision event, and tracing is not configured by default in
`cli/extract.py`, `cli/discover.py`, `scripts/targeted_land.py`,
`scripts/extract_figure_variants.py`, or `scripts/recall_audit/*` (they inherit
tracing only when an outer `gvf-run` exported it, or when `GVF_LLM_TRACE_DIR` is
set). The route-coverage panel in the HTML report lists any such stage that made
model calls without a registered decision event, so a run says so itself.

Gene and PMID come from the **caller**, not from prompt text, wherever the
caller has them; `utils.llm_trace.infer_trace_context` remains only as a
last-resort fallback. Scope is thread-local (a `ContextVar`), which is what
keeps ~20 concurrent Tier-2 filter workers sharing one filter instance from
cross-attributing one paper's evidence to another.

Coverage matching is **stage-local**: a decision is credited to the stage it
describes (via `EXPECTED_DECISION_EVENTS`), and a stage's expected event must
appear in that stage's own event list. Comparing against a repo-wide set of every
event type seen anywhere let a stage read as satisfied while its own row still
showed `decisions: 0`.

## Locked paper evaluation

New `benchmarks/codex_paper_eval` runs require four trace stages for every paper:

1. `representation_route` — exact routing request and response;
2. `representation_route_decision` — requested versus available/selected route;
3. `paper_curation` — exact extraction request and response, including retries;
4. `paper_curation_decision` — the response accepted as the paper prediction.

`lock` hashes `trace_manifest.json`, `trace_index.jsonl` and
`llm_trace_report.html` alongside `selection.json` and `predictions.json`, then
makes the prediction, report, manifest, index, and individual trace records
read-only. `score` recomputes their digests before it opens gold. Because the
lock-time manifest rebuild now cross-checks write-time digests, a record forged
between extraction and lock fails the lock instead of being re-blessed.

Historical evaluation runs created before this recorder retain final
predictions, evidence, route rationale, source hashes, timing, and aggregate
token telemetry, but not exact raw request/response envelopes. Those missing
envelopes cannot be reconstructed honestly. Re-run the fixed manifest to create
fully traced model comparisons.

## Adjudication workflow

For a disputed paper:

1. Open the trace report, select the gene/PMID, and follow the chronological
   route → call → decision timeline. Check the integrity level in the banner
   before treating a record as evidence.
2. Read the **route coverage** panel: a stage with calls and no decision event
   means the "why" for that step is not recorded, only the raw exchange.
3. Open `selection.json` (benchmark) or `source_completeness.json` (production)
   and inspect the source choice — every candidate rendering, its
   hash/richness metrics, and why one source was chosen.
4. Open the paper's `representation_route` trace and compare the exact catalog
   with the model's returned `tool_rationale`.
5. Open `representation_route_decision` to see whether availability logic
   overrode the requested tool.
6. Open the **accepted** call (the "Open accepted model call ↑" jump in the
   decision card) and compare its evidence, `inclusion_rationale`, and
   `count_rationale` with the stored prediction. Use `attempt_trace_links` to
   see the retries and the repair without confusing a discarded response for the
   accepted one.
7. Record human adjudication separately; never edit a locked trace.
