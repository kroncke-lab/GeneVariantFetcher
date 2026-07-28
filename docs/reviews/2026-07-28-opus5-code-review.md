# GVF full read-only audit — 2026-07-28 (independent, Opus 5)

Status: dated audit snapshot. For implementation status and the remaining
forward checklist, see `2026-07-28-opus5-implementation-verification.md` and
`TASKS.md`; do not treat the original line references below as current code.

Scope: HEAD `037a609` + full dirty diff (22 modified, 11 untracked). No network, no live APIs, no live model calls. Probes run with `PYTHONDONTWRITEBYTECODE=1 -B -p no:cacheprovider`. **No files created, modified, deleted, formatted, staged, or committed.**

Model provenance: Claude Code `2.1.218`, explicit model ID `claude-opus-5`, effort `max`, canonical model reported as `claude-opus-5`. Session `b94a088d-1ae2-42c1-a537-a8510d915822`.

## Verdict

Tracing is architecturally right and honestly *worded*, but not yet trustworthy as an audit artifact. **Zero routes bypass `capture_llm_call`** — every provider boundary is wrapped. All failures are downstream of capture: attribution (wrong PMID under concurrency), linkage (no accepted-response pointer in production), integrity (manifest cannot detect pre-manifest tampering; the one write-time digest source is unused), packaging (the HTML generator is unreachable in an installed wheel), scale (one 150–400 MB HTML per gene).

Count recovery has a **new P0 the prior audit did not reach**: the *table* branch of the role gate has no negative-context filter at all, so an explicitly labelled `allele count`, `gnomAD AC`, `Age at diagnosis`, or `families screened` value is accepted as a per-variant carrier count.

Focused suite `128 passed in 103.38s`; `ruff` clean; `git diff --check` clean.

---

## Direct answers

### 1. Where exactly is the code-loading / import error?

`pyproject.toml:60-62` — `[tool.setuptools.packages.find] include` omits `scripts*` (and `benchmarks*`, `tests*`). `scripts/__init__.py` exists (0 bytes) so it is a real package *locally*, resolving off `sys.path[0]`; it is simply absent from a wheel. Seven shipped-package sites import it:

| Site | Guarded? | Effect in installed `gvf` |
|---|---|---|
| `cli/gvf_run.py:1275` `scripts.recover_counts` | caller `try` at `:1796` | Step 3.55 degrades to a `stage_warning`, exit 0 |
| `cli/gvf_run.py:2089` `scripts.build_llm_trace_html` | `try` `:2088`/`except` `:2127` | **no `llm_trace_report.html` for any gvf-run** |
| `cli/automated_workflow.py:1131` same | `try` `:1130`/`except` `:1155` | **no HTML for any workflow run** |
| `cli/compare_variants.py:2336` | **NO** | `ImportError` propagates |
| `cli/compare_variants.py:2367` | **NO** | propagates — kills `load_adjudication_overlay_db` |
| `cli/compare_variants.py:2470` | **NO** | propagates |
| `cli/gvf_run.py:1142` `scripts.backfill_paper_metadata` | `try` `:1141` | stage warning |

Broader than the prior framing in both directions: the three `cli/compare_variants.py` sites are **unguarded hard crashes** in the scorer's adjudication-overlay path, and the two HTML sites make the tracing feature's headline artifact **silently dead** in exactly the deployment an operator would use. `benchmarks/codex_paper_eval/run_eval.py:34` is also a hard module-level import, mitigated only by its own `sys.path.insert` at `:31`. Masked locally: the project is not installed in `.venv` (no dist-info), and `Makefile:17` / `.github/workflows/ci.yml:40,62,86` use `pip install -e .`.

### 2. What is non-functional, dead, misleading, or incomplete?

**Non-functional in an installed package:** automatic HTML report (both paths), count recovery, adjudication-overlay scoring.

**Dead:**

- `config/settings.py:348` `count_recovery_enabled` — nothing reads it; `cli/gvf_run.py:1262` reads the env var directly.
- `scripts/build_llm_trace_html.py:1032-1038` — the "Open accepted model call ↑" jump button fires only on `accepted_response_trace_id`, emitted **only** by `gene_literature/llm_relevance.py:288,442` and `run_eval.py:1094`. In production extraction/filters/triage/claim-verification/final-check/count-recovery it never appears, so the button never renders.
- `collect_trace_report_data` returns `run_id` (`:235`); the template never renders it.
- `utils/llm_trace.py:146` `current_trace_context` and `:237` `infer_trace_context` have no external callers.
- `config/settings.py:452,468` final-adjudicator/arbiter efforts are read only by `scripts/recall_audit/`.

**Misleading:**

- HTML integrity banner (`:956-960`) says "Trace integrity verified" for a manifest generated moments earlier from the files it validates — survives a forged prompt (proved below).
- `reasoning_capture.provider_exposed_reasoning_available` is `true` whenever the provider returns `reasoning_tokens`, **including `0`**.
- `--dry-run` reports counts it will not write.
- `config/settings.py:751-755` still says `vision_reasoning_effort` "is inert until those sites are wired" — four sites are wired.
- `utils/llm_trace.py:486` calls the manifest "bounded"; it is unbounded.

**Incomplete:** per-route scope (see matrix), accepted-response linkage, retry/attempt indexing, decision events for table routing / count recovery / final check / pedigree / figure reads, and **zero `docs/` coverage of count recovery**.

### 3. What must change before tracing can honestly claim all observable reasoning steps and produce a useful HTML artifact for every run path?

(a) Make trace context thread-safe and stop clobbering the per-thread scope with shared instance state. (b) Cross-check `trace_index.jsonl` write-time digests inside `build_trace_manifest`/`validate_trace_manifest`, lock the index, and stop asserting "verified" for a just-generated manifest. (c) Move `build_llm_trace_html` into a packaged module, and configure tracing for the whole `gvf-run` lifetime plus every script entry point. (d) Emit `accepted_response_trace_id` + `attempt` on every decision event. (e) Shard the report per gene/PMID (or lazy-load bodies) with an explicit, logged truncation policy. (f) Narrow `reasoning_capture` to reasoning *content*, never token counts.

---

## Findings

### P0-1 — Count-recovery role gate accepts denominators, cohort totals, and explicitly non-carrier table columns as carrier counts

`pipeline/count_recovery.py:98-113` (`_FIELD_CONTEXT_RE["carriers"]` accepts `patients|individuals|subjects|probands|cases|identified|observed|seen|found`); `:114-118` (`_NON_COUNT_CONTEXT_RE` has no `X of Y`, `X/Y`, `among Y`, `N cases of <phenotype>`); `:518-539` (`:526` short-circuits `return True` when `field is None` **or** the quote is one piped line); `:689-698` `_labeled_table_value`; `:738-757` single-number segment; `:834-845` routing. **Every call from the table branch passes `field=None`** (`:694`, `:696`, `:730`, `:755`), so the table path applies no role check and no negative-context check at all.

Live probes, all **ACCEPTED as `carriers`**:

| Verbatim quote | Accepted |
|---|---|
| `The p.Leu552Ser variant was found in 12 of 812 subjects.` | `12` **and** `812` |
| `p.Leu552Ser was present in 5/120 patients` | `5` **and** `120` |
| `Among 44 cases of Long QT, the p.Leu552Ser variant was identified once.` | `44` **and** `1` |
| `The p.Leu552Ser variant was identified in 3 of 200 probands (1.5%).` | `3` **and** `200` |
| `n=7 among 913 individuals carried p.Leu552Ser` | `7` **and** `913` |
| `\| p.Leu552Ser \| allele count \| 12 \|` | `12` |
| `\| p.Leu552Ser \| gnomAD allele count 37 \|` | `37` |
| `\| Variant \| gnomAD AC \|` / `\| p.Leu552Ser \| 12 \|` | `12` |
| `\| Variant \| Age at diagnosis \|` / `\| p.Leu552Ser \| 45 \|` | `45` |
| `\| p.Leu552Ser \| families screened \| 3 \|` | `3` |

The last five are new. They are exactly the PMID 33013630 annotation-table failure class that `TASKS.md:173-184` tracks as a known precision defect — recovery would reintroduce it *additively*, into previously-NULL slots.

**Impact:** corrupts carrier MAE and penetrance denominators, and reads as a recall improvement. **Confidence:** high (reproduced). **Exposure:** flag defaults off, but `scripts/recover_counts.py` is callable and an obsolete pre-hardening scratch run had already exercised it against copies of the production DBs; that scratch run was removed during main-branch consolidation. **Fix:** structured `count_role` + `evidence_locator` from the model; reject denominator grammars unless the integer carries an explicit local carrier/variant-positive label; fail closed on ≥2 candidate integers in prose (matching what `:753-757` already does for tables); apply a widened `_NON_COUNT_CONTEXT_RE` (`allele count|AC|AN|AF|age|onset|years|follow-up|QTc|ms|bpm`) to the table branch and `_labeled_table_value`; never let `field=None` mean accept. Land all ten quotes as `validate_paper_response` fixtures.

### P0-2 — `.env.swp` still present with recoverable secrets

`.env.swp`, 16384 B, mtime 2026-07-21 17:48, mode 0600. `git check-ignore -v` → `.gitignore:42:*.swp` (untracked). The 2026-07-27 review recorded that both Grok and AGY inspected it inside a review snapshot. **It has not been removed and the represented keys are not visibly rotated.** **Fix:** close the owning session, delete, rotate, add `*.swp` to snapshot/archive exclusion lists (not just `.gitignore`). Left in place — read-only task.

### P1-3 — Shared-instance trace context misfiles Tier-2 traces and decisions under the wrong PMID

`utils/llm_utils.py:491-495` `set_llm_trace_context` overwrites **instance** state; `:497-502` `_trace_scope` merges it *on top of* the ContextVar scope, so a stale value **overrides** a correct per-thread scope; `:504-514` `record_llm_decision` reads the same dict. `pipeline/steps.py:693-698` creates **one** `InternFilter`/`ClinicalDataTriageFilter`; `:938` submits `_classify_pmid` to a `ThreadPoolExecutor` (default `filter_max_workers` approximately 20 for Anthropic); `_classify_pmid` calls the shared instance (`:823`, `:849`) → `pipeline/filters.py:358`, `:647`.

**Evidence:** an 8-thread probe over one shared caller observed **7/8 calls bound to another thread's PMID**.

**Impact:** Tier-2 prompt/response records and every `tier2_relevance_decision` / `clinical_data_triage_decision` event are attributed to the wrong paper — a curator asking "why was PMID X dropped" reads another paper's evidence. Extraction escapes only because `pipeline/steps.py:1795-1805` uses a `threading.local()` extractor. **Confidence:** high (reproduced). **Fix:** delete the instance dict; replace all five sites with `with llm_trace_scope(...)` as `pipeline/extraction.py:6606` already does. Test: N threads on one shared `InternFilter`, assert each record's `context.pmid` matches its own paper.

### P1-4 — Manifest integrity is tautological; the write-time digest source is unused and unlocked

`utils/llm_trace.py:361` records a write-time SHA-256 into `trace_index.jsonl`. `:480-536` `build_trace_manifest` **re-hashes at manifest time** and never consults the index. `:539-577` validates files against that re-hash. `scripts/build_llm_trace_html.py:124-133` *generates* a manifest when absent and returns `[]` errors; `:956-960` then prints "✓ Trace integrity verified." `run_eval.py:531-541` rebuilds the manifest at **lock** time, overwriting the extraction-time one.

**Evidence:** forging a prompt inside a record *before* manifest generation → `validate_trace_manifest` returned **no errors**; the index's write-time SHA (`a070304e…`) demonstrably differed from the manifest SHA (`6dcc990b…`) with nothing comparing them; HTML reported `integrity.valid = True` and embedded the forged text. `trace_index.jsonl` is neither digested in `LOCK.json` (`:547-556`) nor chmod-ed read-only (`:568-572`).

**Impact:** `docs/LLM_TRACING.md:73-78` claims the traces sit "inside the same pre-gold integrity boundary as the predictions". The real guarantee is only "unchanged since lock", and the curator-facing artifact overstates it. **Confidence:** high (reproduced). **Fix:** compare write-time vs current digests in `build_trace_manifest`, emit a `write_time_digest_mismatch` error, include the index in manifest + `LOCK.json` + chmod set, and make the banner distinguish "verified against write-time digests" from "manifest generated now".

### P1-5 — Installed-package `scripts.*` imports

See Direct answer 1. **Fix:** move `default_source_roots`/`make_source_resolver` into `pipeline/count_recovery.py`, `build_trace_html_report`+`TRACE_REPORT_NAME` into `utils/llm_trace_html.py`, and `_variant_key`/`VERDICT_TO_ACTION`/`gold_tier_includes_gene` into a packaged module; leave the `scripts/` files as thin CLIs. Add an installed-wheel smoke test (`gvf gvf-run --help`, import `load_adjudication_overlay_db`, import the HTML builder). Make explicit `COUNT_RECOVERY_ENABLED=true` + import failure a stage *failure*.

### P1-6 — Recovered counts carry no role provenance and inherit stale provenance

`pipeline/count_recovery.py:965-1033` writes only the count column plus a `count_recovery_log` row — never `variant_papers.count_provenance`, `penetrance_data.trust_sources`, `trust_tier`, or `field_trust`. `pipeline/trust_gate.py:240-247` gates on `carriers_count_type`/`carriers_column_label` read from `count_provenance` (`:364`, `:390-406`); `:453-457` only *writes* `trust_sources`. Nothing reads `count_recovery_log`. Two symmetric failures: NULL provenance (linkage rows, and the fresh `INSERT` at `:997-1001`) → no rule fires, value stays `trusted` (prior probe: `evaluate_fact({"carriers": 812,…}, {}, …)` → `[]`); stale `carriers_count_type='cohort_total'` → a correctly recovered per-variant count is spuriously stamped `count_is_total`.

Compounding: `harvesting/migrate_to_sqlite.py:667` gives `trust_tier` DEFAULT `'trusted'`, and `cli/gvf_run.py:1794` gates Step 3.55 **independently of `full_coverage`**, so `--skip trust-gate` with recovery on leaves recovered counts `trusted` with zero evaluation. **Fix:** persist role + locator into the provenance the gate reads (or a recovery trust source); insert as `quarantine` and let 3.7 promote; refuse 3.55 when `trust-gate` is skipped.

### P1-7 — Production extraction has no accepted-response linkage; retries are indistinguishable

`utils/retry_utils.py:135-143` allows 8 attempts; `call_llm_json` re-calls `_make_call()` on empty content (`:719`); `_attempt_json_repair` (`:516-579`) issues a further call whose result can *become* the accepted data. Each writes a separate `llm_call` with an **identical** `context` and **no `attempt` index**. `pipeline/extraction.py:6597-6628` `paper_extraction_selection` records `selected_model` but **no `accepted_response_trace_id`**. `_attempt_continuation` (`:4645`) and `_call_adjudicator` (`:5487`) inherit `stage="paper_variant_extraction"`.

`docs/LLM_TRACING.md:14` claims "attempt when known"; `:100-102` tells the curator to "follow retries without confusing a discarded response for the accepted response" — a mechanism that exists **only** in the benchmark and relevance paths. For production, a curator cannot tell which of up to approximately 16 records produced the stored extraction. **Fix:** thread `attempt` through `_trace_scope`; return/stash the trace summary from `call_llm_json*`; put `accepted_response_trace_id` + `repaired` + `discarded_trace_ids` on `paper_extraction_selection`; give continuation, adjudication, and routing their own `stage`.

### P2-8 — HTML report does not scale

`scripts/build_llm_trace_html.py:136-267` loads every record into memory; `:272-282` serialises one `__DATA__` blob; `:801`/`:1167` parse it synchronously on load. No pagination, per-paper split, size cap, or truncation. **Measured** (40 papers × 4 calls, 120 KB source approximately 30k tokens): traces 19.8 MB → **HTML 19.8 MB**, peak builder RSS **345 MB**, **0.49 MB HTML/paper** → **49 MB @ 100, 148 MB @ 300, 396 MB @ 800 papers**. That undercounts: Tier-2 filter calls fire on every candidate PMID (1000–3000/gene) and land in the same report. **Fix:** shard per gene/PMID (or index page + `fetch()`-on-demand), cap embedded bodies with a "truncated — see `llm_traces/<path>`" marker, and **log** what was dropped.

### P2-9 — Redaction misses the header spellings this repo uses

`utils/llm_trace.py:39-49` + `:150-158` normalise case but not hyphens, and have no token rule.

| Key | Verdict |
|---|---|
| `api_key`, `apikey`, `authorization`, `cookie`, `secret`, `client_secret`, `openai_api_key` | REDACTED |
| `api-key` (Azure), `x-api-key` (**Anthropic**), `insttoken`, `X-ELS-Insttoken` (**Elsevier**), `set-cookie`, `token`, `refresh_token`, `bearer`, `azure-ai-api-key` | **RETAINED** |

`json_safe({"headers": {"api-key": "sk-SECRET", …}})` retains `api-key` and `X-ELS-Insttoken` verbatim. Currently latent — the three Azure wrappers pass `request={"endpoint","body"}` and omit headers (`figure_text_extractor.py:198-204`, `figure_variant_reader.py:353-359`, `pedigree_extractor.py:116-122`), and LiteLLM's injected `api_key` (`utils/llm_utils.py:199-201`) *is* redacted. One `extra_headers=`/`default_headers=` on any traced request writes a live key. Secondary: `json_safe:218-222` falls back to `vars()` and `:223` to `repr()` for unknown SDK objects. **Fix:** normalise `-`↔`_`; add token/cookie/bearer keys; treat any `headers`/`extra_headers`/`default_headers` map as redact-by-default with a small allow-list; drop the bare `repr()`.

### P2-10 — `reasoning_capture` claims exposed reasoning from a token count

`utils/llm_trace.py:263-274` flags any key matching `reasoning|thought|thinking` whose value is not in `(None,"",[],{})`; `:443-446` derives the boolean. **Evidence:** `{"usage":{"completion_tokens_details":{"reasoning_tokens":0}}}` → `['$.usage.completion_tokens_details.reasoning_tokens']`, flag `true`. So every Azure GPT-5-family call is recorded as having exposed reasoning. This is the one claim the feature must never overstate, and it contradicts `docs/LLM_TRACING.md:30-33` and `utils/llm_trace.py:350-354`. **Fix:** scan only reasoning *content* shapes (`output[].type=="reasoning"`, `reasoning.summary`, `reasoning_content`, thinking blocks); report counts separately as `reasoning_token_usage`.

### P2-11 — Tracing unconfigured on several run paths; skip-extract rewrites the previous run's manifest

The only production `configure_llm_tracing` is `cli/automated_workflow.py:179`, inside extraction. `cli/gvf_run.py` never calls it. `gvf-run --skip extract` takes `:1658` (`run_dir = runs[-1]`), so Step 3.55 and Step 3.8 run **untraced**; then `:2085-2128` finds the *previous* run's `llm_traces/`, rebuilds its manifest under the *current* `run_id` (`:2100-2113`), and regenerates the HTML — destroying the earlier manifest and presenting stale traces as this run's artifact. `scripts/recover_counts.py` never configures tracing. Also unconfigured: `cli/extract.py`, `cli/discover.py`, `scripts/refresh_run_db.py`, `scripts/targeted_land.py`, `scripts/extract_figure_variants.py`, `scripts/smoke_azure_models.py`, `scripts/recall_audit/*`. `benchmarks/cold_start_eval/run_cold_start.py:195` pops every `GVF_*_DIR` except `GVF_CORPUS_DIR`, silently discarding a `GVF_LLM_TRACE_DIR` override.

### P2-12 — Per-route scope coverage gaps

`harvesting/figure_text_extractor.py:122-126` — no gene, no pmid (caller `paywall_context_enrichment.py:483` has the PMID, sets no scope) → `_unscoped/`. `pipeline/pedigree_extractor.py:291-295` — same; `harvesting/orchestrator.py:656-689` has the PMID and zero trace calls. `harvesting/figure_variant_reader.py:273-278` — gene present, **pmid dropped** though available at `read_images` (`:173-186`); since `_trace_parent` (`utils/llm_trace.py:253-260`) needs both, every figure read for every paper collapses into one group. `pipeline/table_router.py:1408`, `:1479-1505` — no stage/component, so router calls are indistinguishable from primary extraction. `pipeline/count_recovery.py:1036-1141` — zero trace calls; gene and `gap.pmid` are in hand at `:1083`; attribution survives only because `infer_trace_context` (`:52`) matches the prompt's `(PMID 12345678)`, and with no `stage` the file is named `..._llm_call_llm_call.json` and the UI labels it generically (`build_llm_trace_html.py:838`). `pipeline/paper_final_check.py:1714-1760` — the expensive `gpt-5.6-sol@xhigh` check has **no** tracing wiring and no decision event.

### P2-13 — `--dry-run` overstates what recovery will write

`pipeline/count_recovery.py:988-990` increments `written` for every accepted count **without** consulting existing values, and skips `_ensure_log_table`. **Evidence:** one accepted count against a populated slot → dry-run `{'counts_written': 1, 'already_populated_skipped': 0}`; real run `{'counts_written': 0, 'already_populated_skipped': 1}`; stored value unchanged (44). `scripts/recover_counts.py:12-13` mandates dry-run first, so the pre-flight an operator gates on is inflated.

### P2-14 — No pre-mutation backup; no incremental durability

`cli/gvf_run.py:1265-1310` mutates `penetrance_data` in place with no backup, unlike `scripts/recall_recovery/run_all_layers.py:310-316` which copies to `*.before_layers.db` and warns when absent (the practice `TASKS.md:244` records as landed; `CLAUDE.md:86-87` also forbids direct row patching). Separately `pipeline/count_recovery.py:1070-1114` calls `write_recovered_counts` **once** at `:1112`, with a single implicit transaction and `con.close()` in `finally` (`:1024-1027`) — a crash or interrupt after N papers of paid model calls writes **nothing**. No `PRAGMA busy_timeout`.

### P2-15 — Recovery failure invisible to exit code, severity, and the report

`cli/gvf_run.py:1794-1809` appends to `stage_warnings` only; `:2130-2131` sets `exit_code` from `stage_failures`; `_write_run_status` (`:1433-1482`) sets `severity="ok"`. Recovery can fail on **every** paper and the run exits `0` with `status="completed"`, `severity="ok"` — asserted as intended by `tests/unit/test_gvf_run_pipeline_wiring.py:219-261`. Stats are logged at `:1303-1309` then discarded — they reach neither `RUN_REPORT.md` nor `RUN_STATUS.json`. `scripts/recover_counts.py:164` always returns 0.

### P2-16 — Count-recovery routing bypasses the provider abstraction; provenance omits it

`config/settings.py:360-369` hardcodes `azure_ai/gpt-5.6-sol` while `model_provider` defaults to `anthropic` (`:133-134`). There is **no** `get_count_recovery_model()` beside `get_tier2_model()` (`:1145`), `get_vision_model()` (`:1192`), `get_paper_final_check_model()` (`:1250`), so `cli/gvf_run.py:1296` reads the field raw: an Anthropic-only deployment silently routes 3.55 to Azure and fails every call, with no key preflight. `utils/provenance.py:128-155` omits `count_recovery_model` and `count_recovery_reasoning_effort`, so a run's provenance does not record which model filled its counts. `settings.count_recovery_enabled` is dead.

### P2-17 — Benchmark and production disagree on reasoning-effort routing

`benchmarks/codex_paper_eval/run_eval.py:733-746` gates on `"grok" in model.lower()` and otherwise always sends `{"reasoning":{"effort":…}}`; `utils/llm_utils.py:138`+`:240-286` gate on `gpt-5|gpt5|o1|o3|o4-mini`. A deployment that is neither receives effort in the benchmark and **not** in production, so benchmark results do not describe the production configuration for the same model string.

### P2-18 — The only count-recovery measurement bypassed the harness's validation and lock

The now-removed pre-hardening recovery scratch run used `predictions.json`
schema version 1, had no `llm_trace_refs`, no `curation_rationale`, and no
`llm_traces/` directory. Its bespoke `LOCK.json` was not produced by
`command_lock`, so `validate_predictions` never ran. It was a directional
experiment, not an enablement artifact; `TASKS.md` now requires its replacement
to use the current validation, lock, and trace contract.

### P3 (batch)

| # | Finding | Location |
|---|---|---|
| 19 | `_trace_scope` raises `TypeError: got multiple values for keyword argument 'component'` if `set_llm_trace_context` was given `component=`/`operation=` — latent crash inside `call_llm_json`. Reproduced. | `utils/llm_utils.py:491-502` |
| 20 | Extraction triage gets **no** `TIER2_REASONING_EFFORT` and `max_tokens=1200`, below the 8192 reasoning-model floor `filters.py:260-265` applies to the sibling Tier-2 classifier → truncation risk. | `pipeline/extraction_triage.py:385-389` |
| 21 | `_call_adjudicator` swaps model/tokens/temperature but keeps the primary model's `reasoning_effort`; no `TIER3_ADJUDICATOR_REASONING_EFFORT` exists. | `pipeline/extraction.py:5475-5491` |
| 22 | `_attempt_json_repair` omits `build_reasoning_effort_kwargs` unlike its three siblings — not a correctness defect, but the trace records the accepted-JSON path at a different effort than its parent. | `utils/llm_utils.py:543-558` |
| 23 | An Anthropic API failure is recorded as `decision_source="fail_open_unparseable_response"`. | `gene_literature/llm_relevance.py:139-141`, `:270` |
| 24 | `DATA.run_id` never rendered; no per-record run separation, so a resumed run's mixed-run dir reads as one run. | `build_llm_trace_html.py:235` + template |
| 25 | `GVF_LLM_TRACE_DIR`/`_ENABLED` absent from `.env.example` (0 hits); `GVF_LLM_TRACE_RUN_ID` documented nowhere. | `utils/llm_trace.py:116,124,129` |
| 26 | Count recovery has **zero** `docs/` coverage — not in `ARCHITECTURE.md`, `PROTOCOL_CHANGELOG.md`, or `EXTRACTION_CONTRACT.md`, contrary to `CLAUDE.md:25-27`. | `docs/` |
| 27 | Bare `next(...)` → `StopIteration` if `choose_source` returns a path outside `usable`. | `run_eval.py:174-176` |
| 28 | `penetrance_data` has **no index at all** on `(variant_id, pmid)` → the write path's per-count SELECT is a full scan; `find_count_gaps` fetches the whole gene and filters PMIDs in Python (`:255-259`). | `migrate_to_sqlite.py:657-676`; `count_recovery.py:991-995` |
| 29 | Retained full prompts + 2000-char publisher quotes (`recovery_*.json`, `llm_traces/`, embedded HTML) create copyright/retention obligations; 1.3 MB already sits untracked in-tree. | `count_recovery.py:60`; `recover_counts.py:148-163` |
| 30 | `_paper_metadata` reads only benchmark files, so production reports never render title/source-choice/variant chips — `docs/LLM_TRACING.md:94-96` step 2 does not apply to production. | `build_llm_trace_html.py:59-91` |
| 31 | `stats["results"]` returns live dataclasses; `gvf_run` pops them (`:1302`) and the script hand-serialises them. Any new `json.dumps` consumer raises. | `count_recovery.py:1121` |
| 32 | `recover_counts` never plumbs `paper_derived_only`, so the documented inspection mode (`:224-227`) is unreachable. | `count_recovery.py:1057` |

---

## Route-coverage matrix

**Raw** = wrapped in `capture_llm_call`; **Scope** = gene/PMID/stage/attempt; **Decision** = normalized `decision_event`; **Accepted** = `accepted_response_trace_id`.

| Route | File:line | Provider | Raw | Scope | Decision | Accepted | Retry/parse-fail visible | Config'd by default |
|---|---|---|---|---|---|---|---|---|
| LiteLLM choke point | `utils/llm_utils.py:220-226` | LiteLLM | ✅ | caller | n/a | n/a | ✅ per attempt | via workflow |
| Tier-2 relevance | `filters.py:358-363` | LiteLLM | ✅ | ⚠️ **race → wrong PMID** | ✅ ×3 incl. fail-open | ❌ | ✅ | ✅ |
| Clinical-data triage | `filters.py:647-652` | LiteLLM | ✅ | ⚠️ same race | ✅ ×2 | ❌ | ✅ | ✅ |
| Extraction-priority triage | `extraction_triage.py:295-326` | LiteLLM | ✅ | ✅ (serial) | ✅ | ❌ | ✅ | ✅ |
| Paper extraction | `extraction.py:6491`, scope `:6597-6628` | LiteLLM | ✅ | ✅ no `attempt` | ✅ `paper_extraction_selection` | ❌ | ⚠️ up to approximately 16 identical records | ✅ |
| Continuation | `extraction.py:4645` | LiteLLM | ✅ | ⚠️ inherits parent stage | ❌ | ❌ | ⚠️ | ✅ |
| Tier-3 adjudication | `extraction.py:5487` | LiteLLM | ✅ | ⚠️ inherits parent stage | ❌ | ❌ | ⚠️ | ✅ |
| JSON repair | `llm_utils.py:546` | LiteLLM | ✅ | `operation="json_repair"` | ❌ | ❌ | ⚠️ repair-accepted unmarked | ✅ |
| Table routing | `table_router.py:1481-1505` | LiteLLM | ✅ | ❌ no stage/component | ❌ (extraction metadata only) | ❌ | ⚠️ | ✅ |
| Claim verification | `claim_verifier.py:305-314` | LiteLLM | ✅ | ✅ + `variant` | ✅ | ❌ | ✅ | ✅ |
| Claim debate | `run_claim_debate_pilot.py:171` | LiteLLM | ✅ | ❌ | ❌ | ❌ | ⚠️ | ❌ |
| Paper final check 3.8 | `paper_final_check.py:1733` | LiteLLM | ✅ | ❌ none | ❌ | ❌ | ⚠️ | ⚠️ not on skip-extract |
| Source-grounded summary | `paper_final_check.py:1759` | LiteLLM | ✅ | ❌ none | ❌ | ❌ | ⚠️ | ⚠️ same |
| **Count recovery 3.55** | `count_recovery.py:1140` | LiteLLM | ✅ (prod caller) | ❌ pmid inferred only, no stage | ❌ | ❌ | ⚠️ batch failure logged only | ❌ standalone / skip-extract |
| Figure OCR (chat) | `figure_text_extractor.py:130` | LiteLLM | ✅ | ❌ no gene/pmid | ❌ | ❌ | ⚠️ | ✅ |
| Figure OCR (Responses HTTP) | `figure_text_extractor.py:198` | Azure `requests` | ✅ | ❌ no gene/pmid | ❌ | ❌ | ✅ | ✅ |
| Figure variant read (chat) | `figure_variant_reader.py:291` | LiteLLM | ✅ | ⚠️ gene only, pmid dropped | ❌ | ❌ | ⚠️ | ✅ |
| Figure variant read (Responses HTTP) | `figure_variant_reader.py:353` | Azure `requests` | ✅ | ⚠️ gene only | ❌ | ❌ | ✅ | ✅ |
| Pedigree vision (chat) | `pedigree_extractor.py:306` | LiteLLM | ✅ | ❌ no gene/pmid | ❌ | ❌ | ⚠️ | ✅ |
| Pedigree vision (Responses HTTP) | `pedigree_extractor.py:116` | Azure `requests` | ✅ | ❌ no gene/pmid | ❌ | ❌ | ✅ | ✅ |
| Paper relevance | `llm_relevance.py:122`, `:245-291` | Anthropic SDK | ✅ | ✅ | ✅ | ✅ | ⚠️ API error mislabelled | ✅ |
| Synonym relevance | `llm_relevance.py:397-443` | Anthropic SDK | ✅ | gene+synonym, no pmid | ✅ | ✅ | ⚠️ same | ✅ |
| Benchmark route | `run_eval.py:908` | Azure Responses SDK | ✅ | ✅ + `attempt=1` | ✅ + fallback reason | ✅ | ✅ | ✅ forced |
| Benchmark curation | `run_eval.py:1032` | Azure Responses SDK | ✅ | ✅ + `attempt` + representation | ✅ | ✅ | ✅ both attempts | ✅ forced |
| Azure smoke | `smoke_azure_models.py:68` | LiteLLM | ✅ | ❌ | ❌ | ❌ | ⚠️ | ❌ |

Only the two benchmark routes satisfy the full contract. Event types in use: `representation_route_decision`, `paper_curation_decision`, `paper_relevance_decision`, `synonym_relevance_decision`, `paper_extraction_selection`, `variant_claim_verification_decision`, `tier2_relevance_decision`, `clinical_data_triage_decision`, `extraction_priority_triage_decision`.

**Hidden chain-of-thought:** correctly disclaimed at `utils/llm_trace.py:9-12`, `:350-354`, `:526-530`, `docs/LLM_TRACING.md:30-33`. No route requests, stores, or claims private chain-of-thought. The single honesty defect is P2-10; provider reasoning *summaries*, `tool_rationale`, `inclusion_rationale`, `count_rationale`, and `selection_policy` are legitimate explicit rationales, appropriately retained.

---

## Rejected / superseded prior claims

| Prior claim | Verdict |
|---|---|
| P3 "one pipe character routes prose through table validation → under-recovery" | **REJECTED as stated.** Routing is real (`:830-831`, `:841-842`), but all four probe shapes — including multi-number prose with a stray pipe (`Of the 3 probands, p.Leu552Ser was found in 12 carriers \| Table 2.`) — were still **ACCEPTED**, rescued by `_labeled_table_value` (`:689-698`). The impact claim is unsupported; the real defect runs the opposite way (P0-1). |
| P2 "installed-package recovery import fails" | **CONFIRMED but too narrow** → superseded by P1-5. |
| P2 "redaction misses `api-key`, `insttoken`, token keys" | **CONFIRMED and widened** — `x-api-key` (Anthropic's own header), `X-ELS-Insttoken`, `set-cookie`, `refresh_token`, `bearer`, `azure-ai-api-key`, plus `vars()`/`repr()` fallbacks. |
| P2 "count-recovery calls untraced or weakly scoped" | **CONFIRMED and widened** to the skip-extract manifest-rewrite defect and the missing `stage`. |
| P1 "recovered values lack role provenance" | **CONFIRMED and widened** — stale-provenance spurious quarantine, plus `trust_tier` DEFAULT `'trusted'` + `--skip trust-gate`. |
| Open question: "`penetrance_data` has no uniqueness constraint on `(variant_id,pmid)`" | **Refined.** `variant_papers` *does* have `PRIMARY KEY (variant_id, pmid)` (`migrate_to_sqlite.py:648`), so the fan-out comes only from multiple `penetrance_data` rows. `penetrance_data` has neither a constraint nor **any index** on the pair. Agreed: do not add `UNIQUE` before the cohort-identity question is settled — but add a non-unique index now (P3-28). |
| AGY "Step 3.55 can write 500,000" | Correctly rejected. Verified `:909-914`; the `Total alleles → 251384` probe rejected on the ceiling. |
| AGY "`.gitignore` does not ignore `.env.swp`" | Correctly rejected (`git check-ignore` → `.gitignore:42`). The physical file remains an open incident (P0-2). |
| AGY "`_same_variant` must run the position fallback after normalizer errors" | Agreed not a defect (`:327-334` deliberately `continue`s). |
| AGY "direct Azure vision calls bypass tracing" | Correctly rejected — all three wrappers use `capture_llm_call`. But the *replacement* finding is new: those entry points establish scopes with **no gene and no PMID** (P2-12). |
| AGY "JSON repair omitting effort is a correctness defect" | Agreed: not a correctness defect. Downgraded to P3-22. |
| Grok "full prompt retention is itself a defect" | Agreed: policy, not a defect. Retained as P3-29 because it now also inflates a multi-hundred-MB HTML and sits untracked in-tree. |
| `validate_position` true for unknown gene lengths | Corroborated; already tracked at `TASKS.md:86-88`. Not counted as new. |

---

## Checks performed

Full reads: `CLAUDE.md`, `TASKS.md`, the prior local Grok/AGY audit snapshot
(removed after its claims were accepted, rejected, or superseded here),
`docs/LLM_TRACING.md`, `utils/llm_trace.py`,
`scripts/build_llm_trace_manifest.py`, `scripts/build_llm_trace_html.py`,
`pipeline/count_recovery.py`, `scripts/recover_counts.py`,
`utils/llm_utils.py`, and `benchmarks/codex_paper_eval/README.md`. Full dirty
diff (all 22 modified) + inventory of all 11 untracked paths. Targeted reads:
`pipeline/steps.py` (thread pools, instance lifetimes),
`pipeline/trust_gate.py`, `harvesting/migrate_to_sqlite.py` schema,
`utils/retry_utils.py`, `pyproject.toml`, `config/settings.py`, `.env.example`.
Repo-wide enumeration of every model-invocation site and every
`utils.llm_trace` API use (three parallel read-only sweeps).

Offline probes: prose/table role grounding (10 quotes); redaction coverage (19 keys) + nested header maps; `_reasoning_paths` on `reasoning_tokens` 0 and 1024; `_trace_scope` duplicate-kwarg crash; 8-thread shared-caller attribution race; forged-trace vs `validate_trace_manifest` + `trace_index.jsonl` digest comparison + HTML banner; HTML size/RSS scaling with extrapolation; `write_recovered_counts` dry-run vs real fidelity; stray-pipe prose routing.

`pytest tests/unit/{test_llm_trace,test_llm_trace_html,test_count_recovery,test_settings,test_gvf_run_pipeline_wiring,test_codex_paper_eval}.py` → **128 passed in 103.38s**. `ruff check .` → clean. `git diff --check` → clean. `git check-ignore -v .env.swp`; `ls -la .env*`.

## Limitations

No network, no live provider calls, no live biomedical APIs — real envelopes were not observed; `reasoning_capture`, usage shapes, and redaction exposure were probed with synthetic envelopes matching the shapes the code already parses. No wheel was built and no `pip install .` performed; the packaging finding is static (`pyproject.toml` include list vs. import graph), corroborated by `genevariantfetcher.egg-info/top_level.txt` and the absence of any `genevariantfetcher*` dist-info in `.venv`. The full offline suite was not re-run (six focused files were); the prior review recorded `1519 passed, 1 skipped`. The HTML was not opened in a browser — the scalability finding is measured file size, per-record embedding, and the synchronous `JSON.parse` at `:801`, not an observed hang. Count-recovery accuracy was not re-measured.

## Prioritized punch list (for a second exact Opus 5 task)

1. **P0-1** Rebuild the role gate (structured role + locator; reject denominator grammars; fail closed on at least two prose candidates; widened negative context applied to the *table* branch and `_labeled_table_value`; remove the `field is None → True` short-circuit). Land all ten probe quotes as fixtures. Keep the flag off until remeasured.
2. **P0-2** Treat `.env.swp` as a credential incident: close session, delete, rotate, add `*.swp` to snapshot/archive exclusions.
3. **P1-3** Delete `BaseLLMCaller.llm_trace_context`; convert all five sites to `with llm_trace_scope(...)`; add a concurrent shared-caller attribution test.
4. **P1-4** Cross-check `trace_index.jsonl` write-time digests; include the index in manifest + `LOCK.json` + chmod set; make the banner state which guarantee it has.
5. **P1-5** Relocate the three helper groups into packaged modules; add an installed-wheel smoke test; make explicit-enablement import failure a stage failure.
6. **P1-6** Persist recovery role/evidence into the provenance the trust gate reads; insert as `quarantine`; refuse 3.55 when `trust-gate` is skipped.
7. **P1-7 / P2-12** Add `attempt`; emit `accepted_response_trace_id` (+ `repaired`, `discarded_trace_ids`); give continuation, adjudication, routing, count recovery, and the final check their own stage/component; pass gene+PMID into figure and pedigree scopes; add `count_recovery_decision` and `paper_final_check_decision` events.
8. **P2-11** Configure tracing once for the whole `run_gvf_pipeline` lifetime and in `scripts/recover_counts.py`; refuse to rebuild a manifest whose records carry a different `run_id`.
9. **P2-8** Shard the HTML per gene/PMID (or lazy-load bodies); cap embedded payloads with an explicit marker; log what was dropped.
10. **P2-9 / P2-10** Normalize `-` ↔ `_`, add token/header redaction with an allow-list, drop the `repr()` fallback; narrow `reasoning_capture` to reasoning content.
11. **P2-13/14/15** Make dry-run faithful; back up the DB before 3.55 and commit per paper under a `SAVEPOINT`; surface stats in `RUN_STATUS.json` + `RUN_REPORT.md`; make total failure a stage failure and the standalone script exit non-zero.
12. **P2-16/17** Provider-aware count-recovery resolver + key preflight; add both settings to `utils/provenance.py`; consume or delete `settings.count_recovery_enabled`; unify the benchmark's effort gate with `utils/llm_utils`.
13. **P2-18** Promote `db_to_predictions.py`/`score_and_compare.py` into `run_eval.py` so the enablement measurement goes through `validate_predictions` + `command_lock`.
14. **P3 batch** `_trace_scope` crash; triage effort + token floor; `TIER3_ADJUDICATOR_REASONING_EFFORT`; effort through JSON repair; `llm_relevance` API-error label; render `run_id`; document `GVF_LLM_TRACE_*`; add count recovery to `ARCHITECTURE.md` + `PROTOCOL_CHANGELOG.md`; guard the bare `next()`; index `penetrance_data(variant_id, pmid)`; plumb `paper_derived_only`; make `stats` JSON-serializable; fix the stale `vision_reasoning_effort` comment at `config/settings.py:751-755`.
