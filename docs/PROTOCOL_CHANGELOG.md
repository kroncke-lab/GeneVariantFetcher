# GVF Protocol Changelog

**Per-iteration record of what changed in the extraction/scoring protocol.**
One row per protocol-affecting change (normally one merged PR). This is the
change-anchored ledger; two sibling docs give the other views:

- `docs/RECALL_HISTORY.md` — the metric-anchored trajectory (what moved recall/MAE).
- `docs/ARCHITECTURE.md` — the *current* end-to-end protocol (models, steps, guards).

**Discipline (do this every iteration):** when a PR changes the protocol
(prompts, extraction/table logic, guards, trust gate, scoring, discovery,
acquisition, or model routing), add a row at the **top of the table** below with
the date, PR, category, a one-line summary of the change, and its measured or
expected recall/MAE effect. If it moved a headline number, also add the detail to
`RECALL_HISTORY.md`; if it changed the current protocol shape, update
`ARCHITECTURE.md`. Never rewrite past rows — supersede with a new one.

Categories: `extract` (prompt/extraction/table logic), `trust` (guards, trust
gate, final check), `gold` (gold sync/scoring scope), `acquire` (source/
supplement acquisition), `discover` (search/priority), `infra` (harness, CI,
tooling), `docs`.

## Iterations (newest first)

| Date | PR | Category | What changed | Recall/MAE effect |
|------|----|----------|--------------|-------------------|
| 2026-07-26 | — | extract | **Extraction prompts cut ~51% / ~33%.** Instruction overhead is paid per paper and competes with the paper for context: the pair spent ~10.3K tokens of boilerplate per extraction (`EXTRACTION_PROMPT` alone ~6.8K), which on the measured median call exceeded the source text actually sent. Rules are now stated ONCE in a shared `_CORE_RULES` block (the two prompts previously duplicated ~200 lines that had to be re-synced by hand — the original reason `TABLE_ATTRIBUTION_GUIDANCE` was extracted). Every load-bearing invariant is retained and still pinned: `count_provenance` keys, the full 9-value `count_type` enum, the cohort_total/screened_N/unknown→NULL rule, `fact_provenance` keys, per-observation locators, study-design metadata, and the full output schemas. Overhead 41,058→22,532 chars (~10.3K→~5.6K tok). | **UNMEASURED** — prompt text changed, so recall/MAE effect requires a benchmark *extraction* run (costs tokens), not just a rescore. Do `benchmarks/curated_extraction_eval` before any headline claim. |
| 2026-07-26 | — | trust | **Per-paper final check PARKED (Steps 3.8 + 3.9 default off).** One `@xhigh` reasoning call per paper cost more time and money than its measured effect justified for a step that only *records* findings — the ablation the 2026-07-20 three-model review said to do before optimizing its cost (`docs/PROTOCOL_COST_EVAL.md`, `project_final_check_triage_2026_07_20`). Parked as a **pair**: the composer alone can only refuse *stale* stored findings, and stale findings raise a stage failure demanding a "reviewer replay" that can never happen with the reviewer disabled (`cli/gvf_run.py:1854`) — every run would fail acceptance forever. Code, prompts, reason codes and tests retained and green; revive with `PAPER_FINAL_CHECK_ENABLED=1 PAPER_FINAL_CHECK_GATE_ENABLED=1`. New `tests/unit/test_final_check_parked.py` prevents silent re-enabling. | Removes the dominant per-paper LLM cost. Structural gold-free trust (Step 3.7) is unaffected; field trust already composed into existing DBs is preserved by Step 3.7. Model-found contradictions (incl. the new `attributed_to_other_study`) go dormant. |
| 2026-07-26 | — | trust / docs | **Attribution becomes enforceable.** Counts a paper credits to *another* publication (reference/citation column, footnote marker, "previously reported (ref 12)", compilation-table caption) had prompt guidance (`TABLE_ATTRIBUTION_GUIDANCE`, pinned by a test that only checks the prompt *text*) and no downstream check — the nearest code, `unsupported_count`, is deliberately advisory, so a compilation row could reach the trusted projection. Adds reason code `attributed_to_other_study` to the final-check vocabulary and to `ENFORCEABLE_REASON_CODES`, with per-code guidance drawing the boundary against `unsupported_count` (thin/unlocated support stays advisory). Reviewer prompts bumped `pfc7`→`pfc8`, `pfs12`→`pfs13`; compose policy `pfcg2`→`pfcg3`. New doc `docs/EXTRACTION_CONTRACT.md` (meta prompt) maps every extraction rule to the code that backs it, or records that none does. | No four-gene change expected (additive, source-grounded, high-severity-only). **Version bump invalidates stored final-check results** — papers re-check before enforcement resumes. Unmeasured until a rescore. |
| 2026-07-25 | — | acquire | **Linked-supplement recovery.** Supplement links the markup advertised were write-only (rendered as `_link_:` in FULL_CONTEXT, read by nothing) and the href was captured host-less, so they could never be fetched. Three parts: (1) `figure_extractor` takes a `base_url` and absolutises every href at capture — all 5 call sites pass the real page URL; (2) `harvesting/supplement_link_resolver.py` + `supplement_links` / `supplement_links_unfetched` in `artifacts.json`, feeding the existing `download_supplement()`; (3) `scripts/fetch_linked_supplements.py`, wired default-on into `gvf-run` Step 4 source recovery (self-gating: no-ops unless a paper has recorded links and an empty supplements dir). Route note: NCBI per-file `bin/` URLs now **403** unattended clients — the working route is the Europe PMC archive endpoint (one ZIP per article), with per-file as fallback. | Backlog sweep: **4,442 supplement files across 407 papers, 249 FULL_CONTEXT refolded** (the other 158 held only figure images). Recall effect not yet scored — needs a rescore. |
| 2026-07-20 | #167 | trust / infra | **Final-check triage prototype** (decision engine + no-LLM offline shadow report; three-way reviewed codex/grok/agy). Does NOT change the production auditor. Offline: routing zero-count-with-source papers to a cheap completeness lane could avoid ~80% of `@xhigh` calls (vs ~29% conservative) on canonical KCNH2 — pending 101-paper calibration. See `docs/PROTOCOL_COST_EVAL.md`. | No pipeline change (analysis + tooling only). |
| 2026-07-20 | #165 | trust / extract / discover | Study-design-aware counts + provenance honesty (coworker BRCA critique): trust_gate **tg3** (`negative_count`, full-partition arithmetic, field-scoped `implied_unaffected_zero`); `variant_papers.source_notation`; prompt decouples count from phenotype; penetrance/segregation discovery lane + priority signal; vertical parser stops asserting `pathogenic`. | No four-gene change *by design* (gold-free/additive; 0 new KCNH2 quarantine). Improves BRCA honesty. |
| 2026-07-20 | #164 | trust | Enforce source-grounded per-paper final checks (Step 3.8/3.9): returned quotes verified against the source; only quote-verified exact fact/field findings become enforceable; weak unsupported-count findings stay advisory. | No headline change; raises trusted-tier precision. |
| 2026-07-20 | #163 | gold | Versioned Azure gold snapshots + scoring tiers (cardiac/all/noncardiac); harden tier filtering and bulk exclusions; required-sync scoring reads the selected tier. | Scoring-scope change; headline unchanged. |
| 2026-07-20 | #162 | extract / trust | Generalize table-role validation (reject row-ID / population-frequency / cohort-denominator / clinical-measure columns by class; record selected column + count type). Includes multi-gene / borderless table regression fixes. | Precision/robustness; headline unchanged. |
| 2026-07-20 | #161 | gold | Live gold sync from the Variant_Browser Azure review API into a versioned SQLite cache (immutable snapshots, per-sync change log, reviewer/approver identity, reversible exclusions; native JSON normalized; fail-closed on raw/disputed/withheld/stale/checksum-invalid). | Gold-integrity; headline unchanged. |
| 2026-07-12 | #150 | infra / docs | Dual dashboards + no-gold run-health / worklist / delta cards + trust badges (`cli/dashboard.py`); status dashboard on GH Pages. | No extraction change. |
| 2026-07-12 | — | acquire | Four-gene idempotent supplement reconciliation + gated SCN5A land (per-`mmc` reconciliation; fold gap 289→0 papers). | uniqV 86.1%→**86.2%**, rows 80.8%→**81.2%**, carriers MAE **0.614**. |
| 2026-07-10 | #142 | trust | Trust gate v1 — gold-free soft-quarantine two-tier DB (`arith_inconsistent` / `count_is_total` / `population_count` / `paper_outlier`), default-on `gvf-run` Step 3.7. | Precision/honesty; no recall regression. |
| 2026-07-10 | #140 | trust / infra | Trust/validation hardening: SAVEPOINT rollback, honest metrics, fail-closed regression gate, non-zero exit + `RUN_STATUS.json`, cold-start eval, provenance widened. | No recall regression. |
| 2026-07-08 | #136 | infra | Turnkey on-ramp (`--email` sets NCBI_EMAIL; unified `pytest tests/unit`; Makefile; anthropic default provider). | No extraction change. |
| 2026-06-12 | — | extract / acquire | PDF-linearized table reconstruction + iter-2 quality gate + targeted four-gene lands. | uniqV 85.4%→**86.1%**, MAE 0.614–0.635. |
| 2026-06-05 | — | acquire | Elsevier `mmc` supplement acquisition + on-disk fold (the primary recall lever; figures gave 0). | uniqV 83.8%, rows 78.3%, carriers MAE 0.882. |
| 2026-05-29 | — | extract / trust | cDNA↔protein matcher bridge + count-outlier guard + caption-scoped tables. | uniqV 82.2%, MAE ~0.90. |
| 2026-05-26 | — | acquire / extract | Acceptance-gated source-replay sweep + per-PMID rollback + ClinVar/PubTator recovery (the single biggest jump). | uniqV 62.6%→**81.4%**. |
| 2026-05-21 | — | acquire | Vanderbilt Elsevier `insttoken` unblock (foundational; 242/246 paywalled Elsevier full texts recovered). | Enables the cardiac baselines. |
| 2025-11-18 → 2026-05 | — | infra | Pipeline build-out: discovery → harvest → Tier1/2 → Tier3 extraction → migrate → recall scoring; publisher routes; gold builders. | Pre-systematic-benchmark. |

For metric detail on any recall-moving row, see the matching dated entry in
`docs/RECALL_HISTORY.md`. For the current sample cost/quality of the protocol, see
`docs/PROTOCOL_COST_EVAL.md`.
