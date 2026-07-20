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
