# Protocol cost & quality

## Experimental reason-class routing (2026-08-10; default OFF)

The highest-leverage fixed-48 cost lane was compact claim verification: 146
calls / 394,420 tokens in the 2026-08-08 replay. The A1 experiment changed only
its opening condition. Count-semantic/precision risks still open the verifier;
completeness-only, missing-count, and source blockers are classified for their
appropriate recovery lanes. Models, prompts, evidence cards, ranking, cap,
primary extraction, trust, and recovery behavior stayed fixed.

| Route | 2026-08-08 baseline | A1 | Change |
|---|---:|---:|---:|
| Sol claim verification | 146 calls / 394,420 tokens | 106 / 285,472 | **-27.6% tokens** |
| All fixed-48 LLM work | 279 calls / 1,372,842 tokens | 244 / 1,270,531 | **-7.5% tokens** |
| Summed provider-call time | 3,670.5s | 3,522.4s | -4.0% |

Variant quality stayed inside the predeclared gates: all-layer precision,
recall, and F1 moved -0.18pp, -0.20pp, and -0.20pp; paper-only F1 moved -0.28pp.
The raw carrier-MAE gate failed, however (0.723→0.902). A causal audit did not
locate that regression in the rerouted papers: the eight changed-decision papers
had identical TP/FP/FN and slightly lower total absolute carrier error. Because
the gate was declared on the aggregate run, the switch remains default off
pending a paired same-primary-output ablation or independent locked replicate.

The full 56-paper run added eight BRCA2 diagnostic entries and used 303 calls /
1,538,140 tokens in total. BRCA2 results are directional only because the
curator/derived `gold_overrides` are not fully manual headline gold. Reproduction
scripts, source/prediction hashes, exact route telemetry, and the full decision
record are in
`benchmarks/codex_paper_eval/runs/20260810_failure_routing_a1_56/`.

## Current traced fixed-48 measurement (2026-08-08)

The source-snapshot production replay supplied exact provider usage rather than
the proxies used in the older sample below. Across 48 cardiac papers it made
279 calls: 978,972 input tokens, 314,432 output tokens, and 1,372,842 total
tokens. The four gene jobs ran in parallel and completed in about 17 minutes
wall-clock; summed provider-call duration was 3,670.5 seconds.

| Route | Calls | Input | Output | Total | Summed call time |
|---|---:|---:|---:|---:|---:|
| Kimi table routing | 13 | 15,916 | 30,463 | 46,379 | 121.6s |
| Grok paper extraction | 43 | 567,021 | 196,948 | 843,407 | 2,030.9s |
| Sol claim verification | 146 | 325,510 | 68,910 | 394,420 | 1,029.5s |
| Sol figure reading | 76 | 68,519 | 17,948 | 86,467 | 485.6s |
| Sol paper adjudication | 1 | 2,006 | 163 | 2,169 | 2.9s |

Because this was an explicit `--pmid-file` benchmark, discovery plus Tier 1/2
relevance were skipped. Luna therefore had no legitimate Tier 2 work in this
run. The largest Sol lanes were source-grounded claim verification and figure
vision, both judgment-sensitive; replacing them merely because Luna is cheaper
would confound cost and quality. The clean first Luna A/B is Tier 2 relevance or
extraction-priority triage on a separate labeled discovery set, with all other
routes fixed. Terra is a possible midpoint. Pricing was not embedded in the
traces, so token totals are exact but dollar cost remains deployment-contract
dependent.

Quality results and lock hashes are recorded in
`benchmarks/codex_paper_eval/runs/20260808_fixed48_snapshot_replay/SUMMARY.md`.

## Historical small sample (2026-07-20)

> Historical configuration. This sample deliberately included Steps 3.8/3.9.
> The pair has been parked/default-off since 2026-07-26, so its timings and cost
> are evidence for that decision, not a description of the 2026-08-08 default
> pipeline. See `docs/ARCHITECTURE.md` for current routing.

**Purpose.** Measure that extraction protocol's **cost (time + money)**
and **quality** on a small sample, so we can decide whether to spend a full
cardiac re-extraction before running over all ~1,500 gold papers. This is a
**sample, not the full suite** — do not read these as headline recall (that lives
in `docs/RECALL_STATUS.md`). The protocol shape is in `docs/ARCHITECTURE.md`; the
change history is in `docs/PROTOCOL_CHANGELOG.md`.

## What was run

3 papers through the current default `gvf-run` protocol, **core LLM path only**:
extract (Kimi table-routing → grok-4.3 extraction → gpt-5.4/DeepSeek/Kimi debate)
→ migrate → trust gate (Step 3.7, tg3) → source-grounded final check (Step 3.8,
`gpt-5.6-sol@xhigh`) → composer (3.9). Run with `--pmid-file --no-source-recovery
--skip doctor --skip layers --skip source-qc --skip report` to **isolate LLM
cost** — i.e. no network acquisition, no recovery layers, no supplement folding.
Papers were read from the on-disk `corpus/`. Runs were **sequential** (one paper
at a time); production uses `MAX_WORKERS` / `AZURE_MAX_WORKERS` concurrency.

| Gene | PMID(s) | Kind | Wall-clock | LLM calls | Variants |
|------|---------|------|-----------:|----------:|---------:|
| BRCA2 | 10398279, 15365993 | typical cohort papers | 950s (**475s/paper**) | ~50 (**~25/paper**) | 118 |
| KCNH2 | 10973849 | large variant catalogue (247 variants) | **1322s** | 20 | 247 |

## Cost (time)

- **Typical paper ≈ 8 min** (475s), **~25 LLM calls**, full ensemble incl. the
  `gpt-5.6-sol@xhigh` final check.
- **Large paper ≈ 22 min** (1322s) — the cost tail. Big variant catalogues use
  compact extraction (fewer, larger calls: 20) but the `@xhigh` final check over
  many variants dominates.
- **Full-suite extrapolation** (~1,502 cardiac gold PMIDs), core LLM path only,
  at ~475s/paper: **≈ 198 hours sequential**, i.e. **≈ 198 / C hours** at
  concurrency `C` (e.g. ~25 h at C=8, ~12 h at C=16). Large-paper tail adds to
  this. The full **turnkey** run additionally spends non-LLM time on source
  acquisition/recovery, which this measurement excluded.

## Cost (money) — order of magnitude, needs real token accounting

The pipeline does **not** currently aggregate per-run token usage, so an exact
dollar figure is not yet defensible. Measured proxies: ~25 LLM calls/paper;
per-paper input ≈ 14k tokens (56k-char full text) fanned across grok-4.3
extraction, a 3-model debate, Kimi routing, and an `@xhigh` reasoning final
check. Order of magnitude is **~$0.5–2 / typical paper**, so the full ~1,500-paper
cardiac suite is roughly **~$1k–3k** — dominated by the `@xhigh` final check.
**Recommended before a full run:** add token/cost logging (LiteLLM exposes
`.usage`) so the next sample reports exact spend, and consider a cheaper or
sampled final-check tier to cut the dominant cost.

## Quality — does it improve things?

**The #165 honesty/provenance changes work, live:**

- **Traceability (criticism 7):** `source_notation` populated **365/365**
  (BRCA2 118/118 + KCNH2 247/247), including **legacy BIC notation** (e.g.
  `6174delT`) preserved beside the normalized cDNA — a curator can now trace and
  audit every normalized variant.
- **Study design captured:** BRCA2 papers tagged `cohort_population` /
  `population_screening` and `proband_referral`. (The large KCNH2 catalogue used
  compact mode and left design null — so the design-gated trust rules correctly
  stay dormant there.)
- **No fabricated penetrance (criticisms 1, 2, 6):** BRCA2 cohort penetrance rows
  report `total_carriers` with **NULL** affected/unaffected splits — the new
  prompt stopped the old "all carriers affected, unaffected=0" 100%-penetrance
  fabrication *at the source*. Trust gate: 8/8 trusted, 0 quarantined (no
  fabrication left to quarantine).

**Recall (sample, cardiac):** KCNH2 10973849 fresh **core-only** extraction hit
**93% unique-variant recall (55/59)** vs the canonical **100% (59/59)**. The
4-variant gap is the **supplement acquisition + recovery layers that were
deliberately skipped** to isolate LLM cost — **not** a regression from the
protocol changes. Independent confirmation the changes don't hurt cardiac:
`apply_trust_gate` on the canonical KCNH2 DB adds **0** new quarantine, and the
curated four-gene benchmark score is unchanged.

## Verdict / recommendation

- The current protocol is **~8 min and ~$0.5–2 per typical paper** (large papers
  several× that), with the `@xhigh` final check as the dominant cost.
- The #165 changes **improve honesty and traceability with no measured cardiac
  regression**, and demonstrably fix the BRCA fabrication the coworker flagged.
- **Before a full cardiac re-extraction:** (1) add real token/cost logging;
  (2) decide whether the `@xhigh` final check is worth its cost on every paper or
  should be sampled/downgraded; (3) then run the full suite with source recovery
  on. Until then, the headline in `RECALL_STATUS.md` stands unchanged.

_Sample artifacts (scratch, not committed): out_brca2 / out_kcnh2 DBs, per-gene
logs, and `timing.txt` under the session scratchpad._

---

## Final-check triage prototype (2026-07-20)

The `@xhigh` final check is the dominant cost and runs on **every** paper, so we
prototyped triaging it. Design was reviewed three ways — **codex (gpt-5.6-sol,
high), grok-4.5 (high), and agy/Antigravity (Gemini 3.1 Pro, high)** — which
converged strongly.

**Should it be interleaved per-paper (right after each paper's extraction)?
No** (all three). The check only *records* findings (it doesn't trigger
re-extraction), so moving it earlier changes timing, not efficacy; and it would
run before source recovery / recovery layers settle, forcing a second pass. (The
one concern that turned out **not** to apply: `trust_gate.paper_outlier` is
computed per-PMID, not cross-paper, so interleaving wouldn't corrupt it.) The
real wall-clock lever is **bounded final-check concurrency** (the calls are
currently sequential) + a **same-run targeted re-extraction loop** — both larger,
separate features, not a reorder. Keep batch `3.7 → 3.8 → 3.9`.

**The real lever is risk-triage, but never "clean → skip."** The check does two
jobs: a count-trust audit (mostly duplicated by the free trust gate) and a
**completeness** check (which reported carrier groups we missed). Completeness has
no free substitute and **does not correlate with extraction difficulty**, so a
clean paper demotes to a **cheap, completeness-preserving pass** (same prompt,
cheaper model/effort, escalate to `@xhigh` on any flag/miss/low-confidence),
never a skip. A stable ~10% random audit sample of cheap papers still runs full
to measure silent false negatives.

**Built (this iteration):** the pure decision engine `pipeline/paper_final_check_triage.py`
(`decide_final_check_tier` → full/cheap/skip + escalation predicate, unit-tested)
and a **no-LLM offline shadow report** `scripts/final_check_triage_report.py` that
runs the predicate over any migrated DB. Deliberately **no change to the
production auditor yet** — the experts' unanimous guidance is shadow-calibrate
first.

**Offline shadow results** (canonical KCNH2, 1,131 papers; the go/no-go artifact):

| Zero-count-with-source routing | @xhigh (`full`) | `cheap` | `@xhigh` calls avoided (pre-escalation) |
|---|---:|---:|---:|
| `full` (conservative, experts' default) | 71% | 29% | ~29% |
| `cheap` (aggressive) | 20% | 80% | ~80% |

The dominant cost driver is **papers that extracted zero counts but have source
(~55%)**, sent to full for completeness. Routing them to the cheap completeness
lane is the big lever (71%→20% full) — but only safe if the cheap pass reliably
catches missed groups on that bucket, which is the calibration question. On our
3 recent sample papers all went `full` (2 BRCA2 on `thin_provenance`, 1 KCNH2 on
`zero_counts_with_source`) — correct, since they're either sparsely-quoted or a
zero-count catalogue.

**Rollout gate before enforcing (per codex/grok):** shadow-calibrate on the
current curated staging set + all known count-role/missing-carrier failures; require
~100% escalation of known high-severity count errors, no missed quote-verified
carrier gap, and ≥50% projected `@xhigh` cost reduction. **Deferred:** the live
`apply_paper_final_check` wiring (shadow/enforce modes + cheap-tier + separate
screen storage), persisting the paper census to the DB (to sharpen the predicate),
bounded final-check concurrency, and the same-run re-extraction loop. A
pre-existing **completeness honesty bug** (long tables are head/tail sampled but
not marked `truncated`, so `completeness=complete` can survive a partial view;
`pipeline/paper_final_check.py` ~L602) must be fixed before enforce.
