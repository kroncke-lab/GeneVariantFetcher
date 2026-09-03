# Protocol cost & quality

This is a dated measurement record. Model routes, prices, and default-on stages
are those used by each named run; they are not current configuration guidance.
Use `config/settings.py` and `docs/ARCHITECTURE.md` for current defaults and
`TASKS.md` for the next acceptance gate.

## Mixed all-gold tranche budget (2026-09-03; no extraction calls)

The deterministic registry at
`benchmarks/evaluation_tiers/mixed_gold/registry.json` assigns all **1,422
source-available** gene-paper attempts from the 1,534-attempt named-variant
inventory to 29 article-disjoint tranches. The remaining 111 source-unavailable
attempts and one quarantined wrong-paper attempt stay explicit in
`inventory.tsv`; they are not scored as misses.

Observed per-gene token usage gives a **$4.98–$5.30 one-arm proxy per tranche**
and **$9.96–$10.61 for a paired frozen-baseline/candidate comparison**. With 25%
retry/variance headroom, budget **$12.45–$13.26 per paired tranche**. Running
every tranche once would be $148.83 before headroom; paired execution of the
entire suite would be $297.66, or **$372.08 with headroom**. These are planning
proxies from the dated 2026-08-24 public price card, not invoices or tail bounds;
source acquisition, human review, and consultation are excluded. Compare the
first paid tranche's actual charge with the proxy before relying on it.

The primary endpoint is paper-derived micro variant-identity recall, with
paper-derived precision as a non-inferiority guardrail and gene/gold-provenance
breakouts always reported. Attempts are clustered by PMID. ClinVar/PubTator are
secondary linkage diagnostics; ambiguous origins are unscored. Inspecting a
score burns that tranche for confirmation, so consume tranches in registry
order and confirm changes on a still-unopened tranche. Setup and score enforce
the paired arm/order contract through the append-only `consumption_log.jsonl`.

The paired decision is preregistered before any arm is opened. Deltas are
candidate minus baseline on the same tranche. Discovery passes only if observed
paper-derived recall improves by at least 1 percentage point, its one-sided 95%
PMID-cluster-bootstrap lower bound is at least -1 point, and the corresponding
precision lower bound is at least -2 points. The bootstrap uses 10,000
deterministic resamples. A discovery pass only advances to confirmation: the
same rule must pass on the next unopened tranche before accepting the change.
These are engineering regression thresholds, not clinical-effect thresholds.
Use `run_eval.py compare`; do not choose a rule after viewing scores.

Gold-row denominators vary from 110 to 593 per tranche, and rare non-cardiac
provenance strata are sparse. The paired comparison remains valid within a
tranche, but per-provenance results must show their sample size and remain
diagnostic. Absolute pooled precision is not an exhaustive scientific estimate
for non-exhaustive collaborator gold. Until the ledger supports an explicit
abandonment event, do not open a baseline unless its candidate arm will be
completed.

## Current measured costs and $100 improvement budget (2026-08-24)

The accepted immutable gold-118 baseline is
`benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118`: 579
calls (569 successful), 1,951,910 input tokens, 2,730,460
total-minus-input-inclusive tokens, and a **$11.241 public-list-price proxy**
($0.0953 per attempt). By role, the proxy is $0.236 Kimi routing, $2.296 Grok
4.3 extraction, and $8.708 GPT-5.6 Sol vision/verification. Because the
provider's visible output and reasoning accounting differ for Grok, the
conservative alternate total is about $11.85. The four gene jobs completed
concurrently in 43.4 minutes;
summed provider duration is not wall time.

The active improvement envelope excludes that pre-existing baseline and charges
the following attributable work:

| Item | Calls (successful) | Tokens | Charged proxy |
| --- | ---: | ---: | ---: |
| Targeted exploration with incomplete refresh telemetry | — | — | $10.00000000 |
| `20260824_grouped_structural_gold118` | 580 (557) | 2,692,803 | $11.86999755 |
| `20260824_source_recovery_gold118` | 547 (537) | 2,626,745 | $11.19883675 |
| AHA PMID 15466642 production micro replay | 9 (9) | 35,887 | $0.20972125 |
| `20260824_aha_table_sourcebound_gold118` | 579 (569) | 2,730,460 | $11.24050705 |
| **Used** |  |  | **$44.51906260** |
| **Remaining from $100** |  |  | **$55.48093740** |

The source-recovery split was Kimi K2.6 26/26 calls and $0.26288925,
GPT-5.6 Sol 400/396 and $8.12796000, and Grok 4.3 121/115 and $2.80798750.
The promoted AHA lock split was Kimi 25/25 and $0.23640455, GPT-5.6 Sol
432/429 and $8.70826500, and Grok 4.3 122/115 and $2.29583750.
Packet generation, scoring, source inspection, and tests were deterministic and
made no model calls. Grok/AGY consultation CLIs exposed no billing telemetry;
their account cost is excluded, not asserted to be zero.

The completed historical 146-paper BRCA/BMPR2 run cost $23.664 by the same
method. Reconstructing the exact later 50/50/50 candidate from its available
traces gives an **audit estimate** of 1,026 calls, about 5.08M provider tokens,
and **$25.62** ($4.10 BMPR2, $10.28 BRCA1, $11.24 BRCA2). That estimate is not
an immutable cost manifest, and the cohort has no exhaustive gold standard;
therefore it is a workload/cost measurement, not a quality/cost benchmark.

The first new blind lock improved count error but regressed recall. The second
restored accepted recall and lowered carrier MAE, but failed carrier coverage
and count-bearing precision floors. A separately source-backed AHA repair then
received one authorized final lock and passed every gate: 554/283/78, 87.66%
recall, 66.19% raw precision, carrier coverage 207, and carrier MAE 0.193.

No further gold-118 best-replicate run is authorized. The highest-leverage
remaining FN source, KCNH2 PMID 29650123, has no mutation roster in its only
publisher supplement, and another stochastic pass cannot repair that. The next
$0-model-call step is human curation of the preregistered 90-paper
BRCA1/BRCA2/BMPR2 calibration set. After calibration changes freeze, release
the 60-paper holdout and score it once. Spend the remaining $55.48 only on a
calibration-informed, pre-frozen candidate or on one holdout extraction if the
frozen existing DB is not the candidate.

> **Superseded as a next step on 2026-08-25.** `TASKS.md` descoped the
> BRCA1/BRCA2/BMPR2 150-paper curation: it is not being run, the cardiac
> **gold-120** set is the human-curated standard, and exact-150 metrics stay
> undefined. The dollar envelope above is still current; the paragraph's
> *next-step assignment* is not. `TASKS.md` owns the next gate.

## Canonical rollout sizes (as measured 2026-08-13)

> These are the 2026-08-13 tier sizes, retained because the cost figures in this
> document are derived from them. They are **not** the current tier list —
> `benchmarks/evaluation_tiers/README.md` and its `registry.json` are the count
> authority, and the counts have since changed.

The active progression at that date was 50 gold-scored gene–paper attempts, then
120 manually curated cardiac-gold attempts (116 unique PMIDs), then the
546-attempt full reviewer backlog (507 unique PMIDs). Cross-gene papers are
intentionally processed once per gene–disease workspace, so provider cost scales
with attempts, not unique PMIDs.

The patched gold-120 revalidation supplies the current exact trace-derived
measurement. The 120 attempts used 527 calls, 2,351,247 tokens, and 7,903.6
summed provider-seconds (131.7 minutes). The four gene jobs ran concurrently and
completed in 34.4–38.9 minutes each. Public list-price arithmetic gives **$9.77
total / $0.0814 per attempt**.
This is a cost proxy, not an Azure invoice: Grok and Kimi use their public xAI
and Fireworks prices because the deployment-specific Azure rates are not in the
trace, cached-input discounts are ignored, and provider-returned total-minus-
input is conservatively charged as output/reasoning.

| Gene (30 attempts each) | Summed provider time | Calls (successful) | Tokens | Cost proxy | Cost / attempt |
| --- | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 33.21 min | 108 (106) | 515,098 | $2.365 | $0.0788 |
| KCNQ1 | 33.02 min | 166 (162) | 587,438 | $2.255 | $0.0752 |
| RYR2 | 30.25 min | 132 (132) | 604,475 | $2.671 | $0.0890 |
| SCN5A | 35.24 min | 121 (121) | 644,236 | $2.482 | $0.0827 |
| **Total / mean** | **131.73 min** | **527 (521)** | **2,351,247** | **$9.774** | **$0.0814** |

| Model role | Calls (successful) | Input tokens | Output/reasoning tokens | Provider time | Cost proxy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kimi K2.6 table routing | 19 (19) | 34,671 | 33,799 | 199.2 s | $0.168 |
| Grok 4.3 primary extraction | 115 (114) | 1,165,102 | 452,388 | 4,571.4 s | $2.587 |
| GPT-5.6 Sol verification/figure vision | 393 (388) | 517,604 | 147,683 | 3,133.0 s | $7.019 |

The proxy uses $5/M input and $30/M output for
[GPT-5.6 Sol](https://azure.microsoft.com/en-us/blog/gpt-5-6-now-available-in-microsoft-foundry/),
$1.25/M and $2.50/M for
[Grok 4.3](https://docs.x.ai/developers/models/grok-4.3), and $0.95/M and $4/M
for [Kimi K2.6](https://docs.fireworks.ai/serverless/pricing). Actual Foundry
cost depends on deployment, region, caching, and commercial terms.

The fixed 146-attempt experiment is now complete, so the earlier cardiac-only
projection can be compared with an actual BRCA-heavy run. It used 972 calls,
4,261,341 tokens, 13,362.1 summed provider-seconds, and a **$23.664** cost proxy.
The three genes ran concurrently and completed in 70.7--97.8 minutes. The large
BRCA tables make this **$0.1621 per attempt**, about twice the cardiac-gold
mean; gene/paper density is therefore a load-bearing planning variable.

| Experimental gene | Attempts | Calls (successful) | Tokens | Provider time | Cost proxy |
| --- | ---: | ---: | ---: | ---: | ---: |
| BMPR2 | 50 | 195 (195) | 1,023,662 | 0.730 h | $4.275 |
| BRCA1 | 50 | 407 (405) | 1,559,111 | 1.439 h | $9.112 |
| BRCA2 | 46 | 370 (370) | 1,678,568 | 1.543 h | $10.278 |
| **Total** | **146** | **972 (970)** | **4,261,341** | **3.712 h** | **$23.664** |

The completed 120-paper manual-gold comparison already lies inside the
requested 100--150 comparison range; no experimental gene was widened.

Straight-line calibration bounds are:

| Scope | Attempts | Cardiac-like proxy | BRCA-heavy proxy | Planning allowance |
| --- | ---: | ---: | ---: | ---: |
| Completed experiment | 146 | $11.89 / 2.67 h | **$23.66 / 3.71 h actual** | observed 97.8 min concurrent wall |
| Exact full reviewer tier | 546 | $44.47 / 9.99 h | $88.50 / 13.88 h | roster estimate $56–$72; budget $70–$90 and 10–20 elapsed h |
| Generic 500–600 attempts | 500–600 | $40.72–$48.87 / 9.15–10.98 h | $81.04–$97.25 / 12.71–15.25 h | roughly $50–$120 before human review |

The exact 546 roster has seven measured or directly scaled strata totaling
about $39.95; applying the two observed rates to its four unmeasured 50-paper
genes gives the narrower $56–$72 roster estimate. The allowance covers retry
tails and source/QC variability; it does not include human curation time.
Source recovery can add substantial elapsed time without adding much LLM spend.
The patched Gate 2 revalidation passed its like-for-like counted-extra precision
criterion (95.70% raw / 95.87% trusted versus the 77.3% floor). Detailed
experimental cost, trust, source, and collaborator-gold QC is locked in
`runs/20260813_experimental_146/EXPERIMENTAL_COST_AND_QC.md`.

## Collaborator-gold scope correction (2026-08-11; no new calls)

Active comparison/review scope is now the cardiac 48 plus the two BRCA2 papers
with lead-approved Variant Browser adjudications by Nate. The six internally
derived BRCA2 papers were excluded by PMID from active publishing and metrics.
The 50-paper metrics are a deterministic projection of the locked 56-paper
predictions, so this correction cost **zero additional model tokens and zero
provider calls**. Historical production and Luna costs below remain attributed
to the full 56-paper diagnostic run and must not be presented as the marginal
cost of the active 50-paper cohort.

## Luna count-semantics failure route (2026-08-10; shadow)

The exact A1 56-paper predictions were held locked while the scorer and answer
key were audited. Compact `azure_ai/gpt-5.6-luna@xhigh` claim cards targeted the
few multi-cohort rows dominating carrier MAE.

The underlying A1 production trace used the then-current failure-routing shape, not
Luna: 16 Kimi table-routing calls, 55 Grok 4.3 extraction calls, 90 GPT-5.6 Sol
figure-reading calls, and 142 GPT-5.6 Sol risk-ranked claim-verification calls
(303 provider calls total across 56 papers). The seven Luna calls below were an
additional shadow adjudication of selected ambiguities. This separation matters
when attributing either quality or cost to the production protocol.

| Measure | Broad missing-count probe | Compact semantics cards |
|---|---:|---:|
| Attempted calls | 6 | 7 |
| Completed/grounded gap checks | 162 / 0 | n/a |
| Total tokens | 153,010 | 27,682 |
| Wall time in traces | 525.7s | 80.9s |
| Useful source decisions | 0 | 4 dominant count-scope decisions |

The compact route used about 82% fewer tokens than the already-stopped broad
probe, but these are diagnostic workloads rather than an apples-to-apples model
efficiency benchmark. Source-backed answer-key/scorer repair reduced exact
56-paper carrier MAE from 0.8148 to 0.0794 after a subsequent blind Grok/AGY/
Claude source audit confirmed the headline rows and six control corrections.
All 378 predicted carrier observations remain; count recall is 34.05%. This is
a benchmark-quality result, not yet a production routing change. No Anthropic
model participated in the original Luna run; Claude Fable 5 Max was used only
in the later independent review. Full metrics and locked digests are in
`benchmarks/count_semantics_eval/runs/20260810_luna_xhigh_56/`.

## Historical small sample (2026-07-20)

**Purpose.** Measure the then-current extraction protocol's **cost (time + money)**
and **quality** on a small sample, so we can decide whether to spend a full
cardiac re-extraction before running over all ~1,500 gold papers. This is a
**sample, not the full suite** — do not read these as headline recall (that lives
in `docs/RECALL_STATUS.md`). The protocol shape is in `docs/ARCHITECTURE.md`; the
change history is in `docs/PROTOCOL_CHANGELOG.md`.

## What was run

3 papers through the then-current Azure `gvf-run` route, **core LLM path only**:
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

The pipeline did **not** aggregate per-run token usage in this sample, so an exact
dollar figure is not yet defensible. Measured proxies: ~25 LLM calls/paper;
per-paper input ≈ 14k tokens (56k-char full text) fanned across grok-4.3
extraction, a 3-model debate, Kimi routing, and an `@xhigh` reasoning final
check. Order of magnitude is **~$0.5–2 / typical paper**, so the full ~1,500-paper
cardiac suite is roughly **~$1k–3k** — dominated by the `@xhigh` final check.
**Recommendation recorded with the sample:** add token/cost logging (LiteLLM
exposes `.usage`) before budgeting a full run and consider a cheaper or sampled
final-check tier.

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

- The measured 2026-07-20 route was **~8 min and ~$0.5–2 per typical paper** (large papers
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

The `@xhigh` final check was the dominant measured cost and ran on **every** paper, so we
prototyped triaging it. Design was reviewed three ways — **codex (gpt-5.6-sol,
high), grok-4.5 (high), and agy/Antigravity (Gemini 3.1 Pro, high)** — which
converged strongly.

**Should it be interleaved per-paper (right after each paper's extraction)?
No** (all three). The check only *records* findings (it doesn't trigger
re-extraction), so moving it earlier changes timing, not efficacy; and it would
run before source recovery / recovery layers settle, forcing a second pass. (The
one concern that turned out **not** to apply: `trust_gate.paper_outlier` is
computed per-PMID, not cross-paper, so interleaving wouldn't corrupt it.) The
real wall-clock lever was **bounded final-check concurrency** (the prototype's
calls were sequential) plus a **same-run targeted re-extraction loop** — both larger,
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

**Historical rollout gate before enforcing (per codex/grok):** shadow-calibrate
on the canonical first-tier manifest plus all known count-role/missing-carrier failures; require
~100% escalation of known high-severity count errors, no missed quote-verified
carrier gap, and ≥50% projected `@xhigh` cost reduction. **Deferred:** the live
`apply_paper_final_check` wiring (shadow/enforce modes + cheap-tier + separate
screen storage), persisting the paper census to the DB (to sharpen the predicate),
bounded final-check concurrency, and the same-run re-extraction loop. A
pre-existing **completeness honesty bug** (long tables are head/tail sampled but
not marked `truncated`, so `completeness=complete` can survive a partial view;
`pipeline/paper_final_check.py` ~L602) must be fixed before enforce.
