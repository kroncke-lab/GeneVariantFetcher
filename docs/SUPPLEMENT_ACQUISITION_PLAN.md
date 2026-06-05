# Supplementary-Mutation-Table Acquisition Plan

Generated 2026-06-05 from a 6-way recon + synthesis workflow that re-verified
every load-bearing claim against the live repo and the 4 gold cardiac genes
(KCNH2/KCNQ1/SCN5A/RYR2). This file is the durable source of truth for the
supplement-acquisition effort; it is safe to compact the session against it.

## MEASURED RESULTS (2026-06-05 execution)

Executed Phase 0 (T1/T2/T3) + T5. **Measured recovery: +175 unique gold-missing
variants** (~10% of the 1,673-variant gap), per-paper (fold → Tier-3 re-extract →
match vs the gold-missing list; disjoint paper sets):

| Source | Papers | Recovered / missing |
|---|---|---|
| T4 no-network fold (on-disk unparsed supplement) | SCN5A 22840528 | 8 / 8 |
| T5 Elsevier mmc fetch — single biggest blocker | SCN5A 29325976 | 64 / 87 |
| T5 Elsevier mmc fetch — batch | 16 papers | 103 / 126 |
| **Total** | **18 papers** | **175 / 221 (79%)** |

Elsevier mmc batch hit rate: 16/59 addressable papers actually carry mmc
supplements (KCNH2 1/14, KCNQ1 3/7, SCN5A 8/30, RYR2 4/8); the rest have no mmc
referenced in the authenticated XML. Where a supplement exists, ~79% of its
missing variants are recovered.

### OFFICIAL re-score (corpus→flat bridge, `scripts/corpus_to_harvest.py`)

Built the bridge and re-scored via `refresh_run_db` (re-extract the 18 forced
supplement papers, acceptance-gated) → `run_recall_suite`. To isolate the
supplement effect, BOTH DBs were rebuilt identically from extractions with
`--skip-recovery` (no ClinVar/PubTator/figure layers), so the delta is purely
supplements. Absolute numbers are therefore below the layered published baseline
(82.2% uniqV); the supplement lift STACKS on top of the full pipeline.

| Metric | Before (rebuilt) | After (+supplements) | Δ |
|---|---|---|---|
| unique_variants | 75.1% (2262) | 77.2% (2325) | **+63 (+2.1pp)** |
| variant_rows | 63.4% (4331) | 66.6% (4551) | **+220 (+3.2pp)** |
| patients | 66.8% | 69.8% | +576 |
| affected | 65.1% | 67.6% | +306 |
| carriers MAE | 0.915 | 0.887 | better |

Per-gene unique_variants: SCN5A 75.1→78.1 (+3.0pp, rows +4.9pp), RYR2 76.9→79.0
(+2.1pp), KCNQ1 77.0→78.9 (+1.9pp; PMID 19716085 correctly regression-gated:
re-extract gave 0 vs prior 412), KCNH2 70.8→71.1 (+0.3pp). The official
gene-deduped uniqV gain (+63) is below the per-paper proxy (+175) because uniqV
dedups across papers and the stricter official matcher is used; variant_rows
(+220) tracks the per-paper recovery more closely.

**Status:** official supplement lift confirmed (+63 uniqV / +220 rows from 18
papers, isolated). Done: T1, T2, T3, T5, corpus→flat bridge + re-score.
Remaining: T6 (Wiley/Springer API supplements, ~11% of gap), T7 (max_supplements
cap), T8/T9 (Karger/Sage, deferred ~0 value); and a full-pipeline refresh (with
recovery layers) to land the supplement variants in the published DBs.

## TL;DR — the premise was wrong; here is where the value actually is

The task was framed as "get Karger/Sage/Cloudflare-blocked supplements." The
data says **those publishers carry essentially none of the recall gap**:

- **Karger [10.1159] = 5 missing variants (0.3%). Sage [10.1177] = 0 (0.0%).**
  Chasing the Cloudflare-blocked publishers recovers ~5 variants. It is a
  correctness/observability concern, not a recall lever.

The real supplement value is in two places:

1. **A structural wiring bug (no network, free recall):** the re-fold layer
   that turns on-disk `{pmid}_supplements/` into LLM-visible text exists
   (`harvesting/supplement_fold.py`) but is wired ONLY into the standalone
   `scripts/fold_supplements.py` — `refresh_run_db.py`/`gvf-run` never re-read
   supplements. So every refresh silently drops them.
2. **Publisher full-text APIs fetch ZERO supplements:** Elsevier/Wiley/Springer
   API clients return body markdown only; `fetch_paywalled.py` Wiley/Springer
   fallbacks literally hardcode `supplements_downloaded=0`. A paper recovered
   via API lands body-only, no supplement table.

## Verified numbers

### Supplement-reachable recall ceiling = 132 / 1673 missing gold variants (7.9%)
(gold-missing variants whose bearing-PMID is `advertised_not_converted`)

| Gene | reachable | of missing | PMIDs |
|---|---|---|---|
| KCNQ1 | 121 | 389 (31.1%) | 8 — 30758498(57!), 22949429(19), 21956039(11), 24667783(10), 17470695(10), 19841300(9), 19716085(4), 27114410(1) |
| RYR2 | 11 | 224 (4.9%) | 2 — 28404607(10), 19216760(1) |
| KCNH2 | 0 | 172 | — |
| SCN5A | 0 | 888 | — |
| **TOTAL** | **132** | **1673 (7.9%)** | |

Caveats: theoretical ceiling (the Tier-3 LLM still has to extract from the
supplement, so realized ≤ 132); fragile — 57/121 of the KCNQ1 ceiling is one
paper (PMID 30758498, `10.1001/jamacardio.2018.4925`).

### No-network on-disk reparse opportunity = small (the old "~431 rows" MEMORY is STALE)
Post-2026-06-04 corpus, **1170/1183 (98.9%)** gold papers with a table-bearing
supplement on disk are ALREADY folded into FULL_CONTEXT. Only **13 genuinely
unparsed** (KCNH2 2, KCNQ1 6, SCN5A 4, RYR2 1) — ~857 pre-dedup variant-like
tokens, dominated by KCNQ1 38489124 (634). Deterministic refold = 0 new
variants; realizing the surface needs an LLM re-extraction pass.

### Publisher distribution of the FULL recall gap (1673 variants, 82.4% DOI-resolved)
Elsevier 515 (30.8%, already insttoken-unblocked → a REFRESH problem) ·
AHA/LWW [10.1161] 167 (10.0%, largest NEW paywall) · Wiley 119 (7.1%, TDM works
but API drops supplements) · OUP 76 (4.5%) · Springer 60 (3.6%) · NPG 56 (3.3%)
· **Karger 5 (0.3%) · Sage 0 (0.0%)**. Stub/0-byte acquisition gaps = 145
(8.7%, 127 SCN5A) — that's missing *full text*, not supplements.

## Task plan (sequenced highest-value/lowest-risk first)

### PHASE 0 — no network, free recall, do first
- **T1** (low): QC counter — per gene, count PMIDs with a `{pmid}_supplements/`
  whose converted markdown is NOT in FULL_CONTEXT. Extend
  `scripts/recall_audit/summarize_acquisition_outcome.py`; reuse
  `supplement_fold.build_supplement_markdown` + the `GVF_FOLDED_SUPPLEMENTS`
  sentinel. Verified baseline: 13 gold papers.
- **T2** (low): Fix `supplement_fold.build_supplement_markdown` (lines ~88-92)
  top-level `iterdir()` → recursive `rglob` over convertible suffixes (still
  exclude `.zip`), so files inside extracted-zip subdirs fold. Recovers 2
  nested-zip gold papers (KCNQ1 23098067, SCN5A 26173111). Add a unit test.
- **T3** (med, deps T2): Wire `fold_supplements_into_full_context` into
  `scripts/refresh_run_db.py` (before the `*_FULL_CONTEXT.md` discovery glob at
  ~349-363) and the `gvf-run` replay; for any PMID whose FULL_CONTEXT grew,
  force DATA_ZONES regen (extraction prefers DATA_ZONES when `use_condensed`).
  Non-destructive (`.pre_fold_bak`, sentinel-delimited, idempotent).
- **T4** (net/LLM, deps T3): Re-extract + re-score the 13 unparsed gold papers.
  Measures realized vs the ~857-token bound. Acceptance-gated refresh, no blind
  sweep.

### PHASE 1 — per-publisher supplement fetch, ordered by gap share (the genuine levers)
- **T5** (med, deps T3): Elsevier API supplement discovery+download — largest
  bucket (30.8%; single biggest blocker SCN5A 29325976 `mmc1.docx` = 87
  variants). `harvesting/elsevier_api.py:xml_to_markdown` (~230) emits body only;
  parse + download mmc refs; make `fetch_paywalled.py:try_elsevier_api_fallback`
  (~412-424) actually fetch mmc, not just pass through. Reuse
  `gene_literature/supplements/elsevier_fetcher.py` mmc pattern.
- **T6** (med, deps T3): Wiley + Springer API supplement download — fix the
  hardcoded `supplements_downloaded=0` in `fetch_paywalled.py` Wiley (~470-471)
  and Springer (~514-515) fallbacks; or register Wiley/Springer fetchers in
  `gene_literature/supplements/unified.py` (only PMC/Elsevier/CrossRef today).
- **T7** (low, deps T3): Raise/parameterize `max_supplements=12`
  (`paywall_context_enrichment.py` ~276,357) → e.g. 40 / `GVF_MAX_SUPPLEMENTS`;
  keep the 25 MB per-file cap.

### PHASE 2 — blocked publishers (<0.4% of gap; low priority)
- **T8** (low): Observability — classify Karger (10.1159) as
  `blocked_karger_cloudflare` in the recall-audit `infer_publisher`; add a Sage
  (10.1177) pacing branch in `fetch_paywalled.py:_domain_for`. No recall lift.
- **T9** (high/net, deps T8): DEFER. Cloudflare unblock for Karger/Sage via
  Vanderbilt EZproxy URL-rewrite (`ezproxy.library.vanderbilt.edu/login?url=...`)
  — the cookie is loaded but never used as a proxy prefix. Only if the free
  EZproxy path clears it; do NOT invest in a CF arms race for 0.4% of recall.

Every task ends with `scripts/refresh_run_db.py` (acceptance-gated) +
`scripts/run_recall_suite.py --score` vs
`/tmp/gvf_fig_enrich/baseline/<GENE>/missing_in_sqlite.csv` to measure realized
vs theoretical delta.
