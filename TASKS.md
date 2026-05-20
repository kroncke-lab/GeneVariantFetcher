# GVF Tasks

## Current Focus
- Improve honest cold-start recall for KCNH2, RYR2, and SCN5A to 90% across PMIDs, variant rows, unique variants, patients, affected, and unaffected counts.
  - Current metrics, highest-yield PMIDs, and next run plan are in `docs/CURRENT_RECALL_STATUS_2026-05-20.md`.
  - Do not duplicate live recall tables here; this file is only a short task checklist.
- **Measurement loop is multi-gene now** — score KCNH2, RYR2, and SCN5A from `gene_variant_fetcher_gold_standard/normalized/*_recall_input.csv`; KCNQ1 scoring remains available when its gold input is in scope.
- **KCNE1 is extraction-only until a gold input exists.** There is no `KCNE1_recall_input.csv` in the current gold-standard package, so KCNE1 recall cannot be claimed yet.

## Completed
- [x] Springer API integration (full-text retrieval working)
- [x] Wire `gene_literature/supplements/` UnifiedSupplementFetcher into orchestrator
- [x] Fix all test failures and collection errors (18 → 0)
- [x] Canonical form comparison for variant matching
- [x] Remove 50-paper download cap (was silently dropping 350+ papers)
- [x] Fail-open Tier 1 filter (min_keywords 2→1, no-abstract papers pass through)
- [x] Add source-completeness report (Step 3.5: `source_completeness.json`)
- [x] Add zero-variant QA flagging and re-extraction (GVF_QA_MODEL env var)
- [x] Improve cohort/screening table carrier count extraction prompts
- [x] Fix pytest environment (unit tests passing locally)
- [x] **Paywall recovery pipeline (2026-05-11)**
  - [x] `harvesting/browser_html/cookie_loader.py` reads Chrome cookies → Playwright format
  - [x] `harvesting/browser_html/authenticated_pool.py` injects cookies + uses `channel="chrome"`
  - [x] `harvesting/browser_html/content_quality.py` quality gate (paywall phrases + body floor + padding)
  - [x] `harvesting/browser_html/dom_extract.py` DOM walker + CF detection + `pick_better_markdown`
  - [x] Modernized strategies — AHA (`#bodymatter` + CF wait + reload), Elsevier (`/abstract→/fulltext` + imprint domains), generic (DOM walker fallback + CF wait), Wiley (CF wait)
  - [x] `scripts/fetch_paywalled.py` CLI with DOI resolution, per-domain pacing, CF retry sweep, **PMC fallback**
  - [x] Verified against 10 audit cases: 10/10 correctly classified by the gate
- [x] **`cli/compare_variants.py` matcher refinements (2026-05-11)**
  - [x] `_positions_compatible(a, b)` with subset rule (handles `p.Pro926AlafsTer14` (926, 14) vs `P926fsX` (926,))
  - [x] `_cdna_indel_protein_positions(cdna)` codon math (c.842dupG → aa 281)
  - [x] `_protein_indel_position(variant)` protein indel parser
  - [x] `find_best_match(..., consumed: Set[str])` for 1-to-1 enforcement at match time
  - [x] cDNA/protein bridge pass wired into greedy 1-to-1 matcher
- [x] Cherry-picked from heuristic-chatelet worktree
  - [x] Scanner: parenthesized HGVS pattern `p.(Arg176Trp)` (commit `f6307da`)
  - [x] Harvest: log silent supplement-scrape failures with PMID context (commit `d3a639f`)
  - [x] Extraction: reject DATA_ZONES that override full text with garbage (commit `e6daa24`)
- [x] **Clinical mutation-list table extraction (2026-05-14)**
  - [x] Deterministically infer one carrier per row for clinical mutation-list tables with variant/proband/family/patient context but no explicit count column
  - [x] Preserve affected vs unaffected counts when the row text indicates asymptomatic/control/unaffected status
  - [x] Keep pure nomenclature/list tables out of extraction unless clinical row context is present
- [x] **Bulk KCNH2 stale-context reharvest (2026-05-14)**
  - [x] Reharvested all 213 stale/missed KCNH2 FULL_CONTEXT files in `results/KCNH2/20260506_102238/pmc_fulltext/`
  - [x] Recovered 17 real contexts; 196 still require paywall/manual/source access
  - [x] Wrote audit artifacts under `results/KCNH2/20260506_102238/reharvest_out_20260514/`
- [x] **KCNH2 v12 manual-recovery score + matcher patch (2026-05-15)**
  - [x] Scored `KCNH2_v12_manual_recovery_20260515.db`: PMIDs 184/262 (70.2%), variant rows 542/991 (54.7%), unique variants 323/530 (60.9%), patients 1758/2674 (65.7%)
  - [x] Added frameshift canonicalization for extraction spellings like `fsTer`, `fs/185`, `fs+*49`, and malformed `AlaX14`
  - [x] Extended cDNA/protein bridge matching to multi-base cDNA indel ranges; recovered `P926fsX` (PMID 26496715) and `P1034fsX` (PMID 29622001)
- [x] **Turnkey closeout and honest enrichment cleanup (2026-05-18)**
  - [x] Added `gvf gvf-run` as the one-command doctor -> extract -> recovery layers -> report path.
  - [x] Scored KCNH2/KCNQ1/RYR2/SCN5A under `validation_runs/closeout_20260518_124343/`.
  - [x] Changed ClinVar/PubTator recovery layers to default to DB-observed PMIDs instead of gold PMIDs; gold-PMID enrichment is now explicit diagnostic mode.
  - [x] Turnkey recovery layers now back up the SQLite DB before mutation.
- [x] **Four-gene turnkey long run and cleanup pass (2026-05-19)**
  - [x] Ran KCNH2, KCNE1, RYR2, and SCN5A end to end under `validation_runs/turnkey_e2e_20260518_213934/`.
  - [x] Confirmed KCNH2, RYR2, and SCN5A remain below 90% on all scored recall metrics except RYR2 unaffected.
  - [x] Confirmed KCNE1 extraction completes but cannot be scored without a gold recall input.
  - [x] Removed default KCNH2 v12 auto-merge, removed gold-PMID leakage from default recovery, added DB backups, broadened LLM provider checks, filtered non-article LinkOuts, added retry/challenge handling, and guarded Data Scout against oversized raw contexts.

## Active Tasks
- [ ] **Close source/acquisition gaps to >90%** using the highest-yield PMIDs in `docs/CURRENT_RECALL_STATUS_2026-05-20.md`.
- [ ] **Investigate count semantics and cohort-table over-counting**. Preserve raw count columns and classify study-wide counts versus per-variant carriers before writing patient/affected/unaffected totals.
- [ ] **Create or import KCNE1 per-PMID gold input** before making KCNE1 recall claims.
- [ ] **Obtain Elsevier INSTTOKEN** from Vanderbilt E-resources (eresources@vanderbilt.edu). Exact missing-PMID ranks live in `docs/CURRENT_RECALL_STATUS_2026-05-20.md`; do not carry older per-PMID gap counts here.
- [ ] **Integrate valid paywall recovery artifacts into canonical runs** before re-extraction. Keep paper-specific recovery lists in validation artifacts or `docs/CURRENT_RECALL_STATUS_2026-05-20.md`, not in this checklist.
- [ ] **Manual PDF download for hard-blocked/stub PMIDs** only after the current-status source audit confirms the artifact is still missing. Feed real downloaded PDFs through `harvesting/format_converters.py`; avoid OneDrive on-demand placeholders.
- [ ] Set weekly recall cadence (Friday re-run + compare) — produces trajectory chart for grant.

## Blocked
- [ ] **Karger access** — Cloudflare interstitial doesn't clear even with 30+ cookies + `--no-headless`. TDM request status unknown.
- [ ] **Wiley TDM full-text API** — `WILEY_API_KEY` in `.env` is revoked. Without it, Human Mutation (10.1002/humu.*) is unreachable.
- [ ] **Sage/Liebert (PMID 23631430)** — Sage's CF instance refuses Playwright fingerprint even with the right cookies.

## Backlog
- [ ] Create a KCNE1 gold test set and keep all non-KCNH2 gold inputs source-reconciled.
- [ ] Re-run extraction with regex disabled to measure regex vs LLM contribution.
- [ ] Expand quality gate test set beyond the 10 audit files.
- [x] Parameterize `scripts/recall_recovery/ingest_clinvar.py` and `ingest_pubtator.py` so cold-start genes can run the same recovery layers KCNH2 uses.

See `docs/CURRENT_RECALL_STATUS_2026-05-20.md` for current status and `CLAUDE.md`
for recovery architecture and handoff details.
