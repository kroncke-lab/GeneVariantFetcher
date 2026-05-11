# GVF Tasks

## Current Focus
- Improve unique variant recall from 59% (2025-12-11 baseline) → 90% by June 2026 grant deadline.
- **Measurement loop is broken — re-run KCNH2 end-to-end is the single most important next action.**

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
- [x] Cherry-picked from heuristic-chatelet worktree
  - [x] Scanner: parenthesized HGVS pattern `p.(Arg176Trp)` (commit `f6307da`)
  - [x] Harvest: log silent supplement-scrape failures with PMID context (commit `d3a639f`)
  - [x] Extraction: reject DATA_ZONES that override full text with garbage (commit `e6daa24`)

## Active Tasks
- [ ] **Re-run KCNH2 extraction end-to-end** with patched pipeline; regenerate `comparison_results/`. This is the highest-leverage measurement task.
- [ ] **Obtain Elsevier INSTTOKEN** from Vanderbilt E-resources (eresources@vanderbilt.edu). Unlocks PMID 15840476 (86 gold variants).
- [ ] **Wire codon-math bridge into the matcher loop** — `_cdna_indel_protein_positions` and `_protein_indel_position` are now in `cli/compare_variants.py` but not yet called by the comparison logic. Need to add a pass that maps cDNA indels to protein positions and bridges remaining unmatched pairs.
- [ ] **Manual PDF download for the 5 stub PMIDs** (15840476 Tester, 26496715 Karger, 16922724 Wiley CG, 12402336 Wiley HM, 23631430 Sage). Feed through `harvesting/format_converters.py`.
- [ ] Set weekly recall cadence (Friday re-run + compare) — produces trajectory chart for grant.

## Blocked
- [ ] **Karger access** — Cloudflare interstitial doesn't clear even with 30+ cookies + `--no-headless`. TDM request status unknown.
- [ ] **Wiley TDM full-text API** — `WILEY_API_KEY` in `.env` is revoked. Without it, Human Mutation (10.1002/humu.*) is unreachable.
- [ ] **Sage/Liebert (PMID 23631430)** — Sage's CF instance refuses Playwright fingerprint even with the right cookies.

## Backlog
- [ ] Expand validation to SCN5A, KCNQ1
- [ ] Create golden test sets for non-KCNH2 genes
- [ ] Re-run extraction with regex disabled to measure regex vs LLM contribution
- [ ] Consolidate duplicate documentation (ARCHITECTURE.md vs GVF_architecture.md)
- [ ] Audit remaining unmerged worktrees:
  - [ ] `youthful-lederberg-029a0d` — `harvesting/tier4_html_fallback.py` (513 lines, possibly worth comparing against current Tier 3.5 stack)
  - [ ] `strange-hopper-e72c59` — `harvesting/full_context_rebuilder.py` + `utils/table_variant_extractor.py` (table-targeted variant extraction)
  - [ ] `busy-mahavira-f60eeb` — `harvesting/browser_sso_fetcher.py` (uncommitted; earlier authenticated-pool prototype, may have ideas not in current `authenticated_pool.py`)
- [ ] Expand quality gate test set beyond the 10 audit files

See CLAUDE.md for current status, recovery architecture, the Elsevier INSTTOKEN unblock path, and blocker details.
