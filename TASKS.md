# GVF Tasks

## Current Focus
- Improve KCNH2 recall from the 2026-05-15 measured baseline to 90% by June 2026 grant deadline.
  - Current KCNH2 unique-variant recall: 323/530 (60.9%).
  - Current KCNH2 variant-row recall: 542/991 (54.7%).
- **Measurement loop is multi-gene now** — score KCNH2, KCNQ1, and SCN5A from `gene_variant_fetcher_gold_standard/normalized/*_recall_input.csv`; RYR2 still needs a per-PMID clinical source before recall can be scored.

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

## Active Tasks
- [ ] **Close KCNH2 source/extraction gap to >90%**. Current post-patch gap is 350 variant rows to 90% against `KCNH2_v12_manual_recovery_20260515.db`; top remaining blockers are PMID 15840476 (86 rows, Elsevier INSTTOKEN), 14661677 (24), 29650123 (21), 24667783 (20), 16922724 (20), 23098067 (16), and 23631430 (12).
- [ ] **Investigate `count_mismatches=117`** from the latest recall summary. Many top mismatches look like SQLite over-counting patient/carrier totals from cohort tables, which is separate from missing-row recall.
- [ ] **Re-run extraction end-to-end for KCNQ1 and SCN5A** with patched table/full-context pipeline; score with `scripts/run_recall_suite.py`.
- [ ] **Create or import RYR2 per-PMID gold input**. Variant_Browser currently exposes only aggregate RYR2 variant counts, so RYR2 is documented but not recall-scorable.
- [ ] **Obtain Elsevier INSTTOKEN** from Vanderbilt E-resources (eresources@vanderbilt.edu). Unlocks PMID 15840476 (86 gold variants). Local OneDrive files named `15840476.pdf` and `26496715.pdf` were checked on 2026-05-14 and are zero-filled placeholders, not usable PDFs.
- [ ] **Run fresh KCNQ1 and SCN5A extraction DBs.** The 2026-05-14 recall runner can score these genes, but current local `results/` has no KCNQ1/SCN5A SQLite DBs yet.
- [ ] **Integrate valid paywall recovery artifacts into canonical runs** before re-extraction; known high-value recovered KCNH2 artifacts are PMIDs 10973849, 11854117, 19038855, 23098067, and 19996378.
- [ ] **Manual PDF download for hard-blocked/stub PMIDs** (15840476 Tester, 26496715 Karger, 16922724 Wiley CG, 12402336 Wiley HM, 23631430 Sage). Feed real downloaded PDFs through `harvesting/format_converters.py`; avoid OneDrive on-demand placeholders.
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
