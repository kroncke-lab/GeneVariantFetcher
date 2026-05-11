# GeneVariantFetcher

## Objective
Extract genetic variants from biomedical literature with **90% recall**.
Part of the Kroncke Lab variant interpretation pipeline.

## Current Status (v2.1.0 — 2026-05-11)

- **Unique Variant Recall (last measured baseline):** 59.1% (289 of 489 variants)
- **Deadline:** June 2026 (R01 grant submission)
- **Critical caveat:** the 59.1% baseline is from a 2025-12-11 SQLite DB. Every change since February 2026 is **unmeasured.** Re-running KCNH2 end-to-end and regenerating `comparison_results/` is the single highest-value next step.

### 2026-05-11 Recovery Pipeline Build

Built and validated a paywall-recovery pipeline against the 10 highest-value paywalled PMIDs in the gap (~270 missing variants between them):

- **Tier 3.5 strategies modernized** — AHA primary selector switched to `#bodymatter` (legacy `.hlFld-Fulltext` is the 8-char eLetters wrapper post-2024-redesign); Elsevier imprint sites (Heart Rhythm, JACC, Cell, Lancet) auto-rewrite `/abstract` → `/fulltext`; all strategies emit `"cloudflare"` in `result.notes` on interstitial so the retry sweep fires; DOM-walker fallback (`harvesting/browser_html/dom_extract.py`) replaces brittle class-name extraction.
- **Authenticated Playwright** — `harvesting/browser_html/cookie_loader.py` reads cookies from local Chrome via `browser_cookie3` (146 cookies across AHA / ScienceDirect / Wiley / Karger / Oxford / Sage / VUMC SSO / Microsoft SSO). `harvesting/browser_html/authenticated_pool.py` injects them into Playwright with `channel="chrome"` for fingerprint match. Does **not** conflict with Chrome being open.
- **Quality gate hardened** — `harvesting/browser_html/content_quality.py` rejects: (a) explicit paywall sentences ("Log in, subscribe or purchase for full access", "You do not currently have access to this content", "Get full access to this article", "Available to Purchase"); (b) bodies with <3 KB of substantive paragraph + table content; (c) author-affiliation padding (same line repeated 5+ times) when body content is also marginal. Verified against 10 audit files: 10/10 correctly classified.
- **PMC fallback** — `scripts/fetch_paywalled.py` looks up the PMCID via Europe PMC after a Tier 3.5 stub/empty, fetches the PMC HTML, and overwrites `FULL_CONTEXT.md` only when the new extraction passes the gate. Recovered PMID 19038855 via PMC2677528 (33.7 KB body).
- **Variant matcher refinements** — `cli/compare_variants.py` gained `_positions_compatible(a, b)` with a subset rule (so `p.Pro926AlafsTer14` (926, 14) matches `P926fsX` (926,)), `_cdna_indel_protein_positions(cdna)` codon math (c.842dupG → aa 281), `_protein_indel_position(variant)` extraction, and a `consumed` set parameter on `find_best_match` for 1-to-1 enforcement at match time.

**Honest recovery scorecard for the top-10 paywalled PMIDs:**

| PMID | Source | Result | Path |
|---|---|---|---|
| 10973849 | Splawski 2000 Circulation (59 gold) | ✅ real body | AHA strategy |
| 11854117 | Moss 2002 Circulation pore (44 gold) | ✅ real body | AHA strategy |
| 19038855 | Johnson 2009 Neurology (28 gold) | ✅ real body | **PMC fallback** (PMC2677528) |
| 19996378 | Horigome 2010 Circ EP | ✅ real body | AHA strategy (--no-headless required for CF) |
| 23098067 | Stattin 2012 BMC OA | ✅ real body | OA (no auth needed) |
| **15840476** | **Tester 2005 Heart Rhythm (86 gold — THE biggest prize)** | ❌ stub | **needs Elsevier INSTTOKEN** |
| 26496715 | Lu 2015 Karger Cardiology (53 gold) | ❌ stub | no OA path; Karger CF + paywall |
| 16922724 | Millat 2006 Wiley Clinical Genetics | ❌ stub | no OA path; Wiley TDM key revoked |
| 12402336 | Jongbloed 2002 Wiley Human Mutation | ❌ stub | no OA path; Wiley CF blocks |
| 23631430 | Hofman 2013 Sage/Liebert GTMB | ❌ stub | Sage CF doesn't clear even non-headless |

**Net of recovery work:** 5/10 of the highest-value paywalled PMIDs are now in `results/KCNH2/20260506_102238/paywall_recovery_v*/` as real full-text. The remaining 5 fail for one of three reasons: no PMC deposit, no Wiley TDM access, or Cloudflare-protected paywall stubs we can't bypass without an institutional IP/token.

## The Elsevier INSTTOKEN situation (highest-ROI unblock)

PMID 15840476 alone is 86 missing variants ≈ 18% of the 489-variant baseline. The Elsevier full-text retrieval API at `https://api.elsevier.com/content/article/doi/{doi}?view=FULL` returns the complete article XML (body + tables + supplements) when called with **two** auth headers:

```
X-ELS-APIKey:    <ELSEVIER_API_KEY>      # already in .env
X-ELS-Insttoken: <ELSEVIER_INSTTOKEN>    # MISSING — request from Vanderbilt IT
```

The API key alone returns metadata + abstract (≈6 KB; `originalText` is length 1). The insttoken is what unlocks subscription content. It's a per-institution token issued by Elsevier to libraries that subscribe to Heart Rhythm and other titles. `.env` already has a commented placeholder (`# ELSEVIER_INSTTOKEN=`) — getting one means contacting Vanderbilt's E-resources librarian (eresources@vanderbilt.edu) and asking for the Elsevier API institutional token.

Once set, `harvesting/elsevier_api.py` already supports it (`X-ELS-Insttoken` header is sent automatically when present). One config line unlocks programmatic access to 15840476 and every other Elsevier subscription paper.

**Workaround until INSTTOKEN is obtained:** Brett can open the article in his logged-in Chrome (Vanderbilt SSO unlocks Heart Rhythm), save the PDF, and feed it through the existing PDF→markdown pipeline (`harvesting/format_converters.py` handles this).

## Health Assessment (2026-05-11 update)

### Recall gap
- 30.9 percentage points to close (59.1% → 90%) with ~6 weeks until June 2026 grant deadline.
- 918 missing variant rows; top-10 PMIDs account for ~389 (42% of the gap). Of those 10, 5 are now recovered as real full text (≈131 variants worth of body content) and 5 remain stubs (≈168 variants behind paywalls we can't programmatically pierce).
- The Tester paper (PMID 15840476) at 86 missing variants is by far the largest single bottleneck.

### What's working
- Authenticated cookie-based Playwright fetching for Elsevier imprint sites (Heart Rhythm), AHA, Wiley OA titles, Karger abstracts, Springer/BMC.
- Quality gate that correctly rejects paywall stubs.
- PMC fallback for NIH-funded subscription papers.
- DOM-walker fallback when legacy publisher selectors rot.
- `cli/compare_variants.py` matcher with positional-digit guard + cDNA↔protein bridge + 1-to-1 assignment.

### What's not working
- Karger and Sage Cloudflare instances refuse our Playwright fingerprint even with cookies + non-headless mode.
- Wiley Human Mutation (`10.1002/humu.*`) hits CF after one back-to-back request.
- Elsevier API metadata-only without INSTTOKEN.
- Comparison baseline is 5 months stale; recall trajectory is invisible.

## Architecture
```
GeneVariantFetcher/
├── cli/                            # Typer CLI (entry point: `gvf`)
│   ├── __init__.py                 # App definition: extract, scout, audit-paywalls
│   ├── automated_workflow.py       # Main extraction orchestration (Step 3.5 source-completeness)
│   ├── compare_variants.py         # Validation against curated data
│   │                               # incl. _positions_compatible, _cdna_indel_protein_positions,
│   │                               # _protein_indel_position, _find_cdna_protein_bridge
│   ├── scout.py                    # Standalone Data Scout command
│   ├── audit_paywalls.py           # Paywall/captcha/OA classification
│   ├── browser_fetch.py            # Interactive browser for paywalled papers
│   ├── reharvest.py                # Reharvest stale FULL_CONTEXT.md from a run dir
│   └── fetch_manager.py            # Manual fetch workflow
├── gene_literature/                # Paper discovery
│   ├── pubmed_client.py            # PubMed search
│   ├── pubmind_fetcher.py          # PubMind variant search
│   ├── synonym_finder.py           # Gene alias lookup (NCBI Gene)
│   ├── discovery.py                # Multi-source PMID discovery
│   ├── collector.py                # Coordinates discovery + relevance
│   └── supplements/                # Tiered supplement fetcher
│       ├── base.py                 # SupplementFile dataclass + abstract base
│       ├── pmc_fetcher.py          # Tier 1: Europe PMC + NCBI OA
│       ├── elsevier_fetcher.py     # Tier 2: Elsevier API (uses INSTTOKEN when present)
│       └── unified.py              # UnifiedSupplementFetcher (orchestrates tiers)
├── harvesting/                     # Paper download & conversion
│   ├── orchestrator.py             # Main pipeline coordinator
│   ├── supplement_scraper.py       # Publisher-specific web scrapers (Tier 3 fallback)
│   ├── format_converters.py        # XML/PDF/Excel/Word → Markdown
│   ├── elsevier_api.py             # Elsevier full-text API (needs INSTTOKEN for full body)
│   ├── springer_api.py             # Springer full-text API
│   ├── wiley_api.py                # Wiley TDM API (key revoked as of 2026-05)
│   ├── pmc_api.py                  # PMC BioC XML API
│   ├── unpaywall_api.py            # OA link resolution
│   ├── doi_resolver.py             # DOI → publisher URL routing
│   ├── migrate_to_sqlite.py        # JSON → SQLite migration
│   ├── oa_recovery.py              # OARecoveryClient: Europe PMC / NCBI PMC OA / Unpaywall
│   ├── html_body_fetcher.py        # Focused body-text fetcher for stale FULL_CONTEXT.md
│   └── browser_html/               # Tier 3.5 — Authenticated Playwright pipeline
│       ├── __init__.py
│       ├── base.py                 # PublisherStrategy ABC + FetchResult + FetchContext
│       ├── fetcher.py              # BrowserHTMLFetcher; takes pool=, bypass_embargo=,
│       │                           # enabled_override= for authenticated callers
│       ├── browser_pool.py         # Default headless Chromium pool
│       ├── authenticated_pool.py   # Cookie-inheriting Chrome-channel pool
│       ├── cookie_loader.py        # browser_cookie3 → Playwright cookie shape
│       ├── content_quality.py      # validate_article_content: hard rejects on
│       │                           # paywall phrases, body floor, padding detection
│       ├── dom_extract.py          # extract_body_markdown, pick_better_markdown,
│       │                           # looks_like_cloudflare_challenge, looks_like_paywall_stub
│       ├── embargo.py              # EmbargoChecker + get_pub_date_from_pmid
│       └── strategies/             # Per-publisher (auto-registered)
│           ├── aha.py              # #bodymatter selector + multi-cycle CF wait
│           ├── elsevier_open.py    # /abstract → /fulltext rewrite for imprints
│           ├── generic.py          # DOM walker + CF wait (Karger, Sage, BMJ)
│           ├── oxford.py
│           └── wiley.py            # CF wait + DOM walker fallback
├── pipeline/                       # LLM extraction pipeline
│   ├── extraction.py               # Variant extraction (rejects DATA_ZONES garbage)
│   ├── steps.py                    # Pipeline stages (Tier 1/2/3)
│   ├── filters.py                  # Keyword + LLM relevance (Tier 1 fail-open)
│   ├── data_scout.py               # Data zone identification
│   ├── prompts.py                  # LLM prompt templates (cohort carrier-sum)
│   ├── aggregation.py              # Cross-paper variant aggregation
│   └── pedigree_extractor.py       # Pedigree image analysis
├── utils/                          # Shared utilities
│   ├── variant_normalizer.py       # HGVS normalization + PROTEIN_LENGTHS for 8 cardiac genes
│   ├── variant_scanner.py          # Regex pre-scanner; incl. parenthesized HGVS p.(Arg176Trp)
│   ├── pubmed_utils.py             # PubMed/Entrez helpers
│   ├── retry_utils.py              # Tenacity retry decorators
│   ├── resilience.py               # Circuit breaker for APIs
│   ├── llm_utils.py                # LiteLLM wrapper (provider-aware rate limit)
│   └── manifest.py                 # Run manifest tracking
├── config/
│   ├── settings.py                 # Pydantic settings; incl. BROWSER_HTML_USE_PROFILE
│   └── constants.py                # Shared constants
├── scripts/                        # Standalone CLI scripts
│   ├── fetch_paywalled.py          # ⭐ Authenticated Playwright recovery + PMC fallback
│   ├── recover_via_pmc_routes.py   # Europe PMC / NCBI PMC OA / Unpaywall recovery
│   ├── recover_abstract_only.py    # Resolve abstract-only PMIDs via OA paths
│   ├── harvest_missing.py          # Harvest-only audit (no LLM spend)
│   └── ... (legacy validation/recall scripts)
├── gui/                            # Web GUI (Gradio)
├── comparison_results/             # Baseline Excel files for validation
├── golden_test_set/                # Golden test specifications
└── tests/                          # pytest suite (unit + integration + recall)
```

## Key Commands
```bash
# Run extraction for a gene
gvf extract KCNH2 --email you@email.com --output /mnt/temp2/kronckbm/gvf_output/

# With synonym discovery and Data Scout
gvf extract KCNH2 --email you@email.com --output ./results --auto-synonyms --scout-first

# Run Data Scout standalone
gvf scout ./results/KCNH2/*/pmc_fulltext --gene KCNH2

# Audit paywall status
gvf audit-paywalls --harvest-dir ./results/KCNH2/*/pmc_fulltext --out-dir ./audit

# Compare against baseline (standalone script)
python -m cli.compare_variants --gene KCNH2 --baseline comparison_results/KCNH2*.xls

# Authenticated paywall recovery — uses local Chrome cookies + PMC fallback
.venv/bin/python scripts/fetch_paywalled.py \
    --pmid 15840476 --pmid 10973849 ... \
    -o results/KCNH2/<run>/paywall_recovery/

# OA-only recovery for abstract-only PMIDs (Unpaywall + Europe PMC + bioRxiv)
.venv/bin/python scripts/recover_abstract_only.py \
    --harvest-dir results/KCNH2/<run>/pmc_fulltext \
    --email brett.kroncke@gmail.com

# PMC-routes recovery
.venv/bin/python scripts/recover_via_pmc_routes.py ...

# Reharvest a stale run
gvf reharvest --run-dir results/KCNH2/<run>/

# Launch GUI
python main.py
```

## Publisher Coverage (2026-05-11)
| Publisher | Status | Notes |
|-----------|--------|-------|
| PMC / Europe PMC | ✅ Working | Free; PMC fallback in `fetch_paywalled.py` auto-rescues NIH-funded subscription papers |
| Springer / BMC | ✅ Working | OA + API |
| Nature | ✅ Working | scraper |
| AHA Journals (Circulation, JAHA, Stroke, Circ EP) | ✅ Working with cookies | `#bodymatter` selector; CF clears after multi-cycle wait; sequential same-host requests escalate CF — use 15-30s pacing |
| Elsevier — Heart Rhythm, JACC, Cell, Lancet (imprint sites) | 🟡 Partial | `/abstract` → `/fulltext` rewrite works only when user has institutional access. Without Elsevier INSTTOKEN, the API returns metadata only |
| Wiley OA titles (Clinical Genetics) | ✅ Working for OA papers | Subscription Wiley + Human Mutation: ❌ blocked by CF |
| Oxford Academic | ✅ Working | Embargo 12 months |
| Karger | ❌ Blocked | Cloudflare interstitial doesn't clear even with 30+ Karger cookies + non-headless mode |
| Sage / Mary Ann Liebert | ❌ Blocked | Same as Karger; Sage CF refuses Playwright fingerprint |

## Current Blockers (in priority order)
1. **Elsevier INSTTOKEN for Vanderbilt** — single config line unlocks PMID 15840476 (86 missing variants ≈ 18% of recall gap). Request from `eresources@vanderbilt.edu`.
2. **No re-run since 2025-12-11.** Comparison baseline DB is 5 months old. Every change since February is unmeasured. **Re-run KCNH2 end-to-end and regenerate `comparison_results/`** — this is the highest-leverage measurement task.
3. **Karger/Sage Cloudflare** — no current workaround; would need either VPN/EZproxy, IP allowlisting, or a different fetcher tier (e.g. Karger TDM API agreement).
4. **Wiley TDM key revoked** — `WILEY_API_KEY` in `.env` no longer authenticates. Without it, Human Mutation papers behind CF can't be unlocked.

## Related Repos
- **VariantFeatures** — Aggregates features for extracted variants
- **BayesianPenetranceEstimator** — Uses extracted phenotype data
- **Variant_Browser** — Displays results clinically

## File Locations
- **Output:** `/mnt/temp2/kronckbm/gvf_output/` (external mount on Linux box)
- **Local results (this machine):** `results/KCNH2/<timestamp>/`
- **Paywall recovery artifacts:** `results/KCNH2/20260506_102238/paywall_recovery{,_v2,_v3,_v4}/`
- **Baseline:** `comparison_results/KCNH2*.xls`
- **API keys:** `.env` (Elsevier, Wiley, Springer, NCBI configured; `ELSEVIER_INSTTOKEN` MISSING)
- **Gene validation:** `utils/variant_normalizer.py` (`PROTEIN_LENGTHS` dict, 8 cardiac genes)

## Next Steps (in execution order)

### Tier 0: measurement (do this first)
1. **Re-run KCNH2 extraction end-to-end** with current main on the dev box, regenerate `comparison_results/`. Probably reveals 10-15 points of already-earned recall.
2. **Pull the top-10 missing-PMID list from the fresh comparison.** Verify which are still missing after recovery work landed.

### Tier 1: highest-ROI recall wins
3. **Obtain Elsevier INSTTOKEN** from Vanderbilt — unlocks PMID 15840476 (86 variants) via the API. One config line.
4. **Brett manually downloads PDFs for 3-4 hardest paywalled papers** that can't be auto-recovered (15840476 until INSTTOKEN lands, Karger/Sage/Wiley HM). Feed through `harvesting/format_converters.py`.
5. **Wire codon-math bridge in `cli/compare_variants.py`** — the new helpers `_cdna_indel_protein_positions` and `_protein_indel_position` are present but not yet called in the matcher loop. Wire them so `c.842dupG` matches `R281fsX`.

### Tier 2: pipeline hygiene
6. **Set weekly recall cadence** (Friday re-run + compare → trajectory chart for grant).
7. **Expand quality gate test set** beyond the 10 audit files (regression coverage).
8. **Audit remaining worktrees for valuable code** — `youthful-lederberg-029a0d` has a `harvesting/tier4_html_fallback.py` (513 lines) and a `cdna_protein_bridge_match` function worth comparing. `strange-hopper-e72c59` has `harvesting/full_context_rebuilder.py` and `utils/table_variant_extractor.py`. `busy-mahavira-f60eeb` has a `harvesting/browser_sso_fetcher.py` (uncommitted) — likely an earlier authenticated-pool prototype.

### Tier 3: scope expansion (after KCNH2 ≥ 85%)
9. Expand validation to SCN5A, KCNQ1.
10. Create golden test sets for non-KCNH2 genes.
11. Re-run extraction with regex disabled to measure regex vs LLM contribution.
12. Consolidate duplicate documentation (ARCHITECTURE.md vs GVF_architecture.md).

## Notes for AI Agents Picking This Up

- The repo uses a `.venv` at the project root. `pytest -q` from project root works.
- Pre-commit hooks (`ruff`, `ruff-format`, trailing-whitespace) will modify staged files on commit; **expect to re-stage and recommit** when the hook reformats. Never use `--no-verify`.
- The agent-scratch files (`HEARTBEAT.md`, `IDENTITY.md`, `SOUL.md`, `TOOLS.md`, `USER.md`, `.openclaw/`) are intentionally untracked — **don't commit them**.
- Worktrees under `.claude/worktrees/` are parallel investigation branches. Many predate current main; their compare_variants/strategies are usually OLDER than what's in main. Diff carefully before porting.
- For paywall recovery, the canonical user-facing CLI is `scripts/fetch_paywalled.py`. It auto-resolves DOIs via NCBI esummary, loads cookies from local Chrome, runs Tier 3.5, applies the quality gate, falls back to PMC, and produces `summary.json` per run.
- Brett's machine: Mac Mini with logged-in Chrome (VUMC SSO active via Microsoft). NOT on VPN. NOT on Vanderbilt IP. Cookie-based recovery only works for publishers that accept session cookies as the auth credential — most Elsevier subscription content does NOT (it gates on institutional IP or INSTTOKEN).
