# GeneVariantFetcher

## Objective
Extract genetic variants from biomedical literature with **90% recall**.
Part of the Kroncke Lab variant interpretation pipeline.

## Current Status (v2.1.0 — 2026-02-12)
- **Unique Variant Recall:** 59.1% (289 of 489 variants)
- **Carrier Recall:** 71.3%
- **Affected Recall:** 67.6%
- **Deadline:** June 2026 (R01 grant submission)

### Recent Architecture Changes
1. **Scanner decoupled** — Now runs on FULL_CONTEXT.md, not condensed DATA_ZONES
2. **Concatenated gene+variant regex** — Detects HERGG604S, KCNH2A561V patterns
3. **Unicode arrow normalization** — Handles →, ➔, ⟶ variants
4. **Canonical form comparison** — Matching uses normalized keys for higher accuracy
5. **Tier 3 model** — gpt-4o (upgraded from gemini-2.0-flash)
6. **Unified supplement fetcher** — `gene_literature/supplements/` wired into orchestrator; tiered PMC + Elsevier API fetch with dedup, DOI scraping as fallback

### Recall Improvement Changes (2026-05-10)
7. **Download cap removed** — `max_papers_to_download` default changed from 50 to None (uncapped). The 50-paper cap was silently discarding 350+ filtered PMIDs.
8. **Fail-open Tier 1 filter** — `TIER1_MIN_KEYWORDS` lowered from 2 to 1; papers with no abstract now PASS through to Tier 2 instead of being dropped.
9. **Source-completeness report** — New Step 3.5 emits `source_completeness.json` after extraction: abstract-only count, stub files, supplement coverage, zero-variant papers, single-carrier-only papers.
10. **Zero-variant QA re-extraction** — Papers that pass Tier 2 but yield 0 variants are flagged. Set `GVF_QA_MODEL=anthropic/claude-opus-4-7` to re-extract them with a stronger model.
11. **Cohort carrier count prompt** — Extraction prompts now explicitly handle cohort/screening tables where the same variant appears across multiple tables/cohorts, with instructions to SUM carrier counts.

## Health Assessment (2026-05-02)

**Status: behind schedule, but architecture is sound.**

### Recall gap
- 30.9 percentage points to close (59.1% → 90%) with ~8 weeks until June 2026 grant deadline.
- 918 missing variant rows in `comparison_results/` cover 505 unique variants across 236 PMIDs.
- Top-10 PMIDs account for **~389 missing entries (42% of the gap)** — large cohort/screening papers. PMID 15840476 alone accounts for 86 missing variants.

### Variant-type breakdown of the gap
- ~708 standard missense (`A123B`) — should be in regex+LLM range
- ~156 frameshift/splice (`fsX`, `fs*`) — extractable with current pipeline
- ~32 structural/CNV (`7q36.1q36.2Del`) — likely out of scope
- The missense bulk is the red flag: if straightforward variants from known PMIDs aren't reaching the DB, the bottleneck is **upstream of extraction** (harvest, supplements, paywall) rather than LLM accuracy.

### Velocity warning
- Feb 12–16: 25+ commits, major architecture overhaul (preprocessor, concurrency, tier upgrade, supplement fetcher)
- Mar–Apr: 5 commits total
- May: 0 commits as of 2026-05-02
- The comparison baseline (`comparison_results/summary.json`) references a SQLite DB from 2025-12-11. **Every Feb–Apr fix is currently unmeasured.**

### Critical risks for June
1. **No re-run since architecture changes.** The Feb work (scanner-on-FULL_CONTEXT, concatenated regex, gpt-4o tier 3, unified supplements) very likely improved recall — but with no re-extraction the trajectory is invisible.
2. **Stale comparison.** Baseline DB is 5 months old. Recall could already be 70%+ and we wouldn't know.
3. **pytest broken locally.** `python -m pytest` fails (no module installed in active env). CI gating exists but local validation loop is broken.
4. **Top-10 PMIDs untouched.** April's "browser fetch from top-10 missing PMIDs" commit landed but the 918-row gap predates it; need to verify those papers are actually harvested + extracted.

### Recommended priority order (next 2 weeks)
1. **Re-run KCNH2 extraction end-to-end** with current main; regenerate `comparison_results/`. This closes the measurement loop and probably reveals 10–15 points of already-earned recall.
2. **Triage top-10 missing PMIDs**: for each, confirm (a) harvested, (b) supplements fetched, (c) Data Scout reached the variant tables. A handful of large cohort papers can close most of the gap.
3. **Fix local pytest env** so changes can be validated without a full extraction.
4. **Set a weekly recall cadence** (Friday re-run + compare) — produces a trajectory chart for the grant.

### What's working
- Architecture is right: tiered extraction, modular supplements, canonical normalization
- Test coverage exists for key modules (just needs an env to run in)
- Tooling is mature: paywall audit, browser fetch workflow, comparison harness
- The "is the bug technical or operational" answer is **operational** — the pipeline likely works better than the headline 59.1% reflects.

## Architecture
```
GeneVariantFetcher/
├── cli/                        # Typer CLI (entry point: `gvf`)
│   ├── __init__.py             # App definition: extract, scout, audit-paywalls
│   ├── automated_workflow.py   # Main extraction orchestration
│   ├── compare_variants.py     # Validation against curated data
│   ├── scout.py                # Standalone Data Scout command
│   ├── audit_paywalls.py       # Paywall/captcha/OA classification
│   ├── browser_fetch.py        # Interactive browser for paywalled papers
│   └── fetch_manager.py        # Manual fetch workflow
├── gene_literature/            # Paper discovery
│   ├── pubmed_client.py        # PubMed search
│   ├── pubmind_fetcher.py      # PubMind variant search
│   ├── synonym_finder.py       # Gene alias lookup (NCBI Gene)
│   ├── discovery.py            # Multi-source PMID discovery
│   ├── collector.py            # Coordinates discovery + relevance
│   └── supplements/            # Tiered supplement fetcher
│       ├── base.py             # SupplementFile dataclass + abstract base
│       ├── pmc_fetcher.py      # Tier 1: Europe PMC + NCBI OA
│       ├── elsevier_fetcher.py # Tier 2: Elsevier API
│       └── unified.py          # UnifiedSupplementFetcher (orchestrates tiers)
├── harvesting/                 # Paper download & conversion
│   ├── orchestrator.py         # Main pipeline coordinator (uses UnifiedSupplementFetcher)
│   ├── supplement_scraper.py   # Publisher-specific web scrapers (Tier 3 fallback)
│   ├── format_converters.py    # XML/PDF/Excel/Word → Markdown
│   ├── elsevier_api.py         # Elsevier full-text API
│   ├── springer_api.py         # Springer full-text API
│   ├── wiley_api.py            # Wiley TDM API
│   ├── pmc_api.py              # PMC BioC XML API
│   ├── unpaywall_api.py        # OA link resolution
│   ├── doi_resolver.py         # DOI → publisher URL routing
│   └── migrate_to_sqlite.py    # JSON → SQLite migration
├── pipeline/                   # LLM extraction pipeline
│   ├── extraction.py           # Variant extraction from text
│   ├── steps.py                # Pipeline stages (Tier 1/2/3)
│   ├── filters.py              # Keyword + LLM relevance filters
│   ├── data_scout.py           # Data zone identification
│   ├── prompts.py              # LLM prompt templates
│   ├── aggregation.py          # Cross-paper variant aggregation
│   └── pedigree_extractor.py   # Pedigree image analysis
├── utils/                      # Shared utilities
│   ├── variant_normalizer.py   # HGVS normalization + PROTEIN_LENGTHS for 8 cardiac genes
│   ├── variant_scanner.py      # Regex pre-scanner for variants in text
│   ├── pubmed_utils.py         # PubMed/Entrez helpers
│   ├── retry_utils.py          # Tenacity retry decorators
│   ├── resilience.py           # Circuit breaker for APIs
│   ├── llm_utils.py            # LiteLLM wrapper
│   └── manifest.py             # Run manifest tracking
├── config/
│   ├── settings.py             # Pydantic settings (all env vars + defaults)
│   └── constants.py            # Shared constants
├── gui/                        # Web GUI (Gradio)
├── comparison_results/         # Baseline Excel files for validation
├── golden_test_set/            # Golden test specifications
└── tests/                      # pytest suite (unit + integration + recall)
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

# Launch GUI
python main.py
```

## Publisher Coverage
| Publisher | Handler | API | Status |
|-----------|---------|-----|--------|
| Nature | scraper | - | Working |
| Elsevier | scraper | full-text + supplements | Working |
| Springer/BMC | scraper | full-text | Working (API active) |
| Oxford Academic | scraper | - | Working |
| Wiley | scraper | TDM full-text | Working |
| PMC | scraper | BioC XML + Europe PMC supplements | Working (free) |
| Karger | - | - | Blocked by Cloudflare |

## Current Blockers
1. **Karger access** — Cloudflare blocking; TDM request drafted

## Related Repos
- **VariantFeatures** — Aggregates features for extracted variants
- **BayesianPenetranceEstimator** — Uses extracted phenotype data
- **Variant_Browser** — Displays results clinically

## File Locations
- **Output:** `/mnt/temp2/kronckbm/gvf_output/`
- **Baseline:** `comparison_results/KCNH2*.xls`
- **API keys:** `.env` (Elsevier, Wiley, Springer configured)
- **Gene validation:** `utils/variant_normalizer.py` (PROTEIN_LENGTHS dict)

## Next Steps
See TASKS.md for detailed task list. Priority:
1. Browser fetch paywalled papers from top-10 missing PMIDs
2. Add Wiley supplement handler
3. Expand to SCN5A, KCNQ1
4. Create golden test sets for non-KCNH2 genes
