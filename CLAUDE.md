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
