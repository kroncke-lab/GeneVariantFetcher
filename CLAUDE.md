# GeneVariantFetcher

## Objective
Extract genetic variants from biomedical literature with **90% recall**.
Part of the Kroncke Lab variant interpretation pipeline.

## Current Status
- **Recall:** 45.2% (on golden test set of 13 papers)
- **Latest run:** 2370 variant mentions from 249 papers
- **Deadline:** June 2026 (R01 grant submission)

## Architecture
```
GeneVariantFetcher/
├── cli/                    # Command-line interface
│   ├── extract.py          # Main extraction pipeline
│   └── compare_variants.py # Validation against curated data
├── harvesting/             # Paper discovery & download
│   ├── orchestrator.py     # Main pipeline coordinator
│   ├── pubmed_client.py    # PubMed search
│   ├── supplement_scraper.py # Publisher handlers
│   └── supplement_reference_parser.py # Reference extraction
├── pipeline/               # LLM extraction pipeline
│   ├── extraction.py       # Variant extraction from text
│   └── steps.py            # Pipeline stages
├── config/
│   └── gene_config.py      # Gene-specific validation (protein/transcript lengths)
└── comparison_results/     # Baseline Excel files for validation
```

## Key Commands
```bash
# Run extraction for a gene
python -m cli.extract --gene KCNH2 --output /mnt/temp2/kronckbm/gvf_output/

# Compare against baseline
python -m cli.compare_variants --gene KCNH2 --baseline comparison_results/KCNH2*.xls

# Full pipeline with synonyms
python -m harvesting.orchestrator --gene KCNH2 --use-synonyms
```

## Publisher Coverage
| Publisher | Handler | API | Status |
|-----------|---------|-----|--------|
| Nature | ✅ | - | Working |
| Elsevier | ✅ | ✅ | Working |
| Springer/BMC | ✅ | ❌ | Handler working, API needs registration |
| Oxford Academic | ✅ | - | Working |
| Wiley | ✅ | ✅ | Working |
| Karger | ❌ | - | Blocked by Cloudflare |

## Current Blockers
1. **Springer API** - Brett needs to register with Vanderbilt credentials
2. **Karger access** - Cloudflare blocking; TDM request drafted

## Related Repos
- **VariantFeatures** - Aggregates features for extracted variants
- **BayesianPenetranceEstimator** - Uses extracted phenotype data
- **Variant_Browser** - Displays results clinically

## File Locations
- **Output:** `/mnt/temp2/kronckbm/gvf_output/`
- **Baseline:** `comparison_results/KCNH2*.xls`
- **API keys:** `.env` (Elsevier, Wiley configured)

## Next Steps
See TASKS.md for detailed task list. Priority:
1. Springer API registration
2. Re-run extraction with regex disabled
3. Browser fetch paywalled papers
4. Expand to SCN5A, KCNQ1
