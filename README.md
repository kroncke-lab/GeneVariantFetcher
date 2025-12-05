# GeneVariantFetcher

Extract human genetic variant carriers from biomedical literature into a SQLite database.

## What It Does

Given a gene name, this tool:
1. Discovers variant-related papers via PubMind and PubMed keyword searches, then merges and de-duplicates PMIDs
2. Downloads full-text articles and supplemental materials (Excel, Word)
3. Extracts patient-level variant data using LLM cascade
4. Aggregates penetrance statistics across papers
5. Writes normalized data to SQLite

**Key Insight:** 70-80% of variant data is in supplemental files. This tool extracts it all.

**Why combine PubMind + PubMed?** PubMind's LLM-driven filtering is highly precise but can miss edge cases; PubMed's recall-oriented keyword search helps recover those. Using both by default balances precision and recall while keeping the merged PMID list clean.

## Quick Start

```bash
# Install
pip install -e .

# Set environment variables
export AI_INTEGRATIONS_OPENAI_API_KEY="your-key"
export NCBI_EMAIL="your@email.com"

# Run
python automated_workflow.py BRCA1 --email your@email.com --output ./results

# Control sourcing (optional)
# Disable PubMed queries
export PUBMIND_ONLY=true           # or: export USE_PUBMED=false

# Force PubMed alongside PubMind (recall-heavy)
export USE_PUBMED=true             # keeps default combined sourcing enabled
```

## CLI Options

```bash
python automated_workflow.py GENE --email EMAIL --output DIR [OPTIONS]

Options:
  --max-pmids N       Limit PMIDs to discover (default: 100)
  --max-downloads N   Limit papers to download (default: 50)
  --tier-threshold N  Model cascade threshold (default: 1)
  --verbose           Enable verbose logging

PMID sources:
  • PubMind-only: set PUBMIND_ONLY=true (skips PubMed keyword search)
  • PubMind + PubMed (default): leave both unset, or set USE_PUBMED=true
```

## Output

```
{OUTPUT_DIR}/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids_pubmind.txt        # PMIDs from PubMind
├── {GENE}_pmids_pubmed.txt         # PMIDs from PubMed keyword search
├── {GENE}_pmids.txt                # Merged + de-duplicated PMIDs used by the workflow
├── pmc_fulltext/                   # Full-text + supplements
├── extractions/                    # Per-paper JSON
├── {GENE}_penetrance_summary.json
└── {GENE}.db                       # SQLite database
```

## Database Schema

- **papers** - Paper metadata (pmid, title, doi)
- **variants** - Variant info (gene, cDNA, protein notation)
- **individual_records** - Patient data (age, sex, affected_status, phenotype)
- **penetrance_data** - Cohort statistics

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `AI_INTEGRATIONS_OPENAI_API_KEY` | Yes | OpenAI API key (or `OPENAI_API_KEY`) |
| `NCBI_EMAIL` | Yes | Email for NCBI E-utilities |
| `NCBI_API_KEY` | No | Increases NCBI rate limits |
| `USE_PUBMED` | No | `true`/`false` toggle to include PubMed keyword search (defaults to `true`) |
| `PUBMIND_ONLY` | No | If `true`, skip PubMed/EuropePMC and rely solely on PubMind |

## Development

See [CLAUDE.md](CLAUDE.md) for architecture details.

```bash
# Run tests
pytest tests/
```
