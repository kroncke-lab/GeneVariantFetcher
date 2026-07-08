# GeneVariantFetcher — Quick Start Guide

Get GVF running in 5 minutes and extract genetic variants from the literature.
This is the canonical local setup and first-run guide; README and agent
handoff files should link here instead of duplicating install or `.env` blocks.

Scope: use this for local installation and a first single-gene run. Use
`NEW_GENE_RUNBOOK.md` for no-gold production runs, `RECALL_REFRESH_RUNBOOK.md`
for existing scored runs, and `END_TO_END_RECALL_RUN.md` for portable recall
runs on another machine.

## Prerequisites

- **Python 3.11+** (required)
- **pip** (Python package manager)
- **LLM provider API key** (Anthropic, OpenAI, or Azure AI; required for extraction)
- **Email address** (required for NCBI API compliance)

Optional but recommended:
- Publisher credentials for expanded paper access, especially
  `ELSEVIER_API_KEY` plus `ELSEVIER_INSTTOKEN`

## Installation

```bash
# Clone the repository
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher

# Create and activate a virtual environment
python3.11 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install GVF in development mode
pip install -e ".[browser,dev]"
python -m playwright install chromium
```

## Configuration

Create a `.env` file in the repository root:

```bash
# Required
NCBI_EMAIL=you@example.com

# Required: at least one LLM provider key
ANTHROPIC_API_KEY=your-anthropic-key
# OPENAI_API_KEY=sk-your-openai-key-here
# AZURE_AI_API_KEY=your-azure-ai-key

# Optional: Publisher APIs for better paper coverage
ELSEVIER_API_KEY=your-elsevier-key
ELSEVIER_INSTTOKEN=your-elsevier-insttoken
SPRINGER_API_KEY=your-springer-key
WILEY_API_KEY=your-wiley-key

# Optional: Higher rate limits with NCBI API key
NCBI_API_KEY=your-ncbi-key
```

See [API_KEYS.md](API_KEYS.md) for instructions on obtaining each key.

## First Run

`gvf gvf-run` is the one command you need. Point it at **your own gene** — any
HGNC symbol works, and no gold standard is required:

```bash
gvf gvf-run <YOUR_GENE> --email you@lab.edu --output ./results [--disease "<phenotype>"]
```

`gvf-run` runs the regular end-to-end path (discovery → triage → full-text and
supplement acquisition → extraction → SQLite → recovery layers), with source
recovery on by default. The `--email` flag is used for NCBI compliance; it is
applied automatically, so you do not need an `NCBI_EMAIL` line in `.env` just to
run (you do still need one LLM provider key).

To reproduce a known result first, KCNH2 (a well-studied cardiac gene) is a good
validation target:

```bash
gvf gvf-run KCNH2 --email you@lab.edu --output ./results --disease "Long QT Syndrome"
```

The lower-level `gvf extract` command exists for debugging individual stages;
prefer `gvf gvf-run` for normal work.

### What This Does

1. **Discovers papers** — Searches PubMind and PubMed for KCNH2-related literature
2. **Downloads full-text** — Retrieves papers from PMC and publisher APIs
3. **Extracts variants** — Uses the configured LLM provider to identify genetic variants and patient data
4. **Aggregates results** — Combines findings across papers
5. **Creates database** — Migrates everything to a queryable SQLite database

### Command Options

```bash
# Limit papers for faster testing
gvf extract KCNH2 --email you@example.com --max-pmids 20 --max-downloads 10

# Verbose output
gvf extract KCNH2 --email you@example.com --verbose

# Custom output directory
gvf extract KCNH2 --email you@example.com --output /path/to/results

# Publish the scored run into the Variant_Browser review DB for collaborator
# adjudication (opt-in, best-effort; needs a sibling Variant_Browser checkout)
gvf gvf-run KCNH2 --email you@example.com --output ./results --publish-review
```

To pull collaborator adjudications back into the gold standard, see
[docs/VARIANT_BROWSER_INTEGRATION.md](VARIANT_BROWSER_INTEGRATION.md)
(`scripts/ingest_review_adjudications.py`).

## What to Expect

### Timing (KCNH2 example)

| Stage | Typical Duration | Notes |
|-------|-----------------|-------|
| Paper discovery | 1-2 min | ~150-300 PMIDs found |
| Full-text download | 5-15 min | ~30-50% success rate (normal) |
| Variant extraction | 10-30 min | Depends on paper count |
| Aggregation + DB | 1-2 min | Fast |
| **Total** | **20-50 min** | Varies by gene |

### Output Files

After completion, you'll find:

```
./output/KCNH2/20260210_143022/
├── KCNH2.db                        # ← SQLite database (query this!)
├── KCNH2_pmids.txt                 # Discovered PMIDs
├── KCNH2_penetrance_summary.json   # Aggregated penetrance data
├── KCNH2_workflow_summary.json     # Run statistics
├── extractions/                    # Per-paper extraction JSONs
│   ├── KCNH2_PMID_12345678.json
│   └── ...
└── pmc_fulltext/                   # Downloaded papers
    ├── 12345678_FULL_CONTEXT.md
    └── ...
```

See [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) for detailed format specifications.

### Sample Output

```
================================================================================
WORKFLOW COMPLETE!
================================================================================
Gene: KCNH2
PMIDs discovered: 245
Papers downloaded: 78
Papers with extractions: 72
Total variants found: 234
Success rate: 32%

SQLite database: ./output/KCNH2/20260210_143022/KCNH2.db
================================================================================
```

## Quick Database Queries

```bash
# Open the database
sqlite3 ./output/KCNH2/20260210_143022/KCNH2.db

# Count variants
sqlite> SELECT COUNT(*) FROM variants;

# List pathogenic variants
sqlite> SELECT protein_notation, clinical_significance
        FROM variants
        WHERE clinical_significance = 'pathogenic';

# Get penetrance data
sqlite> SELECT v.protein_notation, p.total_carriers_observed, p.affected_count
        FROM variants v
        JOIN penetrance_data p ON v.variant_id = p.variant_id
        ORDER BY p.total_carriers_observed DESC
        LIMIT 10;
```

## PMC-Only Mode (No API Keys)

If you don't have publisher API keys, GVF still works using only PubMed Central:

```bash
# Works without publisher keys; skips live paywall recovery
gvf gvf-run KCNH2 --email you@example.com --output ./output --no-source-recovery
```

**Limitations:**
- Only ~30% of papers have PMC full-text available
- Supplemental materials may be limited
- Some high-value papers behind paywalls will be skipped

For comprehensive coverage, obtain at least Elsevier and Springer keys (both free for researchers).

## Troubleshooting

### "No PMIDs found"
- Check your internet connection
- Verify the gene symbol is correct
- Some rare genes have limited literature

### "LLM API error"
- Verify one provider key is present, for example `echo $ANTHROPIC_API_KEY`
- Check the selected provider account has credits or quota
- Rate limits may require waiting

### "Few papers downloaded"
- Normal — only ~30% of PubMed has PMC full-text
- Add publisher API keys for better coverage
- Check `pmc_fulltext/paywalled_missing.csv` for inaccessible papers

### Extraction seems slow
- LLM calls take time; this is normal
- Use `--max-downloads 10` for testing
- Check `--verbose` output for progress

## Next Steps

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Understand the pipeline
- **[API_KEYS.md](API_KEYS.md)** — Get publisher API keys
- **[OUTPUT_FORMAT.md](OUTPUT_FORMAT.md)** — Work with results
- **[RECALL_STATUS.md](RECALL_STATUS.md)** — Current recall metrics and blockers

## Getting Help

1. Check the troubleshooting section above
2. Review error logs in `{gene}_workflow.log`
3. Open an issue on GitHub with:
   - Gene symbol you tried
   - Error message
   - Contents of `{gene}_workflow.log`
