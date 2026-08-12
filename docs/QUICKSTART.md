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
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install GVF in development mode
pip install -e ".[browser,dev]"
python -m playwright install chromium
```

## Configuration

Create a `.env` file in the repository root:

```bash
# Required
NCBI_EMAIL=brett.kroncke@gmail.com

# Required: select a provider and configure its key. Anthropic is the shipped
# default; set MODEL_PROVIDER=azure or openai when using those providers.
MODEL_PROVIDER=anthropic
ANTHROPIC_API_KEY=your-anthropic-key
# MODEL_PROVIDER=azure
# AZURE_AI_API_KEY=your-azure-ai-key
# AZURE_AI_API_BASE=https://your-resource.services.ai.azure.com
# MODEL_PROVIDER=openai
# OPENAI_API_KEY=sk-your-openai-key-here

# Optional: Publisher APIs for better paper coverage
ELSEVIER_API_KEY=your-elsevier-key
ELSEVIER_INSTTOKEN=your-elsevier-insttoken
SPRINGER_API_KEY=your-springer-key
WILEY_API_KEY=your-wiley-key

# Optional: Higher rate limits with NCBI API key
NCBI_API_KEY=your-ncbi-key
```

See [API_KEYS.md](API_KEYS.md) for instructions on obtaining each key.

## Local Corpus Preflight

On Brett's current workstation, the repo's `corpus/` path is an absolute
symlink to
`/Volumes/Ezekers/ResearchData/GeneVariantFetcher/corpus`. Before the first run
or any corpus maintenance job:

```bash
test -L corpus && test -d corpus
test "$(readlink corpus)" = "/Volumes/Ezekers/ResearchData/GeneVariantFetcher/corpus"
```

If a check fails, mount the APFS volume named `Ezekers` at `/Volumes/Ezekers`.
Do not replace the broken symlink or create a local `corpus/`, because that
would send new papers to a second corpus on the internal disk. The workstation
symlink is local-only and untracked; a fresh checkout must recreate it after
verifying that the external target exists.

`gvf-run`'s doctor step enforces this: a **dangling** `corpus/` symlink (the drive
was expected and is not attached) stops the run at Step 1 rather than silently
re-fetching the whole corpus over the network. It names the volume from the link
itself, so the message is correct on any machine. `--skip doctor` overrides.

### On a collaborator's machine

You do **not** need Brett's drive. `corpus/` is only a cache of already-fetched
source, so an **absent** `corpus/` never blocks a run — it just starts cold. Pick
one of:

```bash
mkdir corpus                      # plain local cache (needs the flag below)
ln -s /path/to/your/storage corpus   # or your own external storage
export GVF_CORPUS_DIR=/path/to/corpus   # or point at it without a symlink
```

A plain local `corpus/` on the internal disk is a deliberate choice, so confirm it
with `GVF_ALLOW_LOCAL_CORPUS=1` (otherwise corpus-writing jobs refuse, to protect
against an accidentally-detached drive). What you actually need for a first run is
API keys — see [API_KEYS.md](API_KEYS.md).

## First Run

`gvf gvf-run` is the one command you need. Point it at **your own gene** — any
HGNC symbol works, and no gold standard is required:

```bash
gvf gvf-run <YOUR_GENE> --email brett.kroncke@gmail.com --output ./results [--disease "<phenotype>"]
```

`gvf-run` runs the regular end-to-end path (discovery → triage → full-text and
supplement acquisition → extraction → SQLite → recovery layers), with source
recovery on by default. The `--email` flag is used for NCBI compliance; it is
applied automatically, so you do not need an `NCBI_EMAIL` line in `.env` just to
run (you do still need one LLM provider key).

To reproduce a known result first, KCNH2 (a well-studied cardiac gene) is a good
validation target:

```bash
gvf gvf-run KCNH2 --email brett.kroncke@gmail.com --output ./results --disease "Long QT Syndrome"
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
gvf extract KCNH2 --email brett.kroncke@gmail.com --output ./results --max-pmids 20 --max-downloads 10

# Verbose output
gvf extract KCNH2 --email brett.kroncke@gmail.com --output ./results --verbose

# Custom output directory
gvf extract KCNH2 --email brett.kroncke@gmail.com --output /path/to/results

# Publish the scored run into the Variant_Browser review DB for collaborator
# adjudication (opt-in, best-effort; needs a sibling Variant_Browser checkout)
gvf gvf-run KCNH2 --email brett.kroncke@gmail.com --output ./results --publish-review
```

To pull collaborator adjudications directly from the live Azure review database
into GVF's local gold SQLite cache, see
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
./results/KCNH2/20260210_143022/
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

SQLite database: ./results/KCNH2/20260210_143022/KCNH2.db
================================================================================
```

## Quick Database Queries

```bash
# Open the database
sqlite3 ./results/KCNH2/20260210_143022/KCNH2.db

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
gvf gvf-run KCNH2 --email brett.kroncke@gmail.com --output ./results --no-source-recovery
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
- For a bounded fetch test, run `gvf extract ... --output ./results --max-downloads 10`; `gvf-run` does not expose that limit
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
