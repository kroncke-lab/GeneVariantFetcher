# Simple Usage Guide

## The Complete Workflow (One Command!)

Running the entire pipeline from gene name to SQLite database is now **incredibly simple**:

```bash
python automated_workflow.py BRCA1 --email your@email.com --output /path/to/output
```

That's it! This single command does **everything**:

### What Happens Step-by-Step

#### Step 1: Fetch Papers (PubMind)
```
📚 Discovering relevant papers from PubMind...
✓ Found 152 PMIDs for BRCA1
✓ Saved PMID list to: /path/to/output/BRCA1/20231126_143022/BRCA1_pmids.txt
```

#### Step 2: Download Full-Text (PMC)
```
📥 Downloading full-text papers from PubMed Central...
✓ Successfully downloaded 45 full-text papers
```
- Converts PMID → PMCID automatically
- Downloads complete article XML
- Downloads ALL supplemental files (Excel, PDF, Word)
- Converts everything to Markdown
- Creates unified `{PMID}_FULL_CONTEXT.md` files

#### Step 3: Extract Variant Data (AI)
```
🧬 Extracting variant and patient data using AI...
✓ Extracted 12 variants from PMID 35443093
✓ Extracted 8 variants from PMID 33442691
✓ Completed 45/45 papers
```
- Uses OpenAI GPT-4 to extract structured data
- Finds genetic variants (c.1234G>A, p.Arg412Gln)
- Extracts penetrance data (affected/unaffected carriers)
- Captures individual patient records

#### Step 4: Aggregate Penetrance
```
📊 Aggregating penetrance data across papers...
✓ Aggregated penetrance data for 127 variants
```
- Combines variant data across all papers
- Calculates total carriers per variant
- Validates penetrance statistics

#### Step 5: Migrate to SQLite ⭐ NEW!
```
💾 Migrating data to SQLite database...
✓ Migrated 45/45 extractions to SQLite
✓ Database saved to: /path/to/output/BRCA1/20231126_143022/BRCA1.db
```
- Normalizes JSON data into relational tables
- Creates variant, penetrance, patient, and paper tables
- Compresses data for efficient storage

#### Step 6: Final Summary
```
================================================================================
WORKFLOW COMPLETE!
================================================================================
Gene: BRCA1
PMIDs discovered: 152
Papers downloaded: 45
Papers with extractions: 45
Total variants found: 127
Variants with penetrance data: 89
Total carriers observed: 1,234
Total affected carriers: 456
Success rate: 29.6%

💾 Database migrated: 45/45 extractions

All outputs saved to: /path/to/output/BRCA1/20231126_143022
Summary report: /path/to/output/BRCA1/20231126_143022/BRCA1_workflow_summary.json
Penetrance summary: /path/to/output/BRCA1/20231126_143022/BRCA1_penetrance_summary.json
SQLite database: /path/to/output/BRCA1/20231126_143022/BRCA1.db
================================================================================
```

---

## Output Structure

After running, you'll have:

```
/path/to/output/BRCA1/20231126_143022/
├── BRCA1.db                           ← SQLite database (query this!)
├── BRCA1_pmids.txt                    ← List of PMIDs found
├── BRCA1_workflow_summary.json        ← Overall statistics
├── BRCA1_penetrance_summary.json      ← Penetrance aggregation
├── extractions/                       ← Individual paper extractions
│   ├── BRCA1_PMID_35443093.json
│   ├── BRCA1_PMID_33442691.json
│   └── ... (one per paper)
└── pmc_fulltext/                      ← Full-text papers
    ├── PMID_35443093_FULL_CONTEXT.md
    ├── PMID_35443093_supplements/
    │   ├── Table_S1.xlsx
    │   └── Figure_S2.pdf
    └── ... (one per paper)
```

---

## Querying Your Database

Once you have the SQLite database, query it:

```bash
# Get overall statistics
python query_variants_db.py /path/to/output/BRCA1/20231126_143022/BRCA1.db --stats

# Search for a specific variant
python query_variants_db.py /path/to/output/BRCA1/20231126_143022/BRCA1.db --variant "c.1234G>A"

# Get penetrance for a variant
python query_variants_db.py /path/to/output/BRCA1/20231126_143022/BRCA1.db --penetrance "p.Arg412Gln"
```

---

## Customization Options

### Limit Number of Papers
```bash
# Only process first 10 PMIDs (for testing)
python automated_workflow.py BRCA1 --email your@email.com --max-pmids 10 --max-downloads 5
```

### Different Gene
```bash
python automated_workflow.py SCN5A --email your@email.com
python automated_workflow.py TP53 --email your@email.com
python automated_workflow.py KCNQ1 --email your@email.com
```

### Change Output Directory
```bash
python automated_workflow.py BRCA1 --email your@email.com --output my_results
```

### Verbose Logging
```bash
python automated_workflow.py BRCA1 --email your@email.com --verbose
```

---

## Required Setup

### 1. OpenAI API Key (Required)
```bash
export OPENAI_API_KEY="sk-..."
```

Or create a `.env` file:
```env
OPENAI_API_KEY=sk-...
```

### 2. Email (Required for NCBI)
Must provide via `--email` flag for PubMed API compliance.

---

## What Changed? (Simplification)

**Before**: Multiple confusing scripts, manual steps, duplicated code
```bash
# Old way (manual, confusing):
python pubmind_fetcher.py BRCA1            # Step 1
python harvest_pmc_fulltext.py             # Step 2 (edit script first!)
python pipeline.py                         # Step 3 (complex tiered system)
python migrate_to_sqlite.py --data-dir ... # Step 4 (manual migration)
```

**After**: One command does everything
```bash
# New way (automatic, simple):
python automated_workflow.py BRCA1 --email your@email.com
```

**Code cleanup**:
- ❌ Deleted `extractor.py` (duplicate)
- ❌ Deleted `pipeline/harvesting.py` (duplicate)
- ❌ Deleted `pipeline/utils/llm_utils.py` (duplicate)
- ❌ Archived `pipeline.py` → `pipeline_tiered_old.py` (overly complex)
- ✅ Kept `automated_workflow.py` (clean, simple, complete)
- ✅ Integrated SQLite migration (no manual step needed)

**Results**:
- 956 lines of code deleted
- 68 lines added
- **Net: -888 lines** (much simpler!)

---

## Troubleshooting

### "No PMIDs found"
- Normal for very rare genes
- Try a different gene or check PubMind status

### "No papers downloaded"
- Only ~30% of PubMed articles have full-text in PMC
- This is normal and expected
- The workflow continues with available papers

### "OpenAI API error"
- Check your API key is set: `echo $OPENAI_API_KEY`
- Verify you have credits: https://platform.openai.com/usage
- Check for rate limits in error message

### "Email required"
- Must provide `--email your@email.com` for NCBI compliance

---

## Cost Estimation

**PubMind/PMC**: Free ✓

**OpenAI API**: ~$0.01-0.10 per paper
- Depends on paper length
- Full-text with supplements = longer = more expensive
- Typical cost for 50 papers: $2-5

**Example**:
```
Gene: BRCA1
Papers extracted: 45
Estimated cost: $3.50
```

---

## Advanced: Re-running Extraction Only

If you already have full-text papers downloaded and just want to re-extract:

```bash
python rerun_extraction.py /path/to/output/BRCA1/20231126_143022/pmc_fulltext
```

This skips steps 1-2 and only runs extraction + aggregation + SQLite migration.

---

## Questions?

- ✅ **Simple workflow**: Use `automated_workflow.py`
- ✅ **Query database**: Use `query_variants_db.py`
- ✅ **Re-extract**: Use `rerun_extraction.py`
- ❌ **DON'T use**: `pipeline.py` (archived as `pipeline_tiered_old.py`)
