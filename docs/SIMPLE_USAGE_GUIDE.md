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

Once you have the SQLite database, query it directly with SQLite:

```bash
# Open database in SQLite CLI
sqlite3 /path/to/output/BRCA1/20231126_143022/BRCA1.db

# Get variant counts
sqlite> SELECT COUNT(*) FROM variants;

# List all variants with penetrance data
sqlite> SELECT protein_notation, gene_symbol FROM variants WHERE variant_id IN (SELECT variant_id FROM penetrance_data);

# Get affected/unaffected counts for a variant
sqlite> SELECT affected_status, COUNT(*) FROM individual_records WHERE variant_id = 1 GROUP BY affected_status;
```

See [SQLITE_MIGRATION_GUIDE.md](SQLITE_MIGRATION_GUIDE.md) for more query examples.

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

## Advanced: Re-running Extraction Only

If you already have full-text papers downloaded and just want to re-extract, use the GUI's "Folder Jobs" feature:

1. Launch the GUI: `python main.py`
2. Go to the **New Job** tab
3. Select "Process existing folder" mode
4. Point to your existing output directory
5. Choose to start at the extraction step

This skips discovery and download, running only extraction + aggregation + SQLite migration.

---

## Questions?

- **Simple workflow**: Use `automated_workflow.py` or `python main.py` (GUI)
- **Query database**: Use SQLite directly or see [SQLITE_MIGRATION_GUIDE.md](SQLITE_MIGRATION_GUIDE.md)
- **Re-extract**: Use GUI "Folder Jobs" feature
- **Export to CSV**: Use `python scripts/extract_ttr_to_csv.py`
