# GeneVariantFetcher — Output Format Reference

Detailed specification of all files and formats GVF produces.

## Directory Structure

After running GVF, your output directory contains:

```
{OUTPUT_DIR}/{GENE}/{TIMESTAMP}/
│
├── {GENE}.db                           # SQLite database (primary output)
├── {GENE}_pmids.txt                    # Discovered PMIDs (plain text)
├── {GENE}_pmids_pubmind.txt            # PMIDs from PubMind
├── {GENE}_pmids_pubmed.txt             # PMIDs from PubMed
├── {GENE}_workflow_summary.json        # Run statistics
├── {GENE}_penetrance_summary.json      # Aggregated penetrance data
├── {GENE}_workflow.log                 # Execution log
├── run_manifest.json                   # Execution metadata
│
├── abstract_json/                      # Per-paper metadata & abstracts
│   └── {PMID}.json
│
├── pmid_status/                        # Filter decisions
│   ├── filtered_out.csv                # Papers excluded by filters
│   └── extraction_failures.csv         # Failed extractions
│
├── pmc_fulltext/                       # Downloaded papers
│   ├── {PMID}_FULL_CONTEXT.md          # Unified paper content
│   ├── {PMID}_DATA_ZONES.md            # High-value sections (if scouted)
│   ├── {PMID}_supplements/             # Supplemental files
│   │   ├── Table_S1.xlsx
│   │   └── Figure_S2.pdf
│   ├── successful_downloads.csv        # Download success tracking
│   └── paywalled_missing.csv           # Papers we couldn't access
│
└── extractions/                        # Per-paper extraction JSONs
    └── {GENE}_PMID_{PMID}.json
```

---

## SQLite Database Schema

The `{GENE}.db` database is the primary queryable output. It uses a normalized relational schema.

### Entity Relationship Diagram

```
papers (pmid PK)
  │
  ├──→ variant_papers (variant_id, pmid) ←── variants (variant_id PK)
  │                                              │
  │                                              ├──→ penetrance_data
  │                                              │      └──→ age_dependent_penetrance
  │                                              ├──→ individual_records
  │                                              ├──→ functional_data
  │                                              ├──→ phenotypes
  │                                              └──→ variant_metadata
  │
  ├──→ extraction_metadata
  └──→ tables_processed
```

### Core Tables

#### `papers`
Paper metadata.

```sql
CREATE TABLE papers (
    pmid TEXT PRIMARY KEY,
    title TEXT,
    journal TEXT,
    doi TEXT,
    pmc_id TEXT,
    gene_symbol TEXT,
    extraction_timestamp TEXT
);
```

#### `variants`
Unique variants with HGVS notations.

```sql
CREATE TABLE variants (
    variant_id INTEGER PRIMARY KEY,
    gene_symbol TEXT,
    cdna_notation TEXT,
    protein_notation TEXT,
    genomic_position TEXT,
    clinical_significance TEXT,
    UNIQUE(gene_symbol, protein_notation)
);
```

#### `variant_papers`
Links variants to their source papers (many-to-many).

```sql
CREATE TABLE variant_papers (
    variant_id INTEGER,
    pmid TEXT,
    source_location TEXT,
    key_quotes TEXT,  -- JSON array
    PRIMARY KEY (variant_id, pmid),
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

### Data Tables

#### `penetrance_data`
Cohort-level penetrance statistics.

```sql
CREATE TABLE penetrance_data (
    penetrance_id INTEGER PRIMARY KEY,
    variant_id INTEGER,
    pmid TEXT,
    total_carriers_observed INTEGER,
    affected_count INTEGER,
    unaffected_count INTEGER,
    uncertain_count INTEGER,
    penetrance_percentage REAL,
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

#### `age_dependent_penetrance`
Age-stratified penetrance data.

```sql
CREATE TABLE age_dependent_penetrance (
    age_penetrance_id INTEGER PRIMARY KEY,
    penetrance_id INTEGER,
    age_range TEXT,           -- e.g., "0-20", "20-40"
    penetrance_percentage REAL,
    carriers_in_range INTEGER,
    FOREIGN KEY (penetrance_id) REFERENCES penetrance_data(penetrance_id)
);
```

#### `individual_records`
Individual patient/carrier records.

```sql
CREATE TABLE individual_records (
    record_id INTEGER PRIMARY KEY,
    variant_id INTEGER,
    pmid TEXT,
    age_at_onset INTEGER,
    sex TEXT,
    affected_status TEXT,     -- 'affected', 'unaffected', 'uncertain'
    phenotype_details TEXT,
    family_id TEXT,
    proband_status TEXT,
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

#### `functional_data`
Functional study results.

```sql
CREATE TABLE functional_data (
    functional_id INTEGER PRIMARY KEY,
    variant_id INTEGER,
    pmid TEXT,
    summary TEXT,
    assays TEXT,              -- JSON array
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

#### `phenotypes`
Patient phenotype descriptions.

```sql
CREATE TABLE phenotypes (
    phenotype_id INTEGER PRIMARY KEY,
    variant_id INTEGER,
    pmid TEXT,
    patient_count INTEGER,
    phenotype_description TEXT,
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

### Example Queries

```sql
-- Count variants by clinical significance
SELECT clinical_significance, COUNT(*) as count
FROM variants
GROUP BY clinical_significance;

-- Get all pathogenic variants with penetrance data
SELECT v.protein_notation, v.cdna_notation,
       p.total_carriers_observed, p.affected_count, p.penetrance_percentage
FROM variants v
JOIN penetrance_data p ON v.variant_id = p.variant_id
WHERE v.clinical_significance = 'pathogenic'
ORDER BY p.total_carriers_observed DESC;

-- Find papers with most variant reports
SELECT papers.pmid, papers.title, COUNT(DISTINCT vp.variant_id) as variant_count
FROM papers
JOIN variant_papers vp ON papers.pmid = vp.pmid
GROUP BY papers.pmid
ORDER BY variant_count DESC
LIMIT 10;

-- Age-dependent penetrance for a specific variant
SELECT v.protein_notation, adp.age_range, adp.penetrance_percentage, adp.carriers_in_range
FROM variants v
JOIN penetrance_data pd ON v.variant_id = pd.variant_id
JOIN age_dependent_penetrance adp ON pd.penetrance_id = adp.penetrance_id
WHERE v.protein_notation = 'p.Gly628Ser';
```

---

## JSON Extraction Format

Individual extraction files in `extractions/{GENE}_PMID_{PMID}.json`:

```json
{
  "paper_metadata": {
    "pmid": "12345678",
    "title": "Paper Title",
    "extraction_summary": "Brief summary of what was extracted"
  },
  "variants": [
    {
      "gene_symbol": "KCNH2",
      "cdna_notation": "c.1883G>A",
      "protein_notation": "p.Gly628Ser",
      "genomic_position": null,
      "clinical_significance": "pathogenic",
      "patients": {
        "count": 5,
        "demographics": "3 male, 2 female",
        "phenotype": "long QT syndrome"
      },
      "penetrance_data": {
        "total_carriers_observed": 12,
        "affected_count": 8,
        "unaffected_count": 4,
        "uncertain_count": 0,
        "penetrance_percentage": 66.7,
        "age_dependent_penetrance": [
          {
            "age_range": "0-20",
            "penetrance_percentage": 45.0,
            "carriers_in_range": 4
          },
          {
            "age_range": "20-40",
            "penetrance_percentage": 80.0,
            "carriers_in_range": 8
          }
        ]
      },
      "individual_records": [
        {
          "age_at_onset": 15,
          "sex": "M",
          "affected_status": "affected",
          "phenotype_details": "syncope, QTc 520ms"
        }
      ],
      "functional_data": {
        "summary": "Patch-clamp studies showed 60% reduction in IKr current",
        "assays": ["patch-clamp", "Western blot"]
      },
      "segregation_data": "Segregates with disease in 3 families",
      "population_frequency": "0.0001 (gnomAD)",
      "evidence_level": "strong",
      "source_location": "Results, Table 2",
      "additional_notes": "Founder mutation in Finnish population",
      "key_quotes": [
        "G628S was identified in 12 family members across 3 generations"
      ]
    }
  ],
  "tables_processed": [
    "Table 1: Patient characteristics",
    "Supplementary Table S2: Variant details"
  ],
  "extraction_metadata": {
    "total_variants_found": 1,
    "extraction_confidence": "high",
    "study_type": "cohort",
    "challenges": [],
    "notes": "Clear phenotype-genotype data presented"
  }
}
```

### Field Definitions

| Field | Type | Description |
|-------|------|-------------|
| `cdna_notation` | string | cDNA-level HGVS notation (e.g., c.1883G>A) |
| `protein_notation` | string | Protein-level HGVS notation (e.g., p.Gly628Ser) |
| `clinical_significance` | string | pathogenic, likely_pathogenic, VUS, benign, etc. |
| `penetrance_data` | object | Carrier counts and penetrance statistics |
| `individual_records` | array | Per-patient data when available |
| `functional_data` | object | In vitro/functional study results |
| `evidence_level` | string | strong, moderate, weak, supporting |
| `key_quotes` | array | Verbatim quotes supporting the extraction |

---

## FULL_CONTEXT.md Format

Papers are converted to Markdown for consistent LLM input:

```markdown
# Title of the Paper

## Abstract

The abstract text goes here...

## Introduction

Introduction text...

## Methods

### Study Population

Details about the cohort...

### Genetic Analysis

Sequencing methodology...

## Results

### Table 1: Patient Characteristics

| Patient ID | Age | Sex | Variant | Phenotype |
|------------|-----|-----|---------|-----------|
| P001 | 25 | M | G628S | LQTS |
| P002 | 32 | F | G628S | Asymptomatic |

### Variant Findings

Text describing variants...

## Discussion

Discussion text...

## References

1. First reference...
2. Second reference...

---

## Supplementary Materials

### Supplementary Table S1

[Content from Excel file converted to Markdown table]

### Supplementary Figure S1

[Description or placeholder for figures]
```

---

## Workflow Summary Format

`{GENE}_workflow_summary.json`:

```json
{
  "gene_symbol": "KCNH2",
  "run_timestamp": "2026-02-10T14:30:00",
  "duration_seconds": 2400,
  "pmids_discovered": 245,
  "pmids_filtered_out": 89,
  "pmids_passed_filters": 156,
  "papers_downloaded": 78,
  "papers_extracted": 72,
  "extraction_failures": 6,
  "total_variants_found": 234,
  "variants_with_penetrance": 156,
  "total_carriers": 1245,
  "total_affected": 678,
  "success_rate_percent": 32.0,
  "output_directory": "/path/to/output/KCNH2/20260210_143000"
}
```

---

## Penetrance Summary Format

`{GENE}_penetrance_summary.json`:

```json
{
  "gene_symbol": "KCNH2",
  "generated_at": "2026-02-10T14:45:00",
  "total_unique_variants": 234,
  "variants_with_penetrance_data": 156,
  "aggregated_variants": [
    {
      "protein_notation": "p.Gly628Ser",
      "cdna_notation": "c.1883G>A",
      "total_carriers_across_papers": 45,
      "total_affected_across_papers": 32,
      "total_unaffected_across_papers": 13,
      "aggregated_penetrance_percent": 71.1,
      "source_paper_count": 8,
      "source_pmids": ["12345678", "23456789", "..."]
    }
  ],
  "summary_statistics": {
    "mean_penetrance": 65.4,
    "median_penetrance": 68.2,
    "min_penetrance": 15.0,
    "max_penetrance": 100.0
  }
}
```

---

## Run Manifest Format

`run_manifest.json`:

```json
{
  "run_id": "a1b2c3d4-e5f6-7890-abcd-ef1234567890",
  "gene_symbol": "KCNH2",
  "start_timestamp": "2026-02-10T14:00:00",
  "end_timestamp": "2026-02-10T14:40:00",
  "duration_seconds": 2400,
  "status": "completed",
  "config": {
    "max_pmids": null,
    "max_papers_to_download": null,
    "tier_threshold": 1,
    "auto_synonyms": true,
    "scout_first": false
  },
  "statistics": {
    "pmids_discovered": 245,
    "pmids_filtered_out": 89,
    "papers_downloaded": 78,
    "papers_extracted": 72,
    "total_variants_found": 234
  },
  "output_locations": {
    "sqlite_database": "/path/to/KCNH2.db",
    "penetrance_summary": "/path/to/KCNH2_penetrance_summary.json",
    "workflow_summary": "/path/to/KCNH2_workflow_summary.json"
  },
  "errors": [],
  "warnings": [
    "6 papers failed extraction due to timeout"
  ]
}
```

---

## CSV Tracking Files

### successful_downloads.csv

```csv
pmid,source,download_time,file_size,has_supplements
12345678,pmc,2026-02-10T14:15:00,45678,true
23456789,elsevier,2026-02-10T14:16:30,38291,false
```

### paywalled_missing.csv

```csv
pmid,doi,publisher,reason
34567890,10.1234/example,wiley,no_api_key
45678901,10.5678/another,unknown,no_oa_version
```

### filtered_out.csv

```csv
pmid,filter_stage,reason,confidence
11111111,keyword,insufficient_clinical_terms,0.95
22222222,llm,review_article_not_primary_data,0.87
```

---

## Best Practices for Using Output

### Querying the Database

```python
import sqlite3
import pandas as pd

# Connect to database
conn = sqlite3.connect('KCNH2.db')

# Load all variants with pandas
df = pd.read_sql_query("""
    SELECT v.protein_notation, v.clinical_significance,
           p.total_carriers_observed, p.affected_count
    FROM variants v
    LEFT JOIN penetrance_data p ON v.variant_id = p.variant_id
""", conn)

# Export to CSV
df.to_csv('kcnh2_variants.csv', index=False)

conn.close()
```

### Validating Extractions

Always spot-check critical variants:

1. Find the source PMID in `variant_papers`
2. Review `key_quotes` for the extraction
3. Open `pmc_fulltext/{PMID}_FULL_CONTEXT.md` to verify

### Combining Multiple Runs

If you run GVF multiple times for the same gene:

```bash
# Each run creates a timestamped directory
output/KCNH2/20260210_143000/
output/KCNH2/20260215_091500/

# Use the most recent, or merge databases manually
```

---

## Next Steps

- [VALIDATION.md](VALIDATION.md) — Understand recall metrics
- [ARCHITECTURE.md](ARCHITECTURE.md) — Technical deep-dive
- [QUICKSTART.md](QUICKSTART.md) — Get running
