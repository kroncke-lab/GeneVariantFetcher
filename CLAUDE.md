# CLAUDE.md - AI Assistant Architecture Guide

**GeneVariantFetcher: Unified Biomedical Literature Extraction Pipeline**

**Last Updated:** 2025-12-03
**Document Version:** 2.0
**Project Version:** 0.2.0

---

## Table of Contents

1. [Mission & Objective](#1-mission--objective)
2. [Discovery Layer: Dual Source Strategy](#2-discovery-layer-dual-source-strategy)
3. [Content Assembly: Building Complete Context](#3-content-assembly-building-complete-context)
4. [Tiered Filtering Funnel: Cost Efficiency](#4-tiered-filtering-funnel-cost-efficiency)
5. [Data Storage & Schema (CRITICAL)](#5-data-storage--schema-critical)
6. [Gene-Disease Mapping](#6-gene-disease-mapping)
7. [CLI & Production Usage](#7-cli--production-usage)
8. [Architecture & Data Flow](#8-architecture--data-flow)
9. [Key Modules & Components](#9-key-modules--components)
10. [Development Workflows](#10-development-workflows)
11. [Code Conventions](#11-code-conventions)
12. [Common Tasks](#12-common-tasks)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Mission & Objective

### Core Goal

**Extract ALL human carriers of genetic variants from biomedical literature** (including full-text articles and supplemental materials), creating a unified, normalized database of variant-phenotype associations.

### Success Criteria

**Only capture humans with clinical phenotyping that allows classification as:**
- **Affected** - Individual exhibits disease phenotype
- **Unaffected** - Variant carrier without disease phenotype
- **Ambiguous** - Insufficient information to classify

This classification relies on **gene-disease associations** (see [Section 6](#6-gene-disease-mapping)).

### Output

A **single harmonized SQLite database** containing:
- Individual-level carrier data (case reports)
- Aggregate cohort statistics (consortia studies)
- Penetrance calculations across papers
- Full provenance (PMID, evidence quotes)

---

## 2. Discovery Layer: Dual Source Strategy

### Philosophy: Precision + Recall

This repository implements a **merged dual-source discovery engine** combining:

1. **Source A: PubMind** (High Precision)
   - Variant-focused discovery engine
   - Pre-filtered for genetic content
   - Lower false positive rate
   - **Primary discovery mechanism**

2. **Source B: Exhaustive PubMed Search** (High Recall)
   - Comprehensive safety net
   - Catches papers PubMind missed
   - **Critical for supplemental materials** (variants often hidden in Excel tables)
   - Merged from separate "Exhaustive Search" repository

### Why Both Sources?

**Problem:** PubMind may miss papers where:
- Variants are only mentioned in supplemental files
- Gene names use uncommon aliases
- Papers focus on phenotype with variants as secondary data

**Solution:** Exhaustive PubMed search ensures **complete literature coverage**.

### Implementation

```python
# From automated_workflow.py
from gene_literature.pubmind_fetcher import fetch_pmids_for_gene

# Step 1: PubMind discovery (primary)
pmids = fetch_pmids_for_gene(
    gene_symbol="BRCA1",
    email="user@example.com",
    max_results=100
)

# Step 2: Exhaustive search (merged capability)
# Falls back to PubMed if PubMind insufficient
```

**Key Configuration:**
```bash
USE_PUBMIND=true          # Enable PubMind
PUBMIND_ONLY=false        # Allow fallback to PubMed
USE_PUBMED=true           # Enable exhaustive search
```

---

## 3. Content Assembly: Building Complete Context

### The Single .md Context File Requirement

**Critical:** For each PMID, we must create **one unified Markdown file** containing:

1. **Title** (paper title)
2. **Abstract** (complete abstract text)
3. **Full-Text** (complete article body)
4. **All Supplemental Materials** (Excel, Word, PDF → Markdown)

### Why This Matters

**Real-world example:** A paper on SCN5A may have:
- Abstract: "We studied 50 patients..."
- Full-text: Describes clinical phenotypes
- **Supplement Table S1.xlsx**: Contains actual variant calls (c.1234G>A, p.Arg412Gln)

**Without supplements:** We miss 70-80% of extractable variant data.

### Implementation: PMCHarvester

```python
# From harvesting/orchestrator.py
from harvesting import PMCHarvester

harvester = PMCHarvester(output_dir="pmc_fulltext")
harvester.harvest(pmids, delay=2.0)

# Creates: {PMID}_FULL_CONTEXT.md
#   - Article text (from PMC XML)
#   - Supplement S1 (converted from Excel)
#   - Supplement S2 (converted from Word)
```

### Content Assembly Workflow

```
PMID 12345678
    ↓
┌────────────────────────────────┐
│ 1. PMID → PMCID Conversion     │
│    (NCBI E-utilities API)      │
└────────────────────────────────┘
    ↓
┌────────────────────────────────┐
│ 2. Download Full-Text XML      │
│    (PMC Open Access)           │
└────────────────────────────────┘
    ↓
┌────────────────────────────────┐
│ 3. Scrape Supplements          │
│    - Nature/Elsevier/Generic   │
│    - Download .xlsx, .docx     │
└────────────────────────────────┘
    ↓
┌────────────────────────────────┐
│ 4. Convert to Markdown         │
│    - markitdown library        │
│    - Preserve tables           │
└────────────────────────────────┘
    ↓
┌────────────────────────────────┐
│ 5. Unified Context File        │
│    12345678_FULL_CONTEXT.md    │
│    - Title + Abstract          │
│    - Full text                 │
│    - All supplements           │
└────────────────────────────────┘
```

**Output Structure:**
```
automated_output/BRCA1/20251203_120000/pmc_fulltext/
├── 12345678_FULL_CONTEXT.md         ← Single unified file
├── 12345678_supplements/
│   ├── Table_S1.xlsx                ← Raw files preserved
│   └── Figure_S2.pdf
├── successful_downloads.csv          ← Success tracking
└── paywalled_missing.csv            ← Unavailable papers
```

---

## 4. Tiered Filtering Funnel: Cost Efficiency & Precision

### Philosophy: Progressive Cost Reduction

**Problem:** Running expensive LLM extraction on ALL papers is cost-prohibitive.

**Solution:** Multi-tier filtering progressively eliminates irrelevant papers before expensive extraction.

### The Three-Tier Funnel

```
1000 Papers from Discovery
    ↓
┌─────────────────────────────────────────────┐
│ TIER 1: Deterministic Filter (INSTANT)     │
│ Type: Regex/Keyword Matching               │
│ Cost: $0.00 per paper                       │
│ Speed: ~1000 papers/second                  │
│                                              │
│ Filters:                                    │
│ - Non-human studies (mouse, rat, etc.)     │
│ - Review articles (no original data)       │
│ - Computational predictions only            │
│ - Wrong disease/gene context               │
│                                              │
│ Implementation: keyword_filter.py           │
└─────────────────────────────────────────────┘
    ↓ (500 papers remain - 50% filtered)
┌─────────────────────────────────────────────┐
│ TIER 2: Triage LLM (CHEAP)                 │
│ Type: Cheap LLM Classification              │
│ Model: GPT-4o-mini or Claude Haiku          │
│ Cost: ~$0.0001 per paper                    │
│ Speed: ~100 papers/minute                   │
│                                              │
│ Question: "Does this paper contain:         │
│  1. Human subjects                          │
│  2. Clinical phenotype data                 │
│  3. Original patient-level information?"    │
│                                              │
│ Implementation: clinical_data_triage.py     │
└─────────────────────────────────────────────┘
    ↓ (150 papers remain - 70% filtered)
┌─────────────────────────────────────────────┐
│ TIER 3: Expert Extraction (EXPENSIVE)      │
│ Type: Reasoning LLM with Structured Output │
│ Model: GPT-4o, o1-preview, or Claude Opus  │
│ Cost: ~$0.05-0.10 per paper                 │
│ Speed: ~20 papers/minute                    │
│                                              │
│ Task: Extract ALL variant carriers:         │
│ - Individual IDs (P1, Patient 42, etc.)    │
│ - Variants (HGVS notation)                 │
│ - Phenotypes (clinical features)            │
│ - Status (Affected/Unaffected/Ambiguous)   │
│ - Supporting evidence (exact quotes)       │
│                                              │
│ Implementation: automated_workflow.py       │
└─────────────────────────────────────────────┘
    ↓ (150 papers × $0.05 = $7.50)

TOTAL COST: $0.00 + $0.05 + $7.50 = $7.55
vs. 1000 papers × $0.10 = $100.00

85% COST SAVINGS
```

### Tier 1: Deterministic Filter

**Location:** `clinical_data_triage.py` (keyword-based mode)

**Excludes:**
- Animal studies: "mouse", "murine", "rat", "zebrafish"
- Reviews: "systematic review", "meta-analysis"
- Computational only: "in silico", "prediction"
- Wrong abbreviation: Gene symbol used for non-genetic term

**Example:**
```python
# TTR gene (Transthyretin) vs "time to reimbursement"
abstract = "We analyzed TTR (time to reimbursement) for Medicare..."
# ❌ FILTERED OUT - wrong context
```

### Tier 2: Clinical Triage LLM

**Location:** `clinical_data_triage.py`

**Prompt Strategy:**
```
Does this paper contain ORIGINAL human clinical phenotyping data?

Title: [TITLE]
Abstract: [ABSTRACT]

Return JSON:
{
  "has_human_data": true/false,
  "has_clinical_phenotypes": true/false,
  "is_original_data": true/false,  # not review/meta-analysis
  "confidence": 0.0-1.0,
  "reasoning": "brief explanation"
}

Pass only if ALL three are true AND confidence > 0.7
```

**What PASSES Tier 2:**
- Case reports (individual patients)
- Case series (small cohorts with patient details)
- Clinical cohorts with individual-level data
- Functional studies with patient phenotypes

**What FAILS Tier 2:**
- Review articles
- Meta-analyses (aggregate data only)
- Animal-only studies
- Population genetics (no phenotypes)

### Tier 3: Expert Extraction

**Location:** `automated_workflow.py` (extraction step)

**See Section 5 for extraction data model**

---

## 5. Data Storage & Schema (CRITICAL)

### Primary Storage Policy

**🚨 CRITICAL REQUIREMENT:**

**ALL extracted data MUST be written directly to SQLite as papers are processed.**

### Storage Architecture

```
Primary Store:   SQLite Database (canonical)
                      ↓
Intermediate:    JSON files (debugging/caching only)
                      ↓
Utilities:       Re-ingest JSON → SQLite if needed
```

### Why SQLite is Canonical

1. **Normalized schema** prevents duplicates
2. **Relational integrity** maintains data consistency
3. **Efficient querying** for analysis
4. **Portable** single-file database
5. **No separate database server** required

### Database Schema

**Location:** `migrate_to_sqlite.py`

#### Table A: `individuals` (Case Report Data)

**Purpose:** Store individual human variant carriers from case reports/series

```sql
CREATE TABLE individuals (
    individual_id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    patient_identifier TEXT,  -- "P1", "Patient 42", "III-2"
    gene_symbol TEXT NOT NULL,

    -- Variant information
    cdna_notation TEXT,       -- c.1234G>A
    protein_notation TEXT,    -- p.Arg412Gln
    genomic_position TEXT,    -- chr11:g.123456C>T
    zygosity TEXT,            -- heterozygous/homozygous/compound het

    -- Clinical phenotyping
    affected_status TEXT NOT NULL,  -- "affected", "unaffected", "ambiguous"
    primary_phenotype TEXT,         -- "Hypertrophic cardiomyopathy"
    age_at_diagnosis REAL,          -- 45.0 (years)
    sex TEXT,                       -- "male", "female", "unknown"

    -- Evidence
    evidence_quotes TEXT,     -- JSON array of supporting sentences

    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

**Example Record:**
```json
{
  "pmid": "35443093",
  "patient_identifier": "P1",
  "gene_symbol": "SCN5A",
  "cdna_notation": "c.3734A>G",
  "protein_notation": "p.Asp1245Gly",
  "affected_status": "affected",
  "primary_phenotype": "Long QT syndrome",
  "age_at_diagnosis": 12,
  "sex": "female",
  "evidence_quotes": "[\"Patient P1 presented with syncope at age 12...\"]"
}
```

#### Table B: `aggregates` (Cohort/Consortium Data)

**Purpose:** Store summary statistics from large cohort studies

```sql
CREATE TABLE aggregates (
    aggregate_id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    gene_symbol TEXT NOT NULL,

    -- Variant information
    cdna_notation TEXT,
    protein_notation TEXT,

    -- Cohort statistics
    total_carriers INTEGER,
    affected_count INTEGER,
    unaffected_count INTEGER,
    ambiguous_count INTEGER,

    -- Penetrance calculation
    penetrance_percentage REAL,  -- (affected / total) * 100
    confidence_interval_lower REAL,
    confidence_interval_upper REAL,

    -- Cohort description
    cohort_description TEXT,
    study_design TEXT,           -- "prospective", "retrospective"

    FOREIGN KEY (pmid) REFERENCES papers(pmid)
);
```

**Example Record:**
```json
{
  "pmid": "33442691",
  "gene_symbol": "BRCA1",
  "protein_notation": "p.Arg1699Gln",
  "total_carriers": 156,
  "affected_count": 89,
  "unaffected_count": 67,
  "penetrance_percentage": 57.1,
  "cohort_description": "Multi-center breast cancer screening cohort"
}
```

#### Supporting Tables

**`papers` table:**
```sql
CREATE TABLE papers (
    pmid TEXT PRIMARY KEY,
    title TEXT,
    journal TEXT,
    publication_date TEXT,
    doi TEXT,
    gene_symbol TEXT,
    extraction_timestamp TEXT
);
```

**`variants` table (normalized):**
```sql
CREATE TABLE variants (
    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    cdna_notation TEXT,
    protein_notation TEXT,
    genomic_position TEXT,
    UNIQUE(gene_symbol, cdna_notation, protein_notation)
);
```

### SQLite Integration in Workflow

**Location:** `automated_workflow.py` (lines 200-250)

```python
# After extraction completes
logger.info("\n💾 STEP 5: Migrating to SQLite database...")

from migrate_to_sqlite import create_database_schema, insert_extraction_data

db_path = output_path / f"{gene_symbol}.db"
conn = create_database_schema(str(db_path))

# Insert extracted data directly
for extraction_file in extraction_files:
    with open(extraction_file) as f:
        data = json.load(f)
    insert_extraction_data(conn, data, gene_symbol)

conn.commit()
conn.close()

logger.info(f"✓ SQLite database created: {db_path}")
```

### Utilities for JSON Re-ingestion

**If you have legacy JSON files** that need migration:

```bash
# Standalone migration utility
python migrate_to_sqlite.py \
  --data-dir automated_output/BRCA1/20251203_120000/extractions \
  --db BRCA1_reingested.db
```

**This utility:**
1. Scans directory for `*_PMID_*.json` files
2. Parses extraction format
3. Inserts into SQLite with deduplication
4. Handles both `individuals` and `aggregates` data

---

## 6. Gene-Disease Mapping

### Purpose

**Determine if a phenotype = "Affected" for this gene**

### The Mapping File

**Location:** `gene_disease_mapping.json` (to be created)

**Format:**
```json
{
  "BRCA1": {
    "primary_diseases": [
      "breast cancer",
      "ovarian cancer",
      "breast carcinoma",
      "ovarian carcinoma"
    ],
    "aliases": ["breast and ovarian cancer"],
    "exclusions": ["benign breast disease"]
  },
  "SCN5A": {
    "primary_diseases": [
      "long QT syndrome",
      "brugada syndrome",
      "cardiac arrhythmia",
      "sudden cardiac death"
    ]
  },
  "TTR": {
    "primary_diseases": [
      "amyloidosis",
      "familial amyloid polyneuropathy",
      "transthyretin amyloidosis"
    ]
  }
}
```

### How It's Used

**During extraction:**

```python
# Extraction prompt includes gene-disease context
def build_extraction_prompt(paper, gene_symbol, mapping):
    diseases = mapping[gene_symbol]["primary_diseases"]

    prompt = f"""
Extract variant carriers from this paper about {gene_symbol}.

AFFECTED STATUS DEFINITION:
A patient is "affected" if they have one of these diseases:
{', '.join(diseases)}

A patient is "unaffected" if they are a variant carrier BUT do NOT have these diseases.

[Rest of extraction prompt...]
"""
```

**Example Classification:**

**Gene:** BRCA1
**Patient:** 45-year-old female with BRCA1 c.5266dup and breast cancer
**→ Status:** Affected ✓

**Gene:** BRCA1
**Patient:** 30-year-old male with BRCA1 c.5266dup, no cancer, undergoing surveillance
**→ Status:** Unaffected ✓

**Gene:** BRCA1
**Patient:** Patient with BRCA1 variant (no phenotype mentioned)
**→ Status:** Ambiguous ⚠️

---

## 7. CLI & Production Usage

### ⚠️ CRITICAL: Production Entry Point

**DO NOT use `example_*.py` scripts for production!**

**Primary Entry Point:** `automated_workflow.py`

```bash
# Production usage
python automated_workflow.py BRCA1 --email your@email.com

# With limits
python automated_workflow.py SCN5A \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100
```

### Key Flags

```bash
--gene SYMBOL          # Required: Gene symbol (e.g., BRCA1)
--email EMAIL          # Required: For NCBI E-utilities
--max-pmids N          # Limit PMIDs from discovery (default: 100)
--max-downloads N      # Limit full-text downloads (default: 50)
--output-dir PATH      # Output directory (default: automated_output)
--tier-threshold N     # Variant count to trigger next model (default: 1)
```

### What `example_*.py` Scripts Are

**Purpose:** Demonstration and testing only

**Examples:**
- `example_automated_workflow.py` - Thin wrapper around `automated_workflow.py`
- `example_pubmind_only.py` - Shows PubMind API usage
- `example_harvesting.py` - Shows PMCHarvester usage

**DO:**
- ✅ Use as learning resources
- ✅ Reference for API usage patterns
- ✅ Copy code snippets for custom workflows

**DON'T:**
- ❌ Run in production
- ❌ Build pipelines around them
- ❌ Expect them to be maintained

### Query the Database

```bash
# Get statistics
python query_variants_db.py automated_output/BRCA1/20251203_120000/BRCA1.db --stats

# Export to CSV
python query_variants_db.py automated_output/BRCA1/20251203_120000/BRCA1.db --export variants.csv

# SQL query
python query_variants_db.py automated_output/BRCA1/20251203_120000/BRCA1.db \
  --query "SELECT * FROM individuals WHERE affected_status='affected'"
```

---

## 8. Architecture & Data Flow

### Complete End-to-End Pipeline

```
User: "python automated_workflow.py BRCA1 --email user@example.com"
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 1: DISCOVERY                                           │
│ Module: gene_literature/pubmind_fetcher.py                   │
│                                                                │
│ Action: Fetch PMIDs from PubMind/PubMed                      │
│ Output: BRCA1_pmids.txt (100 PMIDs)                          │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 2a: CONTENT ASSEMBLY                                   │
│ Module: harvesting/orchestrator.py (PMCHarvester)            │
│                                                                │
│ Action: Download full-text + supplements from PMC            │
│ Output: pmc_fulltext/{PMID}_FULL_CONTEXT.md files            │
│ Coverage: ~30% of PMIDs (rest paywalled)                     │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 2b: TIER 1 FILTERING (Optional)                        │
│ Module: clinical_data_triage.py (keyword mode)               │
│                                                                │
│ Action: Fast keyword filtering                               │
│ Filters: Non-human, reviews, wrong context                   │
│ Cost: $0.00                                                    │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 2c: TIER 2 FILTERING (Optional)                        │
│ Module: clinical_data_triage.py (LLM mode)                   │
│                                                                │
│ Action: LLM triage for clinical data                         │
│ Question: "Has original human clinical phenotyping?"         │
│ Cost: ~$0.0001/paper                                          │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 3: TIER 3 EXTRACTION                                   │
│ Module: automated_workflow.py (extraction step)              │
│                                                                │
│ Action: Extract variant carriers with reasoning LLM          │
│ Models: GPT-4o → o1-preview (tier cascade)                   │
│ Output: extractions/{PMID}_extraction.json                   │
│ Cost: ~$0.05-0.10/paper                                       │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 4: AGGREGATION                                         │
│ Module: rerun_aggregate_helper.py                            │
│                                                                │
│ Action: Combine penetrance across papers                     │
│ Output: penetrance_aggregated.json                           │
└──────────────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────────────┐
│ STAGE 5: SQLITE MIGRATION (CANONICAL STORAGE)                │
│ Module: migrate_to_sqlite.py                                 │
│                                                                │
│ Action: Write to normalized SQLite database                  │
│ Output: BRCA1.db (canonical data store)                      │
│                                                                │
│ Tables: papers, variants, individuals, aggregates            │
└──────────────────────────────────────────────────────────────┘
    ↓
✓ COMPLETE: automated_output/BRCA1/20251203_120000/BRCA1.db
```

### Repository Structure

```
GeneVariantFetcher/
├── 📄 PRODUCTION ENTRY POINTS
│   ├── automated_workflow.py ......... **MAIN PRODUCTION ENTRY POINT**
│   ├── main.py ....................... Help message
│   ├── migrate_to_sqlite.py .......... Standalone migration utility
│   ├── query_variants_db.py .......... Database query utility
│   ├── rerun_extraction.py ........... Re-run extraction helper
│   ├── rerun_aggregate_helper.py ..... Aggregation utility
│   │
├── 🧬 GENE LITERATURE PACKAGE (Discovery + Metadata)
│   ├── gene_literature/
│   │   ├── pubmind_fetcher.py ........ PubMind API client
│   │   ├── collector.py .............. Literature orchestration
│   │   ├── pubmed_client.py .......... PubMed E-utilities
│   │   ├── synonym_finder.py ......... Gene alias discovery
│   │   ├── relevance_checker.py ...... LLM relevance filtering
│   │   └── writer.py ................. Export (JSON/CSV/URLs)
│   │
├── 📥 HARVESTING MODULES (Full-text Download)
│   ├── harvesting/
│   │   ├── orchestrator.py ........... PMCHarvester (main class)
│   │   ├── pmc_api.py ................ NCBI PMC API
│   │   ├── doi_resolver.py ........... DOI → URL resolution
│   │   ├── supplement_scraper.py ..... Web scraping (Nature/Elsevier)
│   │   └── format_converters.py ...... PDF/DOCX/XLSX → Markdown
│   │
├── 🔬 EXTRACTION & FILTERING
│   ├── clinical_data_triage.py ....... Tier 1 & 2 filtering
│   ├── sourcer.py .................... Discovery coordination
│   │
├── 📊 DATA MODELS & CONFIGURATION
│   ├── models.py ..................... Pydantic data models
│   ├── config.py ..................... Configuration singleton
│   │
├── ⚠️ EXAMPLES (DO NOT USE IN PRODUCTION)
│   ├── example_automated_workflow.py
│   ├── example_pubmind_only.py
│   ├── example_harvesting.py
│   └── [other example_*.py files]
│   │
├── 📚 DOCUMENTATION
│   ├── README.md ..................... User-facing guide
│   ├── CLAUDE.md ..................... **THIS FILE - Architecture**
│   ├── ARCHITECTURE.md ............... System design
│   ├── SIMPLE_USAGE_GUIDE.md ......... Step-by-step workflows
│   └── SQLITE_MIGRATION_GUIDE.md ..... Database schema docs
│   │
├── 📊 OUTPUT DIRECTORY (Created automatically)
│   └── automated_output/
│       └── {GENE}/{TIMESTAMP}/
│           ├── {GENE}_pmids.txt ...... Discovered PMIDs
│           ├── pmc_fulltext/ ......... Full-text + supplements
│           │   └── {PMID}_FULL_CONTEXT.md
│           ├── extractions/ .......... JSON extraction results
│           │   └── {PMID}_extraction.json
│           └── {GENE}.db ............. **CANONICAL SQLITE DATABASE**
│
└── ⚙️ CONFIGURATION
    ├── pyproject.toml ................ Dependencies & metadata
    ├── .env.example .................. Configuration template
    └── .gitignore .................... Git exclusions
```

---

## 9. Key Modules & Components

### automated_workflow.py - Main Production Entry Point

**Purpose:** Orchestrates complete end-to-end workflow

**Key Function:** `automated_variant_extraction_workflow(gene_symbol, email, max_pmids, max_papers_to_download)`

**What it does:**
1. Discovers PMIDs (PubMind/PubMed)
2. Downloads full-text (PMCHarvester)
3. Extracts variant data (multi-model cascade)
4. Aggregates penetrance
5. Migrates to SQLite

**Tier Cascade Strategy:**
```python
models_to_try = ["gpt-4o", "o1-preview", "gpt-4o-mini"]

for model in models_to_try:
    result = extract_with_model(paper, model)
    variant_count = len(result.get("individuals", []))

    if variant_count >= tier_threshold:
        break  # Success, move to next paper
    else:
        logger.info(f"Only {variant_count} variants found, trying next model...")
```

---

### harvesting/orchestrator.py - PMCHarvester

**Purpose:** Download full-text articles + supplements

**Main Class:** `PMCHarvester`

**Critical Methods:**
- `harvest(pmids, delay=2.0)` - Main entry point
- `_get_pmc_id(pmid)` - PMID → PMCID conversion
- `_download_full_text(pmcid)` - Get article XML
- `_download_supplements(pmcid)` - Get supplement files
- `_convert_to_markdown(file_path)` - Convert to .md

**Fallback Chain:**
1. PMC Open Access API (primary)
2. DOI resolution (doi.org)
3. Web scraping (Nature, Elsevier, generic)
4. Direct publisher URL construction

---

### clinical_data_triage.py - Tier 1 & 2 Filtering

**Purpose:** Filter papers before expensive extraction

**Two Modes:**

**Mode 1: Keyword Filter (Tier 1)**
```python
def keyword_filter(abstract):
    # Fast deterministic filtering
    if any(word in abstract.lower() for word in EXCLUDE_KEYWORDS):
        return "FAIL"
    return "PASS"
```

**Mode 2: LLM Triage (Tier 2)**
```python
def llm_triage(title, abstract):
    # Ask cheap LLM: "Has human clinical phenotyping?"
    response = call_llm_mini(build_triage_prompt(title, abstract))
    return response["has_clinical_data"]
```

---

### migrate_to_sqlite.py - Database Migration

**Purpose:** Convert JSON extractions to SQLite

**Main Functions:**
- `create_database_schema(db_path)` - Create tables
- `insert_extraction_data(conn, data, gene)` - Insert records
- `migrate_to_sqlite(data_dir, db_path)` - Full migration

**Auto-detection Features:**
- Finds JSON files in nested directories
- Detects gene symbol from path or filenames
- Auto-names database (e.g., `TTR.db`)

---

### gene_literature/pubmind_fetcher.py - Discovery

**Purpose:** Fetch PMIDs from PubMind/PubMed

**Main Function:** `fetch_pmids_for_gene(gene_symbol, email, max_results)`

**Returns:** List of PMIDs sorted by relevance

---

## 10. Development Workflows

### Setup

```bash
# Clone repository
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher

# Install dependencies (Python 3.11+ required)
pip install -e .

# Configure environment
cp .env.example .env
# Edit .env with your API keys

# Required:
# AI_INTEGRATIONS_OPENAI_API_KEY=your-key
# NCBI_EMAIL=your@email.com
```

### Running Production Workflow

```bash
# Complete workflow
python automated_workflow.py BRCA1 --email your@email.com

# With custom limits
python automated_workflow.py SCN5A \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100 \
  --output-dir custom_output
```

### Git Workflow

**Branch Naming:**
- `main` - Production code
- `claude/{description}-{session-id}` - Claude branches

**Commit Style:**
- Use imperative mood: "Add feature" (not "Added feature")
- Be specific: "Fix PMID parsing in extraction.py:145"
- No periods at end

---

## 11. Code Conventions

### Python Style

- **Python 3.11+** required
- **PEP 8** compliance
- **Type hints** throughout
- **Docstrings** (Google style)
- **Logging** via `logger = logging.getLogger(__name__)`

### Example Pattern

```python
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def extract_variants(
    paper_path: Path,
    gene_symbol: str,
    model: str = "gpt-4o"
) -> dict:
    """Extract genetic variants from a paper.

    Args:
        paper_path: Path to unified markdown file
        gene_symbol: Gene to extract (e.g., "BRCA1")
        model: LLM model identifier

    Returns:
        Dict with individuals and aggregates data

    Raises:
        ValueError: If paper_path doesn't exist
    """
    logger.info(f"Extracting variants from {paper_path} using {model}")
    # Implementation...
```

---

## 12. Common Tasks

### Task 1: Add Support for New Extraction Field

**Scenario:** Extract "family history" for each individual

**Steps:**
1. Update database schema in `migrate_to_sqlite.py`
2. Update extraction prompt in `automated_workflow.py`
3. Update `models.py` Pydantic models
4. Test extraction on sample paper
5. Re-run migration

### Task 2: Add New Publisher to Supplement Scraper

**Scenario:** Add Wiley journal scraping

**Steps:**
1. Add scraper function to `harvesting/supplement_scraper.py`
2. Add domain router logic
3. Test on known Wiley article
4. Handle edge cases (missing supplements, paywalls)

### Task 3: Debug Failed Extractions

**Steps:**
1. Check logs: `grep ERROR automated_output/BRCA1/*/pipeline.log`
2. Inspect paper: Read `{PMID}_FULL_CONTEXT.md`
3. Test extraction manually with single PMID
4. Adjust model or prompt if needed
5. Re-run: `python rerun_extraction.py automated_output/BRCA1/{timestamp}/`

---

## 13. Troubleshooting

### Issue: No PMIDs Found

**Symptoms:** "No PMIDs found for {GENE}"

**Causes:**
- Gene symbol incorrect (check NCBI Gene database)
- PubMind API down
- No literature exists for gene

**Solutions:**
- Verify gene symbol spelling
- Try fallback: `USE_PUBMIND=false USE_PUBMED=true`
- Check PubMed manually

### Issue: All Papers Paywalled

**Symptoms:** `paywalled_missing.csv` has all PMIDs

**Causes:**
- Papers not in PMC Open Access
- Running from cloud IP (blocked by publishers)

**Solutions:**
- Normal for ~70% of papers
- Run from residential network if possible
- Focus analysis on available papers

### Issue: Empty Extraction Results

**Symptoms:** JSON files empty or minimal data

**Causes:**
- Paper has no individual-level data (only aggregate stats)
- Full-text unavailable (abstract-only)
- Wrong model used (too weak)

**Solutions:**
- Check paper passed Tier 2 filter
- Verify `{PMID}_FULL_CONTEXT.md` has content
- Try stronger model (o1-preview)
- Adjust `tier_threshold` to trigger better models

### Issue: SQLite Migration Fails

**Symptoms:** `IntegrityError` or `FOREIGN KEY constraint failed`

**Causes:**
- Malformed JSON
- Missing required fields
- Duplicate entries

**Solutions:**
```bash
# Validate JSON
python -m json.tool extractions/12345678_extraction.json

# Check schema
sqlite3 BRCA1.db ".schema"

# Re-run migration with fresh DB
rm BRCA1.db
python migrate_to_sqlite.py --data-dir extractions/ --db BRCA1.db
```

---

## Quick Reference

### Environment Variables

```bash
# Required
AI_INTEGRATIONS_OPENAI_API_KEY=your-key
NCBI_EMAIL=your@email.com

# Optional
NCBI_API_KEY=your-ncbi-key       # Increases rate limits
USE_PUBMIND=true                 # Enable PubMind
USE_PUBMED=true                  # Enable PubMed fallback
ANTHROPIC_API_KEY=your-key       # For Claude models
```

### CLI Commands

```bash
# Production workflow
python automated_workflow.py GENE --email EMAIL

# Query database
python query_variants_db.py GENE.db --stats
python query_variants_db.py GENE.db --export output.csv

# Re-run extraction
python rerun_extraction.py automated_output/GENE/{timestamp}/

# Migration
python migrate_to_sqlite.py --data-dir extractions/ --db GENE.db
```

### File Locations

```
automated_output/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt              ← PMIDs discovered
├── pmc_fulltext/
│   └── {PMID}_FULL_CONTEXT.md    ← Complete article
├── extractions/
│   └── {PMID}_extraction.json    ← Extracted data
└── {GENE}.db                     ← CANONICAL DATABASE
```

---

## Summary: What Makes This Repository Unique

1. **Dual Discovery:** PubMind (precision) + Exhaustive PubMed (recall)
2. **Complete Context:** Full-text + supplements in single .md file
3. **Tiered Filtering:** 85% cost savings through progressive filtering
4. **SQLite-First:** Direct database writes, not file-based
5. **Clinical Classification:** Affected/Unaffected/Ambiguous via gene-disease mapping
6. **Production-Ready:** `automated_workflow.py` is the entry point (not example scripts)

---

## Contact & Resources

- **Repository:** https://github.com/kroncke-lab/GeneVariantFetcher
- **Issues:** https://github.com/kroncke-lab/GeneVariantFetcher/issues
- **PubMed E-utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **PubMind:** https://pubmind.ai

---

**Last Updated:** 2025-12-03
**Document Version:** 2.0
**Project Version:** 0.2.0

*This is the authoritative architectural guide for AI assistants working on GeneVariantFetcher.*
