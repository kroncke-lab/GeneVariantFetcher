# GeneVariantFetcher — Technical Architecture

A deep dive into GVF's pipeline architecture, module responsibilities, and extension points.

## Pipeline Overview

GVF is a multi-stage pipeline that transforms a gene symbol into a structured database of variants and patient data extracted from the biomedical literature.

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         GeneVariantFetcher Pipeline                              │
└─────────────────────────────────────────────────────────────────────────────────┘

INPUT: Gene Symbol (e.g., "KCNH2")
  │
  ├──────────────────────────────────────────────────────────────────────────────┐
  │ STEP 0: Synonym Discovery (optional)                                         │
  │   • NCBI Gene Database → gene aliases                                         │
  │   • Example: KCNH2 → ["HERG", "LQT2", "Kv11.1"]                               │
  └──────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1: PMID Discovery                                                           │
│   • PubMind API → {gene}_pmids_pubmind.txt                                       │
│   • PubMed E-Utilities → {gene}_pmids_pubmed.txt                                 │
│   OUTPUT: {gene}_pmids.txt (merged, deduplicated)                                │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1.5: Fetch Abstracts                                                        │
│   • PubMed E-Utilities (efetch)                                                  │
│   OUTPUT: abstract_json/{PMID}.json                                              │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1.6: Filter Papers (Two-Tier)                                               │
│   Tier 1: KeywordFilter (regex, fast) — ~65 clinical/variant keywords            │
│   Tier 2: InternFilter (LLM, gpt-4o-mini) — relevance classification             │
│   OUTPUT: pmid_status/filtered_out.csv                                           │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 2: Download Full-Text (PMCHarvester)                                        │
│   Sources (priority order):                                                      │
│     1. PMC OA (Open Access) → BioC XML                                           │
│     2. Publisher APIs (Elsevier, Springer, Wiley) → XML/PDF                      │
│     3. Unpaywall → OA PDF links                                                  │
│     4. CORE API → aggregated OA content                                          │
│   OUTPUT: pmc_fulltext/{PMID}_FULL_CONTEXT.md                                    │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 2.5: Data Scout (optional)                                                  │
│   • Identifies high-value data zones (tables, methods sections)                  │
│   • Creates condensed context for LLM extraction                                 │
│   OUTPUT: pmc_fulltext/{PMID}_DATA_ZONES.md                                      │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 3: Variant Extraction (ExpertExtractor)                                     │
│   • Input: DATA_ZONES.md > FULL_CONTEXT.md > abstract                            │
│   • Model: gpt-4o-mini → gpt-4o (escalation if needed)                           │
│   • Pre-scan: Regex scanner on FULL_CONTEXT.md (not condensed text)              │
│     - Detects concatenated gene+variant (HERGG604S, KCNH2A561V)                  │
│     - Unicode arrow normalization (→, ➔, ⟶)                                      │
│   OUTPUT: extractions/{gene}_PMID_{pmid}.json                                    │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 4: Aggregation                                                              │
│   • HGVS variant name normalization                                              │
│   • Fuzzy matching for variant deduplication                                     │
│   • Cross-paper penetrance aggregation                                           │
│   OUTPUT: {gene}_penetrance_summary.json                                         │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 5: SQLite Migration                                                         │
│   • Normalized relational schema (8 core tables)                                 │
│   • Indexed for efficient querying                                               │
│   OUTPUT: {gene}.db                                                              │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
FINAL OUTPUTS:
  • {gene}.db                       (SQLite database)
  • {gene}_workflow_summary.json    (run statistics)
  • {gene}_penetrance_summary.json  (aggregated data)
  • run_manifest.json               (execution metadata)
```

## Module Responsibilities

### Core Pipeline Modules

| Module | Location | Responsibility |
|--------|----------|----------------|
| **CLI Entry** | `cli/__init__.py` | Typer app, command registration |
| **Automated Workflow** | `cli/automated_workflow.py` | End-to-end orchestration |
| **Pipeline Steps** | `pipeline/steps.py` | Individual step implementations |

### Discovery & Filtering

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Discovery** | `gene_literature/discovery.py` | PMID collection from multiple sources |
| **PubMind Fetcher** | `gene_literature/pubmind_fetcher.py` | PubMind variant-specific search |
| **PubMed Client** | `gene_literature/pubmed_client.py` | NCBI E-Utilities interface |
| **Synonym Finder** | `gene_literature/synonym_finder.py` | NCBI Gene synonym lookup |
| **Filters** | `pipeline/filters.py` | Keyword + LLM paper filtering |

### Harvesting (Full-Text Acquisition)

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Orchestrator** | `harvesting/orchestrator.py` | Multi-source download coordination |
| **PMC API** | `harvesting/pmc_api.py` | PubMed Central BioC access |
| **Elsevier API** | `harvesting/elsevier_api.py` | ScienceDirect downloads |
| **Springer API** | `harvesting/springer_api.py` | SpringerLink + Nature downloads |
| **Wiley API** | `harvesting/wiley_api.py` | Wiley Online Library access |
| **Format Converters** | `harvesting/format_converters.py` | XML/PDF/DOCX → Markdown |

### Extraction & Analysis

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Extraction** | `pipeline/extraction.py` | LLM-based variant extraction |
| **Prompts** | `pipeline/prompts.py` | LLM prompt templates |
| **Data Scout** | `pipeline/data_scout.py` | High-value zone identification |
| **Aggregation** | `pipeline/aggregation.py` | Cross-paper data combination |

### Utilities

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Variant Normalizer** | `utils/variant_normalizer.py` | HGVS standardization |
| **Variant Scanner** | `utils/variant_scanner.py` | Regex-based variant detection (runs on full text, not condensed) |
| **LLM Utils** | `utils/llm_utils.py` | OpenAI/LiteLLM interface |
| **SQLite Migration** | `harvesting/migrate_to_sqlite.py` | JSON → SQLite conversion |

## Data Flow Detail

### Stage 1: PMID → Text

```
Gene Symbol
    │
    ├─→ PubMind API ──────────→ PMIDs (gene-variant focused)
    │
    ├─→ PubMed E-Utils ───────→ PMIDs (broader search)
    │
    └─→ Merge & Deduplicate ──→ {gene}_pmids.txt
                                    │
                                    ▼
                            Abstract Fetch (E-Utils)
                                    │
                                    ▼
                            abstract_json/{PMID}.json
                                    │
                                    ▼
                            Relevance Filtering
                            (Keyword + LLM tiers)
                                    │
                                    ▼
                            Filtered PMID list
```

### Stage 2: PMID → Full Text

```
Filtered PMIDs
    │
    ├─→ PMC OA Check ─────────→ BioC XML (if available)
    │
    ├─→ DOI Resolution ───────→ Publisher identification
    │       │
    │       ├─→ Elsevier API ─→ ScienceDirect XML
    │       ├─→ Springer API ─→ SpringerLink HTML/PDF
    │       └─→ Wiley API ────→ Wiley XML/PDF
    │
    ├─→ Unpaywall ────────────→ OA PDF location
    │
    └─→ CORE API ─────────────→ Aggregated OA content
            │
            ▼
    Format Conversion (XML/PDF/HTML → Markdown)
            │
            ▼
    Supplement Processing (Excel, Word, PDFs)
            │
            ▼
    {PMID}_FULL_CONTEXT.md (unified document)
```

### Stage 3: Text → Variants

```
{PMID}_FULL_CONTEXT.md
    │
    ├─→ Variant Scanner ──────→ Pre-detected variant patterns
    │   (runs on FULL text)      - Concatenated: HERGG604S, KCNH2A561V
    │                            - Unicode arrows normalized (→, ➔, ⟶)
    │
    ├─→ Data Scout (optional)─→ {PMID}_DATA_ZONES.md
    │                           (condensed high-value sections)
    │
    └─→ ExpertExtractor ──────→ LLM extraction
            │                   (gpt-4o-mini or gpt-4o)
            │
            ▼
    Structured JSON extraction
    {
      "paper_metadata": {...},
      "variants": [
        {
          "protein_notation": "p.Gly628Ser",
          "cdna_notation": "c.1883G>A",
          "clinical_significance": "pathogenic",
          "penetrance_data": {...},
          "functional_data": {...}
        }
      ],
      "extraction_metadata": {...}
    }
```

### Stage 4: Variants → Database

```
extractions/{gene}_PMID_*.json
    │
    ├─→ Variant Normalization ─→ Standardized names
    │   (frameshift, nonsense, protein notation)
    │
    ├─→ Fuzzy Matching ────────→ Variant deduplication
    │   (handles notation variants)
    │
    └─→ Aggregation ───────────→ {gene}_penetrance_summary.json
            │
            ▼
    SQLite Migration
            │
            ▼
    {gene}.db (normalized relational schema)
```

## External API Integration

### API Rate Limits & Authentication

| API | Auth Required | Rate Limit | Notes |
|-----|---------------|------------|-------|
| **PubMed E-Utils** | Email (key optional) | 3/sec (10/sec with key) | NCBI_EMAIL required |
| **PubMind** | None | Courteous use | Variant-focused search |
| **PMC OA** | None | 3/sec | Free full-text |
| **Elsevier** | API key | 5/sec | Requires registration |
| **Springer** | API key | Variable | Free for researchers |
| **Wiley** | API key | Variable | TDM agreement needed |
| **Unpaywall** | Email | 100k/day | OA link resolution |
| **OpenAI** | API key | Per plan | Required for extraction |

### Retry & Circuit Breaker

GVF implements exponential backoff with jitter:

```python
@dataclass
class RetryConfig:
    max_retries: int = 3
    base_delay: float = 1.0
    max_delay: float = 60.0
    exponential_base: float = 2.0
    jitter: bool = True
```

Circuit breakers prevent cascading failures:
- Skip extractions on files < 500 chars
- Skip files with < 30% alphanumeric content
- Timeout extraction after 1200 seconds

## Checkpoint/Resume System

GVF supports resumable jobs via the checkpoint system:

```python
# Pipeline steps (in order)
class PipelineStep(Enum):
    PENDING = "pending"
    DISCOVERING_SYNONYMS = "discovering_synonyms"
    FETCHING_PMIDS = "fetching_pmids"
    FETCHING_ABSTRACTS = "fetching_abstracts"
    FILTERING_PAPERS = "filtering_papers"
    DOWNLOADING_FULLTEXT = "downloading_fulltext"
    SCOUTING_DATA = "scouting_data"
    EXTRACTING_VARIANTS = "extracting_variants"
    AGGREGATING_DATA = "aggregating_data"
    MIGRATING_DATABASE = "migrating_database"
    COMPLETED = "completed"
    FAILED = "failed"
```

Checkpoints are stored in `~/.gvf_jobs/{job_id}/checkpoint.json` with atomic writes and file locking.

## Extension Points

### Adding a New Publisher API

1. Create `harvesting/{publisher}_api.py`:
   ```python
   class PublisherAPI:
       def can_handle(self, doi: str) -> bool:
           """Return True if this publisher can handle the DOI."""
           pass
       
       def download(self, doi: str, output_dir: Path) -> Optional[Path]:
           """Download and return path to content, or None on failure."""
           pass
   ```

2. Register in `harvesting/orchestrator.py`:
   ```python
   self.publishers = [
       ElsevierAPI(),
       SpringerAPI(),
       WileyAPI(),
       YourNewPublisherAPI(),  # Add here
   ]
   ```

### Adding a New Gene Configuration

Add to `config/gene_config.py`:
```python
GENE_CONFIGS = {
    "YOUR_GENE": {
        "aliases": ["ALIAS1", "ALIAS2"],
        "chromosome": "chr1",
        "coding_positions": (start, end),
        "domains": [
            {"name": "Domain1", "start": 1, "end": 100},
        ],
    },
}
```

### Customizing Extraction Prompts

Modify `pipeline/prompts.py`:
```python
EXTRACTION_PROMPT = """
Your custom extraction instructions here.
Focus on specific data types relevant to your use case.
"""
```

### Adding New Output Formats

Extend `harvesting/migrate_to_sqlite.py` or create a new exporter:
```python
def export_to_csv(db_path: Path, output_dir: Path):
    """Export SQLite database to CSV files."""
    pass
```

## Performance Considerations

### Typical Resource Usage

| Gene Size | Papers | RAM | Disk | Time |
|-----------|--------|-----|------|------|
| Small (rare) | 10-50 | 2 GB | 500 MB | 10-20 min |
| Medium | 50-200 | 4 GB | 2 GB | 30-60 min |
| Large (BRCA1) | 200-500 | 8 GB | 5 GB | 1-3 hours |

### Optimization Tips

1. **Parallel downloads**: Orchestrator fetches 3 papers concurrently by default
2. **SSD storage**: Significantly improves SQLite migration speed
3. **API keys**: Publisher keys unlock 2-3x more papers
4. **Data Scout**: `--scout-first` reduces LLM token usage for long papers

## Known Limitations

1. **Gene-specific tuning**: Best performance on cardiac genes (KCNH2, SCN5A, etc.)
2. **Variant aliases**: Only KCNH2 has a comprehensive alias dictionary
3. **PDF extraction**: Tables in scanned PDFs may not extract cleanly
4. **LLM hallucination**: Occasional false positives in extraction (validate critical findings)

## Further Reading

- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) — Database schema and file formats
- [VALIDATION.md](VALIDATION.md) — Recall measurement methodology
- [API_KEYS.md](API_KEYS.md) — Obtaining API credentials
