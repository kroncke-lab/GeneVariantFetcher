# GeneVariantFetcher: Technical Architecture

## Overview

GeneVariantFetcher (GVF) is an LLM-powered pipeline for extracting individual-level genetic variant data from biomedical literature. It addresses a critical gap: most published variant data is in supplementary materials (Excel tables, PDFs) that traditional text mining cannot access.

## Pipeline Architecture

### Stage 1: Discovery & Synonym Resolution
**Module:** `pipeline/steps.py::discover_synonyms()`

- Queries NCBI Gene database for official synonyms
- Example: KCNH2 → HERG, LQT2, ERG1, Kv11.1
- Expands search coverage by 5-10x

### Stage 2: PMID Fetching
**Module:** `pipeline/steps.py::fetch_pmids()`

Three parallel sources:
1. **PubMind** - Semantic search via NCBI
2. **PubMed** - Standard MeSH/keyword search  
3. **Europe PMC** - European literature + preprints

Typical yield: 500-2000 PMIDs per gene

### Stage 3: Tiered Filtering
**Module:** `pipeline/filters.py`

| Tier | Name | Method | Speed | Cost |
|------|------|--------|-------|------|
| 1 | KeywordFilter | Regex matching | <1ms/paper | Free |
| 2 | InternFilter | gpt-4o-mini | ~1s/paper | $0.0001 |
| 3 | ExpertExtractor | gpt-4o | ~10s/paper | $0.01 |

**Tier 1 Keywords:**
- Clinical: patient, carrier, proband, family, pedigree
- Phenotype: affected, unaffected, asymptomatic, symptomatic
- Variant: mutation, variant, pathogenic, penetrance

**Tier 2 LLM Prompt:**
"Does this abstract describe individual patients or carriers with variants in {gene}?"

Typical reduction: 2000 PMIDs → 200-400 relevant papers

### Stage 4: Full-Text Download
**Module:** `harvesting/orchestrator.py`

Sources (in priority order):
1. **PMC Open Access** - Direct API download
2. **Elsevier API** - ScienceDirect (requires API key)
3. **Wiley API** - Wiley Online Library
4. **Browser fallback** - Puppeteer for remaining

**Supplement handling:**
- Detects Excel, Word, PDF supplements
- Downloads and converts all formats to text
- Extracts tables with structure preserved

### Stage 5: Genetic Data Scout 🔑
**Module:** `pipeline/data_scout.py`

**Key Innovation:** Before sending full papers to LLM, the Scout identifies "data zones" - sections containing patient-level data.

Detection patterns:
- Patient identifiers: "Patient 1", "Case-62", "II-1" (pedigree)
- Clinical phrases: "presented with", "diagnosed with"
- Variant density: c.XXX, p.XXX notation counts
- Table structure: rows with patient data

**Output:** `{PMID}_DATA_ZONES.md` - condensed file with only high-value content

**Token reduction:** 60-80% (15,000 tokens → 3,000 tokens)

**Note:** The variant scanner runs on FULL_CONTEXT.md (not DATA_ZONES), ensuring no variants are missed due to condensation.

### Stage 6: LLM Extraction
**Module:** `pipeline/extraction.py`

**Prompt engineering highlights:**

1. **Functional vs Clinical study detection:**
   - Functional: "n cells", "replicates", "trafficking assay"
   - Clinical: "patients", "carriers", "probands"
   - Different extraction logic for each

2. **Penetrance data extraction:**
   - Total carriers observed
   - Affected count (with disease)
   - Unaffected count (carriers without disease)

3. **HGVS notation standardization:**
   - Converts legacy formats: R412H → p.Arg412His
   - Handles frameshift: G24fs+34X → p.Gly24fs*58
   - Unicode arrow normalization (→, ➔, ⟶ → standard)
   - Concatenated gene+variant detection (HERGG604S → HERG G604S)

**Model cascade:**
- First try: gpt-4o-mini (fast, cheap)
- If <threshold variants found: retry with gpt-4o
- Threshold adjustable (default: 1)

### Stage 7: Aggregation
**Module:** `pipeline/aggregation.py`

- Validates HGVS notation syntax
- Deduplicates variants across papers
- Merges carrier counts from multiple sources
- Computes aggregate penetrance estimates

**Output:**
- JSON per paper: `{PMID}_extracted.json`
- Aggregated CSV: `{gene}_aggregated.csv`
- SQLite database: `{gene}_variants.db`

## Key Technical Decisions

### Why LLMs instead of NLP/regex?

1. **Table understanding**: LLMs can interpret unstructured tables where column headers vary
2. **Context integration**: Combines information from multiple sentences/paragraphs
3. **Variant normalization**: Converts diverse notation formats to HGVS
4. **Phenotype extraction**: Understands clinical descriptions

### Why 3-tier filtering?

Cost optimization:
- Tier 1 (free): Eliminates 80% of papers
- Tier 2 ($0.0001/paper): Eliminates 50% of remainder
- Tier 3 ($0.01/paper): Only runs on ~10% of original PMIDs

Total cost: ~$5-20 per gene (vs $200+ if running Tier 3 on all)

### Why Data Scout before extraction?

Token economics:
- Full paper: 15,000-50,000 tokens
- Data zones only: 2,000-5,000 tokens
- Cost reduction: 60-80%
- Quality improvement: Less noise, focused extraction

## Current Performance

**KCNH2 benchmark:**
- Input: 6060 PMIDs discovered
- After Tier 1: 1847 papers
- After Tier 2: 749 papers downloaded
- Extracted: 2529 variants (2138 unique)
- Recall vs gold standard: 39.6% (on papers we have)

**Main bottleneck:** 70% of gold-standard papers are paywalled

## Future Development

1. **Institutional access integration**: VPN/proxy for paywalled content
2. **Multi-gene batching**: Process 5-10 genes in parallel
3. **Active learning**: Use extraction feedback to improve filtering
4. **GUI improvements**: Real-time progress, better error handling
