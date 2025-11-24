# Tiered Biomedical Extraction Pipeline

A cost-efficient, multi-tiered pipeline for extracting genetic variant data from scientific literature using LLMs.

## 🎯 Overview

This pipeline implements a **tiered filtering strategy** to minimize costs while maximizing data quality:

```
┌─────────────────────────────────────────────────────────────────┐
│  SOURCER: Query PubMed/EuropePMC APIs                          │
│  → Returns deduplicated PMIDs for a gene symbol                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  TIER 1: Keyword Filter (Fast & Free)                          │
│  → Checks abstracts for clinical keywords                       │
│  → ❌ DISCARD papers with no clinical relevance                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  TIER 2: Intern Filter (gpt-4o-mini, ~$0.0001/paper)          │
│  → LLM classifies: "Original Data" vs "Review/Irrelevant"      │
│  → ❌ DISCARD reviews and non-clinical papers                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  TIER 3: Expert Extractor (gpt-4o, ~$0.05/paper)              │
│  → Deep extraction of variants, phenotypes, clinical data       │
│  → Handles full text + markdown tables                          │
│  → Returns structured JSON with all variant information         │
└─────────────────────────────────────────────────────────────────┘
```

### 💰 Cost Savings

By filtering papers early, this pipeline can save **80-95% in LLM costs** compared to processing all papers with the expensive Tier 3 extractor.

**Example**: Processing 100 papers for BRCA1
- Without filtering: 100 papers × $0.05 = **$5.00**
- With tiered filtering: 10 papers × $0.05 + 90 papers × $0.0001 = **$0.51** (90% savings!)

## 📦 Installation

```bash
# Clone the repository
git clone <your-repo-url>
cd GeneVariantFetcher

# Install dependencies
pip install -r requirements.txt

# Set up API keys
export OPENAI_API_KEY="your-openai-key"
# OR
export ANTHROPIC_API_KEY="your-anthropic-key"
```

## 🚀 Quick Start

```python
from pipeline import run_pipeline_for_gene

# Run pipeline for a gene symbol
results = run_pipeline_for_gene(
    gene_symbol="BRCA1",
    email="your_email@example.com",
    max_papers=10
)

# Access results
print(f"Variants found: {results['stats'].total_variants_extracted}")
print(f"Cost savings: ${results['stats'].estimated_cost_savings:.2f}")
```

## 📚 Module Documentation

### `models.py`
Pydantic data models for type safety and validation:
- `Paper`: Scientific paper with metadata
- `FilterResult`: Result from a filter tier
- `ExtractionResult`: Result from expert extraction
- `PipelineResult`: Complete pipeline result for a paper

### `sourcer.py`
Paper sourcing from APIs:
- `PaperSourcer`: Queries PubMed and EuropePMC
- `query_papers_for_gene()`: Convenience function
- Returns deduplicated PMIDs
- Includes retry logic with exponential backoff

### `filters.py`
Two-tiered filtering system:

#### `KeywordFilter` (Tier 1)
- Fast, rule-based filtering
- Checks abstracts for clinical keywords
- Configurable keyword list and threshold
- Cost: **$0** (no API calls)

```python
from filters import KeywordFilter

filter = KeywordFilter(min_keyword_matches=3)
result = filter.filter(paper)
```

#### `InternFilter` (Tier 2)
- LLM-based classification using `gpt-4o-mini`
- Distinguishes original clinical data from reviews
- Cost: **~$0.0001 per paper**

```python
from filters import InternFilter

filter = InternFilter(model="gpt-4o-mini")
result = filter.filter(paper)
```

### `extractor.py`
Expert-level variant extraction:

#### `ExpertExtractor` (Tier 3)
- Deep extraction using advanced models (GPT-4o, Claude)
- Handles full text + markdown tables
- Extracts structured variant data with phenotypes
- Cost: **~$0.05 per paper**

```python
from extractor import ExpertExtractor

extractor = ExpertExtractor(model="gpt-4o")
result = extractor.extract(paper)
```

### `pipeline.py`
Main orchestrator that connects everything:

#### `BiomedicalExtractionPipeline`
- Runs all tiers sequentially
- **DISCARDS papers at each tier** to save costs
- Comprehensive logging of drops and reasons
- Tracks statistics and cost savings

```python
from pipeline import BiomedicalExtractionPipeline

pipeline = BiomedicalExtractionPipeline(
    email="your_email@example.com"
)

results = pipeline.run(
    gene_symbol="TP53",
    max_papers=50
)
```

## 🔍 Detailed Logging

The pipeline logs exactly which paper was dropped at which tier:

```
🚀 STARTING PIPELINE FOR GENE: BRCA1
📚 STEP 1: Sourcing papers from PubMed/EuropePMC...
Found 247 unique PMIDs for BRCA1

--- Processing paper 1/247 - PMID 12345678 ---
✓ PASSED - PMID 12345678 through Tier 1: Keyword Filter
✓ PASSED - PMID 12345678 through Tier 2: Intern Filter (LLM)
✓ EXTRACTED - PMID 12345678: 5 variants found

--- Processing paper 2/247 - PMID 87654321 ---
🚫 DROPPED - PMID 87654321 at Tier 1: Keyword Filter: Only 1 keyword match(es), minimum 2 required

--- Processing paper 3/247 - PMID 11223344 ---
✓ PASSED - PMID 11223344 through Tier 1: Keyword Filter
🚫 DROPPED - PMID 11223344 at Tier 2: Intern Filter (LLM): Review article, not original clinical data

📊 PIPELINE SUMMARY FOR BRCA1
================================================================================
Total papers sourced: 247
Passed Tier 1 (Keyword): 89
Passed Tier 2 (Intern): 23
Completed Tier 3 (Extraction): 23

❌ Dropped at Tier 1: 158
❌ Dropped at Tier 2: 66
❌ Failed at Tier 3: 0

✓ Successful extractions: 23
🧬 Total variants extracted: 127
💰 Estimated cost savings: $11.23
⏱️ Total processing time: 245.3 seconds
================================================================================
```

## 📊 Output Contract

The extractor should emit one JSON array **per PMID**. Each element in that array represents an individual who can be linked to a specific variant. Keep the keys consistent so downstream consumers can rely on them:

```json
[
  {
    "individual_id": "12345678_case1",
    "pmid": "12345678",
    "gene": "BRCA1",
    "variants": {
      "hgvs": "NM_007294.4:c.5266dupC",
      "protein": "p.Gln1756Profs*74",
      "genomic": "GRCh38:17:43091434:dupC",
      "rsid": "rs80357906"
    },
    "age": "37",
    "sex": "female",
    "ancestry": "Ashkenazi Jewish",
    "penetrance_exposures": ["prior chest radiation"],
    "phenotypes_hpo": ["HP:0003002"],
    "phenotypes_raw": ["early-onset breast cancer"],
    "affected_status": "affected",
    "family_level": false,
    "data_from": "full_text",
    "evidence": {
      "sentence": "The proband (II-2) carried BRCA1 c.5266dupC and was diagnosed with breast cancer at age 37.",
      "section": "Case description"
    }
  }
]
```

Guidance for generating the JSON:

- **Link-level precision**: Only include an individual when the variant-to-person link is explicit. If a phenotype is only provided at the pedigree level or the variant cannot be tied to a specific person, omit the entry or set `family_level: true`.
- **Variants**: Prefer HGVS; include genomic coordinates and rsIDs when present, but do not guess missing details or re-interpret notation.
- **Phenotypes**: Capture raw phrasing in `phenotypes_raw` and map to HPO codes in `phenotypes_hpo` when possible. Include age, sex, ancestry, and any penetrance-relevant exposures.
- **Affected status**: Use `affected`, `unaffected`, or `uncertain`. Mark carriers as unaffected only when they have clearly passed the risk window without disease; otherwise choose `affected` or `uncertain` and explain in `evidence.sentence`.
- **Provenance**: Each data point should cite the exact sentence (or table cell text) plus the section label. Flag non-narrative sources with `data_from` values such as `table`, `supplement`, or `full_text`.
- **Output shape**: Multi-PMID runs should concatenate multiple arrays sequentially (one per PMID), with stable `PMID_caseX` identifiers per article.

## 🧭 Recommended Execution Flow

To reach full coverage (including supplements), the harvester and extraction pipeline should follow this sequence:

1. Convert each PMID to PMCID when available, then pull the full-text XML plus any linked supplemental files.
2. If the EuropePMC supplemental API fails, resolve the DOI, follow redirects, and try site-specific scraping before doing a generic keyword/link scan.
3. Parse every narrative section, figure legend, table, and supplemental sheet to detect individuals (e.g., "proband", "II-2", "Case 3"), using a consistent `PMID_caseX` naming scheme.
4. For each individual, capture all described variants and associated phenotypes, demographics, exposures, and affected status with the provenance sentence/section.
5. Log failures gracefully (missing supplements, DOI errors, etc.) and continue processing remaining PMIDs so each successful paper yields a fully structured JSON array ready for downstream analysis.

## 🛠️ Customization

### Custom Keyword Filter

```python
from filters import KeywordFilter

custom_keywords = [
    "mutation", "variant", "pathogenic",
    "patient", "clinical", "phenotype"
]

keyword_filter = KeywordFilter(
    keywords=custom_keywords,
    min_keyword_matches=3
)
```

### Different LLM Models

```python
from pipeline import BiomedicalExtractionPipeline
from filters import InternFilter
from extractor import ExpertExtractor

pipeline = BiomedicalExtractionPipeline(
    intern_filter=InternFilter(model="gpt-4o-mini"),
    expert_extractor=ExpertExtractor(model="claude-3-opus-20240229")
)
```

### Disable Tiers for Testing

```python
# Only run filtering, skip expensive extraction
pipeline = BiomedicalExtractionPipeline(
    enable_tier3=False
)

results = pipeline.run("BRCA1", max_papers=100)
print(f"Papers that would be extracted: {results['stats'].passed_tier2_intern}")
```

## 📈 Performance Tips

1. **Start small**: Test with `max_papers=10` first
2. **Monitor costs**: Check `stats.estimated_cost_savings`
3. **Adjust thresholds**: Tune `min_keyword_matches` based on your needs
4. **Use cheaper models**: Try `gpt-4o-mini` for extraction if budget is tight
5. **Cache results**: Save extraction results to avoid re-processing

## 🧪 Testing

```bash
# Run with a small test set
python example_usage.py

# Or in Python:
from pipeline import run_pipeline_for_gene

results = run_pipeline_for_gene(
    gene_symbol="BRCA1",
    max_papers=5  # Small test
)
```

## ⚠️ Important Notes

### API Requirements
- **PubMed**: Requires email address (no API key needed)
- **LiteLLM**: Requires `OPENAI_API_KEY` or `ANTHROPIC_API_KEY`

### Rate Limiting
- PubMed: 3 requests/second without API key, 10/second with
- LiteLLM: Depends on your provider's limits

### Full Text Access
- Current implementation uses abstracts (full text requires additional setup)
- To enable full text: integrate with PMC Open Access or institutional subscriptions

## 🤝 Contributing

Contributions welcome! Areas for improvement:
- [ ] Full-text fetching from PMC
- [ ] PDF parsing for papers without PMC access
- [ ] Additional filter strategies
- [ ] Database integration for caching
- [ ] Web interface

## 📄 License

[Your chosen license]

## 🙏 Acknowledgments

- Built with [LiteLLM](https://github.com/BerriAI/litellm) for model abstraction
- Uses NCBI E-utilities and EuropePMC APIs
- Pydantic for robust data validation

## 📞 Support

For issues or questions:
- Open an issue on GitHub
- Check the example_usage.py for common patterns
- Review logs for debugging (all drops are logged with reasons)

---

**Happy variant hunting! 🧬🔬**
