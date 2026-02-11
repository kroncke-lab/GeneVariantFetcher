# Plan to Achieve Near 100% Variant Extraction Coverage

## Current Status
- **Test PMIDs:** 50
- **Expected variants (Excel):** 665
- **Extracted variants:** 81
- **Current coverage:** 12.2%

---

## Gap Analysis Summary

| Issue Type | PMIDs Affected | Variants Lost | Priority |
|------------|----------------|---------------|----------|
| No full-text downloaded | 16 | ~280 | HIGH |
| Downloaded junk (not real paper) | 9 | ~85 | HIGH |
| Partial extraction (missing variants) | 9 | ~120 | MEDIUM |
| Over-extraction (wrong gene?) | 2 | +8 extra | LOW |

---

## PRIORITY 1: Fix Full-Text Retrieval (16 PMIDs, ~280 variants)

### Problem
16 papers have no FULL_CONTEXT.md file. These are paywalled or failed to download.

**Affected PMIDs:** 15840476, 10973849, 14661677, 29650123, 27871843, 12402336, 10862094, 21216356, 22727609, 14998624, 17210839, 15913580, 19695459, 29759541, 29766883, 9544837, 26823142, 29146210

### Solutions

#### A. Use Browser-Based Fetching (Recommended)
```bash
# Generate paywalled CSV for browser fetch
python cli/browser_fetch.py \
    tests/test_run_output/KCNH2/20260125_113957/pmc_fulltext/paywalled_missing.csv \
    --gene KCNH2 \
    --wait-for-captcha \
    --auto-process
```

**Updates needed in `cli/browser_fetch.py`:**
1. Add retry logic for Cloudflare protection
2. Add support for institutional login (university VPN detection)
3. Improve PDF-to-markdown conversion for scanned PDFs

#### B. Manual Download Fallback
For papers that still fail:
1. Create `tests/manual_downloads/` directory
2. Download PDFs manually from:
   - Sci-Hub (for research purposes)
   - University library access
   - Author preprints on ResearchGate/bioRxiv
3. Add script to convert and integrate manual downloads

#### C. Abstract-Only Extraction
For papers that cannot be obtained:
```python
# Enable abstract-only extraction for paywalled papers
# Update pipeline/extraction.py to:
# 1. Detect when full-text is unavailable
# 2. Fall back to abstract + title extraction
# 3. Flag results as "abstract_only" with lower confidence
```

---

## PRIORITY 2: Fix Junk Downloads (9 PMIDs, ~85 variants)

### Problem
9 papers have FULL_CONTEXT.md files but they contain garbage:
- Generic web pages ("We've detected you are running an older version...")
- Unrelated health information pages
- Database landing pages (not article content)

**Affected PMIDs:** 26496715, 11854117, 16922724, 17905336, 19996378, 15176425, 24103226, 28532774, 19843919

### Solutions

#### A. Improve Content Validation in `harvesting/orchestrator.py`
```python
def _validate_extracted_content(self, content: str, pmid: str) -> bool:
    """Reject content that isn't actual paper text."""

    # Reject common junk patterns
    JUNK_PATTERNS = [
        "We've detected you are running an older version",
        "BrowseHappy.com",
        "PhosphoSitePlus",
        "Uniprot Database Entry",
        "Sign in to access",
        "Subscribe to read",
    ]

    for pattern in JUNK_PATTERNS:
        if pattern.lower() in content.lower():
            return False

    # Require minimum meaningful content
    if len(content) < 2000:
        return False

    # Require paper-like structure
    paper_indicators = ['abstract', 'methods', 'results', 'conclusion', 'references']
    matches = sum(1 for ind in paper_indicators if ind in content.lower())
    if matches < 2:
        return False

    return True
```

#### B. Add Content Quality Scoring
```python
def _score_content_quality(self, content: str) -> float:
    """Score 0-1 indicating likelihood this is real paper content."""
    score = 0.0

    # Check for scientific structure
    if 'abstract' in content.lower(): score += 0.15
    if 'methods' in content.lower(): score += 0.15
    if 'results' in content.lower(): score += 0.15
    if 'discussion' in content.lower(): score += 0.15
    if 'references' in content.lower(): score += 0.1

    # Check for variant mentions
    import re
    variant_pattern = r'[A-Z]\d+[A-Z]|c\.\d+[ACGT]>[ACGT]|p\.[A-Z][a-z]{2}\d+'
    variant_count = len(re.findall(variant_pattern, content))
    score += min(0.3, variant_count * 0.03)

    return score
```

#### C. Re-fetch with Better Sources
For failed content, try alternative sources in order:
1. PubMed Central XML (most reliable)
2. Europe PMC
3. Publisher API (Elsevier, Wiley, Springer)
4. DOI resolution + web scraping
5. Unpaywall API

---

## PRIORITY 3: Improve Extraction Quality (9 PMIDs, ~120 variants)

### Problem
Papers have valid full-text but extraction is incomplete:
- 24606995: Found 2/28 variants (7%)
- 24667783: Found 3/23 variants (13%)
- 23098067: Found 8/22 variants (36%)

### Root Causes Identified

1. **Table parsing failures** - Large tables with variant data not fully extracted
2. **Supplementary data missed** - Variants in supplements not processed
3. **DATA_ZONES too aggressive** - Condensation removes important sections (✅ ADDRESSED: scanner now runs on full text)
4. **Variant nomenclature** - Different formats not matched

### Solutions

#### A. Improve Table Extraction in `pipeline/extraction.py`
```python
# Add dedicated table extraction pass
def extract_from_tables(self, content: str) -> List[Dict]:
    """Specifically target table content for variant extraction."""

    # Find markdown tables
    table_pattern = r'\|[^\n]+\|[\n\r]+\|[-:| ]+\|[\n\r]+((?:\|[^\n]+\|[\n\r]+)+)'
    tables = re.findall(table_pattern, content)

    variants = []
    for table in tables:
        # Parse table rows
        # Look for variant notation patterns
        # Extract associated phenotype data
        pass

    return variants
```

#### B. Process Supplementary Files Better
```python
# In harvesting/orchestrator.py
def _process_supplements(self, pmid: str, supplement_files: List[Path]):
    """Extract variant data from supplementary materials."""

    for supp_file in supplement_files:
        if supp_file.suffix in ['.xlsx', '.xls']:
            # Parse Excel supplements - often contain variant tables
            df = pd.read_excel(supp_file)
            # Look for columns like 'Variant', 'Mutation', 'cDNA', 'Protein'

        elif supp_file.suffix == '.docx':
            # Parse Word supplements
            # Often contain additional variant lists
            pass
```

#### C. Less Aggressive DATA_ZONES Condensation

**✅ ADDRESSED (2026-02-10):** The variant scanner is now decoupled from DATA_ZONES condensation. Scanner runs on FULL_CONTEXT.md to ensure no variants are missed by condensation.

```python
# In pipeline/data_scout.py
# Increase minimum relevance threshold for KEEPING sections
# Currently: min_relevance=0.3
# Change to: min_relevance=0.1 for test papers

# Also: Keep all table content regardless of relevance score
def should_keep_zone(self, zone: Zone) -> bool:
    if zone.type == 'TABLE':
        return True  # Always keep tables
    return zone.relevance_score >= self.min_relevance
```

#### D. Multi-Pass Extraction
```python
def extract_with_fallback(self, paper: Paper) -> ExtractionResult:
    """Try multiple extraction strategies."""

    # Pass 1: Standard extraction from DATA_ZONES
    result1 = self.extract_from_zones(paper)

    # Pass 2: If low yield, try FULL_CONTEXT
    if result1.variant_count < 5:
        result2 = self.extract_from_fulltext(paper)
        result1.merge(result2)

    # Pass 3: Table-specific extraction
    result3 = self.extract_from_tables(paper.full_text)
    result1.merge(result3)

    # Pass 4: Supplement extraction
    if paper.supplements:
        result4 = self.extract_from_supplements(paper.supplements)
        result1.merge(result4)

    return result1
```

---

## PRIORITY 4: Fix Variant Nomenclature Matching

### Problem
Excel uses short notation (A561V), extraction uses HGVS (p.A561V).
Need to normalize for comparison.

### Solution
```python
# In tests/compare_variants.py or utils/variant_utils.py

def normalize_variant(variant: str) -> str:
    """Normalize variant notation for comparison."""

    # Remove p. prefix
    variant = re.sub(r'^p\.', '', variant)

    # Convert 3-letter to 1-letter amino acids
    AA_MAP = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    }

    for three, one in AA_MAP.items():
        variant = variant.replace(three, one)

    # Handle stop codons
    variant = variant.replace('Ter', 'X').replace('*', 'X')

    return variant.upper()

def variants_match(v1: str, v2: str) -> bool:
    """Check if two variant notations refer to same variant."""
    return normalize_variant(v1) == normalize_variant(v2)
```

---

## Implementation Roadmap

### Phase 1: Quick Wins (1-2 days)
1. [ ] Add content validation to reject junk downloads
2. [ ] Re-run download for 9 junk PMIDs with better validation
3. [ ] Add variant normalization to comparison script
4. [x] ~~Reduce DATA_ZONES min_relevance~~ — Scanner now decoupled from DATA_ZONES (runs on full text)

### Phase 2: Browser Fetching (2-3 days)
1. [ ] Run browser_fetch.py for 16 paywalled PMIDs
2. [ ] Add Cloudflare bypass improvements
3. [ ] Integrate manual download workflow

### Phase 3: Extraction Improvements (3-5 days)
1. [ ] Implement table-specific extraction
2. [ ] Add supplement parsing for Excel files
3. [ ] Implement multi-pass extraction
4. [ ] Add extraction quality scoring

### Phase 4: Validation & Tuning (2-3 days)
1. [ ] Re-run full pipeline on test PMIDs
2. [ ] Compare results to Excel
3. [ ] Tune parameters based on results
4. [ ] Document final accuracy

---

## Expected Outcome

After implementing all phases:

| Metric | Current | Target |
|--------|---------|--------|
| PMIDs with extraction | 18/50 | 45/50 |
| Variants extracted | 81 | 600+ |
| Coverage | 12.2% | 90%+ |

---

## Files to Modify

1. `harvesting/orchestrator.py` - Content validation, junk detection
2. `pipeline/extraction.py` - Multi-pass extraction, table parsing
3. `pipeline/data_scout.py` - Less aggressive condensation
4. `cli/browser_fetch.py` - Cloudflare handling, retry logic
5. `utils/variant_utils.py` (new) - Variant normalization
6. `tests/compare_variants.py` - Improved comparison with normalization
