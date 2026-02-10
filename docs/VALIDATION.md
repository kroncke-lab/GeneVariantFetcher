# GeneVariantFetcher — Validation & Recall Metrics

How GVF measures extraction quality, current performance numbers, and how to run your own validation.

## Validation Philosophy

GVF is designed for **high recall** (finding as many variants as possible) rather than perfect precision. In clinical genetics research, missing a variant is often worse than having a false positive that can be filtered later.

**Key trade-off:** We accept some noise in exchange for not missing real variants.

---

## Gold Standard Methodology

### What Is a Gold Standard?

A gold standard is a manually curated set of variants that we know exist in a set of papers. We compare GVF's extractions against this truth to measure recall.

### KCNH2 Gold Standard

Our primary validation uses a curated KCNH2 dataset:

| Metric | Count |
|--------|-------|
| **Unique variants** | 489 |
| **Total carriers** | 780 |
| **Affected individuals** | 413 |
| **Source papers (PMIDs)** | 236 |

This gold standard was created by:
1. Systematic literature review by domain experts
2. Manual extraction from each paper
3. Variant name normalization
4. Cross-validation by multiple reviewers

### Limitations of Our Gold Standard

1. **Single gene**: Only validated on KCNH2 (cardiac)
2. **Historical data**: May miss very recent publications
3. **Interpretation variability**: Some edge cases have ambiguous classification
4. **Incomplete for rare variants**: Low-frequency variants may be underrepresented

---

## Metrics Definitions

### 1. Unique Variant Recall

**What it measures:** Did GVF find each distinct variant?

```
Unique Variant Recall = (Variants found by GVF ∩ Gold standard variants) / Gold standard variants
```

**Example:**
- Gold standard has 489 unique variants
- GVF extracts 289 of them
- Recall = 289/489 = **59.1%**

**Why it matters:** This is the primary metric for variant discovery use cases.

### 2. Carrier Recall

**What it measures:** Did GVF identify the correct number of carriers?

```
Carrier Recall = (Carriers found by GVF) / (Total carriers in gold standard)
```

**Example:**
- Gold standard has 780 total carriers
- GVF identifies 556 carriers
- Recall = 556/780 = **71.3%**

**Why it matters:** Important for penetrance calculations that require accurate counts.

### 3. Affected Individual Recall

**What it measures:** Did GVF correctly identify affected carriers?

```
Affected Recall = (Affected found by GVF) / (Total affected in gold standard)
```

**Example:**
- Gold standard has 413 affected individuals
- GVF identifies 279 affected
- Recall = 279/413 = **67.6%**

**Why it matters:** Critical for understanding phenotype-genotype relationships.

---

## Current Performance (KCNH2)

As of version 2.1.0 (February 2026):

| Metric | Recall | Notes |
|--------|--------|-------|
| **Unique variants** | 59.1% | 289 of 489 variants found |
| **Carriers** | 71.3% | 556 of 780 carriers identified |
| **Affected** | 67.6% | 279 of 413 affected found |

### Historical Progression

| Version | Date | Unique Recall | Key Changes |
|---------|------|---------------|-------------|
| Baseline | 2026-01-28 | 19.4% | Initial extraction |
| v2 | 2026-02-04 | 29.4% | Basic normalization fixes |
| v3 | 2026-02-06 | 54.6% | Springer API, table extraction |
| v6 | 2026-02-10 | 50.3% | Tier 1 + Tier 2 fixes |
| **v8** | **2026-02-10** | **59.1%** | **Complete rebuild + scanner** |

### Recall Ceiling Analysis

Not all gold standard variants can be found:

| Category | Variants | Percentage |
|----------|----------|------------|
| In downloaded papers | 480 | 98.2% |
| In papers we couldn't access | 9 | 1.8% |

**Theoretical maximum recall: 98.2%** (if extraction were perfect)

The ~40% gap between current recall (59.1%) and theoretical maximum (98.2%) represents:
- LLM extraction errors
- Variant notation mismatches
- Complex/unusual variant formats
- Tables that don't parse cleanly

---

## Why Some Variants Are Missed

### 1. Notation Mismatches

Gold standard might have:
```
p.Leu987fs
```

But the paper contains:
```
L987fsX10
p.Leu987PhefsTer10
c.2959delC (p.L987fs)
```

**Mitigation:** GVF includes fuzzy matching and normalization, but edge cases exist.

### 2. Complex Variant Types

Harder to extract:
- Splice site variants: `c.1558-1G>C`
- Large deletions: `exon6-14Del`
- Structural variants: `7q34q36.2Del`

**Mitigation:** Variant scanner pre-detects these patterns to improve LLM prompts.

### 3. Paywalled Papers

Some high-value papers aren't accessible:
- No PMC version
- No publisher API key
- Institutional access required

**Mitigation:** Obtain publisher API keys; consider institutional VPN.

### 4. Table Extraction Failures

Supplemental Excel files sometimes:
- Have merged cells that break parsing
- Use non-standard variant columns
- Are scanned images, not text

**Mitigation:** Manual review of `paywalled_missing.csv` for critical papers.

---

## Running Your Own Validation

### Step 1: Create a Gold Standard File

Create a CSV with known variants:

```csv
variant,pmid,carriers,affected
p.Gly628Ser,12345678,5,3
p.Ala561Val,23456789,12,8
p.Arg534Cys,34567890,3,2
```

### Step 2: Run GVF Extraction

```bash
gvf extract YOUR_GENE --email you@example.com --output ./validation_output
```

### Step 3: Compare Results

Use the provided validation script:

```bash
python scripts/analysis/recall_analysis.py \
    --gold-standard your_gold_standard.csv \
    --extractions ./validation_output/YOUR_GENE/*/extractions/ \
    --output recall_results.json
```

### Step 4: Interpret Results

```python
import json

with open('recall_results.json') as f:
    results = json.load(f)

print(f"Unique Variant Recall: {results['recall_metrics']['unique_variants']['percentage']:.1f}%")
print(f"Carrier Recall: {results['recall_metrics']['carriers']['percentage']:.1f}%")
print(f"Matched Variants: {results['recall_metrics']['unique_variants']['matched']}")
print(f"Missed Variants: {len(results['variants']['missed'])}")

# Review missed variants
for v in results['variants']['missed'][:10]:
    print(f"  - {v}")
```

### Step 5: Gap Analysis

Identify why variants were missed:

```python
# Check if missed variants are in downloaded papers
gap = results.get('gap_analysis', {})
print(f"Missed in downloaded papers: {gap.get('missed_in_downloaded', 'N/A')}")
print(f"Missed in inaccessible papers: {gap.get('missed_in_undownloaded', 'N/A')}")

# Top missed variants
for item in gap.get('top_missed_in_downloaded', [])[:5]:
    print(f"  {item['variant']}: {item['carriers']} carriers in PMIDs {item['pmids_have']}")
```

---

## Improving Recall

### For Your Own Runs

1. **Add publisher API keys** — More papers = more variants found
2. **Use `--scout-first`** — Better context for complex papers
3. **Check paywalled_missing.csv** — Manually obtain critical papers
4. **Review extraction_failures.csv** — Re-run failed papers individually

### Common Issues and Fixes

| Issue | Symptom | Solution |
|-------|---------|----------|
| Low paper coverage | <50 papers downloaded | Add Elsevier/Springer keys |
| Notation mismatches | Variants found but not matching | Check variant normalizer |
| Table extraction | Data in tables missed | Review FULL_CONTEXT.md manually |
| Timeout errors | Papers fail extraction | Increase LLM timeout in settings |

### Contributing Improvements

If you improve recall:
1. Document your changes
2. Run validation before/after
3. Submit a pull request with metrics

---

## Precision Considerations

While GVF optimizes for recall, precision matters too:

### Expected False Positives

| Category | Frequency | Notes |
|----------|-----------|-------|
| LLM hallucination | ~5% | Made-up variants |
| Wrong gene | ~2% | Variant from different gene |
| Misread notation | ~3% | c.123G>A vs c.132G>A |

### Validation Recommendations

For clinical or publication use:
1. **Spot-check high-impact variants** against source papers
2. **Verify penetrance counts** manually for key findings
3. **Cross-reference** with ClinVar/gnomAD
4. **Review key_quotes** field for supporting evidence

---

## Comparison with Other Tools

| Tool | Approach | Typical Recall | Notes |
|------|----------|----------------|-------|
| **GVF** | LLM extraction | 55-60% | Automated, scalable |
| Manual curation | Expert review | 90-95% | Gold standard |
| PubTator | NER-based | 30-40% | Fast but misses context |
| tmVar | ML-based | 40-50% | Good for mentions, not clinical data |

GVF's advantage: Extracts **structured clinical data** (penetrance, affected counts), not just variant mentions.

---

## Interpreting Results for Your Gene

### What to Expect

| Gene Type | Expected Recall | Reasoning |
|-----------|-----------------|-----------|
| Well-studied cardiac (KCNH2, SCN5A) | 55-65% | Best optimization |
| Other cardiac (KCNQ1, RYR2) | 45-55% | Good transfer |
| Non-cardiac (BRCA1, TP53) | 35-50% | Less tuning |
| Rare disease genes | 25-40% | Limited validation |

### Why Cardiac Genes Perform Better

1. **Keyword filters** optimized for cardiac terms
2. **Variant normalizer** tested on cardiac variants
3. **Gold standard** validation only for KCNH2
4. **LLM prompts** include cardiac-specific context

### For Non-Cardiac Genes

Consider:
- Adjusting keyword filters in `config/constants.py`
- Creating gene-specific alias dictionaries
- Building your own gold standard for validation

---

## Next Steps

- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) — Understanding extraction output
- [ARCHITECTURE.md](ARCHITECTURE.md) — How extraction works
- [QUICKSTART.md](QUICKSTART.md) — Run your first extraction

---

## References

For methodology background:
- HGVS variant nomenclature: [varnomen.hgvs.org](https://varnomen.hgvs.org/)
- ClinVar variant interpretation: [ncbi.nlm.nih.gov/clinvar](https://www.ncbi.nlm.nih.gov/clinvar/)
- ACMG variant classification guidelines
