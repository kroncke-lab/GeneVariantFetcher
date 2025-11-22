# Clinical Data Triage Assistant

A specialized tool for triaging scientific papers to identify those containing **original clinical data** about genetic variants.

## Overview

The Clinical Data Triage Assistant helps researchers quickly filter papers to find those with new patient-level clinical data, automatically rejecting review articles, animal studies, and other non-clinical research.

### What It Does

Given a paper's title and abstract, the triage assistant determines:
- **KEEP**: Paper contains original clinical case data (case reports, cohorts, clinical trials)
- **DROP**: Paper is a review, animal study, cell study only, or lacks clinical data

### Use Cases

- **Literature curation**: Filter large sets of papers from PubMed searches
- **Systematic reviews**: Identify primary data sources vs. reviews
- **Database building**: Find papers with extractable patient-level data
- **Research triage**: Quickly identify relevant clinical studies

## Features

- ✅ **Accurate classification** using LLM-based analysis
- ✅ **Detailed reasoning** for each decision
- ✅ **Confidence scores** (0.0-1.0) for each triage result
- ✅ **Batch processing** for multiple papers
- ✅ **Flexible input** (single papers, CSV files, interactive mode)
- ✅ **Cost-efficient** using gpt-4o-mini (~$0.0001 per paper)

## Installation

The triage assistant is already included in the GeneVariantFetcher package. Ensure you have the required dependencies:

```bash
pip install -e .
```

Set your OpenAI API key:

```bash
export AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key"
export AI_INTEGRATIONS_OPENAI_BASE_URL="https://api.openai.com/v1"
```

## Quick Start

### 1. Single Paper Triage

Triage a single paper by providing title and abstract:

```bash
python clinical_data_triage.py \
  --title "A novel SCN5A mutation in Brugada syndrome" \
  --abstract "We report a 28-year-old patient..." \
  --gene SCN5A
```

Output:
```json
{
  "decision": "KEEP",
  "reason": "New case report detected",
  "confidence": 0.95,
  "pmid": null,
  "model": "gpt-4o-mini",
  "tokens_used": 245
}
```

### 2. Batch Processing from CSV

Process multiple papers from a CSV file:

```bash
python clinical_data_triage.py \
  --input papers.csv \
  --gene BRCA1 \
  --output triage_results.json
```

**CSV Format:**
```csv
pmid,title,abstract,gene
12345678,"Title 1","Abstract 1",SCN5A
12345679,"Title 2","Abstract 2",KCNQ1
```

Optional columns: `gene` (overrides command-line gene)

### 3. Interactive Mode

Use interactive mode for manual triage:

```bash
python clinical_data_triage.py --interactive
```

You'll be prompted to enter:
1. Gene symbol
2. PMID (optional)
3. Paper title
4. Abstract (paste and press Enter twice)

## Triage Rules

### Papers to REJECT (DROP)

The triage assistant will **DROP** papers that are:

1. **Review articles** - Summaries of existing literature
2. **Meta-analyses** - Statistical pooling of studies (unless individual patient data is included)
3. **Animal studies** - Mouse, rat, zebrafish, or other model organisms only
4. **Cell studies only** - In vitro experiments without patient data
5. **Variant interpretation guidelines** - Standards without new clinical cases
6. **Computational studies** - Bioinformatics predictions without validation
7. **Basic science** - Molecular mechanisms without clinical phenotypes

### Papers to ACCEPT (KEEP)

The triage assistant will **KEEP** papers that contain:

1. **Case reports** - Even single patient descriptions
2. **Case series** - Multiple patients with clinical details
3. **Clinical cohort studies** - Patient populations with genetic variant data
4. **Functional studies + phenotype** - Lab studies that also describe patient clinical data
5. **Clinical trials** - Studies with patient-level genetic and phenotypic data
6. **Family studies** - Segregation analysis with clinical information

**Key Question**: Does this paper report NEW patient-level clinical data?

## Usage Examples

### Example 1: Quick Triage

```python
from filters import ClinicalDataTriageFilter

triage = ClinicalDataTriageFilter()

result = triage.triage(
    title="Novel BRCA1 mutation in a patient with early-onset breast cancer",
    abstract="A 32-year-old woman presented with breast cancer. Genetic testing revealed...",
    gene="BRCA1",
    pmid="12345678"
)

print(result['decision'])  # "KEEP"
print(result['reason'])    # "New case report detected"
print(result['confidence'])  # 0.92
```

### Example 2: Triage from PubMed Search

```python
from filters import ClinicalDataTriageFilter
from sourcer import PaperSourcer

# Fetch papers
sourcer = PaperSourcer(email="your@email.com")
pmids = sourcer.fetch_papers("SCN5A")

# Triage each paper
triage = ClinicalDataTriageFilter()
keep_papers = []

for pmid in pmids[:50]:  # First 50 papers
    paper = sourcer.fetch_paper_metadata(pmid)
    if paper:
        result = triage.triage_paper(paper, gene="SCN5A")
        if result['decision'] == 'KEEP':
            keep_papers.append(paper)

print(f"Kept {len(keep_papers)} out of 50 papers")
```

### Example 3: Integration with Existing Pipeline

```python
from pipeline import BiomedicalExtractionPipeline
from filters import ClinicalDataTriageFilter

# Add triage as a pre-filter
pipeline = BiomedicalExtractionPipeline(email="your@email.com")
triage = ClinicalDataTriageFilter()

# Custom processing with triage
papers = pipeline.sourcer.fetch_papers("KCNQ1")
triaged_papers = []

for pmid in papers:
    paper = pipeline.sourcer.fetch_paper_metadata(pmid)
    if paper:
        result = triage.triage_paper(paper)
        if result['decision'] == 'KEEP' and result['confidence'] > 0.7:
            triaged_papers.append(paper)

# Process only high-confidence clinical papers
for paper in triaged_papers:
    pipeline.process_paper(paper)
```

## Output Format

All triage results follow this JSON schema:

```json
{
  "decision": "KEEP" | "DROP",
  "reason": "Brief explanation of the decision",
  "confidence": 0.0-1.0,
  "pmid": "12345678" | null,
  "model": "gpt-4o-mini",
  "tokens_used": 245
}
```

**Fields:**
- `decision`: Either "KEEP" (has clinical data) or "DROP" (no clinical data)
- `reason`: Human-readable explanation (e.g., "Review article only", "New case report detected")
- `confidence`: Numerical confidence score (0.0 = no confidence, 1.0 = very confident)
- `pmid`: PubMed ID if provided
- `model`: LLM model used for classification
- `tokens_used`: Number of tokens consumed (for cost tracking)

## Command-Line Reference

```
usage: clinical_data_triage.py [-h] [--interactive | --input INPUT]
                               [--title TITLE] [--abstract ABSTRACT]
                               [--pmid PMID] [--gene GENE]
                               [--output OUTPUT] [--model MODEL]
                               [--verbose]

Clinical Data Triage Tool

optional arguments:
  -h, --help            Show help message
  --interactive, -i     Run in interactive mode
  --input INPUT         Input CSV file (columns: pmid, title, abstract)
  --title TITLE         Paper title (for single paper mode)
  --abstract ABSTRACT   Paper abstract (for single paper mode)
  --pmid PMID          PubMed ID (optional)
  --gene GENE, -g GENE Gene symbol (default: "the gene of interest")
  --output OUTPUT, -o   Output JSON file
  --model MODEL, -m     LLM model (default: gpt-4o-mini)
  --verbose, -v         Enable verbose logging
```

## Testing

Run the test suite to verify functionality:

```bash
python test_triage.py
```

**Warning**: Tests make real API calls and consume credits (~$0.01-0.02 total).

The test suite includes 9 test cases covering:
- Case reports (should KEEP)
- Review articles (should DROP)
- Animal studies (should DROP)
- Cell studies (should DROP)
- Clinical cohorts (should KEEP)
- Functional studies with phenotype (should KEEP)
- Meta-analyses (should DROP)
- Clinical guidelines (should DROP)
- Case series (should KEEP)

## Cost Analysis

Using gpt-4o-mini (default):
- **Per paper**: ~$0.0001 (0.01 cents)
- **100 papers**: ~$0.01 (1 cent)
- **1,000 papers**: ~$0.10 (10 cents)
- **10,000 papers**: ~$1.00 (1 dollar)

Average tokens per triage: ~300 tokens
- Input: ~200 tokens (title + abstract)
- Output: ~100 tokens (JSON response)

## Customization

### Change the LLM Model

Use a different model for better accuracy or lower cost:

```python
# Higher accuracy (more expensive)
triage = ClinicalDataTriageFilter(model="gpt-4o")

# Lower cost (less accurate)
triage = ClinicalDataTriageFilter(model="gpt-3.5-turbo")
```

### Adjust Temperature

Control randomness (lower = more deterministic):

```python
triage = ClinicalDataTriageFilter(
    model="gpt-4o-mini",
    temperature=0.0  # Fully deterministic
)
```

### Filter by Confidence

Only keep high-confidence results:

```python
result = triage.triage(title, abstract, gene)

if result['decision'] == 'KEEP' and result['confidence'] > 0.8:
    # Process high-confidence papers only
    process_paper(paper)
```

## Troubleshooting

### Issue: All papers are being dropped

**Solution**: Check if abstracts are actually being passed. Empty abstracts will be dropped automatically.

### Issue: API errors

**Solution**: Verify your OpenAI API key is set correctly:
```bash
echo $AI_INTEGRATIONS_OPENAI_API_KEY
```

### Issue: Incorrect classifications

**Solution**:
1. Check if the abstract contains enough information
2. Try a more powerful model (gpt-4o instead of gpt-4o-mini)
3. Review the `reason` field to understand the LLM's reasoning

### Issue: Too many borderline cases

**Solution**: Use confidence scores to filter:
```python
# Only trust high-confidence decisions
if result['confidence'] > 0.75:
    use_decision = result['decision']
else:
    # Manual review needed
    flag_for_review(paper)
```

## Integration with Pipeline

The triage filter can be used as part of the existing extraction pipeline:

```python
from pipeline import BiomedicalExtractionPipeline
from filters import ClinicalDataTriageFilter, KeywordFilter

# Create pipeline with triage as an additional filter
pipeline = BiomedicalExtractionPipeline(
    email="your@email.com",
    enable_tier1=True,  # Keyword filter
    enable_tier2=True   # Intern filter
)

# Add triage step
triage = ClinicalDataTriageFilter()

# Custom workflow
papers = pipeline.sourcer.fetch_papers("SCN5A")
for pmid in papers:
    paper = pipeline.sourcer.fetch_paper_metadata(pmid)

    # Triage first (cheapest filter)
    triage_result = triage.triage_paper(paper)
    if triage_result['decision'] == 'DROP':
        continue  # Skip non-clinical papers

    # Then run through the full pipeline
    pipeline.process_paper(paper)
```

## Best Practices

1. **Triage early**: Run triage before expensive extraction to save costs
2. **Use confidence scores**: Filter by confidence to reduce manual review
3. **Batch processing**: Process papers in batches for efficiency
4. **Save results**: Always save triage results for reproducibility
5. **Manual review**: For borderline cases (confidence < 0.6), consider manual review
6. **Gene-specific**: Always provide the gene symbol for context

## API Reference

### `ClinicalDataTriageFilter`

Main class for clinical data triage.

**Constructor:**
```python
ClinicalDataTriageFilter(
    model: str = "gpt-4o-mini",
    temperature: float = 0.1,
    max_tokens: int = 200
)
```

**Methods:**

#### `triage(title, abstract, gene, pmid=None)`
Triage a paper based on title and abstract.

**Parameters:**
- `title` (str): Paper title
- `abstract` (str): Paper abstract or introduction
- `gene` (str): Gene symbol being studied
- `pmid` (str, optional): PubMed ID for logging

**Returns:**
- `dict`: Triage result with decision, reason, confidence, etc.

#### `triage_paper(paper, gene=None)`
Triage a Paper object.

**Parameters:**
- `paper` (Paper): Paper object with title and abstract
- `gene` (str, optional): Gene symbol (uses paper.gene_symbol if not provided)

**Returns:**
- `dict`: Triage result

## Contributing

To add new triage rules or improve classification:

1. Edit the `TRIAGE_PROMPT` in `filters.py`
2. Add test cases in `test_triage.py`
3. Run tests to verify changes
4. Submit a pull request

## License

Part of the GeneVariantFetcher project. See main LICENSE file.

## Citation

If you use this tool in your research, please cite:

```
[Citation information to be added]
```

## Support

For issues or questions:
- Check this README and troubleshooting section
- Review test cases in `test_triage.py`
- File an issue on the GitHub repository
