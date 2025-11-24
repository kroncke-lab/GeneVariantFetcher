<!-- ec8b9a29-af6d-411c-9820-f71bd0498732 640c8514-e34d-4ab0-9ca9-1ec44b8fdc3f -->
# Penetrance Extraction Enhancement Plan

## Goal

Enhance the existing variant extraction pipeline to extract penetrance estimates (e.g., "observed in 10 people, 4 with disease") for calibrating disease prediction models. Support both individual-level and aggregated variant-level penetrance data with age tracking.

## Current State

- Extraction system captures basic variant info in `extractor.py`
- Current schema has `patients.count` but lacks detailed penetrance breakdown
- Output format: `TTR_PMID_FULL.json` per paper
- Missing: individual-level records, affected/unaffected counts, age-at-onset, calculated penetrance percentages

## Implementation Plan

### 1. Enhance Extraction Prompt (`extractor.py`)

**File**: `extractor.py` (lines 24-115)

**Changes**:

- Add explicit instructions to extract individual-level person records
- Request penetrance data: total carriers, affected count, unaffected count
- Extract age-at-onset and age-at-evaluation for age-dependent penetrance
- Extract disease status per individual when linked to variants
- Capture cohort study data (e.g., "10 carriers, 4 affected, 6 unaffected")

**Key additions to prompt**:

- Individual-level extraction: proband, case, patient, family member (II-1, P1, etc.)
- Penetrance tracking: distinguish affected vs unaffected carriers
- Age tracking: age-at-onset, age-at-evaluation, age-at-diagnosis
- Cohort data: extract study-level statistics when provided

### 2. Update JSON Schema

**File**: `extractor.py` (EXTRACTION_PROMPT output format)

**New structure**:

```json
{
  "variants": [
    {
      // ... existing fields ...
      "penetrance_data": {
        "total_carriers_observed": integer,
        "affected_count": integer,
        "unaffected_count": integer,
        "uncertain_count": integer,
        "penetrance_percentage": float,
        "age_dependent_penetrance": [
          {
            "age_range": "string (e.g., '40-50 years')",
            "penetrance": float,
            "carriers_in_range": integer
          }
        ]
      },
      "individual_records": [
        {
          "individual_id": "string",
          "age_at_evaluation": integer,
          "age_at_onset": integer,
          "age_at_diagnosis": integer,
          "sex": "string",
          "affected_status": "affected/unaffected/uncertain",
          "phenotype_details": "string",
          "evidence_sentence": "string"
        }
      ]
    }
  ]
}
```

### 3. Create Post-Processing Aggregation Script

**New file**: `penetrance_aggregator.py`

**Purpose**: Aggregate penetrance data across multiple papers for the same variant

**Features**:

- Merge individual records from multiple papers
- Calculate aggregate penetrance statistics per variant
- Handle age-dependent penetrance aggregation
- Output variant-level summary for model calibration

**Output format**: `{GENE}_penetrance_summary.json`:

```json
{
  "gene_symbol": "TTR",
  "variants": [
    {
      "protein_notation": "p.Val30Met",
      "aggregated_penetrance": {
        "total_carriers": 150,
        "affected": 60,
        "unaffected": 85,
        "uncertain": 5,
        "penetrance_percentage": 40.0,
        "sources": ["PMID1", "PMID2", ...]
      },
      "age_dependent_penetrance": [...],
      "individual_records_count": 150
    }
  ]
}
```

### 4. Update Workflow Integration

**File**: `example_automated_workflow.py`

**Changes**:

- Add post-processing step after extractions complete
- Generate penetrance summary for each gene run
- Include penetrance statistics in workflow summary

### 5. Validation & Quality Checks

**Add validation**:

- Verify penetrance calculations: affected + unaffected + uncertain = total_carriers
- Flag inconsistencies (e.g., affected > total_carriers)
- Validate age data (age-at-onset ≤ age-at-evaluation)
- Check for duplicate individual records across papers

## Files to Modify

1. **`extractor.py`**:

   - Lines 24-115: Enhance EXTRACTION_PROMPT
   - Update JSON schema in prompt

2. **`example_automated_workflow.py`**:

   - Lines 179-226: Add penetrance aggregation step

3. **New file: `penetrance_aggregator.py`**:

   - Individual record merging
   - Variant-level aggregation
   - Age-dependent penetrance calculations

4. **`models.py`** (optional):

   - Add Pydantic models for penetrance data structures

## Testing Strategy

1. Test with existing TTR extractions
2. Validate against known penetrance data from literature
3. Test age-dependent penetrance calculations
4. Test aggregation across multiple papers for same variant

## Backward Compatibility

- Existing extraction outputs remain valid
- New fields are additive (optional fields)
- Can run new extraction alongside old extractions

### To-dos

- [ ] Enhance EXTRACTION_PROMPT in extractor.py to extract individual-level records and penetrance data (affected/unaffected counts, age-at-onset, age-at-evaluation)
- [ ] Update JSON output schema in extraction prompt to include penetrance_data and individual_records fields with both counts and percentages
- [ ] Create penetrance_aggregator.py script to aggregate individual records and calculate variant-level penetrance statistics across multiple papers
- [ ] Integrate penetrance aggregation into example_automated_workflow.py to generate gene-level penetrance summaries
- [ ] Add validation logic to check penetrance data consistency (counts add up, ages are valid, no duplicates)
- [ ] Test enhanced extraction on TTR papers and validate penetrance calculations against known literature values